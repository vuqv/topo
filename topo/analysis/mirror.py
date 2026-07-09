#!/usr/bin/env python3
"""
Mirror-image detection for CG (Cα Gō) trajectories.

A Cα structure-based model can fold into a **global mirror image** of the native
structure: an artifact that preserves most native contacts, so it looks "folded"
by the fraction-of-native-contacts Q yet is chirally wrong. Q cannot distinguish
it; the two reference-based metrics here can, and together give a robust per-frame
mirror flag.

Metric 1 -- mirror-sensitive RMSD
    Proper-rotation (det = +1) Kabsch fit of each frame onto (a) the native Cα
    reference and (b) an explicitly reflected copy of it. A mirror fold fits the
    reflected reference better: ``RMSD_reflected < RMSD_native``. The determinant
    correction is mandatory -- a fit that allows an improper rotation would absorb
    the reflection and defeat the test.

Metric 2 -- chirality fraction-agreement K
    A per-window local chirality (normalized signed triple product) on the
    polyline of secondary-structure segment endpoints. It flips sign under
    mirroring, so ``K`` = fraction of windows whose sign matches native falls to 0
    for a mirror and rises to 1 for a native-handed fold.

Combined per-frame rule (thresholds tunable):
    ``mirror_image = (Q > 0.5) AND (K < 0.3) AND (RMSD_reflected / RMSD_native < 0.9)``

Both metrics are Cα-only, reference-based post-processing, in the same spirit as
``topo.analysis.native_contacts`` (per-frame metric vs the all-atom reference,
CLI + importable library functions). The all-atom ``pdb_file`` is the single
source of truth for the native Cα geometry; secondary-structure endpoints come
from topo's own STRIDE output (or an explicit segment file).

CLI::

    python -m topo.analysis.mirror -r ref.pdb -f cg.dcd -o mirror.csv \\
        [-p cg.psf] [-s stride.out | --segments seg.txt] [--reflect-axis x] \\
        [-b start] [-e end]

Residue numbering: STRIDE / segment files use 1-based residue numbers that map to
0-based bead indices internally (same convention as native_contacts.py).
"""

from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd
import MDAnalysis as mda

# Reuse the native-contacts machinery (single source of truth for Q and for the
# reference Cα geometry) and topo's chain / residue-mapping helpers.
from .native_contacts import (
    load_universe,
    reference_residue_geometry,
    build_native_contacts,
    fraction_native_contacts,
)
from ..utils.nonbonded import _norm_chain, get_residue_mapping


# --------------------------------------------------------------------------- #
# Metric 1 -- mirror-sensitive RMSD
# --------------------------------------------------------------------------- #
def center_coords(coords):
    """Center coordinates at their centroid."""
    coords = np.asarray(coords, dtype=np.float64)
    return coords - coords.mean(axis=0, keepdims=True)


def reflect_coords(coords, axis="x"):
    """Reflect coordinates across one Cartesian axis (a det = -1 improper map)."""
    reflected = np.asarray(coords, dtype=np.float64).copy()
    reflected[:, {"x": 0, "y": 1, "z": 2}[axis]] *= -1.0
    return reflected


def kabsch_proper(mobile, target):
    """Rotate centered ``mobile`` onto centered ``target`` (proper rotation only).

    Both inputs must already be centered. The ``diag(1, 1, d)`` determinant
    correction (``d = sign(det(V @ Wt))``) forces ``det(rotation) = +1`` so the
    fit can never silently absorb a reflection -- the crux of mirror detection.

    Returns the rotated ``mobile`` coordinates (not translated).
    """
    cov = mobile.T @ target
    V, _S, Wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    rot = V @ D @ Wt
    return mobile @ rot


def _rmsd(a, b):
    """RMSD between two aligned coordinate sets of shape (N, 3)."""
    diff = a - b
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def aligned_rmsd(mobile, target_centered):
    """RMSD after centering ``mobile`` and a proper-rotation Kabsch fit onto
    ``target_centered`` (which must already be centered)."""
    mobile_c = center_coords(mobile)
    mobile_rot = kabsch_proper(mobile_c, target_centered)
    return _rmsd(mobile_rot, target_centered)


def mirror_rmsd(ca_frame, ca_native, reflect_axis="x"):
    """Per-frame RMSD to the native and to the reflected-native reference.

    Parameters
    ----------
    ca_frame : (N, 3) array
        Cα coordinates of the frame.
    ca_native : (N, 3) array
        Cα coordinates of the native reference (any centering; centered here).
    reflect_axis : {"x", "y", "z"}
        Axis negated to build the reflected reference. Any single axis works --
        the proper-rotation fit re-orients freely, so the choice does not affect
        the result.

    Returns
    -------
    (rmsd_native, rmsd_reflected) : tuple of float
        ``rmsd_reflected < rmsd_native`` indicates a mirror-like fold.
    """
    native_c = center_coords(ca_native)
    reflected_c = reflect_coords(native_c, reflect_axis)
    return aligned_rmsd(ca_frame, native_c), aligned_rmsd(ca_frame, reflected_c)


# --------------------------------------------------------------------------- #
# Metric 2 -- chirality fraction-agreement K
# --------------------------------------------------------------------------- #
def local_chirality(points, eps=1e-8):
    """Normalized local chirality along an ordered point sequence.

    For bond vectors ``v_k = P[k+1] - P[k]`` and each window ``i``::

        chi_i = ((v_i x v_{i+1}) . v_{i+2}) / (|v_i| |v_{i+1}| |v_{i+2}|)

    a signed pseudo-scalar in [-1, 1] (handedness). It flips sign under a
    reflection and is invariant to proper rotation and translation.

    Parameters
    ----------
    points : (M, 3) or (n_frames, M, 3) array
        Ordered points (the SS segment-endpoint polyline). A single (M, 3) frame
        is promoted internally and returned as (M - 3,).
    eps : float
        Added to the denominator to guard zero-length bond vectors (coincident
        endpoints).

    Returns
    -------
    chi : (M - 3,) or (n_frames, M - 3) array
    """
    pts = np.asarray(points, dtype=np.float64)
    single = pts.ndim == 2
    if single:
        pts = pts[None]

    v1 = pts[:, 1:-2] - pts[:, :-3]
    v2 = pts[:, 2:-1] - pts[:, 1:-2]
    v3 = pts[:, 3:] - pts[:, 2:-1]

    cross = np.cross(v1, v2, axis=2)
    dot = np.einsum("ijk,ijk->ij", cross, v3)
    denom = (np.linalg.norm(v1, axis=2)
             * np.linalg.norm(v2, axis=2)
             * np.linalg.norm(v3, axis=2)) + eps
    chi = dot / denom
    return chi[0] if single else chi


def chirality_agreement(chi_frames, chi_ref):
    """Fraction agreement K per frame: fraction of windows whose chirality sign
    matches the native reference. ``K -> 1`` native-handed, ``K -> 0`` mirror.

    ``chi_frames`` may be (n_windows,) for a single frame or (n_frames, n_windows).
    Returns a (n_frames,) array.
    """
    cf = np.asarray(chi_frames)
    if cf.ndim == 1:
        cf = cf[None]
    return np.mean(np.sign(cf) == np.sign(chi_ref), axis=1)


# --------------------------------------------------------------------------- #
# Secondary-structure segment endpoints
# --------------------------------------------------------------------------- #
def run_stride(pdb_file, out_dir=None, timeout=60):
    """Run ``stride -h`` on ``pdb_file`` and cache the output to
    ``{stem}_stride.dat``; return that path.

    Mirrors the STRIDE-runner pattern in ``nonbonded.build_nonbonded_interaction``:
    locate the executable, run it, and validate by output content (the ``LOC``
    secondary-structure records) rather than the exit code -- some STRIDE builds
    return non-zero even on success.
    """
    stride_exe = shutil.which("stride")
    if stride_exe is None:
        raise RuntimeError(
            "No STRIDE output supplied and the 'stride' program was not found. "
            "Provide --segments or -s <stride output>, or install STRIDE on PATH."
        )
    stem = os.path.splitext(os.path.basename(pdb_file))[0]
    out_dir = Path(out_dir) if out_dir is not None else Path(pdb_file).resolve().parent
    stride_path = out_dir / f"{stem}_stride.dat"
    print(f"Running STRIDE (stride -h {pdb_file} -> {stride_path}).")
    result = subprocess.run(
        [stride_exe, "-h", str(pdb_file)],
        capture_output=True, text=True, timeout=timeout,
    )
    if "LOC" not in result.stdout:
        raise RuntimeError(
            f"STRIDE produced no secondary-structure (LOC) records for {pdb_file} "
            f"(exit code {result.returncode}). stderr: {result.stderr}"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    stride_path.write_text(result.stdout)
    return str(stride_path)


def segment_endpoints_from_stride(stride_output_file, key_to_index, min_len=4):
    """0-based Cα bead indices of the head & tail of each STRIDE helix/strand
    element (length >= ``min_len``), flattened across elements and sorted -> the
    ordered point sequence for the chirality polyline.

    Reuses topo's ``(chain, resid) -> 0-based index`` map so the endpoints line up
    with the CG bead order used everywhere else (and multi-chain structures resolve).

    STRIDE LOC record (fixed columns), e.g.::

        LOC  AlphaHelix   MET     1 A      MET     13 A
        LOC  Strand       SER    20 A      PRO     27 A
    """
    endpoints = []
    with open(stride_output_file) as fh:
        for line in fh:
            if not line.startswith("LOC "):
                continue
            sec_name = line[5:17].strip()          # 'AlphaHelix'|'310Helix'|'Strand'|'Turn...'
            if "Helix" not in sec_name and "Strand" not in sec_name:
                continue                           # keep helices + strands only
            chain = _norm_chain(line[28])
            start_resnum = int(line[21:27])
            end_resnum = int(line[39:45])
            if (end_resnum - start_resnum + 1) < min_len:
                continue                           # drop very short elements (e.g. 3-res 3-10)
            for resnum in (start_resnum, end_resnum):
                key = (chain, resnum)
                if key not in key_to_index:
                    raise KeyError(
                        f"STRIDE residue {key} not found in structure "
                        f"(check chain/numbering consistency)"
                    )
                endpoints.append(key_to_index[key])  # 0-based bead index
    if len(endpoints) < 4:
        raise ValueError(
            "fewer than 2 SS segments (need >= 4 endpoints) -- chirality undefined"
        )
    return sorted(endpoints)


def segment_endpoints_from_file(segments_file):
    """0-based Cα bead indices from an explicit ``<label> <start> <end>`` segment
    file (1-based residue numbers, one element per line), flattened and sorted.

    Residue numbers are taken to run 1..n_residues in CG bead order (single-chain
    convention, as in the standalone prototypes); use the STRIDE path for
    multi-chain structures.
    """
    endpoints = []
    with open(segments_file) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split()
            endpoints.append(int(parts[1]) - 1)    # 1-based -> 0-based
            endpoints.append(int(parts[2]) - 1)
    if len(endpoints) < 4:
        raise ValueError(
            f"{segments_file}: fewer than 2 SS segments (need >= 4 endpoints) "
            f"-- chirality undefined"
        )
    return sorted(endpoints)


def resolve_segment_endpoints(reference, u_ref, segments=None, stride_output=None,
                              min_len=4, out_dir=None):
    """Resolve SS segment endpoints with three-way precedence:
    explicit ``segments`` file -> precomputed ``stride_output`` -> auto-run STRIDE
    on the reference PDB. Returns a sorted list of 0-based bead indices.
    """
    if segments is not None:
        return segment_endpoints_from_file(segments)
    if stride_output is None:
        stride_output = run_stride(reference, out_dir=out_dir)
    key_to_index, _idx_to_name, _n = get_residue_mapping(u_ref)
    return segment_endpoints_from_stride(stride_output, key_to_index, min_len=min_len)


# --------------------------------------------------------------------------- #
# Combined classification
# --------------------------------------------------------------------------- #
def classify_mirror(Q, K, rmsd_native, rmsd_reflected,
                    q_fold=0.5, k_thresh=0.3, rmsd_ratio_thresh=0.9, eps=1e-12):
    """Per-frame mirror classification.

    Returns
    -------
    rmsd_ratio : ndarray
        ``rmsd_reflected / rmsd_native`` (denominator floored by ``eps``).
    is_mirror : bool ndarray
        ``(Q > q_fold) AND (K < k_thresh) AND (rmsd_ratio < rmsd_ratio_thresh)``.
    """
    Q = np.asarray(Q, dtype=float)
    K = np.asarray(K, dtype=float)
    rmsd_native = np.asarray(rmsd_native, dtype=float)
    rmsd_reflected = np.asarray(rmsd_reflected, dtype=float)
    rmsd_ratio = rmsd_reflected / np.maximum(rmsd_native, eps)
    is_mirror = (Q > q_fold) & (K < k_thresh) & (rmsd_ratio < rmsd_ratio_thresh)
    return rmsd_ratio, is_mirror


# --------------------------------------------------------------------------- #
# Trajectory topology
# --------------------------------------------------------------------------- #
def build_cg_universe(reference, dcd, psf=None):
    """Load the CG trajectory. If ``psf`` is given, use it as the topology;
    otherwise derive a Cα-only topology straight from the all-atom ``reference``
    (a DCD carries coordinates only, so MDAnalysis needs a matching topology).
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*DCDReader currently makes independent timesteps.*",
            category=DeprecationWarning,
        )
        if psf is not None:
            return mda.Universe(psf, dcd)
        u_ref = load_universe(reference)
        ca = u_ref.select_atoms("protein and name CA")
        return mda.Merge(ca).load_new(dcd)


def resolve_frame_range(start, end, n_frames):
    """Python-style [start, end) with negative-index support (``-1`` end = last
    frame)."""
    if start < 0:
        start = n_frames + start
    start = max(0, start)
    if end == -1:
        end = n_frames
    elif end < 0:
        end = n_frames + end
    end = min(n_frames, end)
    if start >= end:
        raise ValueError(f"Invalid frame range: start={start}, end={end}")
    return start, end


# --------------------------------------------------------------------------- #
# CLI / driver
# --------------------------------------------------------------------------- #
def parse_args():
    """Parse command-line arguments for the mirror-detection driver."""
    p = argparse.ArgumentParser(
        prog="python -m topo.analysis.mirror",
        description="Per-frame mirror-image detection (RMSD-to-reflected + "
                    "chirality K + Q) for a CG trajectory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-r", "--reference", required=True,
                   help="All-atom reference PDB (native Cα + STRIDE source).")
    p.add_argument("-f", "--dcd", required=True, help="CG DCD trajectory.")
    p.add_argument("-o", "--output", required=True, help="Output CSV path.")
    p.add_argument("-p", "--psf", default=None,
                   help="CG PSF topology. If omitted, a Cα-only topology is "
                        "derived from the reference PDB.")
    p.add_argument("-s", "--stride", default=None,
                   help="Precomputed STRIDE output (LOC records). If omitted and "
                        "no --segments, STRIDE is run on the reference.")
    p.add_argument("--segments", default=None,
                   help="Explicit '<label> <start> <end>' segment file (1-based). "
                        "Takes precedence over STRIDE.")
    p.add_argument("--reflect-axis", choices=["x", "y", "z"], default="x",
                   help="Axis negated to build the reflected reference.")
    p.add_argument("--min-seg-len", type=int, default=4,
                   help="Minimum SS element length (residues) kept for chirality.")
    p.add_argument("-b", "--start-frame", type=int, default=0,
                   help="First frame (inclusive; negative indexing supported).")
    p.add_argument("-e", "--end-frame", type=int, default=-1,
                   help="Last frame (exclusive; -1 = last frame).")
    # Native-contact (Q) definition -- kept identical to native_contacts.py.
    p.add_argument("--cutoff", type=float, default=4.5,
                   help="Heavy-atom distance defining a native contact (Å).")
    p.add_argument("--local-separation", type=int, default=3,
                   help="Minimum sequence separation for Q (abs(i - j) > this).")
    p.add_argument("--tolerance", type=float, default=1.2,
                   help="Cα-distance stretch factor for a 'formed' contact.")
    # Classification thresholds.
    p.add_argument("--q-fold", type=float, default=0.5,
                   help="Q above this counts as folded.")
    p.add_argument("--k-thresh", type=float, default=0.3,
                   help="K below this counts as chirally inverted.")
    p.add_argument("--rmsd-ratio", type=float, default=0.9,
                   help="RMSD_reflected/RMSD_native below this counts as mirror-fit.")
    p.add_argument("--summary-last-frames", type=int, default=133,
                   help="Trailing frames over which the per-trajectory mirror "
                        "verdict is taken (the early, unfolded part is excluded). "
                        "Specified in frames because the DCD time step is not "
                        "reliable.")
    return p.parse_args()


def main():
    """Score per-frame mirror metrics (Q, K, RMSD_native, RMSD_reflected) for a
    CG trajectory, classify each frame, and write the combined CSV."""
    args = parse_args()

    for label, path in (("reference", args.reference), ("dcd", args.dcd)):
        if not Path(path).exists():
            raise FileNotFoundError(f"Missing {label} file: {path}")
    if args.psf is not None and not Path(args.psf).exists():
        raise FileNotFoundError(f"Missing psf file: {args.psf}")

    # --- reference Cα geometry (single source of truth) ---
    u_ref = load_universe(args.reference)
    ca_native, resindex_to_pos = reference_residue_geometry(u_ref)
    n_res = ca_native.shape[0]

    # --- native contacts for Q (whole molecule) ---
    heavy = u_ref.select_atoms("protein and not name H*")
    heavy_positions = heavy.positions.copy()
    heavy_res = np.fromiter(
        (resindex_to_pos[ri] for ri in heavy.resindices), dtype=int,
        count=heavy.n_atoms,
    )
    # Whole-molecule Q: every heavy-atom contact (|i-j| > local_separation),
    # NOT restricted to secondary-structure residue pairs. A prior pipeline may
    # define Q only over SS-SS pairs; that yields different absolute Q values.
    pairs, dnat = build_native_contacts(None, None, heavy_positions, heavy_res,
                                        ca_native, args.cutoff, args.local_separation)
    print(f"=== Native contacts (whole protein, not SS-restricted): "
          f"{pairs.shape[0]} (heavy-atom <= {args.cutoff} Å, "
          f"|i-j| > {args.local_separation}) ===")

    # --- SS segment endpoints for chirality ---
    endpoints = resolve_segment_endpoints(
        args.reference, u_ref, segments=args.segments, stride_output=args.stride,
        min_len=args.min_seg_len, out_dir=Path(args.output).resolve().parent,
    )
    endpoints = np.asarray(endpoints, dtype=int)
    if endpoints.max() >= n_res:
        raise ValueError(
            f"Segment endpoint index {int(endpoints.max())} out of range for "
            f"{n_res} residues (check numbering)."
        )
    chi_ref = local_chirality(ca_native[endpoints])
    print(f"=== SS endpoints: {endpoints.size} points "
          f"({chi_ref.size} chirality windows) ===")

    # --- precompute centered native + reflected reference for RMSD ---
    native_c = center_coords(ca_native)
    reflected_c = reflect_coords(native_c, args.reflect_axis)

    # --- CG trajectory ---
    u_cg = build_cg_universe(args.reference, args.dcd, psf=args.psf)
    n_beads = u_cg.atoms.n_atoms
    if n_beads != n_res:
        raise ValueError(
            f"Trajectory has {n_beads} beads but reference has {n_res} residues; "
            f"the CG model must be one Cα bead per residue (atom-count mismatch)."
        )

    n_frames_total = u_cg.trajectory.n_frames
    start, end = resolve_frame_range(args.start_frame, args.end_frame, n_frames_total)
    n_frames = end - start
    print(f"=== Scoring {n_frames} frames (of {n_frames_total}) ===")

    cg_atoms = u_cg.atoms
    rmsd_native = np.empty(n_frames, dtype=float)
    rmsd_reflected = np.empty(n_frames, dtype=float)
    Q = np.empty(n_frames, dtype=float)
    ep_traj = np.empty((n_frames, endpoints.size, 3), dtype=np.float64)

    for k, _ts in enumerate(u_cg.trajectory[start:end]):
        pos = cg_atoms.positions.astype(np.float64)
        rmsd_native[k] = aligned_rmsd(pos, native_c)
        rmsd_reflected[k] = aligned_rmsd(pos, reflected_c)
        Q[k] = fraction_native_contacts(pos, pairs, dnat, args.tolerance)
        ep_traj[k] = pos[endpoints]

    chi_frames = local_chirality(ep_traj)
    K = chirality_agreement(chi_frames, chi_ref)

    rmsd_ratio, is_mirror = classify_mirror(
        Q, K, rmsd_native, rmsd_reflected,
        q_fold=args.q_fold, k_thresh=args.k_thresh, rmsd_ratio_thresh=args.rmsd_ratio,
    )

    # --- write combined CSV ---
    df = pd.DataFrame({
        "Frame": np.arange(start, end, dtype=int),
        "Q": Q,
        "K": K,
        "RMSD_native": rmsd_native,
        "RMSD_reflected": rmsd_reflected,
        "RMSD_ratio": rmsd_ratio,
        "is_mirror": is_mirror,
    })
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    # --- summary: per-frame count (all frames) + trailing-window verdict ---
    # Per-frame is_mirror in the CSV is the authoritative output. The
    # per-trajectory verdict is taken over the LAST --summary-last-frames frames,
    # because the early trajectory is typically unfolded. The window is in frames
    # (not ns) because the DCD time step is not trusted. The verdict thresholds
    # the tail MEANS of each signal (matches the QGKR.npz convention).
    n_mirror = int(np.count_nonzero(is_mirror))
    n_win = max(1, min(args.summary_last_frames, n_frames))
    tail = slice(n_frames - n_win, n_frames)
    q_tail = float(np.nanmean(Q[tail]))
    k_tail = float(np.nanmean(K[tail]))
    ratio_tail = float(np.nanmean(rmsd_ratio[tail]))
    n_mirror_tail = int(np.count_nonzero(is_mirror[tail]))
    trajectory_is_mirror = (q_tail > args.q_fold and k_tail < args.k_thresh
                            and ratio_tail < args.rmsd_ratio)
    win_frames = np.arange(start, end)[tail]

    print("\n=== Summary ===")
    print(f"  mirror frames (all)    : {n_mirror}/{n_frames} "
          f"({n_mirror / n_frames:.1%})")
    print(f"\n  --- verdict over last {n_win} frames "
          f"(Frame {int(win_frames[0])}..{int(win_frames[-1])}) ---")
    print(f"  mean Q (tail)          : {q_tail:.4f}")
    print(f"  mean K (tail)          : {k_tail:.4f}")
    print(f"  mean RMSD_ratio (tail) : {ratio_tail:.4f}")
    print(f"  mirror frames (tail)   : {n_mirror_tail}/{n_win} "
          f"({n_mirror_tail / n_win:.1%})")
    print(f"  TRAJECTORY MIRROR      : {trajectory_is_mirror}  "
          f"[(Q>{args.q_fold}) & (K<{args.k_thresh}) & (ratio<{args.rmsd_ratio})]")
    print(f"\nSaved: {out_path.resolve()}")


if __name__ == "__main__":
    main()
