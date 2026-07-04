#!/usr/bin/env python3
"""
Stage 3 of the ribosome-prep pipeline: **orient** the assembled subunit.

Apply a single rigid-body transform (translation + proper rotation, **no scaling,
no reflection**) that puts the exit tunnel on a known axis, so stages 4-5
(coarse-grain + truncate) can use the X-axis directly as the tunnel centre line
and ``d = sqrt(y^2 + z^2)`` as the radial distance.

Landmarks (from the organism config's ``[orient]`` section)
-----------------------------------------------------------
* ``P_PTC``  -- the PTC atom            (E. coli ``23S:2602:N6``)
* ``P_exit`` -- the tunnel-exit atom    (E. coli ``L24:51:N``)
* ``P_A``    -- centroid of the A-site tRNA 3'-acceptor ribose ring (``AtR:last:ribose``)
* ``P_P``    -- centroid of the P-site tRNA 3'-acceptor ribose ring (``PtR:last:ribose``)

A landmark is written ``segid:resid:atom`` for a single atom, or
``segid:resid:ribose`` for the centroid of the ribose carbons ``C1'..C5'`` (the
same atom set ``cg_ribosome.py`` uses for its ``R`` bead, so the
landmark and the CG bead agree). ``resid`` may be the literal ``last`` -- the
tRNA's 3'-terminal acceptor -- so the tail landmark is numbering-robust (it is
residue 76 for a standard 76-nt tRNA but 77 for a 77-nt tRNA). ``gen_subunit``
additionally renumbers tRNA acceptors to 76, so ``last`` and ``76`` coincide in
the assembled models.

Transform, applied in order
---------------------------
1. **Translate** so ``P_PTC -> (0, 0, 0)``.
2. **Rotate** so the tunnel vector ``v_x = P_exit - P_PTC`` lies along **+x**.
3. **Rotate about the now-fixed x-axis** so the tRNA-tail vector
   ``v_y = P_P - P_A`` has its yz-projection pointing along **+y** (the
   x-component of ``v_y`` is untouched, so step 2 is preserved). z follows from
   right-handedness.

The rotation is built explicitly and checked (``det(R) = +1``, ``R^T R = I``)
before it is applied; after transforming, the code **asserts** ``|P_PTC| ~ 0``,
``v_x`` is ``(+, 0, 0)`` and ``v_y`` is ``(~0, +, ~0)`` within tolerance and
fails loudly otherwise.

Usage
-----
    python fix_orient.py -c configs/<organism>.ini
"""
from __future__ import annotations

import argparse
import configparser
import os
import sys
from collections import OrderedDict

import numpy as np

RIBOSE_RING = ["C1'", "C2'", "C3'", "C4'", "C5'"]   # matches cg_ribosome's R bead


# ---------------------------------------------------------------------------
# PDB I/O (minimal, fixed-column; preserves segID)
# ---------------------------------------------------------------------------
def read_pdb_atoms(path):
    """Read ATOM/HETATM records. Returns ``(lines, coords)``.

    ``lines`` is the raw record list (order preserved); ``coords`` is an
    ``(N, 3)`` float array aligned with the ATOM records only. Non-ATOM lines
    are kept in ``lines`` with a ``None`` coordinate slot so the file can be
    rewritten verbatim.
    """
    records = []          # (line, atom_index or None)
    coords = []
    idx = {}              # (segid, resid, atomname) -> coord row (first hit)
    for line in open(path):
        tag = line[:6].strip()
        if tag in ("ATOM", "HETATM"):
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            row = len(coords)
            coords.append([x, y, z])
            records.append((line, row))
            seg = line[72:76].strip()
            name = line[12:16].strip()
            resid = line[22:26].strip()
            idx.setdefault((seg, resid, name), row)
        else:
            records.append((line, None))
    return records, np.asarray(coords, float), idx


def write_pdb(path, records, coords):
    """Rewrite ``records`` to ``path`` with ATOM coordinates taken from ``coords``."""
    with open(path, "w") as fh:
        fh.write("REMARK  TOPO stage-3 oriented model "
                 "(PTC at origin; tunnel on +x; tRNA tails on +y).\n")
        for line, row in records:
            if row is None:
                fh.write(line)
            else:
                x, y, z = coords[row]
                fh.write("%s%8.3f%8.3f%8.3f%s" % (line[:30], x, y, z, line[54:]))


# ---------------------------------------------------------------------------
# Landmarks
# ---------------------------------------------------------------------------
def landmark_point(spec, idx, coords):
    """Resolve a ``segid:resid:atom`` (or ``:ribose``) landmark to an xyz point.

    ``resid`` may be the literal ``last``, meaning the highest-numbered nucleotide
    of ``segid`` that carries a full ribose ring -- i.e. the tRNA's 3'-terminal
    acceptor (CCA ``A``). This is numbering-robust: it is residue 76 for a
    standard 76-nt tRNA but 77 for the human P-site tRNA (a 77-nt tRNA), so the
    tail landmark always sits on the true acceptor rather than a hardcoded index.
    """
    seg, resid, atom = spec.split(":")
    if resid == "last":
        cands = [int(r) for (s, r, a) in idx
                 if s == seg and a == "C1'" and r.lstrip("-").isdigit()]
        if not cands:
            raise KeyError(f"landmark {spec}: no ribose nucleotides found for segid {seg}")
        resid = str(max(cands))
    if atom == "ribose":
        pts = [coords[idx[(seg, resid, a)]] for a in RIBOSE_RING
               if (seg, resid, a) in idx]
        if len(pts) < 3:
            raise KeyError(f"landmark {spec}: only {len(pts)}/5 ribose atoms found")
        return np.mean(pts, axis=0)
    key = (seg, resid, atom)
    if key not in idx:
        raise KeyError(f"landmark {spec}: atom not found in structure")
    return coords[idx[key]].copy()


# ---------------------------------------------------------------------------
# Rigid-body transform
# ---------------------------------------------------------------------------
def _rot_align(a, b):
    """Proper rotation mapping unit vector ``a`` onto unit vector ``b`` (Rodrigues)."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-12:                 # parallel or anti-parallel
        if c > 0:
            return np.eye(3)
        # 180 deg: rotate about any axis perpendicular to a
        perp = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, perp); axis /= np.linalg.norm(axis)
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[1], axis[0], 0]])
        return np.eye(3) + 2 * K @ K
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * (1.0 / (1.0 + c))


def build_transform(p_ptc, p_exit, p_a, p_p):
    """Return ``(R, t)`` implementing the ordered orient transform.

    Final position is ``x' = R @ (x - p_ptc)``; ``t = -R @ p_ptc``.
    """
    v_x = p_exit - p_ptc
    v_y = p_p - p_a

    # step 2: align tunnel vector to +x
    R1 = _rot_align(v_x, np.array([1.0, 0.0, 0.0]))

    # step 3: rotate about x so the yz-projection of v_y points to +y
    vy1 = R1 @ v_y
    phi = -np.arctan2(vy1[2], vy1[1])             # bring (y,z) onto (+, 0)
    c, s = np.cos(phi), np.sin(phi)
    R2 = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

    R = R2 @ R1
    t = -R @ p_ptc
    return R, t


def orient(config_path, verbose=True, tol=1e-4):
    """Run stage 3 from an INI config; write the oriented all-atom PDB."""
    cfgdir = os.path.dirname(os.path.abspath(config_path))
    cp = configparser.ConfigParser()
    cp.optionxform = str
    if not cp.read(config_path):
        raise FileNotFoundError(config_path)

    def rel(p):
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(cfgdir, p))

    asm = cp["assembly"]
    orient_cfg = cp["orient"]
    in_pdb = rel(asm["out"])                       # stage-2 output is stage-3 input
    out_pdb = rel(orient_cfg["oriented_out"])

    records, coords, idx = read_pdb_atoms(in_pdb)
    p_ptc = landmark_point(orient_cfg["ptc"], idx, coords)
    p_exit = landmark_point(orient_cfg["exit"], idx, coords)
    p_a = landmark_point(orient_cfg["a_site"], idx, coords)
    p_p = landmark_point(orient_cfg["p_site"], idx, coords)

    R, t = build_transform(p_ptc, p_exit, p_a, p_p)

    # --- verify the rotation is proper and orthonormal ---
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-8), "R is not orthonormal"
    assert abs(np.linalg.det(R) - 1.0) < 1e-8, \
        f"det(R) = {np.linalg.det(R):.6f} != +1 (reflection!)"

    new = (coords @ R.T) + t
    coords[:] = new

    # --- orientation asserts (§5) ---
    q_ptc = R @ p_ptc + t
    q_vx = R @ (p_exit - p_ptc)
    q_vy = R @ (p_p - p_a)
    len_vx = np.linalg.norm(p_exit - p_ptc)
    len_vy = np.linalg.norm(p_p - p_a)

    checks = {
        "|P_PTC| ~ 0": np.linalg.norm(q_ptc) < tol,
        "v_x on +x (y,z~0)": abs(q_vx[1]) < tol * len_vx + 1e-6
                             and abs(q_vx[2]) < tol * len_vx + 1e-6,
        "v_x points +x": q_vx[0] > 0,
        "v_y in xy-plane (z~0)": abs(q_vy[2]) < tol * len_vy + 1e-6,
        "v_y points +y": q_vy[1] > 0,
    }
    if verbose:
        print(f"Orient {os.path.basename(in_pdb)}:")
        print(f"  landmarks: PTC={orient_cfg['ptc']}  exit={orient_cfg['exit']}  "
              f"A={orient_cfg['a_site']}  P={orient_cfg['p_site']}")
        print(f"  det(R) = {np.linalg.det(R):+.6f}  (proper rotation)")
        print(f"  transformed P_PTC  = {q_ptc}")
        print(f"  transformed v_x    = {q_vx}   (|v_x|={len_vx:.2f})")
        print(f"  transformed v_y    = {q_vy}   (|v_y|={len_vy:.2f})")
        for name, ok in checks.items():
            print(f"    [{'OK' if ok else 'FAIL'}] {name}")
    if not all(checks.values()):
        failed = [n for n, ok in checks.items() if not ok]
        raise AssertionError(f"orientation asserts failed: {failed}")

    os.makedirs(os.path.dirname(out_pdb), exist_ok=True)
    write_pdb(out_pdb, records, coords)
    if verbose:
        print(f"  wrote {out_pdb}")
    return out_pdb


def main(argv=None):
    """CLI entry point: ``fix_orient.py -c configs/<organism>.ini``."""
    p = argparse.ArgumentParser(
        description="Stage 3: rigid-body orient (PTC->origin, tunnel->+x, "
                    "tRNA tails->+y).")
    p.add_argument("-c", "--config", required=True, help="organism INI config")
    args = p.parse_args(argv)
    orient(args.config)


if __name__ == "__main__":
    main()
