"""Stitch nascent-chain synthesis trajectories into one VMD-playable movie.

A synthesis run writes a **separate** trajectory per growth step, and each step has
a **different number of beads** (length ``L`` has ``L`` CA atoms), so a single VMD
molecule (which needs a constant atom count) cannot just concatenate them. This
module pads every frame up to the final length -- parking the not-yet-synthesized
beads at a far sentinel coordinate -- and writes the segments back-to-back into one
fixed-width movie plus a ready-to-run ``movie.tcl`` that hides the parked beads, so
the chain appears to grow N->C.

It handles **two output layouts** with a shared stitching core
(:func:`stitch_segments`):

* **CSP (per-stage)** -- ``<out>/L_<L>/stage_<1,2,3>/`` from the Continuous Synthesis
  Protocol (:mod:`topo.csp.protocol`), plus ``ejection/`` / ``dissociation/``.
  Use :func:`stitch_movie` / :func:`find_stage_segments`. The movie plays the chain
  growing **stage by stage** (new residue at the A-site, settle, translocate to P).
* **Flat (per-length)** -- ``<out>/L_<L>/traj.dcd`` from a fixed-rate per-length loop
  (e.g. the Tutorial-9 cylinder runner), plus ``ejection/`` / ``dissociation/``.
  Use :func:`stitch_length_movie` / :func:`find_lengths`.

The :func:`main` CLI (``topo-csp-movie``) **auto-detects** the layout (per-stage if any
``L_<L>/stage_<s>/`` folders exist, else per-length).

Usage::

    topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
    python -m topo.csp.movie -o synth_out
    vmd -e synth_out/movie.tcl
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import sys
from typing import List, Optional, Tuple

import numpy as np

# Sentinel coordinate (angstrom) for parked (not-yet-synthesized) beads. Far
# enough that the VMD selection ``x > 9000`` cleanly isolates them; the view is
# fit to the full-length frame (where nothing is parked), so this does not zoom
# the camera out (see the generated movie.tcl).
SENTINEL_A = 99999.0

# Per-residue elongation sub-stages, in playback order (CSP layout).
STAGES = (1, 2, 3)
# Post-synthesis phase folders, in the order they should appear after growth.
# Both the CSP protocol and the flat per-length (cylinder) loop run the same two
# free-run phases: ejection then dissociation.
POST_PHASES_CSP = ("ejection", "dissociation")
POST_PHASES_FLAT = ("ejection", "dissociation")


# ==========================================================================
# Shared stitching core
# ==========================================================================
def stitch_segments(out_root: str, segments: List[Tuple[str, int, str, str]],
                    out_prefix: str = "movie", park: str = "sentinel",
                    ribosome_pdb: Optional[str] = None,
                    verbose: bool = True) -> Tuple[str, str, str]:
    """Pad and concatenate an ordered list of trajectory ``segments`` into a movie.

    The reusable core shared by the per-stage CSP movie (:func:`stitch_movie`) and
    the per-length flat movie (:func:`stitch_length_movie`). ``segments`` is an
    ordered list of ``(label, n_atoms, psf, dcd)`` -- each a standalone trajectory
    with ``n_atoms`` beads; every frame is padded up to the widest segment
    (``N`` = max ``n_atoms``) by parking the extra (not-yet-synthesized) beads, and
    the segments are written back-to-back in the given order. Also writes a
    ready-to-run ``<out_prefix>.tcl``.

    Parameters
    ----------
    out_root : str
        Directory the stitched movie files are written into.
    segments : list of tuple of (str, int, str, str)
        Ordered playback segments, each ``(label, n_atoms, psf, dcd)``: a human
        label, the segment's bead count, and its topology and trajectory paths.
    out_prefix : str
        Basename for the movie files (default ``movie``).
    park : {'sentinel', 'cterm'}
        Where to put not-yet-synthesized beads in each frame. ``sentinel`` (far
        away, hidden by the VMD script -- cleanest) or ``cterm`` (stacked on the
        C-terminus -- no VMD selection needed, but leaves a small bead cluster at
        the growing tip).
    ribosome_pdb : str, optional
        Path to a static CG ribosome PDB to copy next to the movie and load as
        fixed scenery in the generated ``.tcl``; ``None`` to skip.
    verbose : bool
        When true, print progress messages while stitching.

    Returns
    -------
    tuple of (str, str, str)
        The ``(psf_path, dcd_path, tcl_path)`` of the written movie files.

    Raises
    ------
    ValueError
        If ``park`` is not ``'sentinel'`` or ``'cterm'``.
    SystemExit
        If ``segments`` is empty.
    """
    import MDAnalysis as mda  # heavy; import only when actually stitching

    def log(msg: str) -> None:
        """Print ``msg`` to stdout when this stitch is running verbosely.

        Parameters
        ----------
        msg : str
            The progress message to emit (suppressed when ``verbose`` is false).
        """
        if verbose:
            print(msg)

    if park not in ("sentinel", "cterm"):
        raise ValueError(f"park must be 'sentinel' or 'cterm', got {park!r}.")
    if not segments:
        raise SystemExit("no trajectory segments to stitch.")

    out_psf = os.path.join(out_root, f"{out_prefix}.psf")
    out_dcd = os.path.join(out_root, f"{out_prefix}.dcd")
    out_tcl = os.path.join(out_root, f"{out_prefix}.tcl")

    # Movie topology = a segment with the most beads (the final length); every
    # frame is padded up to N atoms.
    N = max(n for _, n, _, _ in segments)
    psf_full = next(psf for _, n, psf, _ in segments if n == N)
    shutil.copyfile(psf_full, out_psf)
    full = mda.Universe(psf_full)
    # Give the topology-only universe a single in-memory coordinate frame so we
    # can assign per-frame positions into it before writing.
    full.load_new(np.zeros((1, N, 3), dtype=np.float32))
    log(f"Movie topology: {out_psf}  ({N} beads = final length)")
    log(f"Parking not-yet-synthesized beads: {park}")

    total_frames = 0
    skipped = 0
    with mda.Writer(out_dcd, n_atoms=N) as writer:
        for label, n, psf, dcd in segments:
            # A segment can be unreadable if its DCD is still being written (an
            # in-progress run) or was truncated by a crash -- skip it (with a warning)
            # rather than aborting the whole movie. Open + read inside the guard so a
            # premature-EOF header error or a zero-frame DCD just drops that segment.
            try:
                u = mda.Universe(psf, dcd)
                nfr = 0
                for _ in u.trajectory:
                    coords = np.empty((N, 3), dtype=np.float32)
                    coords[:n] = u.atoms.positions
                    if n < N:
                        if park == "sentinel":
                            coords[n:] = SENTINEL_A
                        else:  # stack on the current C-terminus (last real bead)
                            coords[n:] = u.atoms.positions[n - 1]
                    full.atoms.positions = coords
                    writer.write(full.atoms)
                    nfr += 1
            except Exception as exc:
                skipped += 1
                log(f"  {label:>12}: SKIPPED (unreadable/empty DCD: "
                    f"{type(exc).__name__}: {exc})")
                continue
            total_frames += nfr
            log(f"  {label:>12}: {nfr} frames")
    if skipped:
        log(f"  ({skipped} segment(s) skipped -- truncated or still being written)")

    log(f"Movie trajectory: {out_dcd}  ({total_frames} frames total)")

    # Optional static ribosome reference (v2): copy it next to the movie so the
    # generated tcl can load it as fixed scenery the chain grows inside.
    ribo_name = None
    if ribosome_pdb is not None:
        ribo_name = f"{out_prefix}_ribosome.pdb"
        shutil.copyfile(ribosome_pdb, os.path.join(out_root, ribo_name))
        log(f"Static ribosome reference: {os.path.join(out_root, ribo_name)}")

    _write_tcl(out_tcl, os.path.basename(out_psf), os.path.basename(out_dcd),
               park=park, ribosome_name=ribo_name)
    log(f"VMD script: {out_tcl}")
    log("")
    log(f"View it with:  vmd -e {out_tcl}")
    return out_psf, out_dcd, out_tcl


def _write_tcl(path: str, psf_name: str, dcd_name: str, park: str,
               ribosome_name: Optional[str] = None) -> None:
    """Write a VMD script that loads the movie and grows the chain N->C.

    Parameters
    ----------
    path : str
        Destination path for the generated ``.tcl`` script.
    psf_name : str
        Basename of the movie topology file the script should load.
    dcd_name : str
        Basename of the movie trajectory file the script should load.
    park : {'sentinel', 'cterm'}
        Parking scheme used when stitching; ``sentinel`` adds a per-frame
        selection that hides the far-parked future beads, ``cterm`` shows all
        beads (no hiding needed).
    ribosome_name : str, optional
        Basename of a static ribosome PDB to load as a separate scenery
        molecule; ``None`` to omit the ribosome block.

    Returns
    -------
    None
    """
    # The hiding selection is only needed for the 'sentinel' parking scheme.
    if park == "sentinel":
        sel = "not (x > 9000)"
        hide_note = ("# Not-yet-synthesized beads are parked far away and hidden "
                     "each frame\n# (selection re-evaluated per frame via selupdate).")
    else:
        sel = "all"
        hide_note = ("# 'cterm' parking: future beads are stacked on the "
                     "C-terminus (no hiding).")

    # Optional: load the static ribosome (v2) as a separate molecule for context.
    ribo_block = ""
    if ribosome_name is not None:
        ribo_block = f"""
# Static ribosome scenery (v2): a separate molecule the chain grows inside.
mol new {ribosome_name} type pdb waitfor all
mol delrep 0 top
mol representation Points 1.0
mol color ColorID 6
mol selection {{all}}
mol material Transparent
mol addrep top
mol top 0
"""

    tcl = f"""# VMD visualization of the co-translational synthesis movie.
# Generated by topo.csp.movie.
#   vmd -e {os.path.basename(path)}      (run from this folder)

mol new {psf_name} type psf waitfor all
mol addfile {dcd_name} type dcd waitfor all

mol delrep 0 top

{hide_note}
# Beads (van der Waals); colored by residue id so the growing chain is a gradient.
mol representation VDW 1.5 16.0
mol color ResID
mol selection {{{sel}}}
mol material Opaque
mol addrep top
mol selupdate 0 top on

# Backbone trace (bonds between consecutive synthesized beads).
mol representation Licorice 0.4 16.0
mol color ResID
mol selection {{{sel}}}
mol addrep top
mol selupdate 1 top on
{ribo_block}
# Fit the camera to the LAST frame (full length -> nothing parked), then rewind.
set nf [molinfo top get numframes]
if {{$nf > 0}} {{
    animate goto [expr {{$nf - 1}}]
    display resetview
    animate goto 0
}}

axes location off
display projection Orthographic
color Display Background white
animate speed 0.85

puts "Loaded $nf frames. Press Play (or run: animate forward) to watch the chain grow."
"""
    with open(path, "w") as fh:
        fh.write(tcl)


# ==========================================================================
# CSP layout: <out>/L_<L>/stage_<1,2,3>/  (+ ejection/ dissociation/)
# ==========================================================================
def _pick_traj(phase_dir: str, outname: str = "traj") -> Optional[str]:
    """Pick the trajectory file to read for one stage/phase, or ``None`` if absent.

    Prefer ``<outname>.dcd`` when it actually holds frames, but a coarse ``nstout``
    (relative to the steps run per stage) can leave a stage's DCD **empty** (0 bytes)
    -- every frame interval landed past the end of the short run. In that case fall
    back to ``<outname>_final.pdb``, the last-frame snapshot the runner always writes,
    so the stage still contributes its (single) final conformation instead of being
    silently dropped. Without this, a coarse-output CSP run yields a movie that skips
    most lengths.

    Parameters
    ----------
    phase_dir : str
        Directory of a single stage or post-synthesis phase (e.g.
        ``<out_root>/L_<L>/stage_<s>`` or ``<out_root>/ejection``).
    outname : str, optional
        Per-stage output basename used by the runner (default ``"traj"``); the
        candidate files are ``<outname>.dcd`` and ``<outname>_final.pdb``.

    Returns
    -------
    str or None
        Path to ``<outname>.dcd`` if it exists and is non-empty, else
        ``<outname>_final.pdb`` if present, else the empty ``<outname>.dcd`` if it
        exists, else ``None`` when no trajectory file is found.
    """
    dcd = os.path.join(phase_dir, f"{outname}.dcd")
    if os.path.isfile(dcd) and os.path.getsize(dcd) > 0:
        return dcd
    final_pdb = os.path.join(phase_dir, f"{outname}_final.pdb")
    if os.path.isfile(final_pdb):
        return final_pdb
    if os.path.isfile(dcd):  # empty DCD and no final-frame fallback -- read it anyway
        return dcd
    return None


def find_stage_segments(out_root: str, outname: str = "traj"
                        ) -> List[Tuple[str, int, str, str]]:
    """Return ordered ``[(label, n_atoms, psf, traj), ...]`` for a CSP run.

    Walks ``<out_root>/L_<L>/stage_<1,2,3>/`` in increasing ``L`` then stage order
    (only stages with a ``.psf`` and a readable trajectory are kept), then appends any
    ``ejection/`` / ``dissociation/`` phase. ``n_atoms`` is ``L`` for a stage of
    length ``L`` (the nascent-only output) and the post-phase psf's atom count for
    the post phases. The trajectory is the stage's ``.dcd`` when it has frames, else
    its ``_final.pdb`` snapshot (see :func:`_pick_traj`).

    Parameters
    ----------
    out_root : str
        CSP run output root, containing the ``L_<L>/`` length folders and any
        ``ejection/`` / ``dissociation/`` post-synthesis phases.
    outname : str, optional
        Per-stage output basename used by the runner (default ``"traj"``); the
        per-stage psf is ``<outname>.psf`` and the trajectory is chosen by
        :func:`_pick_traj`.

    Returns
    -------
    list of tuple of (str, int, str, str)
        Segments in synthesis order, each ``(label, n_atoms, psf, traj)``:
        ``label`` is ``"L=<L> s<s>"`` for a stage or the phase name for a post
        phase; ``n_atoms`` is the nascent bead count (``L`` for a stage, the psf's
        atom count for a post phase); ``psf`` and ``traj`` are the topology and
        trajectory paths. Stages or phases lacking a ``.psf`` or a readable
        trajectory are skipped.
    """
    lengths = []
    for d in glob.glob(os.path.join(out_root, "L_*")):
        m = re.search(r"L_(\d+)$", os.path.basename(d))
        if m:
            lengths.append((int(m.group(1)), d))
    lengths.sort(key=lambda t: t[0])

    segments: List[Tuple[str, int, str, str]] = []
    for L, d in lengths:
        for s in STAGES:
            sd = os.path.join(d, f"stage_{s}")
            psf = os.path.join(sd, f"{outname}.psf")
            traj = _pick_traj(sd, outname)
            if os.path.isfile(psf) and traj is not None:
                segments.append((f"L={L} s{s}", L, psf, traj))

    # Post-synthesis phases (at full length). Read the bead count from the psf.
    import MDAnalysis as mda  # heavy; only when post phases are present
    for name in POST_PHASES_CSP:
        pd = os.path.join(out_root, name)
        psf = os.path.join(pd, f"{outname}.psf")
        traj = _pick_traj(pd, outname)
        if os.path.isfile(psf) and traj is not None:
            n = len(mda.Universe(psf).atoms)
            segments.append((name, n, psf, traj))

    return segments


def stitch_movie(out_root: str, out_prefix: str = "movie",
                 park: str = "sentinel", outname: str = "traj",
                 ribosome_pdb: Optional[str] = None,
                 verbose: bool = True) -> Tuple[str, str, str]:
    """Discover the CSP per-stage segments and stitch them into a movie.

    Thin CSP-aware wrapper around :func:`stitch_segments` (the shared core). Calls
    :func:`find_stage_segments` to enumerate the per-stage and post-synthesis
    segments, then hands them to the shared padding/stitching core, which writes
    ``<out_prefix>.psf`` / ``<out_prefix>.dcd`` / ``<out_prefix>.tcl``.

    Parameters
    ----------
    out_root : str
        CSP run output root (contains the ``L_<L>/stage_<s>/`` folders).
    out_prefix : str, optional
        Basename for the stitched movie files (default ``"movie"``).
    park : str, optional
        Where to put not-yet-synthesized beads in each frame (default
        ``"sentinel"``); passed through to the stitching core.
    outname : str, optional
        Per-stage output basename used by the runner (default ``"traj"``).
    ribosome_pdb : str or None, optional
        Optional coarse-grained ribosome PDB to load as static scenery in the
        generated ``<out_prefix>.tcl`` (default ``None``).
    verbose : bool, optional
        If ``True`` (default), print the number of discovered segments and let the
        stitching core report progress.

    Returns
    -------
    tuple of (str, str, str)
        Paths to the written ``psf``, ``dcd`` and ``tcl`` movie files, as returned
        by :func:`stitch_segments`.

    Raises
    ------
    SystemExit
        If no per-stage trajectories are found under ``out_root``.
    """
    segments = find_stage_segments(out_root, outname=outname)
    if not segments:
        raise SystemExit(
            f"no per-stage trajectories found under {out_root!r} "
            f"(expected {out_root}/L_<L>/stage_<1,2,3>/{outname}.dcd + .psf). "
            f"Did you run `topo-csp -f csp.ini` first?")
    if verbose:
        print(f"Found {len(segments)} CSP segments under {out_root}/ "
              f"(per residue x 3 stages + post-synthesis).")
    return stitch_segments(out_root, segments, out_prefix=out_prefix, park=park,
                           ribosome_pdb=ribosome_pdb, verbose=verbose)


# ==========================================================================
# Flat layout: <out>/L_<L>/traj.dcd  (+ ejection/ dissociation/)
# ==========================================================================
def find_lengths(out_root: str,
                 outname: str = "traj") -> List[Tuple[int, str, str]]:
    """Return ``[(L, psf, dcd), ...]`` for each ``<out_root>/L_<L>/`` length.

    Scans ``<out_root>`` for ``L_<L>/`` length folders, sorted by ``L``; only
    lengths with both a ``.psf`` and a ``.dcd`` are kept. This is the **flat**
    per-length layout (one trajectory per length, e.g. the Tutorial-9 cylinder
    runner), in contrast to the CSP per-stage layout (:func:`find_stage_segments`).

    Parameters
    ----------
    out_root : str
        The run's output root (contains the ``L_<L>/`` folders).
    outname : str
        Per-length output basename used by the runner (default ``traj``); the
        files looked for are ``<outname>.psf`` and ``<outname>.dcd``.

    Returns
    -------
    list of tuple of (int, str, str)
        One ``(L, psf, dcd)`` per length found, sorted by ascending ``L``;
        ``psf`` and ``dcd`` are the paths to that length's topology and trajectory.
    """
    items = []
    for d in glob.glob(os.path.join(out_root, "L_*")):
        m = re.search(r"L_(\d+)$", os.path.basename(d))
        if not m:
            continue
        L = int(m.group(1))
        psf = os.path.join(d, f"{outname}.psf")
        dcd = os.path.join(d, f"{outname}.dcd")
        if os.path.isfile(psf) and os.path.isfile(dcd):
            items.append((L, psf, dcd))
    items.sort(key=lambda t: t[0])
    return items


def find_post(out_root: str, outname: str = "traj") -> List[Tuple[str, str, str]]:
    """Return ``[(name, psf, dcd), ...]`` for present post-synthesis phases (flat).

    Looks for ``<out_root>/ejection/`` and ``<out_root>/dissociation/`` (in that
    order) -- the optional post-synthesis runs written after the chain reaches its
    final length in the flat per-length layout. Only phases with both a ``.psf`` and
    a ``.dcd`` are returned.

    Parameters
    ----------
    out_root : str
        The run's output root (may contain the ``ejection/`` and/or ``dissociation/``
        folders).
    outname : str
        Per-length output basename used by the runner (default ``traj``); the
        files looked for are ``<outname>.psf`` and ``<outname>.dcd``.

    Returns
    -------
    list of tuple of (str, str, str)
        One ``(name, psf, dcd)`` per present post-synthesis phase, in
        ``POST_PHASES_FLAT`` order.
    """
    found = []
    for name in POST_PHASES_FLAT:
        d = os.path.join(out_root, name)
        psf = os.path.join(d, f"{outname}.psf")
        dcd = os.path.join(d, f"{outname}.dcd")
        if os.path.isfile(psf) and os.path.isfile(dcd):
            found.append((name, psf, dcd))
    return found


def stitch_length_movie(out_root: str, out_prefix: str = "movie",
                        park: str = "sentinel", outname: str = "traj",
                        ribosome_pdb: Optional[str] = None,
                        verbose: bool = True) -> Tuple[str, str, str]:
    """Discover the flat per-length segments and stitch them into a movie.

    The flat-layout counterpart of :func:`stitch_movie`: enumerates
    ``<out_root>/L_<L>/`` per-length trajectories (:func:`find_lengths`), appends any
    ``ejection/`` / ``dissociation/`` phase (:func:`find_post`), and hands them to the
    shared :func:`stitch_segments` core. Used by the Tutorial-9 cylinder runner.

    Parameters
    ----------
    out_root : str
        The run's output root (contains the ``L_<L>/`` folders).
    out_prefix : str
        Basename for the movie files (default ``movie``).
    park : {'sentinel', 'cterm'}
        Where to put not-yet-synthesized beads in each frame (passed through to the
        stitching core).
    outname : str
        Per-length output basename used by the runner (default ``traj``).
    ribosome_pdb : str, optional
        Optional CG ribosome PDB to load as static scenery in the generated
        ``.tcl`` (default ``None``).
    verbose : bool
        When true, let the stitching core report progress.

    Returns
    -------
    tuple of (str, str, str)
        The ``(psf_path, dcd_path, tcl_path)`` of the written movie files.

    Raises
    ------
    SystemExit
        If no per-length trajectories are found under ``out_root``.
    """
    import MDAnalysis as mda  # heavy; import only when actually stitching

    items = find_lengths(out_root, outname=outname)
    if not items:
        raise SystemExit(
            f"no per-length trajectories found under {out_root!r} "
            f"(expected {out_root}/L_<L>/{outname}.dcd + .psf).")

    # Ordered playback segments: each growth length, then the post-synthesis phases
    # (ejection / dissociation) at full length. Each segment is (label, n_atoms,
    # psf, dcd); the post phases run at the final length, so n_atoms = max L.
    segments = [(f"L={L}", L, psf, dcd) for L, psf, dcd in items]
    for name, psf, dcd in find_post(out_root, outname=outname):
        n = len(mda.Universe(psf).atoms)
        segments.append((name, n, psf, dcd))

    return stitch_segments(out_root, segments, out_prefix=out_prefix, park=park,
                           ribosome_pdb=ribosome_pdb, verbose=verbose)


# ==========================================================================
# CLI
# ==========================================================================
def _is_csp_layout(out_root: str) -> bool:
    """Return True if ``out_root`` has the CSP per-stage layout.

    Parameters
    ----------
    out_root : str
        Output root to inspect.

    Returns
    -------
    bool
        True if any ``<out_root>/L_<L>/stage_<s>/`` folder exists (CSP), else False
        (treat as the flat per-length layout).
    """
    return bool(glob.glob(os.path.join(out_root, "L_*", "stage_*")))


def main(argv: Optional[List[str]] = None) -> None:
    """Command-line entry point for ``topo-csp-movie`` / ``python -m topo.csp.movie``.

    Parses arguments and stitches the synthesis trajectories under ``--out-root``
    into one VMD-playable movie, **auto-detecting** the layout: the per-stage CSP
    layout (:func:`stitch_movie`) if any ``L_<L>/stage_<s>/`` folders are present,
    otherwise the flat per-length layout (:func:`stitch_length_movie`). The parser
    exposes ``-o/--out-root`` (required), ``--prefix``, ``--park``
    (``sentinel``/``cterm``), ``--outname`` and ``--ribosome``. With no arguments it
    prints help and exits ``0``.

    Parameters
    ----------
    argv : list of str or None, optional
        Argument vector to parse. If ``None`` (default), arguments are taken from
        ``sys.argv``.

    Returns
    -------
    None
    """
    p = argparse.ArgumentParser(
        prog="topo-csp-movie",
        description="Stitch nascent-chain synthesis trajectories into one "
                    "VMD-playable movie that grows the chain, and write a movie.tcl "
                    "to view it. Auto-detects the CSP per-stage layout "
                    "(<out>/L_<L>/stage_<1,2,3>/) or the flat per-length layout "
                    "(<out>/L_<L>/), plus any ejection/ dissociation/ phase.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--out-root", required=True,
                   help="synthesis run output root (contains the L_<L>/ folders).")
    p.add_argument("--prefix", default="movie",
                   help="basename for the stitched movie files.")
    p.add_argument("--park", default="sentinel", choices=["sentinel", "cterm"],
                   help="where to put not-yet-synthesized beads each frame.")
    p.add_argument("--outname", default="traj",
                   help="per-segment output basename used by the runner.")
    p.add_argument("--ribosome", default=None,
                   help="optional CG ribosome PDB to load as static scenery in the "
                        "generated movie.tcl.")
    if argv is None and len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)
    args = p.parse_args(argv)
    stitch = stitch_movie if _is_csp_layout(args.out_root) else stitch_length_movie
    stitch(args.out_root, out_prefix=args.prefix, park=args.park,
           outname=args.outname, ribosome_pdb=args.ribosome)


if __name__ == "__main__":
    main()
