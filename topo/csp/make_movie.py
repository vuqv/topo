"""Stitch the CSP per-residue/-stage trajectories into one VMD-playable movie.

The continuous-synthesis runner (:mod:`topo.csp.csp`) writes a standalone trajectory
for **every sub-stage of every residue**: ``<out>/L_<L>/stage_<1,2,3>/traj.dcd`` (and,
after the chain reaches full length, ``ejection/`` and optional ``dissociation/``).
Length ``L`` has ``L`` CA beads, so the per-stage DCDs cannot just be concatenated
(VMD needs a constant atom count).

This tool discovers those segments **in synthesis order** -- ``L=5`` stage 1, 2, 3,
then ``L=6`` stage 1, 2, 3, ... then ``ejection`` / ``dissociation`` -- and hands them
to the shared padding/stitching core (:func:`topo.translation.make_movie.stitch_segments`),
which pads every frame up to the final length (parking the not-yet-made beads) and
writes ``movie.psf`` / ``movie.dcd`` / ``movie.tcl``. The movie therefore plays the
chain growing **stage by stage**, so you can watch the new residue appear at the A-site
(stage 1), settle (stage 2) and translocate to the P-site (stage 3) for each codon.

Usage::

    topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
    vmd -e synth_out/movie.tcl
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import List, Optional, Tuple

from topo.translation.make_movie import stitch_segments

# Per-residue elongation sub-stages, in playback order.
STAGES = (1, 2, 3)
# Post-synthesis phase folders, in the order they should appear after growth.
POST_PHASES = ("ejection", "dissociation")


def _pick_traj(phase_dir: str, outname: str = "traj") -> Optional[str]:
    """Pick the trajectory file to read for one stage/phase, or ``None`` if absent.

    Prefer ``<outname>.dcd`` when it actually holds frames, but a coarse ``nstout``
    (relative to the steps run per stage) can leave a stage's DCD **empty** (0 bytes)
    -- every frame interval landed past the end of the short run. In that case fall
    back to ``<outname>_final.pdb``, the last-frame snapshot the runner always writes,
    so the stage still contributes its (single) final conformation instead of being
    silently dropped. Without this, a coarse-output CSP run yields a movie that skips
    most lengths.
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
    for name in POST_PHASES:
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
    """Discover the CSP stage segments and stitch them into a movie.

    Thin CSP-aware wrapper around
    :func:`topo.translation.make_movie.stitch_segments` (the shared core). See the
    module docstring for the output layout it expects.
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


def main(argv: Optional[List[str]] = None) -> None:
    """CLI: ``topo-csp-movie -o <out_root>``."""
    p = argparse.ArgumentParser(
        prog="topo-csp-movie",
        description="Stitch the CSP per-residue/-stage trajectories "
                    "(<out_root>/L_<L>/stage_<1,2,3>/traj.dcd) -- plus any "
                    "ejection/ or dissociation/ phase -- into one VMD-playable movie "
                    "that grows the nascent chain stage by stage, and write a "
                    "movie.tcl to view it.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--out-root", required=True,
                   help="CSP run output root (contains the L_<L>/ folders).")
    p.add_argument("--prefix", default="movie",
                   help="basename for the stitched movie files.")
    p.add_argument("--park", default="sentinel", choices=["sentinel", "cterm"],
                   help="where to put not-yet-synthesized beads each frame.")
    p.add_argument("--outname", default="traj",
                   help="per-stage output basename used by the runner.")
    p.add_argument("--ribosome", default=None,
                   help="optional CG ribosome PDB to load as static scenery in the "
                        "generated movie.tcl.")
    if argv is None and len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)
    args = p.parse_args(argv)
    stitch_movie(args.out_root, out_prefix=args.prefix, park=args.park,
                 outname=args.outname, ribosome_pdb=args.ribosome)


if __name__ == "__main__":
    main()
