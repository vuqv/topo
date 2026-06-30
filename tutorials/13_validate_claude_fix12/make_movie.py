#!/usr/bin/env python
"""Stitch this validation run's CSP trajectories into one VMD movie.

The continuous-synthesis runner (``topo-csp -f csp.ini``) writes a standalone
trajectory for every sub-stage of every residue under ``synth_out_debug/``::

    synth_out_debug/L_<L>/stage_<1,2,3>/traj.{psf,dcd}   # nascent chain at length L
    synth_out_debug/ejection/traj.{psf,dcd}              # full-length chain leaving

Each length ``L`` has a different bead count, so the per-stage DCDs cannot be
concatenated directly (VMD needs a constant atom count). This script discovers
the segments in synthesis order and hands them to the shared stitching core
(:func:`topo.csp.movie.stitch_movie`), which pads every frame up to the
final length -- parking not-yet-made beads -- and writes ``movie.psf`` /
``movie.dcd`` / ``movie.tcl``. The truncated ribosome is loaded as static
scenery so you can watch the nascent chain grow stage by stage out of the exit
tunnel.

This is just a thin, path-pinned wrapper around the ``topo-csp-movie`` console
tool; equivalently you could run::

    topo-csp-movie -o synth_out_debug --ribosome ribosome_trunc.pdb

Usage::

    python make_movie.py
    vmd -e synth_out_debug/movie.tcl
"""
from __future__ import annotations

import os

from topo.csp.movie import stitch_movie

# Run next to this script regardless of the caller's cwd.
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.path.join(HERE, "synth_out")   # matches `outdir` in csp_val.ini
RIBOSOME = os.path.join(HERE, "ribosome_trunc.pdb")  # static scenery


def main() -> None:
    ribosome = RIBOSOME if os.path.isfile(RIBOSOME) else None
    psf, dcd, tcl = stitch_movie(
        OUT_ROOT,
        out_prefix="movie",
        park="sentinel",         # park not-yet-synthesized beads off to the side
        ribosome_pdb=ribosome,
        verbose=True,
    )
    print("\nWrote:")
    print(f"  {psf}")
    print(f"  {dcd}")
    print(f"  {tcl}")
    print(f"\nView with:\n  vmd -e {tcl}")


if __name__ == "__main__":
    main()
