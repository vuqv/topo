#!/usr/bin/env python3
"""Collapse sentinel-parked beads onto the growing tip so a Tube renders cleanly.

The synthesis-movie stitcher (``topo.csp.movie`` / ``make_movie_cylinder.py``)
pads every frame to the full residue count and parks the **not-yet-synthesized**
beads at a far sentinel (x ≈ 99999 Å). VMD's *VDW/Licorice* reps hide those with a
per-frame ``x < 9000`` selection, but a *Tube* rep still splines a long streak out
to the parked bead (the selection hides the sphere, not the spline vertex).

This rewrites the movie DCD so that, in each frame, every parked bead is moved
onto the **highest-index synthesized (visible) bead** — the growing tip. The Tube
then stays on the real chain (the collapsed future beads pile invisibly at the
tip), so no ``--select``/``--selupdate`` is needed and there are no streaks.

    python collapse_sentinels.py -p movie.psf -d movie.dcd -o movie_clean.dcd
"""
from __future__ import annotations

import argparse

import mdtraj as md
import numpy as np

# Parked-bead threshold in nm (movie sentinel is ~99999 Å = 9999.9 nm).
SENTINEL_NM = 900.0


def collapse(psf: str, dcd: str, out: str, verbose: bool = True) -> None:
    t = md.load(dcd, top=psf)
    xyz = t.xyz  # (n_frames, n_atoms, 3), nm
    for i in range(t.n_frames):
        vis = xyz[i, :, 0] < SENTINEL_NM
        if vis.all() or not vis.any():
            continue
        tip = np.where(vis)[0].max()          # growing tip = highest visible index
        xyz[i, ~vis] = xyz[i, tip]            # pile parked beads on the tip
    t.xyz = xyz
    t.save_dcd(out)
    if verbose:
        print(f"collapsed sentinels: {dcd} -> {out} ({t.n_frames} frames)")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-p", "--psf", required=True)
    ap.add_argument("-d", "--dcd", required=True)
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args(argv)
    collapse(args.psf, args.dcd, args.out)


if __name__ == "__main__":
    main()
