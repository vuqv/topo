#!/usr/bin/env python3
"""Emit a VMD scenery script that draws Tutorial 7's analytic exit tunnel.

The cylinder synthesis (`topo-cylinder`) has **no ribosome beads** -- the tunnel
is a pure force. To *show* it in the process GIF we draw it as VMD ``graphics``
on their own molecule: the cylindrical **bore**, the closed **PTC** end cap, and
the semi-transparent **exit-face wall** (an annulus whose hole is the bore). The
geometry is read from the same ``cylinder.ini`` that drove the run, so the drawn
tunnel always matches the forces the chain actually felt.

The output is meant to be passed to ``render_cg.py --extra-tcl <file>``.

    python cylinder_tunnel_tcl.py -f cylinder.ini -o tunnel.tcl

This is a headless-render companion to ``make_movie_cylinder.py`` (which writes
an *interactive* tunnel-aware ``movie.tcl``); both read the same geometry.
"""
from __future__ import annotations

import argparse

# Single source of truth for the tunnel geometry -> VMD graphics (also used by
# `topo-csp-movie --tunnel`); this headless companion just writes it to a file.
from topo.csp.cylinder import tunnel_tcl


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--ini", default="cylinder.ini", help="the run's cylinder.ini.")
    ap.add_argument("-o", "--out", default="tunnel.tcl", help="output TCL path.")
    ap.add_argument("--wall-outer", type=float, default=None,
                    help="outer radius of the drawn exit wall (nm; default r+3).")
    args = ap.parse_args(argv)
    with open(args.out, "w") as fh:
        fh.write(tunnel_tcl(args.ini, wall_outer_nm=args.wall_outer))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
