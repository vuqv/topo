#!/usr/bin/env python3
"""
Drive the full 5-stage ribosome-prep pipeline for one organism.

    python run.py -c configs/<organism>.ini [--keep-segids PtR,AtR] \
                  [--r-cyl 30 --x-lo -8 --x-exit 58]

Stages
------
1. (acquire)  -- the raw cif(s) are downloaded separately (see fetch.sh / PROVENANCE)
2. gen        -- gen_subunit.py   : assemble large subunit + tRNAs  -> *_model.pdb
3. orient     -- fix_orient.py    : PTC->origin, tunnel->+x, tails->+y -> *_oriented.pdb
4. coarse-grain -- cg_ribosome.py                                    -> *_cg.pdb
5. truncate   -- truncate_ribosome.py (cylinder + exit cap)          -> *_cg_trunc.pdb

All four stage scripts live in this bundle.
"""
from __future__ import annotations

import argparse
import configparser
import os

from gen_subunit import gen_subunit
from fix_orient import orient
from cg_ribosome import coarse_grain
from truncate_ribosome import truncate


def _cfg(config_path):
    cp = configparser.ConfigParser()
    cp.optionxform = str
    cp.read(config_path)
    return cp


def run(config_path, r_cyl=30.0, x_lo=-8.0, x_exit=58.0, keep_segids=None):
    """Run stages 2-5 for the organism described by ``config_path``."""
    cfgdir = os.path.dirname(os.path.abspath(config_path))
    cp = _cfg(config_path)
    tag = cp["assembly"]["tag"]
    oriented = os.path.normpath(os.path.join(cfgdir, cp["orient"]["oriented_out"]))
    outdir = os.path.dirname(oriented)
    cg = os.path.join(outdir, f"{tag}_model_cg.pdb")
    trunc = os.path.join(outdir, f"{tag}_model_cg_trunc.pdb")

    print("=" * 70, f"\nSTAGE 2  gen_subunit  ({tag})\n", "=" * 70, sep="")
    gen_subunit(config_path)
    print("=" * 70, "\nSTAGE 3  fix_orient\n", "=" * 70, sep="")
    orient(config_path)
    print("=" * 70, "\nSTAGE 4  cg_ribosome\n", "=" * 70, sep="")
    coarse_grain(oriented, cg)
    print("=" * 70, "\nSTAGE 5  truncate_ribosome\n", "=" * 70, sep="")
    truncate(cg, trunc, r_cyl=r_cyl, x_lo=x_lo, x_exit=x_exit,
             keep_segids=tuple(keep_segids or ()))
    print(f"\nDELIVERABLE: {trunc}")
    return trunc


def main(argv=None):
    p = argparse.ArgumentParser(description="Run the 5-stage ribosome-prep pipeline.")
    p.add_argument("-c", "--config", required=True)
    p.add_argument("--r-cyl", type=float, default=30.0)
    p.add_argument("--x-lo", type=float, default=-8.0)
    p.add_argument("--x-exit", type=float, default=58.0)
    p.add_argument("--keep-segids", default=None,
                   help="comma-separated segIDs to always keep whole (e.g. PtR,AtR)")
    args = p.parse_args(argv)
    keep = args.keep_segids.split(",") if args.keep_segids else None
    run(args.config, r_cyl=args.r_cyl, x_lo=args.x_lo, x_exit=args.x_exit,
        keep_segids=keep)


if __name__ == "__main__":
    main()
