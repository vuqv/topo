# Split a multi-copy trajectory into per-chain single-chain trajectories.
#
# run_simulation.py writes one DCD with `n_copies` non-interacting chains, where
# copy k occupies the contiguous atom block [k*n : (k+1)*n], plus two topologies:
#   {protein_code}.psf        -- single-chain (the topology of each copy)
#   {protein_code}_multi.psf  -- multi-chain (matches the combined .dcd)
# This script writes one single-chain DCD per copy; load each with the
# single-chain PSF ({protein_code}.psf) for ordinary analysis.
import argparse

import MDAnalysis as mda

import topo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-input', '-f', type=str, help='the same md.ini used for the run')
    args = parser.parse_args()

    cfg = topo.read_simulation_config(args.input, verbose=False)
    prefix = cfg.protein_code
    n_copies = cfg.n_copies

    u = mda.Universe(f'{prefix}_multi.psf', f'{prefix}.dcd')   # multi-chain topology + trajectory
    n_total = u.atoms.n_atoms
    if n_total % n_copies != 0:
        raise SystemExit(f"{n_total} atoms is not divisible by n_copies={n_copies}; "
                         f"check that md.ini matches the run.")
    n_single = n_total // n_copies
    print(f"{n_total} atoms, {n_copies} copies -> {n_single} atoms/chain, "
          f"{len(u.trajectory)} frames")

    for k in range(n_copies):
        ag = u.atoms[k * n_single:(k + 1) * n_single]
        with mda.Writer(f'{prefix}_chain{k}.dcd', ag.n_atoms) as w:
            for _ in u.trajectory:
                w.write(ag)
        print(f"  wrote {prefix}_chain{k}.dcd")

    print(f"Done: {n_copies} single-chain trajectories. Load each with the "
          f"single-chain topology, e.g.\n"
          f"    mda.Universe('{prefix}.psf', '{prefix}_chain0.dcd')")


if __name__ == '__main__':
    main()
