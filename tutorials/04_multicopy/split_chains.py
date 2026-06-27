# Split a multi-copy trajectory into per-chain single-chain trajectories.
#
# run_simulation.py writes one DCD with `n_copies` non-interacting chains, where
# copy k occupies the contiguous atom block [k*n : (k+1)*n], plus the single-chain
# topology <outname>.psf. This thin wrapper just calls the package routine
# topo.split_chains, which streams the combined DCD in chunks (so it
# handles trajectories too large to fit in memory) and writes one single-chain DCD
# per copy: <outname>_<k>.dcd. Each chain is recentred per frame by default (pass
# center=False to preserve raw coordinates). Load each with the single-chain PSF
# (<outname>.psf) for ordinary analysis.
import argparse

import topo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-input', '-f', type=str, help='the same md.ini used for the run')
    args = parser.parse_args()

    cfg = topo.read_simulation_config(args.input, verbose=False)
    out_paths = [cfg.output_path(f'_{k}.dcd') for k in range(cfg.n_copies)]

    n_copies, atoms_per_copy, n_frames = topo.split_chains(
        cfg.output_path('.dcd'), out_paths)

    print(f"{n_frames} frames, {n_copies} copies -> {atoms_per_copy} atoms/chain")
    for path in out_paths:
        print(f"  wrote {path}")
    print(f"Load each with the single-chain topology, e.g.\n"
          f"    mda.Universe('{cfg.output_path('.psf')}', '{cfg.output_path('_0.dcd')}')")


if __name__ == '__main__':
    main()
