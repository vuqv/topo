#!/usr/bin/env python3
"""
remd_openmmtools.py -- true parallel-tempering REMD via openmmtools.multistate.

This is the "batteries-included" replica-exchange implementation: openmmtools
manages the N replicas, the Gibbs replica-mixing (all-pairs) exchange scheme, the
online free-energy estimation, and NetCDF storage/checkpointing/restart. You give
it the system, a temperature ladder, and an MCMC move; it runs the rest.

Requires openmmtools (NOT installed in this env by default)::

    conda install -c conda-forge openmmtools

Built on TOPO's public building blocks: it reuses the ``md.ini`` system via
``topo.engine.build_system``, then hands the ``openmm.System`` to openmmtools.

Outputs (under the md.ini ``output_dir``)
    <outname>_remd.nc          MultiStateReporter NetCDF (all replicas, all data)
    <outname>_remd_checkpoint.nc   periodic checkpoint for restart

Analyze afterwards with openmmtools::

    from openmmtools.multistate import MultiStateReporter, ParallelTemperingAnalyzer
    reporter = MultiStateReporter('<outname>_remd.nc', open_mode='r')
    analyzer = ParallelTemperingAnalyzer(reporter)
    dF, dF_err = analyzer.get_free_energy()             # MBAR free energies
    mixing, n_states, tau = analyzer.generate_mixing_statistics()   # exchange mixing

Run::

    python examples/parallel_tempering_openmmtools/remd_openmmtools.py -f md.ini \
        --tmin 300 --tmax 500 --nreplicas 8 --exchange-interval 1000 --cycles 100

``ref_t`` in the md.ini is ignored (the ladder sets the temperatures). Total steps
per replica = cycles * exchange-interval.
"""
import argparse
import logging
import sys

import openmm as mm
from openmm import unit

from topo import engine
from topo.utils.config import read_simulation_config

try:
    from openmmtools import states, mcmc, cache
    from openmmtools.multistate import ParallelTemperingSampler, MultiStateReporter
except ImportError:
    sys.exit("openmmtools is required: conda install -c conda-forge openmmtools")


def parse_args():
    parser = argparse.ArgumentParser(
        prog="remd_openmmtools.py",
        description="Parallel-tempering REMD on a TOPO CG model via openmmtools.")
    parser.add_argument('-input', '-f', required=True,
                        help='simulation config file (md.ini) -- supplies the system')
    parser.add_argument('--tmin', type=float, required=True,
                        help='lowest replica temperature (Kelvin)')
    parser.add_argument('--tmax', type=float, required=True,
                        help='highest replica temperature (Kelvin)')
    parser.add_argument('--nreplicas', type=int, required=True,
                        help='number of replicas / temperature rungs')
    parser.add_argument('--exchange-interval', type=int, default=1000,
                        help='MD steps per iteration between exchange attempts (default 1000)')
    parser.add_argument('--cycles', type=int, default=None,
                        help='number of exchange iterations (default: md_steps // '
                             'exchange-interval from the md.ini)')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"OpenMM version: {mm.__version__}")

    # --- Base system from the md.ini (ref_t overridden by the ladder) --------
    cfg = read_simulation_config(args.input)
    cfg.prepare_output_dir()

    # openmmtools reports iteration progress, per-replica energies, exchange
    # mixing and timing through the stdlib `logging` module (not stdout). Route
    # it to <output_dir>/<outname>.log so there is a readable simulation log
    # alongside the NetCDF.
    log_path = cfg.output_path('.log')
    handler = logging.FileHandler(log_path, mode='w')
    handler.setFormatter(logging.Formatter('%(asctime)s  %(levelname)-7s  %(message)s',
                                           datefmt='%H:%M:%S'))
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.addHandler(handler)
    print(f"[log] openmmtools iteration log -> {log_path}")

    tau_t = cfg.tau_t if cfg.tau_t is not None else 0.05 / unit.picoseconds
    cycles = args.cycles if args.cycles is not None else max(1, cfg.md_steps // args.exchange_interval)

    built = engine.build_system(cfg)

    # --- Honour the md.ini device: point openmmtools' context cache at it ----
    # openmmtools creates/destroys contexts through a global ContextCache; setting
    # its platform is how you pin CUDA vs CPU (there is no per-move platform arg).
    platform, _ = cfg.make_platform()
    cache.global_context_cache.platform = platform

    # --- Reference thermodynamic + sampler state -----------------------------
    # The reference state carries the System; ParallelTempering clones it across
    # the ladder, changing only the temperature. All replicas start from the same
    # (built) coordinates. box_vectors only matter under PBC.
    box_vectors = built.system.getDefaultPeriodicBoxVectors() if cfg.pbc else None
    ref_state = states.ThermodynamicState(system=built.system,
                                          temperature=args.tmin * unit.kelvin)
    sampler_state = states.SamplerState(positions=built.positions,
                                        box_vectors=box_vectors)

    # --- MCMC move: Langevin dynamics for `exchange-interval` steps -----------
    move = mcmc.LangevinDynamicsMove(timestep=cfg.dt,
                                     collision_rate=tau_t,
                                     n_steps=args.exchange_interval,
                                     reassign_velocities=False)

    # --- Build and run the parallel-tempering sampler ------------------------
    sampler = ParallelTemperingSampler(mcmc_moves=move, number_of_iterations=cycles)
    reporter = MultiStateReporter(cfg.output_path('_remd.nc'),
                                  checkpoint_interval=max(1, cycles // 10))
    sampler.create(thermodynamic_state=ref_state,
                   sampler_states=sampler_state,
                   storage=reporter,
                   min_temperature=args.tmin * unit.kelvin,
                   max_temperature=args.tmax * unit.kelvin,
                   n_temperatures=args.nreplicas)

    if cfg.minimize:
        print("Minimizing all replicas...")
        sampler.minimize()

    print(f"Parallel tempering: {args.nreplicas} replicas {args.tmin:.1f}-{args.tmax:.1f} K, "
          f"{cycles} iterations x {args.exchange_interval} steps")
    print("REMD started")
    sampler.run()

    # Exchange/mixing statistics and per-state free energies are written to the
    # NetCDF; read them back with MultiStateReporter + MultiStateSamplerAnalyzer
    # (see the module docstring).
    print(f"Done. NetCDF: {cfg.output_path('_remd.nc')}")


if __name__ == '__main__':
    main()
