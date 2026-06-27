"""
Run a TOPO coarse-grained simulation from a control file (md.ini).

This is the canonical runner for the package. Use it as a CLI::

    python -m topo.mdrun -f md.ini

or call :func:`mdrun` from your own script. Control-file parsing lives in
``topo.read_simulation_config``. After building the single-chain model, if
``md.ini`` sets ``n_copies > 1`` the model is replicated into that many
non-interacting copies with ``topo.make_noninteracting_copies`` (default
``n_copies = 1`` = single chain).

The temperature *protocol* is selected by the control file:

    * ``anneal = no``  (default) -- constant-temperature equilibrium at ``ref_t``.
    * ``anneal = yes``           -- hold at ``t_high`` then quench/cool to ``ref_t``
                                    (``ref_t`` is the low / refold temperature).

Both share the same build / setup / finalize machinery (:mod:`topo.engine`);
only the temperature schedule differs (:mod:`topo.mdrun.protocol`).
"""
import argparse
import sys
import time

import openmm as mm

from topo import engine
from topo.utils.config import read_simulation_config
from .protocol import temperature_schedule, run_protocol, describe_schedule


def mdrun():
    """
    Run a simulation using the TOPO library and parameters specified in a config file.

    Usage: python -m topo.mdrun -f md.ini
    """
    parser = argparse.ArgumentParser(
        prog="topo-mdrun",
        description="Run a TOPO coarse-grained simulation from a control file "
                    "(md.ini). Replicates into non-interacting copies when "
                    "n_copies > 1. Supports constant-temperature equilibrium "
                    "(default) and temperature annealing/quenching (anneal = yes).")
    parser.add_argument('-input', '-f', type=str, help='simulation config file')
    # A bare `topo-mdrun` (no arguments) prints help, like `-h`.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    print(f"OpenMM version: {mm.__version__}")

    cfg = read_simulation_config(args.input)
    # All generated files go into a single run folder (default traj/), named
    # <outname>.* (default traj.dcd, traj.log, traj.psf, ...).
    cfg.prepare_output_dir()

    # 1. Build the (optionally multi-copy) coarse-grain system.
    built = engine.build_system(cfg)

    print('Simulation started')
    start_time = time.time()

    # 2. Set up the OpenMM Simulation (integrator, platform, restart, reporters).
    ctx = engine.setup_simulation(cfg, built, control_file=args.input)

    # 3. Run the temperature protocol. Equilibrium is a single constant-ref_t
    #    stage; annealing holds at t_high then quenches/cools to ref_t. On a
    #    restart, completed stages are skipped via ctx.done_steps.
    schedule = temperature_schedule(cfg)
    print(f"Temperature protocol: {describe_schedule(schedule)}")
    run_protocol(ctx.simulation, schedule, done_steps=ctx.done_steps)

    # 4. Save the final checkpoint + conformation and close out run metadata.
    engine.finalize_simulation(cfg, ctx, built.topology, start_time)


if __name__ == '__main__':
    mdrun()
