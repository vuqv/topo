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

    * ``anneal = no``  (default) -- constant-temperature equilibrium at ``ref_t``,
      written to ``<outname>.dcd`` / ``.log``.
    * ``anneal = yes``           -- two phases. A **quench** phase holds at
      ``t_high`` (and, for a linear ramp, cools to ``ref_t``) and is written to
      ``<outname>_quench.dcd`` / ``.log``; a **production** phase then runs
      ``md_steps`` at ``ref_t`` and is written to the usual ``<outname>.dcd`` /
      ``.log``. ``anneal_steps`` is separate from ``md_steps`` (grand total =
      ``quench_steps + md_steps``); ``ref_t`` is the low / refold temperature.

Both share the same build / setup / finalize machinery (:mod:`topo.engine`);
only the temperature schedule and output files differ (:mod:`topo.mdrun.protocol`).
"""
import argparse
import sys
import time

import openmm as mm

from topo import engine
from topo.utils.config import read_simulation_config
from .protocol import (quench_schedule, production_schedule, run_protocol,
                       describe_schedule)


def mdrun():
    """
    Run a simulation using the TOPO library and parameters specified in a config file.

    Usage: python -m topo.mdrun -f md.ini
    """
    # Silence MDAnalysis's cosmetic UserWarnings (missing CRYST1 unit cell, absent
    # formalcharges) for the CLI run -- a process-local filter that only fires when
    # someone actually invokes topo-mdrun.
    import warnings
    warnings.filterwarnings("ignore", category=Warning, module=r"MDAnalysis")

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

    # 2. Set up the OpenMM Simulation (integrator, platform, restart). Reporters
    #    are attached per phase below so the quench and production phases write
    #    separate files.
    ctx = engine.setup_simulation(cfg, built, control_file=args.input)
    sim = ctx.simulation

    # 3. Run the temperature protocol.
    prod = production_schedule(cfg)
    if cfg.anneal:
        qsched = quench_schedule(cfg)
        print(f"Temperature protocol [quench -> _quench.*]: {describe_schedule(qsched)}")
        print(f"Temperature protocol [production -> .*]:    {describe_schedule(prod)}")
        if ctx.restart_active:
            # A restart resumes PRODUCTION only. The quench is a one-time prep
            # whose step clock was reset, so the checkpoint is production-relative.
            print(f"Restart: skipping quench, resuming production from step {ctx.done_steps}")
            prod_done = ctx.done_steps
        else:
            # Fresh run: run the quench phase (writes <outname>_quench.*, no
            # checkpoint -- it is short and never restarted), then reset the
            # step/time clock so production is a clean run from step 0. Positions
            # and velocities carry over in the same context, so a `jump` is simply
            # the hot structure suddenly thermostatted at ref_t.
            engine.attach_reporters(cfg, sim, suffix='_quench', append=False,
                                    total_steps=cfg.quench_steps(), checkpoint=False)
            run_protocol(sim, qsched, done_steps=0)
            sim.context.setStepCount(0)
            sim.context.setTime(0.0)
            sim.currentStep = 0
            prod_done = 0

        # --- Production phase (writes <outname>.*) at ref_t, step clock from 0 ---
        engine.attach_reporters(cfg, sim, suffix='', append=ctx.restart_active,
                                total_steps=cfg.md_steps, checkpoint=True)
        run_protocol(sim, prod, done_steps=prod_done)
    else:
        # Plain equilibrium: a single constant-ref_t production run.
        print(f"Temperature protocol: {describe_schedule(prod)}")
        engine.attach_reporters(cfg, sim, suffix='', append=ctx.restart_active,
                                total_steps=cfg.md_steps)
        run_protocol(sim, prod, done_steps=ctx.done_steps)

    # 4. Save the final checkpoint + conformation and close out run metadata.
    engine.finalize_simulation(cfg, ctx, built.topology, start_time)


if __name__ == '__main__':
    mdrun()
