#!/usr/bin/env python3
"""
custom_md.py -- a hackable folded-protein MD workflow.

This is the *open* version of the ``topo-mdrun`` console command. It reproduces
the exact build -> setup -> run -> finalize flow of ``topo.mdrun.mdrun`` using
TOPO's public building blocks (``topo.engine``), but spells every step out so you
can insert your own logic between them. Copy this file, edit the sections marked

    # === EDIT: ... ===

and run it like the real runner::

    python examples/custom_md.py -f md.ini

It reads the same ``md.ini`` control file, so all the usual options
(``pdb_file``, ``ref_t``, ``md_steps``, ``dt``, ``device``, ``n_copies``, output
naming, ...) still apply. Nothing here is private API -- every call is exported
from ``topo`` or ``topo.engine``.

Three common customizations are demonstrated below:
  1. Add your own OpenMM force to the built system (e.g. a pulling/restraint force).
  2. Define a custom temperature schedule instead of plain constant-T production.
  3. Run in segments with an analysis/callback hook between them.

For the fully-featured anneal/quench and restart logic, see the canonical runner
``topo/mdrun/mdrun.py`` -- this example intentionally keeps the constant-T path
simple so the customization points stand out.
"""
import argparse
import sys
import time

import openmm as mm
from openmm import unit

from topo import engine
from topo.utils.config import read_simulation_config
from topo.mdrun.protocol import run_protocol, describe_schedule


def main():
    parser = argparse.ArgumentParser(
        prog="custom_md.py",
        description="Hackable TOPO folded-protein MD workflow (edit me).")
    parser.add_argument('-input', '-f', type=str, required=True,
                        help='simulation config file (md.ini)')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    print(f"OpenMM version: {mm.__version__}")

    # --- Read the control file and prepare the output folder -----------------
    cfg = read_simulation_config(args.input)
    cfg.prepare_output_dir()

    # === Step 1: build the coarse-grained system =============================
    # Returns a BuiltSystem with .cgModel, .system, .topology, .positions.
    built = engine.build_system(cfg)

    # === EDIT: modify the system before it is handed to OpenMM ===============
    # `built.system` is a plain openmm.System, so you can add custom forces here
    # (they must be added BEFORE setup_simulation creates the Context). Example:
    # a harmonic positional restraint on the first bead --
    #
    #   restraint = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    #   restraint.addGlobalParameter("k", 1000.0)  # kJ/mol/nm^2
    #   for p in ("x0", "y0", "z0"):
    #       restraint.addPerParticleParameter(p)
    #   x0, y0, z0 = built.positions[0].value_in_unit(unit.nanometer)
    #   restraint.addParticle(0, [x0, y0, z0])
    #   built.system.addForce(restraint)
    #
    # (No-op by default.)

    # === Step 2: set up the OpenMM Simulation ===============================
    # Integrator, platform/device, starting coordinates & velocities, restart.
    ctx = engine.setup_simulation(cfg, built, control_file=args.input)
    sim = ctx.simulation

    # === Step 3: attach reporters (DCD/log/checkpoint) ======================
    engine.attach_reporters(cfg, sim, suffix='', append=ctx.restart_active,
                            total_steps=cfg.md_steps)

    # === EDIT: define your temperature protocol =============================
    # A schedule is just a list of (temperature, n_steps) stages. The default
    # below is a single constant-ref_t production run -- identical to a plain
    # `topo-mdrun` equilibrium run. Replace it with your own stages, e.g. a
    # step-cooling ramp:
    #
    #   schedule = [(500 * unit.kelvin, 200000),
    #               (400 * unit.kelvin, 200000),
    #               (cfg.ref_t,         cfg.md_steps)]
    schedule = [(cfg.ref_t, cfg.md_steps)]
    print(f"Temperature protocol: {describe_schedule(schedule)}")

    print('Simulation started')
    start_time = time.time()

    # === Step 4: run ========================================================
    # Simplest form -- run the whole schedule in one call (handles restart):
    run_protocol(sim, schedule, done_steps=ctx.done_steps)

    # === EDIT (alternative): run in segments with a callback ================
    # Comment out the run_protocol call above and use this instead to do
    # something every `segment` steps (on-the-fly analysis, adaptive control,
    # extra logging). This ignores restart for brevity.
    #
    #   segment = 10000
    #   sim.integrator.setTemperature(cfg.ref_t)
    #   done = 0
    #   while done < cfg.md_steps:
    #       n = min(segment, cfg.md_steps - done)
    #       sim.step(n)
    #       done += n
    #       state = sim.context.getState(getEnergy=True)
    #       pe = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    #       print(f"  step {done:>10d}  PE = {pe:12.2f} kJ/mol")

    # === Step 5: finalize ===================================================
    # Save the final checkpoint + conformation, close out run metadata.
    engine.finalize_simulation(cfg, ctx, built.topology, start_time)


if __name__ == '__main__':
    main()
