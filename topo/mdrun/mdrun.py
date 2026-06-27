"""
Run a TOPO coarse-grained simulation from a control file (md.ini).

This is the canonical runner for the package. Use it as a CLI::

    python -m topo.mdrun -f md.ini

or call :func:`mdrun` from your own script. Control-file parsing lives in
``topo.read_simulation_config``. After building the single-chain model, if
``md.ini`` sets ``n_copies > 1`` the model is replicated into that many
non-interacting copies with ``topo.make_noninteracting_copies`` (default
``n_copies = 1`` = single chain).
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit

import parmed as pmd
from parmed.exceptions import OpenMMWarning

import topo

# Suppress OpenMM warnings
warnings.filterwarnings("ignore", category=OpenMMWarning)


def mdrun():
    """
    Run a simulation using the TOPO library and parameters specified in a config file.

    Usage: python -m topo.mdrun -f md.ini
    """
    parser = argparse.ArgumentParser(
        prog="topo-mdrun",
        description="Run a TOPO coarse-grained simulation from a control file "
                    "(md.ini). Replicates into non-interacting copies when "
                    "n_copies > 1.")
    parser.add_argument('-input', '-f', type=str, help='simulation config file')
    # A bare `topo-mdrun` (no arguments) prints help, like `-h`.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    print(f"OpenMM version: {mm.__version__}")

    cfg = topo.read_simulation_config(args.input)
    # All generated files go into a single run folder (default traj/), named
    # <outname>.* (default traj.dcd, traj.log, traj.psf, ...).
    cfg.prepare_output_dir()

    # Build the single-chain coarse-grained model.
    cgModel = topo.models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())
    print("Model built successfully...")

    # Single chain by default; replicate into non-interacting copies if requested.
    if cfg.n_copies > 1:
        system, topology, positions = topo.make_noninteracting_copies(
            cgModel.system, cgModel.topology, cgModel.positions,
            n_copies=cfg.n_copies, shift=cfg.copy_shift * unit.nanometer)
    else:
        system, topology, positions = cgModel.system, cgModel.topology, cgModel.positions

    # Write the topology as PSF (we keep only the DCD trajectory + PSF topology;
    # no initial/final PDB). The single-chain topology is always <outname>.psf
    # (used to load the per-chain trajectories from topo.split_chains). A multi-copy
    # run additionally writes <outname>_multi.psf for the combined DCD.
    cgModel.dumpTopology(cfg.output_path('.psf'))               # single-chain
    if cfg.n_copies > 1:
        pmd.openmm.load_topology(topology, system=system, xyz=positions).save(
            cfg.output_path('_multi.psf'), overwrite=True)      # multi-chain

    # Remove center-of-mass motion only when nstcomm is set. COM removal is
    # appropriate for a single chain but couples the drift of multiple
    # independent chains, so it is off unless explicitly requested in the config.
    if cfg.nstcomm is not None:
        system.addForce(mm.CMMotionRemover(cfg.nstcomm))

    # Select the compute platform (CPU/GPU) from the config.
    platform, properties = cfg.make_platform()

    print('Simulation started')
    start_time = time.time()

    integrator = mm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)
    # Constraint tolerance applies to rigid (AllBonds) bonds; harmless otherwise.
    integrator.setConstraintTolerance(cfg.constraint_tolerance)
    simulation = mm.app.Simulation(topology, system, integrator, platform, properties)

    # Resolve the checkpoint path (defaults to <output_dir>/<outname>.chk).
    checkpoint = cfg.checkpoint_path()

    # Coordinate / velocity sourcing:
    #   restart=yes + checkpoint present -> positions AND velocities from checkpoint.
    #   otherwise (fresh start, including restart with a missing checkpoint):
    #       coordinates from init_position if set, else the built model (pdb_file);
    #       velocities drawn from the Boltzmann distribution at ref_t.
    restart_active = cfg.restart and Path(checkpoint).is_file()
    if cfg.restart and not restart_active:
        src = f"init_position '{cfg.init_position}'" if cfg.init_position else f"pdb_file '{cfg.pdb_file}'"
        print(f"[warning] restart=yes but checkpoint not found ({checkpoint}); "
              f"starting fresh from {src} with Boltzmann velocities at {cfg.ref_t}.")

    if restart_active:
        simulation.loadCheckpoint(checkpoint)
        print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
        # Report the energy of the loaded checkpoint state (what will actually be
        # simulated), not the input structure. The build-time check was skipped
        # on restart (cfg.build_kwargs() sets check_forces=False).
        cgModel.reportEnergy(simulation, header='Restart-state potential energy')
        nsteps_remain = cfg.md_steps - simulation.context.getState().getStepCount()
        coord_source = f"checkpoint ({checkpoint})"
        vel_source = f"checkpoint ({checkpoint})"
    else:
        # Choose the starting coordinates: explicit init_position, else the
        # coordinates of the structure used to build the system (pdb_file).
        if cfg.init_position:
            init_pos = mm.app.PDBFile(cfg.init_position).getPositions()
            if len(init_pos) != system.getNumParticles():
                raise SystemExit(
                    f"init_position '{cfg.init_position}' has {len(init_pos)} atoms but the "
                    f"system has {system.getNumParticles()}; they must match.")
            coord_source = f"init_position ({cfg.init_position})"
        else:
            init_pos = positions
            coord_source = f"pdb_file ({cfg.pdb_file})"

        # Shift coordinates into the positive octant.
        xyz = np.array(init_pos / unit.nanometer)
        xyz[:, 0] -= np.amin(xyz[:, 0])
        xyz[:, 1] -= np.amin(xyz[:, 1])
        xyz[:, 2] -= np.amin(xyz[:, 2])
        simulation.context.setPositions(xyz * unit.nanometer)
        simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        nsteps_remain = cfg.md_steps
        vel_source = f"Boltzmann distribution at {cfg.ref_t}"

    simulation.reporters = []
    # Checkpoint frequency (nstchk) is independent of the trajectory frequency
    # (nstxout); nstchk defaults to nstxout when not set in the config.
    simulation.reporters.append(mm.app.CheckpointReporter(checkpoint, cfg.nstchk))
    simulation.reporters.append(
        mm.app.DCDReporter(cfg.output_path('.dcd'), cfg.nstxout,
                           enforcePeriodicBox=bool(cfg.pbc), append=restart_active))
    # topo.topoReporter writes a clean, fixed-width log: each float column uses
    # log_precision decimals and every column is padded to log_width characters so
    # the columns line up. Columns are separated by two spaces (aligned and still
    # machine-parsable via topo.reporter.readOpenMMReporterFile).
    simulation.reporters.append(
        topo.topoReporter(cfg.output_path('.log'), cfg.nstlog,
                          precision=cfg.log_precision, width=cfg.log_width,
                          step=True, time=True,
                          potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                          temperature=True, remainingTime=True, speed=True,
                          totalSteps=cfg.md_steps, separator='  ', append=restart_active))

    # Record run provenance (package versions, hardware, GPU, timing). This is a
    # side channel only and does not affect the simulation.
    runinfo_path = cfg.output_path('_runinfo.log')
    topo.runinfo.write_run_start(
        runinfo_path,
        control_file=args.input,
        checkpoint_file=checkpoint,
        restart=restart_active,
        steps_planned=nsteps_remain,
        simulation=simulation,
        use_gpu=(cfg.device == 'GPU'),
        ppn=cfg.ppn,
        coord_source=coord_source,
        vel_source=vel_source,
    )
    print(f"[tracking] writing run metadata to {runinfo_path}")

    simulation.step(nsteps_remain)
    simulation.saveCheckpoint(checkpoint)

    # Write the last conformation as a PDB so it can seed a later run via
    # init_position. Uses the (possibly multi-copy) simulation topology.
    final_pdb = cfg.output_path('_final.pdb')
    last_positions = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=bool(cfg.pbc)).getPositions()
    with open(final_pdb, 'w') as fh:
        mm.app.PDBFile.writeFile(topology, last_positions, fh)
    print(f"Wrote last conformation to {final_pdb}")

    topo.runinfo.write_run_end(runinfo_path, simulation=simulation,
                               start_epoch=start_time, final_structure=final_pdb)
    print("--- Finished in %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    mdrun()
