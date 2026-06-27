"""Shared simulation engine: build the system and set up an OpenMM ``Simulation``.

This is the layer common to every run mode. It does the parts that do not depend
on the *temperature protocol*:

    * :func:`build_system`      -- coarse-grain model, multi-copy replication,
                                   topology (PSF) dump, COM-motion removal.
    * :func:`setup_simulation`  -- integrator, platform, restart/checkpoint,
                                   starting coordinates/velocities, reporters,
                                   run-metadata header.
    * :func:`finalize_simulation` -- save checkpoint, write final structure,
                                   close out run metadata, print wall time.

The temperature *protocol* (constant-temperature equilibrium vs. annealing /
quenching) lives with each runner -- see :mod:`topo.mdrun.protocol`. Equilibrium
is just the degenerate one-stage schedule, so both protocols share everything in
this module. A future nascent-chain (translation) runner can reuse the same
build/setup/finalize helpers and supply its own driver loop.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional
import time
import warnings

import numpy as np
import openmm as mm
from openmm import unit

import parmed as pmd
from parmed.exceptions import OpenMMWarning

import topo

# Suppress OpenMM warnings (same as the original runner).
warnings.filterwarnings("ignore", category=OpenMMWarning)


@dataclass
class BuiltSystem:
    """Outputs of :func:`build_system`."""
    cgModel: Any
    system: Any
    topology: Any
    positions: Any


@dataclass
class RunContext:
    """Outputs of :func:`setup_simulation` needed to run and finalize."""
    simulation: Any
    checkpoint: str
    restart_active: bool
    done_steps: int            # steps already completed (0 unless restart)
    runinfo_path: str


def build_system(cfg) -> BuiltSystem:
    """Build the coarse-grain model and (optionally) replicate it into copies.

    Mirrors the original ``mdrun`` build block: build the single-chain model,
    replicate into ``n_copies`` non-interacting copies when requested, dump the
    PSF topology (and a ``_multi.psf`` for multi-copy runs), and add a
    ``CMMotionRemover`` only when ``nstcomm`` is set.
    """
    cgModel = topo.models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())
    print("Model built successfully...")

    # Single chain by default; replicate into non-interacting copies if requested.
    if cfg.n_copies > 1:
        system, topology, positions = topo.make_noninteracting_copies(
            cgModel.system, cgModel.topology, cgModel.positions,
            n_copies=cfg.n_copies, shift=cfg.copy_shift * unit.nanometer)
    else:
        system, topology, positions = cgModel.system, cgModel.topology, cgModel.positions

    # Write the topology as PSF. The single-chain topology is always <outname>.psf;
    # a multi-copy run additionally writes <outname>_multi.psf for the combined DCD.
    cgModel.dumpTopology(cfg.output_path('.psf'))               # single-chain
    if cfg.n_copies > 1:
        pmd.openmm.load_topology(topology, system=system, xyz=positions).save(
            cfg.output_path('_multi.psf'), overwrite=True)      # multi-chain

    # Remove center-of-mass motion only when nstcomm is set. COM removal suits a
    # single chain but couples the drift of multiple independent chains, so it is
    # off unless explicitly requested in the config.
    if cfg.nstcomm is not None:
        system.addForce(mm.CMMotionRemover(cfg.nstcomm))

    return BuiltSystem(cgModel=cgModel, system=system,
                       topology=topology, positions=positions)


def setup_simulation(cfg, built: BuiltSystem,
                     control_file: Optional[str] = None) -> RunContext:
    """Create the OpenMM ``Simulation`` and bring it to a ready-to-step state.

    Handles platform selection, the Langevin integrator (created at ``ref_t``;
    the protocol driver sets the per-stage temperature), restart/checkpoint
    loading, starting coordinates/velocities for a fresh run, output reporters,
    and the run-metadata header. The integrator's initial temperature is
    irrelevant to the run because the protocol calls ``setTemperature`` before
    the first step, but ``ref_t`` keeps equilibrium runs byte-for-byte identical.
    """
    integrator = mm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)
    # Constraint tolerance applies to rigid (AllBonds) bonds; harmless otherwise.
    integrator.setConstraintTolerance(cfg.constraint_tolerance)

    platform, properties = cfg.make_platform()
    simulation = mm.app.Simulation(built.topology, built.system, integrator,
                                   platform, properties)

    # Resolve the checkpoint path (defaults to <output_dir>/<outname>.chk).
    checkpoint = cfg.checkpoint_path()

    # Coordinate / velocity sourcing:
    #   restart=yes + checkpoint present -> positions AND velocities from checkpoint.
    #   otherwise (fresh start, including restart with a missing checkpoint):
    #       coordinates from init_position if set, else the built model (pdb_file);
    #       velocities drawn from the Boltzmann distribution at ref_t.
    restart_active = cfg.restart and Path(checkpoint).is_file()
    if cfg.restart and not restart_active:
        src = (f"init_position '{cfg.init_position}'" if cfg.init_position
               else f"pdb_file '{cfg.pdb_file}'")
        print(f"[warning] restart=yes but checkpoint not found ({checkpoint}); "
              f"starting fresh from {src} with Boltzmann velocities at {cfg.ref_t}.")

    if restart_active:
        simulation.loadCheckpoint(checkpoint)
        done_steps = simulation.context.getState().getStepCount()
        print(f"Restart simulation from step: {done_steps}")
        # Report the energy of the loaded checkpoint state (what will actually be
        # simulated), not the input structure. The build-time check was skipped
        # on restart (cfg.build_kwargs() sets check_forces=False).
        built.cgModel.reportEnergy(simulation, header='Restart-state potential energy')
        coord_source = f"checkpoint ({checkpoint})"
        vel_source = f"checkpoint ({checkpoint})"
    else:
        # Choose the starting coordinates: explicit init_position, else the
        # coordinates of the structure used to build the system (pdb_file).
        if cfg.init_position:
            init_pos = mm.app.PDBFile(cfg.init_position).getPositions()
            if len(init_pos) != built.system.getNumParticles():
                raise SystemExit(
                    f"init_position '{cfg.init_position}' has {len(init_pos)} atoms but the "
                    f"system has {built.system.getNumParticles()}; they must match.")
            coord_source = f"init_position ({cfg.init_position})"
        else:
            init_pos = built.positions
            coord_source = f"pdb_file ({cfg.pdb_file})"

        # Shift coordinates into the positive octant.
        xyz = np.array(init_pos / unit.nanometer)
        xyz[:, 0] -= np.amin(xyz[:, 0])
        xyz[:, 1] -= np.amin(xyz[:, 1])
        xyz[:, 2] -= np.amin(xyz[:, 2])
        simulation.context.setPositions(xyz * unit.nanometer)
        simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        done_steps = 0
        vel_source = f"Boltzmann distribution at {cfg.ref_t}"

    # Reporters are attached per phase by the runner (see attach_reporters), so a
    # quench phase can write <outname>_quench.* while production writes
    # <outname>.*. Here we only record run provenance (package versions, hardware,
    # GPU, timing) -- a side channel that does not affect the simulation.
    #
    # Planned steps: a fresh run does the whole protocol (quench + production =
    # total_steps); a restart skips the quench entirely and only finishes the
    # remaining production (md_steps - done_steps, since done_steps is now
    # production-relative).
    if restart_active:
        steps_planned = cfg.md_steps - done_steps
    else:
        steps_planned = cfg.total_steps()
    runinfo_path = cfg.output_path('_runinfo.log')
    topo.runinfo.write_run_start(
        runinfo_path,
        control_file=control_file,
        checkpoint_file=checkpoint,
        restart=restart_active,
        steps_planned=steps_planned,
        simulation=simulation,
        use_gpu=(cfg.device == 'GPU'),
        ppn=cfg.ppn,
        coord_source=coord_source,
        vel_source=vel_source,
    )
    print(f"[tracking] writing run metadata to {runinfo_path}")

    return RunContext(simulation=simulation, checkpoint=checkpoint,
                      restart_active=restart_active, done_steps=done_steps,
                      runinfo_path=runinfo_path)


def attach_reporters(cfg, simulation, suffix='', append=False, total_steps=None,
                     checkpoint=True):
    """(Re)attach the checkpoint / trajectory / log reporters for one phase.

    Replaces ``simulation.reporters`` so the runner can switch output file sets
    between the quench phase (``suffix='_quench'``) and production (``suffix=''``).

    Parameters
    ----------
    suffix : str
        Inserted before the extension: ``''`` -> ``<outname>.dcd`` / ``.log``;
        ``'_quench'`` -> ``<outname>_quench.dcd`` / ``.log``.
    append : bool
        Append to (rather than truncate) the DCD/log -- used when a restart
        resumes this phase, so the log header is not re-written.
    total_steps : int, optional
        Value reported as the run length for the log's remaining-time/speed
        columns; defaults to ``cfg.md_steps``. Each phase keeps its own step
        clock (the quench clock is reset before production), so this is the
        phase length, not a grand total.
    checkpoint : bool
        Whether to write the ``<outname>.chk`` checkpoint. The quench phase sets
        this ``False`` -- it is short and never restarted, so the checkpoint only
        ever holds production state and a restart resumes production cleanly.
    """
    if total_steps is None:
        total_steps = cfg.md_steps

    # Checkpoint frequency (nstchk) is independent of the trajectory frequency
    # (nstxout); nstchk defaults to nstxout when not set in the config.
    simulation.reporters = []
    if checkpoint:
        simulation.reporters.append(mm.app.CheckpointReporter(cfg.checkpoint_path(), cfg.nstchk))
    simulation.reporters.append(
        mm.app.DCDReporter(cfg.output_path(suffix + '.dcd'), cfg.nstxout,
                           enforcePeriodicBox=bool(cfg.pbc), append=append))
    # topo.topoReporter writes a clean, fixed-width log: each float column uses
    # log_precision decimals and every column is padded to log_width characters so
    # the columns line up. Columns are separated by two spaces (aligned and still
    # machine-parsable via topo.reporter.readOpenMMReporterFile).
    simulation.reporters.append(
        topo.topoReporter(cfg.output_path(suffix + '.log'), cfg.nstlog,
                          precision=cfg.log_precision, width=cfg.log_width,
                          step=True, time=True,
                          potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                          temperature=True, remainingTime=True, speed=True,
                          totalSteps=total_steps, separator='  ', append=append))


def finalize_simulation(cfg, ctx: RunContext, topology, start_epoch: float) -> None:
    """Save the final checkpoint, write the last conformation, close metadata."""
    simulation = ctx.simulation
    simulation.saveCheckpoint(ctx.checkpoint)

    # Write the last conformation as a PDB so it can seed a later run via
    # init_position. Uses the (possibly multi-copy) simulation topology.
    final_pdb = cfg.output_path('_final.pdb')
    last_positions = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=bool(cfg.pbc)).getPositions()
    with open(final_pdb, 'w') as fh:
        mm.app.PDBFile.writeFile(topology, last_positions, fh)
    print(f"Wrote last conformation to {final_pdb}")

    topo.runinfo.write_run_end(ctx.runinfo_path, simulation=simulation,
                               start_epoch=start_epoch, final_structure=final_pdb)
    print("--- Finished in %s seconds ---" % (time.time() - start_epoch))
