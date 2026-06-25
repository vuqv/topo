# Run a TOPO coarse-grained simulation from a control file (md.ini).
#
# Control-file parsing lives in the package (topo.read_simulation_config). After
# building the single-chain model, if md.ini sets n_copies > 1 the model is
# replicated into that many non-interacting copies with
# topo.make_noninteracting_copies (default n_copies = 1 = single chain).
import argparse
import time
import warnings

import numpy as np
import openmm as mm
from openmm import unit

import parmed as pmd
from parmed.exceptions import OpenMMWarning

import topo

# Suppress OpenMM warnings
warnings.filterwarnings("ignore", category=OpenMMWarning)


def main():
    """
    Run a simulation using the TOPO library and parameters specified in a config file.

    Usage: python run_simulation.py -f md.ini
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-input', '-f', type=str, help='simulation config file')
    args = parser.parse_args()

    cfg = topo.read_simulation_config(args.input)

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
    # no initial/final PDB). The single-chain topology is always {protein_code}.psf
    # (used to load the per-chain trajectories from split_chains.py). A multi-copy
    # run additionally writes {protein_code}_multi.psf for the combined DCD.
    cgModel.dumpTopology(f'{cfg.protein_code}.psf')              # single-chain
    if cfg.n_copies > 1:
        pmd.openmm.load_topology(topology, system=system, xyz=positions).save(
            f'{cfg.protein_code}_multi.psf', overwrite=True)     # multi-chain

    # Remove center-of-mass motion.
    system.addForce(mm.CMMotionRemover(cfg.nstcomm))

    # Select the compute platform (CPU/GPU) from the config.
    platform, properties = cfg.make_platform()

    print('Simulation started')
    start_time = time.time()

    integrator = mm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)
    simulation = mm.app.Simulation(topology, system, integrator, platform, properties)

    if cfg.restart:
        simulation.loadCheckpoint(cfg.checkpoint)
        print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
        nsteps_remain = cfg.md_steps - simulation.context.getState().getStepCount()
    else:
        # Shift coordinates into the positive octant.
        xyz = np.array(positions / unit.nanometer)
        xyz[:, 0] -= np.amin(xyz[:, 0])
        xyz[:, 1] -= np.amin(xyz[:, 1])
        xyz[:, 2] -= np.amin(xyz[:, 2])
        positions = xyz * unit.nanometer
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        nsteps_remain = cfg.md_steps

    simulation.reporters = []
    simulation.reporters.append(mm.app.CheckpointReporter(cfg.checkpoint, cfg.nstxout))
    simulation.reporters.append(
        mm.app.DCDReporter(f'{cfg.protein_code}.dcd', cfg.nstxout,
                           enforcePeriodicBox=bool(cfg.pbc), append=cfg.restart))
    simulation.reporters.append(
        mm.app.StateDataReporter(f'{cfg.protein_code}.log', cfg.nstlog, step=True, time=True,
                                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                 temperature=True, remainingTime=True, speed=True,
                                 totalSteps=cfg.md_steps, separator='\t', append=cfg.restart))
    simulation.step(nsteps_remain)
    simulation.saveCheckpoint(cfg.checkpoint)

    print("--- Finished in %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    print(f"OpenMM version: {mm.__version__}")
    main()
