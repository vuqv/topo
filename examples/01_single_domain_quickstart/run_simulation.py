# Run a TOPO coarse-grained simulation from a control file (md.ini).
#
# All control-file parsing now lives in the package:
#     cfg = topo.read_simulation_config("md.ini")   -> topo.SimulationConfig
# so this script only contains the simulation logic itself.
import argparse
import time
import warnings

import numpy as np
import openmm as mm
from openmm import unit

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

    # Read every simulation parameter from the control file. See
    # topo.read_simulation_config / topo.SimulationConfig for the full list.
    cfg = topo.read_simulation_config(args.input)

    # Build the coarse-grained model from the input structure.
    cgModel = topo.models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())
    print("Model built successfully...")

    # Remove center-of-mass motion.
    cgModel.system.addForce(mm.CMMotionRemover(cfg.nstcomm))

    # Dump the CA structure + topology for visualization / analysis.
    cgModel.dumpStructure(f'{cfg.protein_code}_init.pdb')
    cgModel.dumpTopology(f'{cfg.protein_code}.psf')

    # Select the compute platform.
    if cfg.device == 'GPU':
        print("Running simulation on GPU CUDA")
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
    else:  # CPU
        print(f"Running simulation on CPU using {cfg.ppn} cores")
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {'Threads': str(cfg.ppn)}

    print('Simulation started')
    start_time = time.time()

    integrator = mm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)
    simulation = mm.app.Simulation(cgModel.topology, cgModel.system, integrator,
                                   platform, properties)

    if cfg.restart:
        simulation.loadCheckpoint(cfg.checkpoint)
        print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
        nsteps_remain = cfg.md_steps - simulation.context.getState().getStepCount()
    else:
        # Shift coordinates so the structure sits in the positive octant.
        xyz = np.array(cgModel.positions / unit.nanometer)
        xyz[:, 0] -= np.amin(xyz[:, 0])
        xyz[:, 1] -= np.amin(xyz[:, 1])
        xyz[:, 2] -= np.amin(xyz[:, 2])
        cgModel.positions = xyz * unit.nanometer
        simulation.context.setPositions(cgModel.positions)
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

    # Write the last frame.
    last_frame = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=bool(cfg.pbc)).getPositions()
    mm.app.PDBFile.writeFile(cgModel.topology, last_frame, open(f'{cfg.protein_code}_final.pdb', 'w'))
    simulation.saveCheckpoint(cfg.checkpoint)

    print("--- Finished in %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    print(f"OpenMM version: {mm.__version__}")
    main()
