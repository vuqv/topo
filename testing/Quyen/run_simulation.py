# Import OpenMM library
# Import hpsOpenMM library
import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import numpy as np
import openmm as mm
from openmm import unit

from parmed.exceptions import OpenMMWarning

import topo
from topo.reporter import topoReporter

# Suppress OpenMM warnings
warnings.filterwarnings("ignore", category=OpenMMWarning)


def main():
    """
        Run a simulation using the hpsOpenMM library and parameters specified in a config file.

        Usage: python run_simulation.py -f md.ini
        or cosmo-simulation -f md.ini
        """

    # Default values (aligned with md.ini):
    md_steps = 1000
    dt = 0.01
    nstxout = 10
    nstlog = 10
    nstcomm = 100
    model = 'topo'
    tcoupl = True
    ref_t = 300.0
    pcoupl = False
    ref_p = 1.0
    frequency_p = 25
    pbc = False
    device = 'CPU'
    ppn = 1
    restart = False
    minimize = True

    # Other attributes
    tau_t = None
    box_dimension = None
    protein_code = None
    checkpoint = None
    pdb_file = None
    domain_def = None
    stride_output_file = None

    # Parse config file
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-input', '-f', type=str, help='simulation config file')
    args = parser.parse_args()
    config_file = args.input

    # Read config file
    print(f"Reading simulation parameters from {config_file} file...")
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(config_file)
    params = config['OPTIONS']

    # Update simulation parameters
    md_steps = int(str(params.get('md_steps', md_steps)).replace('_', ''))
    print(f'Setting number of simulation steps to: {md_steps}')
    dt = float(params.get('dt', dt)) * unit.picoseconds
    print(f'Setting timestep for integration of equations of motion to: {dt}')
    nstxout = int(params.get('nstxout', nstxout))
    print(f'Setting number of steps to write checkpoint and coordinate: {nstxout}')
    nstlog = int(params.get('nstlog', nstlog))
    print(f'Setting number of steps to write logfile: {nstlog}')
    nstcomm = int(params.get('nstcomm', nstcomm))
    print(f'Setting frequency of center of mass motion removal to every {nstcomm} steps')
    model = params.get('model', model)
    print(f'Setting model: {model}')
    tcoupl = bool(strtobool(str(params.get('tcoupl', tcoupl))))
    if tcoupl:
        ref_t = float(params.get('ref_t', ref_t)) * unit.kelvin
        tau_t = float(params.get('tau_t', 0.01)) / unit.picoseconds
        print(
            f'Turning on temperature coupling with reference temperature: {ref_t} and time constant: {tau_t}')
    else:
        print("Temperature coupling is off")
    pbc = bool(strtobool(str(params.get('pbc', pbc))))
    if pbc:
        box_val = params.get('box_dimension', '').strip()
        if box_val.startswith('['):
            box_dimension = loads(box_val)
        else:
            try:
                L = float(box_val)
                box_dimension = [L, L, L]
            except (ValueError, TypeError):
                box_dimension = None
        if box_dimension is not None:
            print(f'Turning on periodic boundary conditions with box dimension: {box_dimension} nm')
        else:
            print('Periodic boundary conditions are off (invalid box_dimension)')
            box_dimension = None
    else:
        box_dimension = None
        print('Periodic boundary conditions are off')
    pcoupl = bool(strtobool(str(params.get('pcoupl', pcoupl))))
    if pcoupl:
        assert pbc, "Pressure coupling requires box dimensions and periodic boundary condition is on"
        ref_p = float(params.get('ref_p', ref_p)) * unit.bar
        frequency_p = int(params.get('frequency_p', frequency_p))
        print(f'Pressure is set to reference of {ref_p} with frequency of coupling {frequency_p}')
    else:
        print("Pressure coupling is off")
    pdb_file = params.get('pdb_file', pdb_file)
    print(f'Input structure: {pdb_file}')
    protein_code = params.get('protein_code', protein_code)
    print(f'Prefix use to write file: {protein_code}')
    domain_def = params.get('domain_def', domain_def)
    if domain_def:
        print(f'Domain definition file: {domain_def}')
    stride_output_file = params.get('stride_output_file', stride_output_file)
    if stride_output_file:
        print(f'STRIDE output file: {stride_output_file}')
    checkpoint = params.get('checkpoint', checkpoint)
    device = params.get('device', device)
    print(f'Running simulation on {device}')
    if device == "CPU":
        ppn = int(params.get('ppn', ppn))
        print(f'Using {ppn} threads')
    restart = bool(strtobool(str(params.get('restart', restart))))
    print(f'Restart simulation: {restart}')
    if restart:
        minimize = False
    else:
        minimize = bool(strtobool(str(params.get('minimize', minimize))))
    print(f'Perform Energy minimization of input structure: {minimize}')

    """
    End of reading parameters
    """

    build_kwargs = dict(minimize=minimize, model=model, box_dimension=box_dimension)
    if domain_def is not None:
        build_kwargs['domain_def'] = domain_def
    if stride_output_file is not None:
        build_kwargs['stride_output_file'] = stride_output_file
    cgModel = topo.models.buildCoarseGrainModel(pdb_file, **build_kwargs)
    print("Model built successfully...")


    # Remove center of mass motion
    cgModel.system.addForce(mm.CMMotionRemover(nstcomm))

    # dump Forcefield File
    # if model in ['topo']:
    #     """ current dumpForceFieldData function can only write the standard format of forcefield which require
    #      sigma, epsilon for each residue.
    #     """
    #     cgModel.dumpForceFieldData(f'{protein_code}_ff.dat')
    # dump Structure into PDB file for visualize
    cgModel.dumpStructure(f'{protein_code}_init.pdb')
    cgModel.dumpTopology(f'{protein_code}.psf')

    if device == 'GPU':
        # Run simulation on CUDA
        print(f"Running simulation on GPU CUDA")
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
        # in case of many GPUs present, we can select which one to use

    elif device == 'CPU':
        print(f"Running simulation on CPU using {ppn} cores")
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {'Threads': str(ppn)}

    print('Simulation started')
    start_time = time.time()

    integrator = mm.LangevinIntegrator(ref_t, tau_t, dt)
    simulation = mm.app.Simulation(cgModel.topology, cgModel.system, integrator, platform,
                                   properties)
    if restart:
        simulation.loadCheckpoint(checkpoint)
        print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
        nsteps_remain = md_steps - simulation.context.getState().getStepCount()
    else:
        xyz = np.array(cgModel.positions / unit.nanometer)
        xyz[:, 0] -= np.amin(xyz[:, 0])
        xyz[:, 1] -= np.amin(xyz[:, 1])
        xyz[:, 2] -= np.amin(xyz[:, 2])
        cgModel.positions = xyz * unit.nanometer
        simulation.context.setPositions(cgModel.positions)
        simulation.context.setVelocitiesToTemperature(ref_t)
        nsteps_remain = md_steps

    simulation.reporters = []
    simulation.reporters.append(mm.app.CheckpointReporter(checkpoint, nstxout))
    simulation.reporters.append(
        mm.app.DCDReporter(f'{protein_code}.dcd', nstxout, enforcePeriodicBox=bool(pbc),
                           append=restart))
    simulation.reporters.append(
        topoReporter(f'{protein_code}.log', nstlog, sbmObject=cgModel, step=True, time=True,
                     potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                     temperature=True, remainingTime=True, speed=True,
                     totalSteps=md_steps, separator='\t', append=restart))
    simulation.step(nsteps_remain)

    # write the last frame
    last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(pbc)).getPositions()
    mm.app.PDBFile.writeFile(cgModel.topology, last_frame, open(f'{protein_code}_final.pdb', 'w'))
    simulation.saveCheckpoint(checkpoint)

    print("--- Finished in %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    print(f"OpenMM version: {mm.__version__}")
    main()
