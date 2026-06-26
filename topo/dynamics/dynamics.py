"""
The Dynamics class is used to perform molecular dynamics simulations of a protein system using the
TOPO (topology-based coarse-grained) force field. The class has a constructor that takes a configuration file as an input,
which contains the parameters for the simulation. The class has several attributes such as the number of steps,
the time step, the frequency for writing coordinates and energies, the coarse-grained model, the temperature and
pressure coupling, the reference temperature and pressure, the frequency for pressure coupling, the PBC, the device,
the number of threads, whether to restart the simulation or not, and whether to minimize energy or not.
The class also has several methods such as reading the config file, setting up the system, running the simulation, and
writing the output files.
The class also uses the parmed library to handle the system and the openmm library to perform the simulations.
"""
import configparser
import time
import warnings
from pathlib import Path
from distutils.util import strtobool
from json import loads

import numpy as np
import openmm as mm
import openmm.unit as unit
from parmed.exceptions import OpenMMWarning

from ..core import models
from ..reporter import topoReporter

warnings.filterwarnings("ignore", category=OpenMMWarning)


class Dynamics:
    """
    The Dynamics class is used to perform molecular dynamics simulations.
    It contains two main functions: reading a configuration file and running a simulation.
    To use the class, a user needs to provide a configuration file (e.g md.ini)
    and specify the parameters that control the simulation.


    Parameters
    ----------
    config_file : str
        The path to the configuration file that contains the parameters for the simulation.

    Attributes
    ----------
    md_steps : int, optional (default: 1000)
        The number of steps to perform in the molecular dynamics simulation.
    dt : float, optional (default: 0.01) [ps]
        The time step for integration in picoseconds.
    nstxout : int, optional (default: 10)
        The number of steps between writing coordinates to the output trajectory file. The last coordinates are always written.
    nstlog : int, optional (default: 10)
        The number of steps between writing energies to the log file. The last energies are always written.
    nstcomm : int, optional (default: 100)
        The frequency for center of mass motion removal.
    model : str, optional (default: 'topo')
        Coarse-grained model used in the simulation. Currently only 'topo' is supported.
    tcoupl : bool, optional (default: False)
        Indicates whether temperature coupling is used in the simulation.
    ref_t : float, optional (default: 300.0) [K]
        The reference temperature for coupling in Kelvin.
    tau_t : float, optional (default: 0.1) [ps]
        The time constant for temperature coupling in picoseconds.
    pcoupl : bool, optional (default: False)
        Indicates whether pressure coupling is used in the simulation.
    ref_p : float, optional (default: 1.0) [bar]
        The reference pressure for coupling in bar.
    frequency_p : int, optional (default: 25)
        The frequency for coupling the pressure.
    pbc : bool, optional (default: True)
        Indicates whether periodic boundary conditions are used in the simulation.
    box_dimension : float or list of float, optional (default: None)
        The dimension of the box used in the simulation. It is better to use rectangular for simplicity.
    checkpoint : str, optional (default: None)
        The name of the checkpoint file.
    pdb_file : str, optional (default: None)
        The input structure used to generate the model.
    device : str, optional (default: 'GPU')
        The device used to perform the simulation. Options are 'GPU' or 'CPU'.
    ppn : int, optional (default: 1)
        In case the simulation is run on a CPU, this parameter controls the number of threads used to run the simulation.
    restart : bool, optional (default: False)
        Indicates whether the simulation should be run from the beginning or restarted from a checkpoint.
    minimize : bool, optional (default: True)
        Indicates whether energy minimization should be performed at the start of the simulation


    Returns
    -------

    """

    def __init__(self, config_file):
        # Set default values
        self.md_steps = 1000
        self.dt = 0.01
        self.nstxout = 10
        self.nstlog = 10
        self.nstcomm = 100
        self.log_precision = 4
        self.log_width = 14
        self.constraint_tolerance = 1e-5
        self.model = 'topo'
        self.tcoupl = True
        self.ref_t = 300.0
        self.pcoupl = False
        self.ref_p = 1.0
        self.frequency_p = 25
        self.pbc = False
        self.device = 'CPU'
        self.ppn = 1
        self.restart = False
        self.minimize = True

        # Other attributes
        self.tau_t = None
        self.box_dimension = None
        self.checkpoint = None
        self.pdb_file = None
        self.init_position = None   # optional PDB of starting coordinates
        # Output: all generated files go to <output_dir>/<outname>.* (default traj/).
        self.output_dir = 'traj'
        self.outname = 'traj'

        # Initialize attributes with config file
        self.read_config(config_file)

        self.cgModel = None

    def read_config(self, config_file):
        """
           Read simulation control parameters from config file *.ini into class attributes.

           Parameters
           ----------
           config_file: string, default=None
               Control file store simulation parameters.

           TODO: check parameters in control file more carefully.
                   Raise error and exit immediately if something wrong.
           """
        config = configparser.ConfigParser()
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
        config.read(config_file)
        params = config['OPTIONS']

        self.md_steps = int(params.get('md_steps', self.md_steps))
        self.dt = float(params.get('dt', self.dt)) * unit.picoseconds
        self.nstxout = int(params['nstxout'])
        self.nstlog = int(params['nstlog'])
        self.nstcomm = int(params.get('nstcomm', self.nstcomm))
        prec_val = params.get('log_precision', None)
        if prec_val is not None:
            self.log_precision = None if str(prec_val).strip().lower() in ('none', '') else int(prec_val)
        width_val = params.get('log_width', None)
        if width_val is not None:
            self.log_width = None if str(width_val).strip().lower() in ('none', '') else int(width_val)
        self.constraint_tolerance = float(params.get('constraint_tolerance', self.constraint_tolerance))
        self.model = params.get('model', self.model)
        self.tcoupl = bool(strtobool(params.get('tcoupl', self.tcoupl)))
        if self.tcoupl:
            self.ref_t = float(params['ref_t']) * unit.kelvin
            self.tau_t = float(params['tau_t']) / unit.picoseconds
        self.pbc = bool(strtobool(params.get('pbc', self.pbc)))
        if self.pbc:
            self.box_dimension = loads(params['box_dimension'])
        else:
            self.box_dimension = None
        self.pcoupl = bool(strtobool(params.get('pcoupl', self.pcoupl)))
        if self.pcoupl:
            assert self.pbc, "Pressure coupling requires box dimensions and periodic boundary condition is on"
            self.ref_p = float(params['ref_p']) * unit.bar
            self.frequency_p = int(params['frequency_p'])
        self.pdb_file = params['pdb_file']
        self.init_position = params.get('init_position', self.init_position)
        self.output_dir = params.get('output_dir', self.output_dir)
        self.outname = params.get('outname', self.outname)
        self.checkpoint = params.get('checkpoint', self.checkpoint)
        self.device = params.get('device', self.device)
        if self.device == "CPU":
            self.ppn = int(params.get('ppn', self.ppn))
        self.restart = bool(strtobool(params.get('restart', self.restart)))
        if self.restart:
            self.minimize = False
        else:
            self.minimize = bool(strtobool(params.get('minimize', self.minimize)))

        # Concise summary of the resolved run configuration.
        device_str = self.device + (f' ({self.ppn} threads)' if self.device == 'CPU' else '')
        tcoupl_str = f'ref_t={self.ref_t}, tau_t={self.tau_t}' if self.tcoupl else 'off'
        pbc_str = f'box={self.box_dimension} nm' if self.pbc else 'off'
        pcoupl_str = f'ref_p={self.ref_p}, freq={self.frequency_p}' if self.pcoupl else 'off'

        print('=' * 66)
        print('TOPO: TOPOlogy-based coarse-grained model for folded prOteins')
        print(f'OpenMM {mm.__version__}  |  config: {config_file}')
        print('=' * 66)
        print('[ Simulation parameters ]')
        print(f'  model={self.model}  steps={self.md_steps}  dt={self.dt}')
        print(f'  output: coord/checkpoint every {self.nstxout}, log every {self.nstlog}, '
              f'com-removal every {self.nstcomm}')
        print(f'  T-coupling: {tcoupl_str}')
        print(f'  PBC: {pbc_str}   P-coupling: {pcoupl_str}')
        print(f'  input={self.pdb_file}  output={self._output_path("")}.*')
        print(f'  device={device_str}  restart={self.restart}  minimize={self.minimize}')
        print('=' * 66)
        """
        End of reading parameters
        """

    def _output_path(self, suffix=''):
        """Path for a generated output file: <output_dir>/<outname><suffix>."""
        return str(Path(self.output_dir) / f'{self.outname}{suffix}')

    def run(self):
        """
        Run simulation
        Returns
        -------

        """

        # Ensure the output folder exists (mkdir is a no-op if it already does).
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        checkpoint = self.checkpoint or self._output_path('.chk')

        # Initialize model
        self.cgModel = models.buildCoarseGrainModel(self.pdb_file, minimize=self.minimize, model=self.model,
                                              box_dimension=self.box_dimension)
        # dump topology PSF file and initial coordinate pdb file
        if self.model in ['topo']:
            """ current dumpForceFieldData function can only write the standard format of forcefield which require
             sigma, epsilon for each residue.
            """
            self.cgModel.dumpForceFieldData(self._output_path('_ff.dat'))
        self.cgModel.dumpStructure(self._output_path('_init.pdb'))
        self.cgModel.dumpTopology(self._output_path('.psf'))

        # In case we want to serialize the system to inspect, use the follow functionality of openmm:
        #  mm.XmlSerializer.serialize(system)-it may change in new version

        # Remove center of mass motion
        self.cgModel.system.addForce(mm.CMMotionRemover(self.nstcomm))

        # setup integrator and simulation object
        integrator = mm.LangevinIntegrator(self.ref_t, self.tau_t, self.dt)
        integrator.setConstraintTolerance(self.constraint_tolerance)
        if self.pcoupl:
            # if pressure coupling is on, add barostat force to the system.
            barostat = mm.MonteCarloBarostat(self.ref_p, self.ref_t, self.frequency_p)
            self.cgModel.system.addForce(barostat)

        # Setup platform to run simulation
        if self.device == 'GPU':
            # Run simulation on CUDA
            print(f"Running simulation on GPU CUDA")
            platform = mm.Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'mixed', "DeviceIndex": "0"}
            # in case of many GPUs present, we can select which one to use

        elif self.device == 'CPU':
            print(f"Running simulation on CPU using {self.ppn} cores")
            platform = mm.Platform.getPlatformByName('CPU')
            properties = {'Threads': str(self.ppn)}

        simulation = mm.app.Simulation(self.cgModel.topology, self.cgModel.system, integrator, platform,
                                       properties)
        start_time = time.time()
        if self.restart:
            simulation.loadCheckpoint(checkpoint)
            print(f"Restart simulation from step: {simulation.context.getState().getStepCount()}")
            nsteps_remain = self.md_steps - simulation.context.getState().getStepCount()
        else:
            # Starting coordinates: explicit init_position, else the built model.
            if self.init_position:
                init_pos = mm.app.PDBFile(self.init_position).getPositions()
                if len(init_pos) != self.cgModel.system.getNumParticles():
                    raise ValueError(
                        f"init_position '{self.init_position}' has {len(init_pos)} atoms but the "
                        f"system has {self.cgModel.system.getNumParticles()}; they must match.")
            else:
                init_pos = self.cgModel.positions
            xyz = np.array(init_pos / unit.nanometer)
            xyz[:, 0] -= np.amin(xyz[:, 0])
            xyz[:, 1] -= np.amin(xyz[:, 1])
            xyz[:, 2] -= np.amin(xyz[:, 2])
            self.cgModel.positions = xyz * unit.nanometer
            simulation.context.setPositions(self.cgModel.positions)
            simulation.context.setVelocitiesToTemperature(self.ref_t)
            nsteps_remain = self.md_steps

        simulation.reporters = []
        simulation.reporters.append(mm.app.CheckpointReporter(checkpoint, self.nstxout))
        simulation.reporters.append(
            mm.app.DCDReporter(self._output_path('.dcd'), self.nstxout, enforcePeriodicBox=bool(self.pbc),
                               append=self.restart))
        simulation.reporters.append(
            topoReporter(self._output_path('.log'), self.nstlog,
                         precision=self.log_precision, width=self.log_width,
                         step=True, time=True,
                         potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                         temperature=True, remainingTime=True, speed=True,
                         totalSteps=self.md_steps, separator='  ', append=self.restart))
        simulation.step(nsteps_remain)

        # write the last frame
        last_frame = simulation.context.getState(getPositions=True, enforcePeriodicBox=bool(self.pbc)).getPositions()
        with open(self._output_path('_final.pdb'), 'w') as f:
            mm.app.PDBFile.writeFile(self.cgModel.topology, last_frame, f)
        simulation.saveCheckpoint(checkpoint)
        print("--- Finished in %s seconds ---" % (time.time() - start_time))
