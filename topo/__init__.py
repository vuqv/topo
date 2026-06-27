"""
TOPO: TOPOlogy-based coarse-grained model for folded prOteins

TOPO is a Python library to run coarse-grained simulations of globular proteins using the OpenMM toolkit.
The library offers flexibility for creating CG models that can be customised to implement different potential models.

Considering an input structure, the library automatizes the creation of forces to specify it.

The library offers methods to tailor forcefield parameters for each force term.

TOPO is divided in three main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures.

The library is open-source and offers flexibility.
"""
__all__ = ['system', 'models', 'runinfo',
           'topoReporter',
           'read_simulation_config', 'SimulationConfig',
           'make_noninteracting_copies', 'replicate_system_intra_only',
           'replicate_topology', 'replicate_positions', 'split_indices',
           'split_chains',
           'analysis',
           'load_domains', 'reference_residue_geometry', 'build_native_contacts',
           'fraction_native_contacts']

from .core import geometry
from .core import models
from .core import system
from .parameters import model_parameters
from .reporter import topoReporter
from .utils import runinfo
from .utils.config import read_simulation_config, SimulationConfig
from .utils.multichain import (make_noninteracting_copies,
                               replicate_system_intra_only, replicate_topology,
                               replicate_positions, split_indices,
                               split_chains)
from . import analysis
from .analysis.native_contacts import (load_domains, reference_residue_geometry,
                                       build_native_contacts,
                                       fraction_native_contacts)
