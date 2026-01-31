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
__all__ = ['system', 'models', 'dynamics', 'build_structure']

from .core import geometry
from .core import models
from .core import system
from .dynamics import dynamics
from .parameters import model_parameters
from .reporter import topoReporter
from .utils import build_structure
