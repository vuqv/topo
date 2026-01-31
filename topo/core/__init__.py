"""
core package of the TOPO package that contains the main TOPO classes.

The topo.core package contains the three main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures.

The second class, models, allows to easily set up predefined CG models.
The final class, system, is the main class that holds all the methods to define, modify and create CG to be simulated
with OpenMM.

"""
__all__ = ['system', 'models']
from .geometry import geometry
from .models import models
from .system import system
