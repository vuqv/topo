System
======

A class containing methods and parameters for generating TOPO coarse-grained systems to be simulated with OpenMM. It is typically constructed via :func:`topo.models.buildCoarseGrainModel`, which sets bonds, angles, torsions, Yukawa electrostatics, and structure-based non-bonded forces.

.. seealso::

   :doc:`../usage/model_theory` for the physics of each force these methods add,
   and :doc:`../usage/python_api` for how to build and inspect a system from
   Python.

.. autoclass:: topo.core.system

   .. automethod:: __init__
   .. automethod:: getCAlphaOnly
   .. automethod:: getAtoms
   .. automethod:: getBonds
   .. automethod:: getAngles
   .. automethod:: getTorsions
   .. automethod:: setBondForceConstants
   .. automethod:: setParticlesMass
   .. automethod:: setParticlesRadii
   .. automethod:: setParticlesCharge
   .. automethod:: addHarmonicBondForces
   .. automethod:: addGaussianAngleForces
   .. automethod:: addPeriodicTorsionForce
   .. automethod:: addYukawaForces
   .. automethod:: addCustomNonBondedForce

   .. automethod:: addIDRNonBondedForce
   .. automethod:: createSystemObject
   .. automethod:: checkBondDistances
   .. automethod:: checkLargeForces
   .. automethod:: addParticles
   .. automethod:: addSystemForces
   .. automethod:: dumpStructure
   .. automethod:: dumpTopology
   .. automethod:: dumpForceFieldData
   .. automethod:: setCAMassPerResidueType
   .. automethod:: setCARadiusPerResidueType
   .. automethod:: setCAChargePerResidueType
   .. automethod:: _setParameters