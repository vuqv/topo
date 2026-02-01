System
======

A class containing methods and parameters for generating TOPO coarse-grained systems to be simulated with OpenMM. It is typically constructed via :func:`topo.models.buildCoarseGrainModel`, which sets bonds, angles, torsions, Yukawa electrostatics, and structure-based non-bonded forces.

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
   .. automethod:: setCAHPSPerResidueType
   .. automethod:: _setParameters