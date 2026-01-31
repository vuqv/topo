System 
========================================================= 

A class containing methods and parameters for generating CG systems to be simulated using the OpenMM interface.
    It offers flexibility to create default and custom CG systems and to easily modify their parameters.

.. autoclass:: topo.core.system

        .. automethod:: __init__
        .. automethod:: topo.core.system.getAtoms
        .. automethod:: topo.core.system.getBonds
        .. automethod:: topo.core.system.setBondForceConstants
        .. automethod:: topo.core.system.setParticlesMasses
        .. automethod:: topo.core.system.setParticlesRadii
        .. automethod:: topo.core.system.setParticlesCharge
        .. automethod:: topo.core.system.addYukawaForces
        .. automethod:: topo.core.system.addAshbaughHatchForces
        .. automethod:: topo.core.system.createSystemObject
        .. automethod:: topo.core.system.checkBondDistances
        .. automethod:: topo.core.system.checkLargeForces
        .. automethod:: topo.core.system.addParticles
        .. automethod:: topo.core.system.addSystemForces
        .. automethod:: topo.core.system.dumpStructure
        .. automethod:: topo.core.system.dumpTopology
        .. automethod:: topo.core.system.dumpForceFieldData
        .. automethod:: topo.core.system.setCAMassPerResidueType
        .. automethod:: topo.core.system.setCARadiusPerResidueType
        .. automethod:: topo.core.system.setCAChargePerResidueType
        .. automethod:: topo.core.system.setCAHPSPerResidueType
        .. automethod:: topo.core.system._setParameters