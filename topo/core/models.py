#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system
from ..utils.build_nonbonded_interaction import build_nonbonded_interaction

class models:
    """
    A class to hold functions for the automated generation of default TOPO models.

    Methods
    -------
    """

    @staticmethod
    def buildCoarseGrainModel(structure_file: str,
                     minimize: bool = False,
                     model: str = 'topo',
                     domain_def: str = None,
                     stride_output_file: str = None,
                     box_dimension: Any = None):
        """
        Build a topology-based coarse-grained model for a folded protein system.

        Creates an alpha-carbon only system with bonds, angles, periodic torsions,
        Yukawa electrostatics, and structure-based (contact) non-bonded interactions.
        Optionally uses domain definitions and STRIDE output for contact-based potentials.

        Parameters
        ----------
        structure_file : str
            Path to the input structure file (PDB/CIF).
        minimize : bool, optional (default: False)
            If True, run energy minimization on the initial structure.
        model : str, optional (default: 'topo')
            Model name; currently only 'topo' is supported.
        domain_def : str, optional
            Path to domain definition file (e.g. domain.yaml) for contact-based non-bonded.
        stride_output_file : str, optional
            Path to STRIDE output file for secondary structure.
        box_dimension : float or array, optional
            If set, use PBC (cubic if float, rectangular if [x,y,z]).

        Returns
        -------
        topo_model : topo.core.system
            Initialized coarse-grained system ready for simulation.
        """

        # common for all model:
        print(f'Generating CA coarse-grained model for structure from file {structure_file}')
        print('')
        topo_model = system(structure_file, model)
        print("Checking input structure file ...")
        print("Be sure that you do not have missing residues in the initial structure. At the moment, I will not take "
              "care of that")

        # Set up geometric parameters of the model
        print('Setting up geometrical parameters ...')
        print('__________________________________________________________________')
        print('Keeping only alpha carbon atoms in topology')
        topo_model.getCAlphaOnly()

        print(f'There are {topo_model.n_chains} chain(s) in the input file.')

        # set particle's properties
        # Common for all
        topo_model.getAtoms()
        print('Added ' + str(topo_model.n_atoms) + ' CA atoms')

        topo_model.getBonds()
        print('Added ' + str(topo_model.n_bonds) + ' bonds')
        # Add constraints to all bonds
        for bond in topo_model.bonds:
            topo_model.system.addConstraint(bond[0].index, bond[1].index, topo_model.bonds[bond][0])

        print("Setting alpha-carbon masses to their average residue mass.")
        topo_model.setCAMassPerResidueType()

        print("Setting alpha-carbon charge to their residue charge.")
        topo_model.setCAChargePerResidueType()

        
        
        # set particle interactions
        # add forces to system
        print('Adding default bond force constant:', end=' ')
        topo_model.setBondForceConstants()
        print('')
        print('__________________________________________________________________')

        # all models have bonded interactions
        print('Adding Forces:')
        topo_model.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')
        print("---")


        # this model has angle bonded potential.
        # angle
        topo_model.getAngles()
        print(f'Added {topo_model.n_angles} angles ')
        topo_model.addGaussianAngleForces()
        print('Added Gaussian Angle Forces')
        print("---")

        # add Periodic Torsion angle for topo model
        topo_model.getTorsions()
        topo_model.addPeriodicTorsionForce()

        if box_dimension:
            use_pbc = True
            if isinstance(box_dimension, list):
                """
                OpenMM use this to write dimension in PDB and dcd file. Require one-argument, so zip box dimension into 
                one variable.
                Rectangular box, given parameter is array of three number
                """
                topo_model.topology.setPeriodicBoxVectors(
                    ((box_dimension[0], 0, 0), (0, box_dimension[1], 0), (0, 0, box_dimension[2])))
            else:
                # cubic box, given parameter is single float
                topo_model.topology.setPeriodicBoxVectors(
                    ((box_dimension, 0, 0), (0, box_dimension, 0), (0, 0, box_dimension)))

            unit_cell = topo_model.topology.getPeriodicBoxVectors()
            # use this to write coordinate in PBC box. requires 3 numbers, unzip to 3
            topo_model.system.setDefaultPeriodicBoxVectors(*unit_cell)

        else:
            use_pbc = False

        topo_model.addYukawaForces(use_pbc)
        print('Added Yukawa Force')
        print("---")

        # non-bonded interaction
        print("Building non-bonded interactions for TOPO model...")
        try:
            distance_matrix, energy_matrix = build_nonbonded_interaction(
                structure_file,
                domain_def,
                stride_output_file # TODO: If Stride output is None, then run stride on structure_file to generate the stride output file
            )
            print(f"Built non-bonded interaction matrices: {distance_matrix.shape}, {energy_matrix.shape}")
            
            # Store the matrices in the topo_model object for later use
            topo_model.distance_matrix = distance_matrix
            topo_model.energy_matrix = energy_matrix
            
            # You can also add them to the system if needed
            topo_model.addCustomNonBondedForce(distance_matrix, energy_matrix, use_pbc)
            
        except Exception as e:
            print(f"Warning: Could not build non-bonded interactions: {e}")
            print("Continuing with default non-bonded interactions...")
            distance_matrix = None
            energy_matrix = None




        print('')
        print('__________________________________________________________________')

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        topo_model.createSystemObject(minimize=minimize, check_bond_distances=True)
        print('topo system Object created')
        print('')

        return topo_model
