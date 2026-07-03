#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system
from ..utils.nonbonded import build_nonbonded_interaction

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
                     box_dimension: Any = None,
                     constraints: Any = 'AllBonds',
                     check_forces: bool = True):
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
        constraints : str or None, optional (default: 'AllBonds')
            Controls the treatment of covalent bonds. These two modes are mutually
            exclusive (a bond is never both constrained and harmonic):

            - 'AllBonds' : rigid bonds. A distance constraint is added at each
              bond's equilibrium length and no harmonic bond force is created.
            - None (or the string 'None'/'none') : flexible bonds. A harmonic
              bond force is created and no constraints are added.
        check_forces : bool, optional (default: True)
            If True, run the build-time energy/large-force check on the initial
            input structure. Set False when restarting from a checkpoint, where
            the input-structure energy is irrelevant (the loaded state, not the
            PDB geometry, is what gets simulated).

        Returns
        -------
        topo_model : topo.core.system
            Initialized coarse-grained system ready for simulation.
        """

        print('')
        print('=' * 66)
        print('[ System build ]')
        print('=' * 66)
        print(f'Building CA coarse-grained model (model={model}) from {structure_file}')

        topo_model = system(structure_file, model)

        # Build alpha-carbon topology (atoms, bonds, angles).
        topo_model.getCAlphaOnly()
        topo_model.getAtoms()
        topo_model.getBonds()
        topo_model.getAngles()

        # Resolve the bond-constraint mode. Accepted values: 'AllBonds' (rigid,
        # default) or None / 'None' / 'none' (flexible). Rigid and flexible are
        # mutually exclusive, so a bond is never both constrained and harmonic.
        if constraints is None or str(constraints).strip().lower() == 'none':
            use_constraints = False
        elif str(constraints).strip().lower() == 'allbonds':
            use_constraints = True
        else:
            raise ValueError(
                f"Invalid constraints option: {constraints!r}. "
                f"Expected 'AllBonds' or None.")

        bond_mode = 'rigid (AllBonds)' if use_constraints else 'flexible (harmonic)'
        print(f'  chains={topo_model.n_chains}  CA atoms={topo_model.n_atoms}  '
              f'bonds={topo_model.n_bonds}  angles={topo_model.n_angles}')
        print(f'  bonds: {bond_mode}')

        # Rigid bonds: constrain every bond at its equilibrium length (no harmonic
        # bond force is added later in that case).
        if use_constraints:
            for bond in topo_model.bonds:
                topo_model.system.addConstraint(bond[0].index, bond[1].index, topo_model.bonds[bond][0])

        # Per-residue particle properties. Mass/charge are per-AA (model_parameters).
        # The excluded-volume radius is the per-residue K-B Rmin/2 (structure-derived),
        # set from the contact build below -- NOT the fixed per-AA model_parameters value
        # (that is the rigid ribosome scenery radius). See setParticlesRadii after the
        # build_nonbonded_interaction call.
        topo_model.setCAMassPerResidueType()
        topo_model.setCAChargePerResidueType()

        # set particle interactions and add forces to system
        topo_model.setBondForceConstants()

        # Only add the harmonic bond force when bonds are flexible. With rigid
        # bonds (constraints=AllBonds) the distance is pinned, so a harmonic term
        # would be redundant (it contributes ~0 energy/force) and would also show
        # a misleading non-zero bond energy on the unconstrained input geometry.
        if not use_constraints:
            topo_model.addHarmonicBondForces()

        topo_model.addGaussianAngleForces()

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

        # non-bonded interaction
        # The structure-based (contact) non-bonded term is a core part of the TOPO
        # model: without it there are no native/Go contacts, no H-bond/sidechain
        # energies, and no domain scaling. A failure here must be fatal rather than
        # silently swallowed, otherwise the simulation runs an incomplete force field.
        try:
            rmin_matrix, energy_matrix, rmin_2 = build_nonbonded_interaction(
                structure_file,
                domain_def,
                stride_output_file,
                return_rmin_2=True,
            )
        except Exception as e:
            raise RuntimeError(
                "Failed to build the TOPO structure-based non-bonded interactions "
                f"(domain_def={domain_def!r}, stride_output_file={stride_output_file!r}): {e}"
            ) from e

        print(f'  contact matrices: {rmin_matrix.shape}')

        # Store the matrices on the model for later use
        topo_model.rmin_matrix = rmin_matrix
        topo_model.energy_matrix = energy_matrix

        # Excluded-volume radius = per-residue K-B Rmin/2 (structure-derived), matching the
        # non-native contacts. particle_rmin_2 feeds only dumpForceFieldData; this keeps that dump
        # consistent with the actual per-residue radii rather than a fixed per-AA lookup.
        topo_model.setParticlesRadii(list(rmin_2))

        # Add the custom non-bonded (contact) force to the system
        topo_model.addCustomNonBondedForce(rmin_matrix, energy_matrix, use_pbc)

        # Generate the system object and add previously generated forces. The
        # bond-distance check always runs (it validates the built geometry); the
        # large-force / initial-energy check is skipped on restart (check_forces
        # is False there) since the loaded checkpoint state is what matters.
        topo_model.createSystemObject(minimize=minimize, check_bond_distances=True,
                                      check_large_forces=check_forces)

        return topo_model
