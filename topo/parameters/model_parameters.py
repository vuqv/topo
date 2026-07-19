"""
Dictionary contains parameters for TOPO model (topology-based coarse-grained model for folded proteins).

First-level key is the model name. Currently only the "topo" model is defined:
  - topo: topology-based / structure-based model for globular (folded) proteins,
    with residue parameters (mass, Rmin_2, charge) and non-bonded
    interaction matrices.

``Rmin_2`` is the collision radius Rmin/2 (nm): half the pair distance of minimum
non-bonded energy. Pairs combine by the sum rule R_ij = Rmin_2_i + Rmin_2_j
(O'Brien's convention). (Field was named ``radii`` before; renamed for clarity --
it was always an Rmin/2, never a diameter or a Lennard-Jones sigma.)

IMPORTANT -- what these ``Rmin_2`` values are used for. They are the **fixed per-type**
radii for **rigid scenery** beads: the ribosome's RNA (P/R/BR) and ribosomal-protein
(per-AA) beads, read by :func:`topo.csp.ribosome.load_ribosome`. A **mobile protein
chain** (the nascent chain, or a folded-protein simulation) does **not** use these:
its excluded-volume Rmin/2 is **per-residue and structure-derived** (Karanicolas-Brooks,
from ``build_nonbonded_interaction``), so the same residue name (e.g. ``ALA``) has a
K-B value in a nascent chain but this fixed table value in a ribosomal protein. This
mirrors O'Brien, who separates the two by atom type (nascent ``A_i`` vs ribosomal
``S<aa>``) while keeping standard residue names. The protein ``Rmin_2`` below are
O'Brien's per-AA sidechain values (his ribosome .prm ``S<aa>`` types == the topo
``OBRIEN_SC_RMIN_2_NM`` table).

Attributes
----------
parameters : dict
    Dictionary of model parameters, keyed by model name (e.g. "topo").
"""
import numpy as np

from .dihedral import load_dihedral_params

protein_list = ["MET", "GLY", "LYS", "THR", "ARG", "ALA", "ASP", "GLU", "TYR", "VAL", "LEU", "GLN", "TRP", "PHE", "SER",
                "HIS", "ASN", "PRO", "CYS", "ILE"]
nucleic_list = ["A", "C", "G", "U"]
parameters = {
    "topo": {
        "bond_length_protein": 0.381,
        "bond_length_nucleic": 0.5,  # RNA
        "bond_force_constant": 41840.0,  # kJ/mol/nm^2 = 50 kcal/mol/A^2 * 4.184 * 100 * 2.
        # The *2 converts CHARMM's Kb(r-r0)^2 to OpenMM HarmonicBondForce's 1/2*k(r-r0)^2
        # (matches O'Brien's parse_cg_prm.py). Was 20920 (missing the *2 -> 2x too soft);
        # moot under AllBonds constraints (bonds are rigid) but corrects flexible-bond runs.
        "bonded_exclusions_index": 2,
        "ALA": {
            "mass": 71.08,
            "Rmin_2": 0.2862278,
            "charge": 0.0,
        },
        "ARG": {
            "mass": 156.19,
            "Rmin_2": 0.3704125,
            "charge": 1.0,
        },
        "ASN": {
            "mass": 114.10,
            "Rmin_2": 0.3199017,
            "charge": 0.0,
        },
        "ASP": {
            "mass": 115.09,
            "Rmin_2": 0.3142894,
            "charge": -1.0,
        },
        "CYS": {
            "mass": 103.14,
            "Rmin_2": 0.3030648,
            "charge": 0.0,
        },
        "GLU": {
            "mass": 129.12,
            "Rmin_2": 0.3367386,
            "charge": -1.0,
        },
        "GLN": {
            "mass": 128.13,
            "Rmin_2": 0.3423509,
            "charge": 0.0,
        },
        "GLY": {
            "mass": 57.05,
            "Rmin_2": 0.252554,
            "charge": 0.0,
        },
        "HIS": {
            "mass": 137.14,
            "Rmin_2": 0.3423509,
            "charge": 0.0,
        },
        "ILE": {
            "mass": 113.16,
            "Rmin_2": 0.3423509,
            "charge": 0.0,
        },
        "LEU": {
            "mass": 113.16,
            "Rmin_2": 0.3423509,
            "charge": 0.0,
        },
        "LYS": {
            "mass": 128.17,
            "Rmin_2": 0.3535755,
            "charge": 1.0,
        },
        "MET": {
            "mass": 131.19,
            "Rmin_2": 0.3423509,
            "charge": 0.0,
        },
        "PHE": {
            "mass": 147.18,
            "Rmin_2": 0.3535755,
            "charge": 0.0,
        },
        "PRO": {
            "mass": 97.12,
            "Rmin_2": 0.3086771,
            "charge": 0.0,
        },
        "SER": {
            "mass": 87.08,
            "Rmin_2": 0.2918401,
            "charge": 0.0,
        },
        "THR": {
            "mass": 101.10,
            "Rmin_2": 0.3142894,
            "charge": 0.0,
        },
        "TRP": {
            "mass": 186.21,
            "Rmin_2": 0.3816371,
            "charge": 0.0,
        },
        "TYR": {
            "mass": 163.18,
            "Rmin_2": 0.3591879,
            "charge": 0.0,
        },
        "VAL": {
            "mass": 99.13,
            "Rmin_2": 0.3311263,
            "charge": 0.0,
        },
        # Parameters for RNA.
        # Nucleotides containing pyrimidines and purines are represented as 3 and 4
        # interaction sites, respectively, with one interaction at the phosphate
        # position (q = -1e), another at the centroid of the ribose ring, and one at
        # the centroid of each conjugated ring in the base.
        #   pyrimidine bases: C, U (one ring)
        #   purine bases:     A, G (two rings) - two BR beads in the CG model.
        # Be careful with sigma_ij (only used for non-native contacts): two
        # conventions exist, sigma_ij = 0.5 * (sigma_i + sigma_j) or R_ij = R_i + R_j.
        # RNA Rmin_2 are O'Brien's per-type Rmin/2 (his ribosome .prm NONBONDED block,
        # combine_ribo_L24_Yang.prm): P 6.44766 A, R 5.231399 A, BR 5.342436 A (all
        # base beads PU1/PU2/PY share the one BR type). These are used as the ribosome
        # per-bead Rmin/2 in the NC<->ribosome 12-10-6 sum-rule EV (append_ribosome).
        # These bead types are ribosome-only (the nascent Ca chain has no P/R/BR), so
        # setting them here does not affect the nascent-chain model. Previously all three
        # were a single 0.710 nm placeholder (why topo used to load O'Brien's .cor/.prm).
        "P": {
            "mass": 95.00,
            "Rmin_2": 0.644766,
            "charge": -1.0,
        },
        "R": {
            "mass": 92.00,
            "Rmin_2": 0.523140,
            "charge": 0.0,
        },
        "BR": {
            "mass": 64.00,
            "Rmin_2": 0.534244,
            "charge": 0.0,
        },
        # Dihedral parameters loaded from data/dihedral_params.csv
        'dihedral_params': load_dihedral_params(),
    },
}
