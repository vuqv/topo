"""
Dictionary contains parameters for TOPO model (topology-based coarse-grained model for folded proteins).

First-level key is the model name. Currently only the "topo" model is defined:
  - topo: topology-based / structure-based model for globular (folded) proteins,
    with residue parameters (mass, radii, charge) and non-bonded
    interaction matrices.

Attributes
----------
parameters : dict
    Dictionary of model parameters, keyed by model name (e.g. "topo").
"""
import numpy as np

from .load_dihedral_params import load_dihedral_params

protein_list = ["MET", "GLY", "LYS", "THR", "ARG", "ALA", "ASP", "GLU", "TYR", "VAL", "LEU", "GLN", "TRP", "PHE", "SER",
                "HIS", "ASN", "PRO", "CYS", "ILE", "ALY", "SEP", "TPO", "PTR"]
nucleic_list = ["A", "C", "G", "U"]
parameters = {
    "topo": {
        "bond_length_protein": 0.381,
        "bond_length_nucleic": 0.5,  # RNA
        "bond_force_constant": 20920.0,  # kj/mol/nm^2, converted from 50kcal/mol/A^2
        "bonded_exclusions_index": 2,
        "ALA": {
            "mass": 71.00,
            "radii": 0.504,
            "charge": 0.0,
        },
        "ARG": {
            "mass": 114.00, # should be around 156.19
            "radii": 0.656,
            "charge": 1.0,
        },
        "ASN": {
            "mass": 114.00,
            "radii": 0.568,
            "charge": 0.0,
        },
        "ASP": {
            "mass": 114.00,
            "radii": 0.558,
            "charge": -1.0,
        },
        "CYS": {
            "mass": 114.00,  # concern, this should be around 103.10
            "radii": 0.548,
            "charge": 0.0,
        },
        "GLU": {
            "mass": 128.00,
            "radii": 0.592,
            "charge": -1.0,
        },
        "GLN": {
            "mass": 128.00,
            "radii": 0.602,
            "charge": 0.0,
        },
        "GLY": {
            "mass": 57.00,
            "radii": 0.450,
            "charge": 0.0,
        },
        "HIS": {
            "mass": 114.00,  # should be around 137
            "radii": 0.608,
            "charge": 0.0,
        },
        "ILE": {
            "mass": 113.00,
            "radii": 0.618,
            "charge": 0.0,
        },
        "LEU": {
            "mass": 113.00,
            "radii": 0.618,
            "charge": 0.0,
        },
        "LYS": {
            "mass": 128.00,
            "radii": 0.636,
            "charge": 1.0,
        },
        "MET": {
            "mass": 131.00,
            "radii": 0.618,
            "charge": 0.0,
        },
        "PHE": {
            "mass": 147.00,
            "radii": 0.636,
            "charge": 0.0,
        },
        "PRO": {
            "mass": 114.00,  # should be about 97.12
            "radii": 0.556,
            "charge": 0.0,
        },
        "SER": {
            "mass": 87.00,
            "radii": 0.518,
            "charge": 0.0,
        },
        "THR": {
            "mass": 101.00,
            "radii": 0.562,
            "charge": 0.0,
        },
        "TRP": {
            "mass": 186.00,
            "radii": 0.678,
            "charge": 0.0,
        },
        "TYR": {
            "mass": 163.00,
            "radii": 0.646,
            "charge": 0.0,
        },
        "VAL": {
            "mass": 99.00,
            "radii": 0.586,
            "charge": 0.0,
        },
        # Parameters for RNA 
        """
        nucleotides containing pyrimidines and purines were represented as 3 and 4 interaction sites, respectively, 
        with one interaction located at the phosphate position and having a q=-1e charge, another at the centroid of the ribose ring, 
        and one at the centroid of each conjugated ring in the base.
        pyrimidine bases: C, T (one ring)
        purine bases: A, G (two rings)- in CG model will have two BR beads.

        Be careful on sigma_ij: (only used for non-native contacts)
        Ed's paper: sigma_ij = 0.5 * (sigma_i + sigma_j)
        Yang's code: R_ij =  (R_i + R_j)
        """
        "P": {
            "mass": 95.00,
            "radii": 0.710,
            "charge": -1.0,
        },
        "R": {
            "mass": 92.00,
            "radii": 0.710,
            "charge": 0.0,
        },
        "BR": {
            "mass": 64.00,
            "radii": 0.710,
            "charge": 0.0,
        },
        # Dihedral parameters loaded from data/dihedral_params.csv
        'dihedral_params': load_dihedral_params(),
    },
}
