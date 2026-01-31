"""
Dictionary contains parameters for TOPO model (topology-based coarse-grained model for folded proteins).

First-level key is the model name. Currently only the "topo" model is defined:
  - topo: topology-based / structure-based model for globular (folded) proteins,
    with residue parameters (mass, radii, charge, hps, eps_di) and non-bonded
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
        "bonded_exclusions_index": 3,
        "ALA": {
            "mass": 71.00,
            "radii": 0.504,
            "charge": 0.0,
            "hps": 0.522942,
            "eps_di": -2.59  # parameter control torsion angle.
        },
        "ARG": {
            "mass": 114.00, # should be around 156.19
            "radii": 0.656,
            "charge": 1.0,
            "hps": 0.478824,
            "eps_di": -1.37
        },
        "ASN": {
            "mass": 114.00,
            "radii": 0.568,
            "charge": 0.0,
            "hps": 0.508236,
            "eps_di": -0.42
        },
        "ASP": {
            "mass": 114.00,
            "radii": 0.558,
            "charge": -1.0,
            "hps": 0.214119,
            "eps_di": -0.80
        },
        "CYS": {
            "mass": 114.00,  # concern, this should be around 103.10
            "radii": 0.548,
            "charge": 0.0,
            "hps": 0.56706,
            "eps_di": -0.15
        },
        "GLU": {
            "mass": 128.00,
            "radii": 0.592,
            "charge": -1.0,
            "hps": -0.08,
            "eps_di": -1.80
        },
        "GLN": {
            "mass": 128.00,
            "radii": 0.602,
            "charge": 0.0,
            "hps": 0.478824,
            "eps_di": -1.25
        },
        "GLY": {
            "mass": 57.00,
            "radii": 0.450,
            "charge": 0.0,
            "hps": 0.49353,
            "eps_di": 0.65
        },
        "HIS": {
            "mass": 114.00,  # should be around 137
            "radii": 0.608,
            "charge": 0.0,
            "hps": 0.684707,
            "eps_di": 0.8
        },
        "ILE": {
            "mass": 113.00,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.625883,
            "eps_di": -1.39
        },
        "LEU": {
            "mass": 113.00,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.640589,
            "eps_di": -2.05
        },
        "LYS": {
            "mass": 128.00,
            "radii": 0.636,
            "charge": 1.0,
            "hps": 0.302354,
            "eps_di": -0.95
        },
        "MET": {
            "mass": 131.00,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.596471,
            "eps_di": -1.60
        },
        "PHE": {
            "mass": 147.00,
            "radii": 0.636,
            "charge": 0.0,
            "hps": 0.74353,
            "eps_di": -0.68
        },
        "PRO": {
            "mass": 114.00,  # should be about 97.12
            "radii": 0.556,
            "charge": 0.0,
            "hps": 0.678824,
            "eps_di": 3.70
        },
        "SER": {
            "mass": 87.00,
            "radii": 0.518,
            "charge": 0.0,
            "hps": 0.508236,
            "eps_di": -0.69
        },
        "THR": {
            "mass": 101.00,
            "radii": 0.562,
            "charge": 0.0,
            "hps": 0.508236,
            "eps_di": -0.30
        },
        "TRP": {
            "mass": 186.00,
            "radii": 0.678,
            "charge": 0.0,
            "hps": 0.92,
            "eps_di": -1.15
        },
        "TYR": {
            "mass": 163.00,
            "radii": 0.646,
            "charge": 0.0,
            "hps": 0.817059,
            "eps_di": -0.68
        },
        "VAL": {
            "mass": 99.00,
            "radii": 0.586,
            "charge": 0.0,
            "hps": 0.584707,
            "eps_di": -0.75
        },
        # Dihedral parameters loaded from data/dihedral_params.csv
        'dihedral_params': load_dihedral_params(),
    },
}
