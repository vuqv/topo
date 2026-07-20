r"""
Module defines model parameters for TOPO (a unified coarse-grained model for globular
and disordered proteins).

The only supported model is "topo" (structure-based for folded residues; residues declared
disordered switch to the transferable IDR potential -- see topo.utils.nonbonded).
"""

from .model_parameters import parameters, protein_list, nucleic_list

__all__ = ["parameters", "protein_list", "nucleic_list"]