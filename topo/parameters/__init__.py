r"""
Module defines model parameters for TOPO (topology-based coarse-grained model for folded proteins).

The only supported model is "topo" (topology/structure-based model for globular proteins).
"""

from .model_parameters import parameters, protein_list, nucleic_list

__all__ = ["parameters", "protein_list", "nucleic_list"]