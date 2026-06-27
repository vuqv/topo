"""Post-simulation trajectory analysis for TOPO CG models.

Currently provides native-contact order parameters (Q), per domain and per
interface; future order parameters (dRMS, RMSD, inter-domain distances) belong
here too.
"""

from .native_contacts import (
    load_domains,
    reference_residue_geometry,
    build_native_contacts,
    fraction_native_contacts,
)

__all__ = [
    'load_domains',
    'reference_residue_geometry',
    'build_native_contacts',
    'fraction_native_contacts',
]
