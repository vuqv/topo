Parameters
==========

The TOPO model parameters are defined in :mod:`topo.parameters.model_parameters`. The ``"topo"`` model includes:

- **Per-residue**: mass, radii (sigma), and charge for each residue type (no hydropathy scale).
- **Bonded**: bond length, bond force constant, and bonded exclusion rule (1–2).
- **Dihedral/torsion**: Karanicolas–Brooks style periodic torsion parameters loaded from ``topo/parameters/data/``.

.. automodule:: topo.parameters.model_parameters
