Parameters
==========

The TOPO force-field constants are defined in
:mod:`topo.parameters.model_parameters`. The ``"topo"`` model provides:

- **Per-residue properties** — ``mass`` (amu), ``radii`` (excluded-volume radius,
  nm), and ``charge`` (e) for each residue type. Only ASP/GLU (−1) and ARG/LYS
  (+1) are charged; all others are neutral. The full table is in
  :doc:`../usage/model_theory`.
- **Bonded constants** — ``bond_length_protein`` = 0.381 nm,
  ``bond_force_constant`` = 20920 kJ mol\ :sup:`-1` nm\ :sup:`-2`, and the
  ``bonded_exclusions_index`` = 2 rule (pairs **two or fewer bonds apart** — i.e.
  1–2 bonded and 1–3 angle neighbours — are excluded from the non-bonded forces).
- **Sequence-dependent dihedral parameters** — periodic-torsion phases and force
  constants keyed by the two central residues of each dihedral, loaded from
  ``topo/parameters/data/dihedral_params.csv`` by
  :func:`topo.parameters.dihedral.load_dihedral_params` (Karanicolas–Brooks-style
  values, scaled by a 0.756 calibration factor).

The sidechain–sidechain contact energies use the Betancourt–Thirumalai (BT)
residue-pair potential in ``topo/parameters/data/bt_potential.csv`` (loaded by
:func:`topo.utils.nonbonded.load_bt_potential`); see :ref:`theory-contacts`.

.. seealso::

   :doc:`../usage/model_theory` for the functional forms that consume these
   constants.


API reference
+++++++++++++

.. automodule:: topo.parameters.model_parameters

.. autofunction:: topo.parameters.dihedral.load_dihedral_params
