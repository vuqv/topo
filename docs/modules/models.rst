Models
======

The ``models`` module provides a single entry point for building the TOPO
coarse-grained system: :meth:`~topo.core.models.models.buildCoarseGrainModel`. It
reads a structure, keeps the alpha-carbon atoms, and assembles the full force
field (bonds, angles, periodic torsions, Yukawa electrostatics, and the
structure-based contact potential), returning a :class:`topo.core.system.system`
object ready to simulate.

.. note::

   The **complete force-field theory** — every energy term, its functional form,
   all constants, and where the parameters come from — lives in
   :doc:`../usage/model_theory`. For programmatic use (arguments, the returned
   object's attributes and methods) see :doc:`../usage/python_api`.


Coarse-grained alpha-carbon (CA) model
++++++++++++++++++++++++++++++++++++++

TOPO represents the protein as one bead per residue, centred on the alpha
carbon. The Hamiltonian is

.. math::

   U_\mathrm{total} = \sum_\mathrm{bonds} U_\mathrm{bond}
       + \sum_\mathrm{angles} U_\mathrm{angle}
       + \sum_\mathrm{torsions} U_\mathrm{torsion}
       + \sum_{i<j} U^\mathrm{el}_{ij}
       + \sum_{i<j} U^\mathrm{nb}_{ij}

with the terms summarised below (full details and constants:
:doc:`../usage/model_theory`).

**Bonds.** Rigid distance constraints at :math:`r_0 = 0.381` nm by default
(``constraints = AllBonds``); or, with ``constraints = None``, a harmonic force
:math:`U_\mathrm{bond} = \tfrac12 k_b (r-r_0)^2`,
:math:`k_b = 20920\ \mathrm{kJ\,mol^{-1}\,nm^{-2}}`.

**Angles (bimodal Gaussian).**

.. math::

   U_\mathrm{angle}(\theta) = -\frac{1}{\gamma}
   \ln\!\Big[ e^{-\gamma[\,k_\alpha (\theta-\theta_\alpha)^2 + \varepsilon_\alpha\,]}
        + e^{-\gamma k_\beta (\theta-\theta_\beta)^2} \Big]

with two basins at :math:`\theta_\alpha = 91.7^\circ` and
:math:`\theta_\beta = 130.0^\circ`.

**Torsions (periodic, sequence-dependent).** A standard periodic torsion with
four periodicities,

.. math::

   U_\mathrm{torsion}(\varphi) =
   \sum_{n=1}^{4} k_{D,n}\,\big[\,1 + \cos(n\varphi - \delta_n)\,\big],

where :math:`k_{D,n}` and :math:`\delta_n` depend on the **two central residues**
of the dihedral (table ``topo/parameters/data/dihedral_params.csv``, scaled by
0.756). This is implemented by
:meth:`~topo.core.system.system.addPeriodicTorsionForce`.

**Electrostatics (Debye–Hückel / Yukawa).**

.. math::

   U^\mathrm{el}_{ij}(r) = f\,\frac{q_i q_j}{\varepsilon_r\, r}\, e^{-r/l_D}

with :math:`f = 138.935458`, :math:`\varepsilon_r = 78.5`,
:math:`l_D = 1.0` nm, a **2.0 nm cutoff** and switching at 1.8 nm. Only ASP/GLU
(−1) and ARG/LYS (+1) carry charge.

**Structure-based contacts (12-10-6).**

.. math::

   U^\mathrm{nb}_{ij}(r) = \varepsilon_{ij}
   \Big[\,13(R_{ij}/r)^{12} - 18(R_{ij}/r)^{10} + 4(R_{ij}/r)^{6}\,\Big]

with the well minimum :math:`-\varepsilon_{ij}` at :math:`r = R_{ij}`. For native
contacts :math:`R_{ij}` is the native Cα–Cα distance and
:math:`\varepsilon_{ij}` sums hydrogen-bond (STRIDE), backbone–sidechain, and
(domain-scaled) sidechain–sidechain energies; non-native pairs get a soft
excluded-volume repulsion. Built by
:func:`topo.utils.nonbonded.build_nonbonded_interaction`; see
:ref:`theory-contacts`.

**Non-bonded exclusions.** Pairs two or fewer bonds apart (1–2 bonded and 1–3
angle neighbours) are excluded from both non-bonded forces
(``bonded_exclusions_index = 2``).


API reference
+++++++++++++

.. autoclass:: topo.core.models.models
   :members:
   :undoc-members:
