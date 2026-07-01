The TOPO model: theory and force field
======================================

This page explains **what TOPO actually simulates** — the coarse-graining choice
and every term in the potential energy function, with its functional form,
constants, and where the parameters come from. It is written for a reader who
knows classical molecular dynamics but has never used TOPO. If you only want to
*run* a simulation, the :doc:`tutorials <../tutorials/index>` and the
:doc:`simulation_control` reference are enough; come here when you want to know
*why* the model behaves the way it does, or to cite the energy terms.

Everything below is implemented in :mod:`topo.core.system` (the force objects)
and :mod:`topo.utils.nonbonded` (the contact energies); the numbers quoted are
the values hard-coded there.


What kind of model is this?
---------------------------

TOPO is a **structure-based (Gō-like) coarse-grained model** for **globular,
folded proteins**.

* **Coarse-graining — one bead per residue.** TOPO reads an all-atom PDB/CIF
  structure and keeps **only the alpha-carbon (Cα) atom of each residue**
  (:meth:`~topo.core.system.system.getCAlphaOnly`). A 106-residue protein
  becomes 106 beads. Each bead carries the mass, excluded-volume radius, and
  charge of its amino-acid type (table below). Consecutive Cα beads *of the same
  chain* are bonded in sequence.

* **Structure-based — the native fold is the energy minimum.** Unlike a transferable
  force field, the energy function is built *from the input structure itself*:
  the pairs of residues that are in contact in your folded PDB are given
  attractive wells centred at their native distances, and everything else is
  repulsive. The crystal/NMR structure you provide therefore **defines** the
  global energy minimum. This is what makes TOPO efficient for studying
  **folding, unfolding, domain motions, and thermal/mechanical stability** — the
  thermodynamic reference state is known by construction.

* **Implicit solvent, implicit ions.** There is no explicit water. Solvent enters
  through the Langevin thermostat (friction + random force) and through screened
  electrostatics (Debye–Hückel, below). A simulation is therefore typically run
  **without a periodic box** (``pbc = no``); the box and barostat options exist
  for the rarer cases where you want them.


The potential energy function
------------------------------

The total potential energy is a sum of **bonded** terms (chain geometry) and
**non-bonded** terms (electrostatics + structure-based contacts):

.. math::

   U_\mathrm{total}
   = \underbrace{\sum_\mathrm{bonds} U_\mathrm{bond}
              + \sum_\mathrm{angles} U_\mathrm{angle}
              + \sum_\mathrm{torsions} U_\mathrm{torsion}}_{\text{bonded (local geometry)}}
   + \underbrace{\sum_{i<j} U^\mathrm{el}_{ij}
              + \sum_{i<j} U^\mathrm{nb}_{ij}}_{\text{non-bonded (long-range)}}

Each term maps to one OpenMM force object, and each appears as its own column in
the run log (see :doc:`outputs`), so you can monitor them separately:

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - Term
     - OpenMM force
     - Log column (force group)
   * - Harmonic bonds *(flexible mode only)*
     - ``HarmonicBondForce``
     - ``Harmonic Bond Energy``
   * - Angles (Gaussian)
     - ``CustomAngleForce``
     - ``Gaussian Angle Energy``
   * - Torsions (periodic)
     - ``PeriodicTorsionForce``
     - ``Periodic Torsion Energy``
   * - Electrostatics (Yukawa)
     - ``CustomNonbondedForce``
     - ``Yukawa Energy``
   * - Structure-based contacts
     - ``CustomNonbondedForce``
     - ``Custom Non-Bonded Energy``

.. note::

   **Units.** TOPO works internally in OpenMM's MD unit system: **nm** for
   length, **kJ/mol** for energy, **ps** for time, **K** for temperature, **bar**
   for pressure, and elementary charge **e**. Many parameters were originally
   expressed in kcal/mol and Å; both the original and converted values are given
   below.


Bonded terms (chain geometry)
-----------------------------

These three terms reproduce the local geometry of the Cα backbone: the
fixed Cα–Cα spacing, the backbone bending preference, and the
sequence-dependent torsional preference.

Bonds
~~~~~

Adjacent Cα beads are held at the canonical virtual-bond length. TOPO offers two
mutually exclusive treatments, selected by ``constraints`` in ``md.ini``:

* **Rigid (default, ``constraints = AllBonds``).** Each bond is a *holonomic
  distance constraint* at the equilibrium length — there is **no harmonic bond
  force**. Removing the stiff bond-stretch vibration is what lets TOPO use a
  large **15 fs (0.015 ps) time step**.

* **Flexible (``constraints = None``).** Each bond is instead a harmonic spring
  and there are no constraints:

  .. math::

     U_\mathrm{bond}(r) = \tfrac{1}{2}\,k_b\,(r - r_0)^2

  with force constant :math:`k_b = 20920\ \mathrm{kJ\,mol^{-1}\,nm^{-2}}`
  (= 50 kcal mol\ :sup:`-1` Å\ :sup:`-2`). Flexible bonds are physically softer
  but require a smaller time step.

In both cases the equilibrium length is :math:`r_0 = 0.381\ \mathrm{nm}` for
protein Cα–Cα bonds (0.5 nm is reserved for nucleic backbones).

Angles — a bimodal Gaussian potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every triplet of consecutive bonded beads (i–j–k) gets a backbone-angle term.
Real protein backbones populate **two** Cα–Cα–Cα angle basins — a tighter one
typical of helices and a wider one typical of extended/β geometry — so a single
harmonic well is inadequate. TOPO uses a **double-Gaussian (log-sum-exp)** angle
potential that interpolates smoothly between the two basins:

.. math::

   U_\mathrm{angle}(\theta) =
   -\frac{1}{\gamma}\,
   \ln\!\Big[
       e^{-\gamma\,[\,k_\alpha (\theta - \theta_\alpha)^2 + \varepsilon_\alpha\,]}
     + e^{-\gamma\,k_\beta (\theta - \theta_\beta)^2}
   \Big]

(the logarithm is natural). The two wells sit at :math:`\theta_\alpha` (the
"α/helical" basin) and :math:`\theta_\beta` (the "β/extended" basin);
:math:`\varepsilon_\alpha` offsets their relative depth and :math:`\gamma`
controls how sharply the lower of the two wells dominates.

.. list-table::
   :header-rows: 1
   :widths: 22 30 30

   * - Parameter
     - Value (MD units)
     - Value (original)
   * - :math:`\theta_\alpha`
     - 1.60047 rad
     - 91.7°
   * - :math:`\theta_\beta`
     - 2.26893 rad
     - 130.0°
   * - :math:`k_\alpha`
     - 445.18 kJ mol\ :sup:`-1` rad\ :sup:`-2`
     - 106.4 kcal mol\ :sup:`-1` rad\ :sup:`-2`
   * - :math:`k_\beta`
     - 110.04 kJ mol\ :sup:`-1` rad\ :sup:`-2`
     - 26.3 kcal mol\ :sup:`-1` rad\ :sup:`-2`
   * - :math:`\varepsilon_\alpha`
     - 17.99 kJ mol\ :sup:`-1`
     - 4.3 kcal mol\ :sup:`-1`
   * - :math:`\gamma`
     - 0.023901 mol kJ\ :sup:`-1`
     - 0.1 mol kcal\ :sup:`-1`

The same parameters are used for every angle in the chain (the potential is not
residue-specific); sequence specificity enters through the torsion and contact
terms instead.

Torsions — sequence-dependent periodic dihedrals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every quadruplet of consecutive bonded beads (i–j–k–l) gets a backbone dihedral
term. TOPO uses a standard **periodic torsion** with **four periodicities**:

.. math::

   U_\mathrm{torsion}(\varphi) =
   \sum_{n=1}^{4} k_{D,n}\,\big[\,1 + \cos(n\,\varphi - \delta_n)\,\big]

The force constants :math:`k_{D,n}` and phases :math:`\delta_n` are
**sequence-dependent**: they are looked up by the **two central residues** of the
dihedral (beads *j* and *k*) from a parameter table
(``topo/parameters/data/dihedral_params.csv``), so an Ala–Gly junction torsions
differently from a Val–Ile junction. These are Karanicolas–Brooks-style
knowledge-based dihedral parameters; TOPO applies a global **0.756 calibration
factor** to every tabulated :math:`k_{D,n}` (see
:func:`topo.parameters.dihedral.load_dihedral_params`).

.. admonition:: Heads-up for readers of older docs
   :class: caution

   Earlier documentation printed a *Gaussian-quartic* dihedral formula
   (a log-sum of exponentials with quartic terms). That described an inherited
   COSMO-style potential — it is **not** what the ``topo`` model uses. The
   implemented torsion is the periodic form above
   (:meth:`~topo.core.system.system.addPeriodicTorsionForce`).


Non-bonded terms (long-range)
-----------------------------

Two pairwise terms act between beads that are **more than two bonds apart**
(see :ref:`exclusions <theory-exclusions>`): screened electrostatics, and the
structure-based contact potential that is the heart of the model.

Electrostatics — Debye–Hückel (Yukawa)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Charged residues interact through a **screened Coulomb (Yukawa) potential**,
which models monovalent-salt screening implicitly:

.. math::

   U^\mathrm{el}_{ij}(r) =
   f\,\frac{q_i\,q_j}{\varepsilon_r\, r}\; e^{-r/l_D}

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Symbol
     - Value
     - Meaning
   * - :math:`f`
     - 138.935458 kJ nm mol\ :sup:`-1` e\ :sup:`-2`
     - Coulomb constant :math:`1/4\pi\varepsilon_0` in MD units
   * - :math:`\varepsilon_r`
     - 78.5
     - Relative dielectric of water (~25 °C)
   * - :math:`l_D`
     - 1.0 nm
     - Debye screening length (≈ 100 mM monovalent salt)
   * - cutoff
     - 2.0 nm
     - Interactions beyond this are ignored
   * - switching
     - 1.8 nm
     - Smooth switch to zero between 1.8 and 2.0 nm

Only four residue types carry charge; everything else is neutral:

* **Negative (−1 e):** ASP, GLU
* **Positive (+1 e):** ARG, LYS
* **Neutral (0 e):** all other residues (including HIS)

So the electrostatic term only matters between acidic/basic residues; it is
weak and short-ranged at physiological salt because of the 1 nm screening length.

.. _theory-contacts:

Structure-based contacts — the heart of the model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every residue pair :math:`(i,j)` that is **not** excluded gets a pairwise
contact potential with a **12-10-6 (Gō-type)** functional form:

.. math::

   U^\mathrm{nb}_{ij}(r) =
   \varepsilon_{ij}\Big[\,
     13\Big(\tfrac{R_{ij}}{r}\Big)^{12}
   - 18\Big(\tfrac{R_{ij}}{r}\Big)^{10}
   +  4\Big(\tfrac{R_{ij}}{r}\Big)^{6}
   \Big]

This well has its **minimum exactly at** :math:`r = R_{ij}`, where
:math:`U^\mathrm{nb}_{ij} = -\varepsilon_{ij}` — so :math:`R_{ij}` is the
preferred distance and :math:`\varepsilon_{ij}` is the **well depth**. The same
cutoff (2.0 nm) and switching (1.8 nm) as the electrostatics apply.

What distinguishes TOPO from a textbook Gō model is **how the well position**
:math:`R_{ij}` **and the well depth** :math:`\varepsilon_{ij}` **are assigned**.
Pairs split into two classes (all built by
:func:`topo.utils.nonbonded.build_nonbonded_interaction`):

**1. Native contacts** — pairs that are genuinely in contact in your input
structure. A pair counts as a native contact if it has at least one of:

* a **backbone hydrogen bond** (from STRIDE — see below),
* a **backbone–sidechain (BS) atomic contact** (any backbone heavy atom of one
  residue within 4.5 Å of a sidechain heavy atom of the other), or
* a **sidechain–sidechain (SS) atomic contact** (sidechain heavy atoms within
  4.5 Å),

*and* the two residues are more than two apart in sequence (``LOCAL_SEPARATION =
2``; this filter is applied per chain, so contacts *between* chains are always
kept). For a native contact:

* :math:`R_{ij}` = the **Cα–Cα distance in the input structure** (the native
  distance), and
* :math:`\varepsilon_{ij}` = the **sum** of three physically distinct
  contributions:

  .. math::

     \varepsilon_{ij}
     = \underbrace{E_\mathrm{HB}}_{\text{H-bonds}}
     + \underbrace{E_\mathrm{BS}}_{\text{backbone–sidechain}}
     + \underbrace{n^{ij}_\mathrm{scale}\; E_\mathrm{SS}}_{\text{scaled sidechain–sidechain}}

  .. list-table::
     :header-rows: 1
     :widths: 20 36 44

     * - Contribution
       - Value
       - Source
     * - :math:`E_\mathrm{HB}`
       - 0.75 kcal/mol per H-bond, **capped at 1.5** (i.e. 0, 0.75, or 1.5)
       - STRIDE backbone H-bond count for the pair
     * - :math:`E_\mathrm{BS}`
       - 0.37 kcal/mol × (number of directional BS contacts: 0, 1, or 2)
       - Heavy-atom distances (≤ 4.5 Å)
     * - :math:`E_\mathrm{SS}`
       - :math:`\lvert\,\mathrm{BT}(t_i,t_j) - 0.6\,\rvert` kcal/mol *(if the pair has an SS contact)*
       - Betancourt–Thirumalai residue-pair potential, by residue **types** :math:`t_i,t_j`
     * - :math:`n^{ij}_\mathrm{scale}`
       - 1.0 by default; set per domain/interface
       - ``domain.yaml`` (see :doc:`domain_definition`)

  All energies are converted to kJ/mol internally (1 kcal/mol = 4.184 kJ/mol).

**2. Non-native pairs** — every other (non-excluded) pair. These get a soft
**excluded-volume repulsion**: a negligible well depth
:math:`\varepsilon_{ij} = 1.32\times10^{-4}` kcal/mol with a distance
:math:`R_{ij} = \tfrac12(\sigma_i + \sigma_j)`, where each residue's
:math:`\sigma_i = 2^{1/6}\times(\text{nearest non-contact Cα distance})`. In
effect, non-native pairs feel almost no attraction but cannot interpenetrate —
they provide chain self-avoidance without biasing toward any non-native fold.

.. admonition:: Why this design matters
   :class: tip

   Splitting the well depth into H-bond + backbone–sidechain + sidechain–sidechain
   parts is what makes the **per-domain scaling** (:doc:`domain_definition`) and
   the **nscale optimizer** (:doc:`../tutorials/05_opt_nscal`) possible: the
   scale factor :math:`n_\mathrm{scale}` multiplies **only the sidechain–sidechain
   part**, leaving the backbone hydrogen-bond and backbone–sidechain energies
   untouched. You can therefore tune the stability of one domain or one interface
   without distorting backbone-driven structure.

The role of STRIDE
~~~~~~~~~~~~~~~~~~~

The hydrogen-bond contribution :math:`E_\mathrm{HB}` requires knowing which
backbone H-bonds exist in the native structure. TOPO obtains these from
`STRIDE <http://webclu.bio.wzw.tum.de/stride/>`_, a standard secondary-structure
and H-bond assignment program:

* If you do **not** set ``stride_output_file``, TOPO runs ``stride -h`` on your
  PDB automatically (STRIDE must be on your ``PATH``) and **caches** the result
  to ``<pdb_prefix>_stride.dat`` next to the structure. Subsequent runs reuse the
  cache; delete it to force regeneration.
* If STRIDE is not installed, precompute the file once
  (``stride -h protein.pdb > stride.dat``) and point ``stride_output_file`` at it.

STRIDE output is parsed for donor/acceptor residue pairs; each physical H-bond
is counted once, and pairs with two or more H-bonds are capped so
:math:`E_\mathrm{HB} \le 1.5` kcal/mol.

.. _theory-exclusions:

Exclusion rule
~~~~~~~~~~~~~~~

Both non-bonded forces (electrostatics and contacts) **skip pairs that are two
or fewer bonds apart** — i.e. **1–2 (bonded) and 1–3 (angle) neighbours** are
excluded, because their geometry is already governed by the bond and angle
terms. This is the ``bonded_exclusions_index = 2`` rule applied via OpenMM's
``createExclusionsFromBonds``. Pairs **1–4 and beyond** *do* feel the non-bonded
terms (subject to the sequence-local contact filter described above).


Per-residue parameters
----------------------

Each Cα bead inherits three properties from its amino-acid type (defined in
:mod:`topo.parameters.model_parameters`): a **mass** (≈ residue molar mass,
amu), a **charge** (e), and a ``radii`` value (nm) reserved for future
inter-chain excluded volume (see the note below the table).

.. list-table::
   :header-rows: 1
   :widths: 14 14 16 12 14 14 16 12

   * - Residue
     - Mass
     - Radii (nm)
     - Charge
     - Residue
     - Mass
     - Radii (nm)
     - Charge
   * - ALA
     - 71.0
     - 0.504
     - 0
     - MET
     - 131.0
     - 0.618
     - 0
   * - ARG
     - 114.0
     - 0.656
     - **+1**
     - PHE
     - 147.0
     - 0.636
     - 0
   * - ASN
     - 114.0
     - 0.568
     - 0
     - PRO
     - 114.0
     - 0.556
     - 0
   * - ASP
     - 114.0
     - 0.558
     - **−1**
     - SER
     - 87.0
     - 0.518
     - 0
   * - CYS
     - 114.0
     - 0.548
     - 0
     - THR
     - 101.0
     - 0.562
     - 0
   * - GLU
     - 128.0
     - 0.592
     - **−1**
     - TRP
     - 186.0
     - 0.678
     - 0
   * - GLN
     - 128.0
     - 0.602
     - 0
     - TYR
     - 163.0
     - 0.646
     - 0
   * - GLY
     - 57.0
     - 0.450
     - 0
     - VAL
     - 99.0
     - 0.586
     - 0
   * - HIS
     - 114.0
     - 0.608
     - 0
     - ILE
     - 113.0
     - 0.618
     - 0
   * - LEU
     - 113.0
     - 0.618
     - 0
     - LYS
     - 128.0
     - 0.636
     - **+1**

Only the **20 standard amino acids** are parameterized; structures containing
non-standard or modified residues (e.g. phosphorylated residues) are not
supported and raise an error at build time. RNA sites P, R, BR are defined for
planned nucleic-acid support.

Of the three properties, the **charge** enters the Yukawa electrostatics and the
**mass** sets the particle dynamics. The **radii are not used by any force in the
current single-chain model**: every contact distance :math:`R_{ij}` — native
*and* non-native — is derived from the input structure (native distances are the
Cα–Cα distances; non-native distances come from the nearest non-contact Cα
distance, see :ref:`theory-contacts`). The ``radii`` values are reserved for a
planned **inter-chain excluded-volume** term (for ribosome–nascent-chain complex
systems), which is not yet implemented.


Temperature, dynamics, and ensembles
------------------------------------

TOPO integrates **Langevin dynamics** (``LangevinIntegrator``): the friction
coefficient ``tau_t`` and reference temperature ``ref_t`` set the thermostat, and
the implicit solvent's viscous drag and random kicks are what make the dynamics
diffusive (as for a protein in water) rather than ballistic.

* **Constant-temperature equilibrium** (default) holds ``ref_t`` for the whole
  run — the standard production protocol.
* **Annealing/quenching** (``anneal = yes``) adds a hot *quench* phase before
  production, to unfold the protein and watch it refold; see
  :doc:`../tutorials/06_anneal`.
* **Pressure coupling** (a Monte-Carlo barostat) and **periodic boundary
  conditions** are available (``pcoupl``/``pbc``) but are rarely needed for a
  single implicit-solvent chain.

Because the native structure is the energy minimum, raising ``ref_t`` toward the
protein's melting temperature breaks native contacts (rising
``Custom Non-Bonded Energy``, falling fraction of native contacts *Q*; see
:doc:`native_contacts`), which is the basis for studying thermal stability.


Calibrating contact nscale
----------------------------

The single most important *adjustable* quantity in the model is
:math:`n_\mathrm{scale}` (the ``nscale`` field in ``domain.yaml``), which
multiplies the sidechain–sidechain well depths. The raw, unscaled model
(:math:`n_\mathrm{scale} = 1`) is usually **under-stabilized** — proteins sit
only marginally folded. Two pages cover how to set it:

* :doc:`domain_definition` — the file format, per-domain and per-interface
  scaling, discontiguous domains, and decoupling.
* :doc:`../tutorials/05_opt_nscal` — the automatic optimizer that searches a
  discrete nscale ladder for the smallest :math:`n_\mathrm{scale}` that keeps
  each domain and interface folded across many trajectories.


Where to go next
----------------

* :doc:`simulation_control` — every ``md.ini`` option that turns these terms on
  and off and sets the run.
* :doc:`domain_definition` — scaling the contact energies per domain.
* :doc:`native_contacts` — measuring how folded the protein is (the *Q* score).
* :doc:`outputs` — the files and the energy/temperature log a run writes.
* :doc:`python_api` — building a model and inspecting its forces in Python.
* :doc:`../tutorials/index` — hands-on, ready-to-run examples.
