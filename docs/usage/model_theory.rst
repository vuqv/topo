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

TOPO is a **unified structure-based (Gō-like) coarse-grained model for globular
and disordered proteins**. *Unified* is the operative word: one Cα force field
covers the whole spectrum below, with **all contact energetics — folded and
disordered alike — tabulated into the same matrices and evaluated by a single**
``CustomNonbondedForce``. No term switches on or off between regimes; a residue's
regime is decided entirely by the *parameter values* it receives:

* **Folded proteins and folded domains** — the classical Gō regime described in
  the rest of this page: native contacts from the input structure, per-domain
  ``nscale`` calibration (:doc:`domain_definition`).
* **Multi-domain proteins with intrinsically disordered regions (IDRs)** —
  ordered and disordered residues coexist in one chain. Residues listed in a
  ``disordered:`` block get no native contacts; instead they interact through a
  transferable, sequence-dependent well (generic cohesion + a scaled
  Betancourt–Thirumalai term) and excluded volume, while the folded domains keep
  their structure-based energetics untouched. Folded–IDR pairs are
  excluded-volume only, so an IDR acts as a genuine flexible linker/tail rather
  than a spurious binder.
* **Fully disordered proteins (IDPs)** — the limiting case where *every* residue
  is disordered. No Gō contacts remain and the model reduces to a transferable
  Cα polymer; the ``eps_gen_kj``/``idr_scale`` defaults are calibrated against
  SAXS :math:`R_g` for a 24-IDP benchmark.

This page describes the shared force field — every term below applies in all
three regimes, and the sections on the :ref:`contact potential
<theory-contacts>` and on :ref:`Rmin_2 <theory-rmin2-transferable>` say exactly
which parameters change when a residue is disordered. The disorder-specific
energy function, the ``domain.yaml`` syntax, and the calibration are documented
in :doc:`disordered_regions`.

Three choices define the model in all three regimes:

* **Coarse-graining — one bead per residue.** TOPO reads an all-atom PDB/CIF
  structure and keeps **only the alpha-carbon (Cα) atom of each residue**
  (:meth:`~topo.core.system.system.getCAlphaOnly`). A 106-residue protein
  becomes 106 beads. Each bead carries the mass and charge of its amino-acid type
  (table below) and an excluded-volume radius that is **derived from the
  structure** for folded residues and taken from the per-type table for
  disordered ones (see :ref:`theory-rmin2-transferable`). Consecutive Cα beads
  *of the same chain* are bonded in sequence.

* **Structure-based — the native fold is the energy minimum.** Unlike a transferable
  force field, the energy function is built *from the input structure itself*:
  the pairs of residues that are in contact in your folded PDB are given
  attractive wells centred at their native distances, and everything else is
  repulsive (residues declared disordered are exempt — they get no native
  contacts at all). The crystal/NMR structure you provide therefore **defines**
  the global energy minimum of the folded part. This is what makes TOPO efficient for studying
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

Non-bonded terms (long-range)
-----------------------------

Two pairwise terms act between beads that are **more than two bonds apart**
(see :ref:`exclusions <theory-exclusions>`): screened electrostatics, and the
structure-based contact potential (:ref:`below <theory-contacts>`).

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

Structure-based contacts
~~~~~~~~~~~~~~~~~~~~~~~~

Every residue pair :math:`(i,j)` that is **not** excluded gets a pairwise
contact potential with a **12-10-6 (Gō-type)** functional form:

.. math::

   U^\mathrm{nb}_{ij}(r) =
   \varepsilon_{ij}\left[\,
     13\left(\frac{R_{ij}}{r}\right)^{12}
   - 18\left(\frac{R_{ij}}{r}\right)^{10}
   +  4\left(\frac{R_{ij}}{r}\right)^{6}
   \,\right]

This well has its **minimum exactly at** :math:`r = R_{ij}`, where
:math:`U^\mathrm{nb}_{ij} = -\varepsilon_{ij}` — so :math:`R_{ij}` is the
preferred distance and :math:`\varepsilon_{ij}` is the **well depth**. The same
cutoff (2.0 nm) and switching (1.8 nm) as the electrostatics apply.

The well position :math:`R_{ij}` and the well depth :math:`\varepsilon_{ij}` are
assigned **per pair class**. Pairs of **ordered** residues split into two classes
(both built by
:func:`topo.utils.nonbonded.build_nonbonded_interaction`); a pair touching a
**disordered** residue forms a third and fourth class, imposed afterwards by
:func:`topo.utils.nonbonded.apply_disorder` (see :doc:`disordered_regions`).

.. important::

   **The same symbol** :math:`R_{ij}` **means two physically different things** —
   a real distance in one class, a collision radius in the others. For a
   **native** contact, :math:`R_{ij}` is the *measured* Cα–Cα **distance of that
   pair** in your input structure, so the attractive well sits exactly at the
   native geometry. For every other class, :math:`R_{ij}` is *not* a distance of
   the pair at all: it is a **sum of two independent per-residue collision
   radii** — how close beads *i* and *j* may approach before they overlap
   (excluded volume). What makes one an attractive **well** and another a soft
   repulsive **wall** is the paired **well depth** :math:`\varepsilon_{ij}`.

   .. list-table::
      :header-rows: 1
      :widths: 20 40 40

      * - Pair class
        - Well position :math:`R_{ij}`
        - Well depth :math:`\varepsilon_{ij}`
      * - **Native** (ordered–ordered, in contact)
        - measured Cα–Cα distance
        - real energy, :math:`E_\mathrm{HB}+E_\mathrm{BS}+n_\mathrm{scale}E_\mathrm{SS}`
      * - **Non-native** (ordered–ordered, not in contact)
        - :math:`(R_\mathrm{min}/2)_i+(R_\mathrm{min}/2)_j`, both Karanicolas–Brooks
        - :math:`1.32\times10^{-4}` kcal/mol (excluded volume only)
      * - **IDR–ordered**
        - :math:`(R_\mathrm{min}/2)_\mathrm{ordered} + r_\mathrm{vdw}^\mathrm{IDR}`
          — the ordered bead keeps its K–B radius, the IDR bead uses the
          transferable :math:`r_\mathrm{vdw}`
        - :math:`1.32\times10^{-4}` kcal/mol (excluded volume only)
      * - **IDR–IDR**
        - :math:`r_\mathrm{vdw}^i + r_\mathrm{vdw}^j` (~11 % tighter than the bare
          :math:`R_\mathrm{min}` sum — O'Brien's generic-bt convention)
        - :math:`\max\big(1.32\times10^{-4}\ \text{kcal/mol},\;
          \varepsilon_\mathrm{gen} + s_\mathrm{IDR}\,\varepsilon_\mathrm{BT}\big)`

   The disordered radius is the **transferable** per-amino-acid value
   :math:`r_\mathrm{vdw} = (R_\mathrm{min}/2)_\text{table}/2^{1/6}`, because a
   Karanicolas–Brooks radius measured from coordinates is meaningless for a
   residue whose coordinates are not its native ones — see
   :ref:`theory-rmin2-transferable`. So a disordered residue always **loses its
   native contact**, and only an IDR–IDR pair can regain attraction, through the
   generic + sequence-dependent term rather than through the structure.

   All four classes are stored in the **same** ``rmin_matrix`` / ``energy_matrix``
   and evaluated by a **single** ``CustomNonbondedForce`` — every distinction
   lives entirely in the *tabulated parameter values*, not in the functional form.
   Merging them into one matrix is an exact computational convenience, not an
   approximation.

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

  .. note::

     **Why** :math:`E_\mathrm{BS}` **is 0, 0.37, or 0.74 — BS contacts are
     directional.** Unlike an SS contact, which is a single symmetric fact about
     the pair (sidechain *i* touches sidechain *j*), a BS contact has a *donor*
     and an *acceptor* role: the backbone belongs to one residue and the
     sidechain to the other. A pair *(i, j)* therefore admits **two independent
     contacts**, and either, both, or neither may be present:

     .. list-table::
        :header-rows: 1
        :widths: 12 44 44

        * - Count
          - Geometry
          - :math:`E_\mathrm{BS}`
        * - **0**
          - Neither backbone reaches the other's sidechain — the pair has no BS
            contact at all (it may still be native via an H-bond or an SS contact).
          - 0
        * - **1**
          - Exactly one direction is satisfied: backbone of *i* is within 4.5 Å of
            a sidechain heavy atom of *j*, **or** backbone of *j* is within 4.5 Å
            of a sidechain heavy atom of *i* — but not both.
          - 0.37 kcal/mol
        * - **2**
          - Both directions are satisfied simultaneously: backbone of *i* touches
            sidechain of *j* **and** backbone of *j* touches sidechain of *i* —
            a mutually interdigitated pair.
          - 0.74 kcal/mol

     Mechanically (:func:`topo.utils.nonbonded.get_bs_contact_matrix`), the two
     directions are found separately and stored in an **asymmetric** matrix whose
     entry :math:`M_{ij}=1` means "backbone of *i* within 4.5 Å of sidechain of
     *j*". The value used by the energy is the symmetrized count
     :math:`M_{ij}+M_{ji} \in \{0,1,2\}`, which is exactly the number of
     directions realized. A count of 2 marks a tighter, more mutually buried
     packing than a count of 1, so it earns twice the well depth. Note that
     glycine, having no sidechain, can only ever be the backbone partner, so any
     BS contact involving it caps at 1.

     This is a genuinely different counting rule from :math:`E_\mathrm{HB}`, which
     is also capped at two units (0, 0.75, 1.5) but for an unrelated reason: there
     the 2 counts *two H-bonds* reported by STRIDE, not two directions.

  All energies are converted to kJ/mol internally (1 kcal/mol = 4.184 kJ/mol).

**2. Non-native pairs** — every other (non-excluded) pair. These get a soft
**excluded-volume repulsion**: a negligible well depth
:math:`\varepsilon_{ij} = 1.32\times10^{-4}` kcal/mol placed at a collision
distance built from a **per-residue collision radius** :math:`(R_\mathrm{min}/2)_i`
combined by the **sum rule**

.. math::

   R_{ij} = \left(\frac{R_\mathrm{min}}{2}\right)_i
          + \left(\frac{R_\mathrm{min}}{2}\right)_j,
   \qquad
   \left(\frac{R_\mathrm{min}}{2}\right)_i
   = \frac{2^{1/6}}{2}\times(\text{nearest non-contact Cα distance}).

Each :math:`(R_\mathrm{min}/2)_i` is residue *i*'s collision **radius** (half the
collision diameter :math:`R_\mathrm{min}`), so the sum rule puts the well minimum
at the sum of the two radii. This is computed by
:func:`topo.utils.nonbonded.calculate_rmin_2_values`, mirroring the per-type
``Rmin_2`` entries in :mod:`topo.parameters.model_parameters`. In effect,
non-native pairs feel almost no attraction but cannot interpenetrate — they
provide chain self-avoidance without biasing toward any non-native fold.

.. admonition:: Naming — ``Rmin/2``, not σ
   :class: note

   TOPO uses the **Rmin/2 (collision-radius) convention** throughout — ``Rmin_2``
   in the parameter tables, :func:`~topo.utils.nonbonded.calculate_rmin_2_values`
   and ``rmin_matrix`` in the builder — to match O'Brien's CHARMM ``.prm``
   NONBONDED blocks. Earlier notes (and textbook Lennard-Jones notation) wrote the
   same quantity with a per-residue **σ** and the arithmetic-mean combining rule
   :math:`R_{ij}=\tfrac12(\sigma_i+\sigma_j)`. **The two are numerically identical**
   — :math:`\sigma_i \equiv R_\mathrm{min} = 2\,(R_\mathrm{min}/2)_i`, and
   :math:`\tfrac12(\sigma_i+\sigma_j) = (R_\mathrm{min}/2)_i + (R_\mathrm{min}/2)_j`
   — only the name changed. Note too that in this **12-10-6** potential the well
   minimum sits at :math:`r = R_{ij}` (i.e. at ``Rmin``), *not* at the potential's
   zero-crossing, so any "σ" here denotes the collision *diameter at the minimum*,
   not the 12-6 σ where :math:`U = 0`.

.. admonition:: Why this design matters
   :class: tip

   Splitting the well depth into H-bond + backbone–sidechain + sidechain–sidechain
   parts is what makes the **per-domain scaling** (:doc:`domain_definition`) and
   the **nscale optimizer** (:doc:`../tutorials/A5_opt_nscal`) possible: the
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
:mod:`topo.parameters.model_parameters`): a **mass**, a **charge** (e), and a
``Rmin_2`` collision radius (nm).

The **mass** is the *residue* molar mass in amu — the free amino acid minus one
water (18.02), because each residue in a chain has lost a water to the peptide
bond. It is a per-type constant and is the only one of the three properties that
every bead uses unconditionally: it sets the particle dynamics.

The ``Rmin_2`` column is a **transferable, per-type** parameter (O'Brien's
per-amino-acid sidechain Rmin/2 values, plus his per-type RNA values below).
"Transferable" is the key point: unlike everything else in a structure-based
model, these numbers depend only on the residue *type*, never on your input
structure — and consequently they are **not** what a folded protein or a nascent
chain uses. See :ref:`theory-rmin2-transferable` immediately after the table.

.. list-table::
   :header-rows: 1
   :widths: 14 14 16 12 14 14 16 12

   * - Residue
     - Mass
     - Rmin_2 (nm)
     - Charge
     - Residue
     - Mass
     - Rmin_2 (nm)
     - Charge
   * - ALA
     - 71.08
     - 0.2862
     - 0
     - MET
     - 131.19
     - 0.3424
     - 0
   * - ARG
     - 156.19
     - 0.3704
     - **+1**
     - PHE
     - 147.18
     - 0.3536
     - 0
   * - ASN
     - 114.10
     - 0.3199
     - 0
     - PRO
     - 97.12
     - 0.3087
     - 0
   * - ASP
     - 115.09
     - 0.3143
     - **−1**
     - SER
     - 87.08
     - 0.2918
     - 0
   * - CYS
     - 103.14
     - 0.3031
     - 0
     - THR
     - 101.10
     - 0.3143
     - 0
   * - GLU
     - 129.12
     - 0.3367
     - **−1**
     - TRP
     - 186.21
     - 0.3816
     - 0
   * - GLN
     - 128.13
     - 0.3424
     - 0
     - TYR
     - 163.18
     - 0.3592
     - 0
   * - GLY
     - 57.05
     - 0.2526
     - 0
     - VAL
     - 99.13
     - 0.3311
     - 0
   * - HIS
     - 137.14
     - 0.3424
     - 0
     - ILE
     - 113.16
     - 0.3424
     - 0
   * - LEU
     - 113.16
     - 0.3424
     - 0
     - LYS
     - 128.17
     - 0.3536
     - **+1**

Only the **20 standard amino acids** are parameterized; structures containing
non-standard or modified residues (e.g. phosphorylated residues) are not
supported and raise an error at build time. RNA sites P, R, BR are defined for
planned nucleic-acid support.

.. _theory-rmin2-transferable:

``Rmin_2`` is transferable — and folded structure overrides it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Of the three properties, the **charge** enters the Yukawa electrostatics and the
**mass** sets the particle dynamics; both apply to every bead. ``Rmin_2`` is
different, and the distinction is the single most common source of confusion:

* ``Rmin_2`` is a **transferable** parameter — one fixed number per residue
  *type*, taken from O'Brien's CHARMM ``.prm`` NONBONDED blocks and independent of
  any structure. It is used where a bead has **no native structure to measure**:

  - **rigid scenery beads** — the ribosome's RNA (P/R/BR) and ribosomal-protein
    beads in continuous synthesis (:func:`topo.csp.ribosome.load_ribosome`);
  - **disordered (IDR) beads** — a residue in a ``disordered:`` region gets the
    transferable radius
    :math:`r_\mathrm{vdw} = (R_\mathrm{min}/2)_{\text{table}}\,/\,2^{1/6}`
    built from its ``Rmin_2`` entry, instead of a structure-derived one (see
    :doc:`disordered_regions`);
  - as the **fallback** per-bead radius for a nascent chain when no
    structure-derived array is supplied ("Option B" in :mod:`topo.csp.ribosome`).

* For a **folded domain or folded protein — the table values are not used.**
  Every radius and every contact distance :math:`R_{ij}` is **derived from your
  input structure**: native wells sit at the measured Cα–Cα distances, and each
  folded residue's collision radius is the per-residue Karanicolas–Brooks value
  :math:`\tfrac12 2^{1/6}\times(\text{nearest non-contact Cα distance})` computed
  by :func:`topo.utils.nonbonded.calculate_rmin_2_values` (see
  :ref:`theory-contacts`). The same residue type therefore carries **different**
  radii in different structures — that is exactly what makes the model
  structure-based rather than transferable.

  A consequence worth stating: in a **multi-domain protein with IDRs**, the two
  conventions coexist inside one chain. Folded beads keep their K-B radii; IDR
  beads take the transferable :math:`r_\mathrm{vdw}`; a folded–IDR cross pair is
  placed at
  :math:`R_{ij} = (R_\mathrm{min}/2)_\mathrm{folded} + r_\mathrm{vdw}^\mathrm{IDR}`.

The transferable values **are** also used by the **inter-chain**
ribosome–nascent-chain excluded-volume term in protein synthesis
(:mod:`topo.csp.ribosome`), which reproduces O'Brien's **12-10-6** form (the
*same* functional form as the intra-chain contacts, **not** a pure
:math:`(\sigma/r)^{12}` repulsion):

.. math::

   U_{ij}(r) = \varepsilon\left[\,
       13\left(\frac{R_{ij}}{r}\right)^{12}
     - 18\left(\frac{R_{ij}}{r}\right)^{10}
     +  4\left(\frac{R_{ij}}{r}\right)^{6}
   \,\right],
   \qquad
   R_{ij} = \left(\frac{R_\mathrm{min}}{2}\right)_i
          + \left(\frac{R_\mathrm{min}}{2}\right)_j ,

with the **sum** combination rule (each bead's collision radius
:math:`R_\mathrm{min}/2` follows the split above: nascent beads use the
per-residue Karanicolas–Brooks radii from the native structure by default, or the
per-amino-acid table values as a fallback; ribosome beads always use the
transferable per-type/per-residue values) and a
very weak well depth :math:`\varepsilon = 1.32\times10^{-4}` kcal/mol, so in
practice the term acts as a **soft excluded volume**. An earlier pure
:math:`(\sigma/r)^{12}` repulsion with the arithmetic-mean combining rule was
~1000× too soft and has been replaced by this O'Brien-consistent form.


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
  :doc:`../tutorials/A6_anneal`.
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
only marginally folded. Three pages cover why and how to set it:

* :doc:`nscale_optimization` — **why** calibration is needed and why it is done
  per domain and per interface (start here).
* :doc:`domain_definition` — the file format, per-domain and per-interface
  scaling, discontiguous domains, and decoupling.
* :doc:`../tutorials/A5_opt_nscal` — the automatic optimizer that searches a
  discrete nscale ladder for the smallest :math:`n_\mathrm{scale}` that keeps
  each domain and interface folded across many trajectories.


Where to go next
----------------

* :doc:`simulation_control` — every ``md.ini`` option that turns these terms on
  and off and sets the run.
* :doc:`synthesis_overview` — how this base force field is **extended for
  protein synthesis** (a growing nascent chain in a rigid ribosome:
  ribosome↔chain excluded volume + electrostatics, the tRNA tether, codon kinetics).
* :doc:`domain_definition` — scaling the contact energies per domain.
* :doc:`native_contacts` — measuring how folded the protein is (the *Q* score).
* :doc:`outputs` — the files and the energy/temperature log a run writes.
* :doc:`python_api` — building a model and inspecting its forces in Python.
* :doc:`../tutorials/index` — hands-on, ready-to-run examples.
