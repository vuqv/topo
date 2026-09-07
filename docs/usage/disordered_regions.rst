Disordered regions and intrinsically disordered proteins
========================================================

TOPO can treat part of a chain as an **intrinsically disordered region (IDR)** or
an entire chain as an **intrinsically disordered protein (IDP)**. An IDR is
declared in the optional ``disordered:`` section of the same domain-definition
YAML file supplied through ``domain_def``. No separate file, argument, or INI key
is required.

This page explains:

- how to declare an IDR;
- what changes in the coarse-grained model;
- why TOPO uses an Ashbaugh–Hatch 12–6 interaction for IDR-involving pairs;
- how ``idr_scale`` and ``eps_ev_kj`` affect the ensemble;
- how IDRs affect native-contact analysis, stability optimization, and
  continuous synthesis.

.. note::
   **Scope.** This treatment applies only to TOPO's Cα model. It is implemented by
   :func:`topo.utils.nonbonded.apply_disorder` at the end of the
   nonbonded build. The same ``disordered:`` declaration is used automatically by
   isolated-protein simulations, native-contact (Q) analysis, the ``nscale``
   optimizer, and continuous synthesis (CSP).

Quick start
-----------

Add a ``disordered:`` block to the domain-definition file:

.. code:: yaml

   n_residues: 283

   intra_domains:
     A: {residues: [25-149], nscale: 1.0}
     B: {residues: [150-283], nscale: 1.0}
   inter_domains:
     A-B: 0.5

   disordered:
     residues: [1-24, 150-165]
     idr_scale: 0.10
     eps_ev_kj: 0.8368

The defaults are:

.. code:: yaml

   idr_scale: 0.10
   eps_ev_kj: 0.8368

Use these defaults unless you have experimental information that supports
recalibration for your system.

Declaring a residue disordered has one central consequence:

   Its native Gō contacts are removed and replaced by weak,
   sequence-dependent, non-native interactions, while excluded volume and the
   transferable backbone remain active.

Bonds, angles, dihedrals, and electrostatics are not changed.

Physical picture
----------------

An IDR model must represent two distinct physical properties:

1. **Excluded volume:** every residue has a finite size, so two beads cannot
   occupy the same space.
2. **Residue-dependent attraction:** different residue pairs can have different
   weak tendencies to associate.

A useful analogy is a solid ball with an adjustable sticky surface. The hardness
and size of the ball determine whether particles can overlap; the surface
stickiness determines how strongly they attract after approaching one another.

For IDRs, these properties should be tunable independently. A polar pair may be
only weakly attractive, but its beads must not become smaller or easier to
overlap. Conversely, strengthening hydrophobic attraction should not
simultaneously inflate the repulsive core.

This separation is particularly important for IDPs. Unlike a folded protein,
an IDP has no single native structure stabilized by a fixed contact network. Its
dimensions emerge from a delicate balance among:

- excluded volume;
- weak residue–residue attraction;
- electrostatic attraction and repulsion;
- backbone conformational entropy;
- solvent-mediated effects.

Small changes in that balance can move an ensemble from expanded to compact.

What marking a region as disordered changes
-------------------------------------------

TOPO already uses a transferable, non-Gō local backbone for every residue:

- 3.81 Å bonds;
- a double-well transferable angle potential;
- the Karanicolas transferable dihedral potential.

Folded structure is encoded through native Gō contacts. Marking a region
disordered therefore removes every native contact involving that region and
replaces it with a weak interaction that does not encode a particular fold.

Every nonlocal residue pair belongs to one of three classes:

+---------------+----------------+-------------+-----------------------------+--------------------+
| Pair class    | Nonbonded      | Native      | Well depth                  | Well position      |
|               | potential      | contacts    |                             |                    |
+===============+================+=============+=============================+====================+
| Folded–folded | Gō 12–10–6     | Retained    | Existing native-contact     | Native Cα          |
|               |                |             | energy, or the non-native   | distance, or the   |
|               |                |             | floor                       | Karanicolas–Brooks |
|               |                |             |                             | sum rule           |
+---------------+----------------+-------------+-----------------------------+--------------------+
| IDR–IDR       | Ashbaugh–Hatch | Removed     | ``max(εNN, sIDR εBT(i,j))`` | Sum of the two     |
|               | 12–6           |             |                             | residue radii      |
+---------------+----------------+-------------+-----------------------------+--------------------+
| Folded–IDR    | Ashbaugh–Hatch | Removed     | Same rule as IDR–IDR        | Sum of the folded  |
|               | 12–6           |             |                             | and IDR residue    |
|               |                |             |                             | radii              |
+---------------+----------------+-------------+-----------------------------+--------------------+

The two potentials act on disjoint OpenMM interaction groups:

- the 12–10–6 force evaluates ``{folded} × {folded}``;
- the Ashbaugh–Hatch force evaluates ``{idr} × {idr}`` and
  ``{idr} × {folded}``.

Consequently, each pair is evaluated exactly once. This is important because
OpenMM takes the union of a force’s interaction groups. If the same pair were
admitted to both forces, OpenMM would add both potentials without raising an
error.

.. _idr-ashbaugh-hatch:

How IDR-involving pairs interact: the current AH–LJ form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define the dimensionless Lennard–Jones shape

.. math::


   L(r)=4\left[
   \left(\frac{\sigma}{r}\right)^{12}
   -
   \left(\frac{\sigma}{r}\right)^6
   \right],

with

.. math::


   \sigma=2^{-1/6}R_{ij}.

The minimum of :math:`L(r)` is therefore at :math:`r=R_{ij}`, where
:math:`L(R_{ij})=-1`. TOPO evaluates every IDR–IDR and folded–IDR pair with

.. math::


   U_{ij}(r)=
   \begin{cases}
   \varepsilon_{\mathrm{EV}}L(r)
   +\left(\varepsilon_{\mathrm{EV}}-\varepsilon_{ij}\right),
   & r\le R_{ij},\\[4pt]
   \varepsilon_{ij}L(r),
   & r>R_{ij}.
   \end{cases}

The three quantities have separate roles:

- :math:`R_{ij}` sets the geometric length scale and the minimum-energy separation;
- :math:`\varepsilon_{\mathrm{EV}}` (``eps_ev_kj``) controls the energetic hardness of
  the repulsive core—that is, the energy penalty for bead overlap;
- :math:`\varepsilon_{ij}` controls the residue-pair-specific attractive well depth,
  or bead “stickiness.”

The additive term
:math:`\left(\varepsilon_{\mathrm{EV}}-\varepsilon_{ij}\right)` in the repulsive
branch makes the potential continuous at :math:`R_{ij}`. Because :math:`L(R_{ij})=-1`,

.. math::


   \begin{aligned}
   U_{ij}(R_{ij}^{-})
   &=-\varepsilon_{\mathrm{EV}}
     +\left(\varepsilon_{\mathrm{EV}}-\varepsilon_{ij}\right)
     =-\varepsilon_{ij},\\
   U_{ij}(R_{ij}^{+})
   &=-\varepsilon_{ij}.
   \end{aligned}

Without this additive term, the two branches would meet at
:math:`-\varepsilon_{\mathrm{EV}}` and :math:`-\varepsilon_{ij}` and would generally be
discontinuous. The force is also continuous: :math:`R_{ij}` is the minimum of
:math:`L(r)`, so :math:`\mathrm{d}L/\mathrm{d}r=0` there, and the additive constant has
zero derivative.

This construction **largely decouples bead size and hardness from bead
attraction**. :math:`R_{ij}` sets where the core is located,
:math:`\varepsilon_{\mathrm{EV}}` sets how energetically difficult it is to push two
beads into that core, and :math:`\varepsilon_{ij}` sets how strongly the beads attract
outside the minimum. Thus, residue-pair attraction can be tuned without using
the same parameter to rescale the repulsive wall. The decoupling is not
mathematically perfect if bead size is defined through an energy-dependent
effective-core criterion, but it is the main practical advantage of the AH
split.

Thus, marking a residue as disordered gives every nonlocal pair involving that
residue a finite excluded-volume core and a sequence-dependent attractive well,
without encoding a native contact.

Why Ashbaugh–Hatch LJ is appropriate for IDRs
---------------------------------------------

The AH construction splits the potential at its minimum. The repulsive branch
uses :math:`\varepsilon_{\mathrm{EV}}`, whereas the attractive branch uses the
pair-specific :math:`\varepsilon_{ij}`. Consequently, TOPO can change residue
stickiness without using the same parameter to rescale the repulsive wall.

In conventional Ashbaugh–Hatch notation, the pair stickiness can be written as

.. math::


   \lambda_{ij}=\frac{\varepsilon_{ij}}
   {\varepsilon_{\mathrm{EV}}}.

Two useful limits follow:

- If :math:`\varepsilon_{ij}=0`, the mathematical expression becomes a WCA-like
  repulsive core with zero energy for :math:`r>R_{ij}`.
- If :math:`\varepsilon_{ij}=\varepsilon_{\mathrm{EV}}`, it becomes ordinary LJ
  12–6.

In the actual TOPO parameterization, the non-native energy floor means that
``idr_scale: 0`` leaves a very small residual well rather than producing the exact
:math:`\varepsilon_{ij}=0` mathematical limit. It is therefore best described as an
**approximately self-avoiding or WCA-like reference**, not a strictly pure WCA
chain.

Under this split, varying :math:`\varepsilon_{ij}` from 0 to 2 kJ/mol moves the
:math:`U(r)=k_BT` core only from approximately :math:`0.846R` to :math:`0.819R`. The attraction
changes strongly while the effective bead size changes by only about 3.2%.

Why IDR pairs do not use TOPO’s Gō 12–10–6 potential
----------------------------------------------------

TOPO’s folded-contact potential is

.. math::


   U_{12-10-6}(r)
   =
   \varepsilon\left[
   13\left(\frac{R}{r}\right)^{12}
   -18\left(\frac{R}{r}\right)^{10}
   +4\left(\frac{R}{r}\right)^6
   \right].

It has a minimum of :math:`-\varepsilon` at :math:`r=R`. This form remains appropriate for
TOPO’s folded native-contact model, but it has two undesirable properties when
applied to every nonlocal IDR pair.

1. The repulsive core and attractive well are coupled
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same parameter :math:`\varepsilon` multiplies the entire potential. Increasing
:math:`\varepsilon` therefore changes both:

- the attractive well;
- the energetic repulsive wall.

Thus, the model cannot request “stronger attraction with the same excluded
volume.” A parameter intended to describe residue stickiness would also alter the
distance at which the repulsive energy becomes comparable with :math:`k_BT`. This makes
it difficult to interpret that parameter as a clean solvent-quality control: an
observed change in chain dimensions could arise from stronger attraction, altered
excluded volume, or both.

2. This particular 12–10–6 form has an outer repulsive barrier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The coefficients :math:`(13,-18,+4)` produce more than a contact minimum. They also
produce a positive bump outside the minimum:

- the potential crosses zero near :math:`1.25R`;
- it reaches a local maximum near :math:`1.45R`;
- the maximum is approximately :math:`+0.143\varepsilon`.

A pair approaching from large separation must cross this bump before entering
the contact well.

Such a barrier can be useful as an effective description within a
native-contact model. However, it is not automatically justified as the same
isotropic barrier for every pair involving an IDR residue. IDR residues still
reorganize and displace solvent when they associate; the narrower claim is that
they do not possess the predefined native-contact event for which TOPO’s Gō
barrier was parameterized.

The barrier is also important thermodynamically because it occupies a larger
shell of space than the short-range contact well. It can therefore strongly
affect the net balance between attraction and repulsion, as described below
using the second virial coefficient.

   **The name “12–10–6” does not itself guarantee a barrier.**

   The barrier arises from the particular signs and coefficients used by TOPO.
   Other potentials described as 12–10–6 may have different shapes.

Why ordinary Lennard–Jones 12–6 is not sufficient
-------------------------------------------------

The ordinary Lennard–Jones potential is

.. math::


   U_{\mathrm{LJ}}(r)
   =
   4\varepsilon
   \left[
   \left(\frac{\sigma}{r}\right)^{12}
   -
   \left(\frac{\sigma}{r}\right)^6
   \right].

It has no outer barrier. Its minimum occurs at

.. math::


   r_{\min}=2^{1/6}\sigma.

Replacing the 12–10–6 potential with ordinary LJ would remove the barrier, but
it would not solve the coupling problem. The same :math:`\varepsilon` scales both the
repulsive term, :math:`4\varepsilon(\sigma/r)^{12}`, and the attractive term,
:math:`-4\varepsilon(\sigma/r)^6`. Therefore, making a residue pair more attractive
by increasing :math:`\varepsilon` simultaneously makes close bead overlap more
energetically costly. Although :math:`\sigma` and the zero-crossing distance remain
fixed, the distance at which the repulsive energy reaches a thermal scale such
as :math:`k_BT` changes. Ordinary LJ therefore does not provide an independent
“stickiness” parameter at fixed energetic core hardness.

AH-LJ supplies this missing control by using two energy scales. The inner branch
is governed primarily by :math:`\varepsilon_{\mathrm{EV}}`, whereas the attractive
well depth is :math:`\varepsilon_{ij}`. Consequently, the residue-specific
:math:`\varepsilon_{ij}` values can tune solvent quality and sequence-dependent
attraction without proportionally rescaling the repulsive wall. Conversely,
:math:`\varepsilon_{\mathrm{EV}}` can be chosen to prevent excessive bead overlap
without forcing every residue pair to have the same attraction strength. This
separation is especially useful for IDRs because their sequence-dependent
interactions must vary among residue pairs while their beads should retain a
comparable excluded-volume core.

The separation is best described as **practical or approximate decoupling**:
:math:`R_{ij}` fixes the geometric contact scale, but an energy-defined effective core
can still move slightly when :math:`\varepsilon_{ij}` changes. As shown above, that
movement is small for the parameter range used here.

The progression is therefore:

.. code:: text

   TOPO Gō 12–10–6
       ├── couples the repulsive core to the well depth
       └── contains an outer repulsive barrier
                       │
                       │ remove the barrier
                       ▼
   Standard LJ 12–6
       └── still couples the repulsive core to the well depth
                       │
                       │ separate the two energy scales
                       ▼
   Ashbaugh–Hatch LJ

Thermodynamic interpretation: the second virial coefficient
-----------------------------------------------------------

Well depth alone does not tell us whether a potential is effectively attractive
or repulsive. The width and position of every attractive or repulsive region
also matter. The **second virial coefficient**, :math:`B_2`, summarizes this net
pairwise effect.

At low concentration, the osmotic pressure can be expanded as

.. math::


   \frac{\Pi}{k_BT}=\rho+B_2\rho^2+B_3\rho^3+\cdots,

where :math:`\rho` is the particle number density. An ideal gas has only the first
term. The :math:`B_2\rho^2` term is the leading correction caused by interactions
between pairs.

Begin with the non-interacting reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For completely non-interacting particles,

.. math::


   U(r)=0
   \qquad\text{at every }r.

Their osmotic pressure is exactly the ideal-gas result,

.. math::


   \frac{\Pi}{k_BT}=\rho,

and there is no pair-interaction correction. Therefore,

.. math::


   B_2=0.

This can also be seen directly from the integral below. If :math:`U(r)=0`, then
:math:`e^{-U(r)/(k_BT)}=1`, so the integrand is zero at every distance.

For a spherically symmetric pair potential,

.. math::


   B_2
   =
   2\pi\int_0^\infty
   \left[1-e^{-U(r)/(k_BT)}\right]r^2\,dr.

The sign follows directly from the integrand:

- If :math:`U(r)>0`, the region is repulsive and contributes positively to :math:`B_2`.
- If :math:`U(r)<0`, the region is attractive and contributes negatively to :math:`B_2`.

Therefore:

+----------------+-----------------+-----------------+-------------+-------------------------+
| Microscopic    | :math:`B_2`     | Solvent regime  | Net pair    | Expected polymer        |
| situation      |                 |                 | behavior    | tendency                |
+================+=================+=================+=============+=========================+
| :math:`U(r)=0` | Exactly zero    | Non-interacting | No pair     | Ideal random-walk       |
| everywhere     |                 | ideal reference | interaction | chain, :math:`\nu=1/2`  |
+----------------+-----------------+-----------------+-------------+-------------------------+
| Repulsion      | Positive        | **Good          | Effective   | Expanded, self-avoiding |
| dominates      |                 | solvent**       | excluded    | chain,                  |
|                |                 |                 | volume      | :math:`\nu\approx0.588` |
+----------------+-----------------+-----------------+-------------+-------------------------+
| Repulsion and  | Approximately   | **Theta         | Net pair    | Ideal-like large-scale  |
| attraction     | zero            | solvent**       | interaction | dimensions,             |
| cancel         |                 |                 | cancels     | :math:`\nu\approx1/2`   |
+----------------+-----------------+-----------------+-------------+-------------------------+
| Attraction     | Negative        | **Poor          | Net         | Compact chain;          |
| dominates      |                 | solvent**       | attraction  | approaching             |
|                |                 |                 |             | :math:`\nu\approx1/3`   |
|                |                 |                 |             | in the dense-globule    |
|                |                 |                 |             | limit                   |
+----------------+-----------------+-----------------+-------------+-------------------------+

The non-interacting and theta cases both have :math:`B_2=0`, but they are not
microscopically identical:

- In a **non-interacting system**, :math:`U(r)=0` and every distance contributes zero.
- At the **theta condition**, repulsive regions contribute positively and
  attractive regions contribute negatively, but their total contributions cancel.

The theta chain therefore looks ideal at sufficiently large length scales even
though its residues still interact locally.

The factor :math:`r^2` counts the volume of the spherical shell at distance :math:`r`.
Consequently, a modest repulsive barrier at a relatively large distance may
contribute strongly. Likewise, two potentials with the same minimum depth need
not have the same :math:`B_2`: a broad shallow attractive well can contribute more
than a narrow deep well.

For TOPO’s current AH 12–6 parameterization, the attractive energy in this
calculation is the already defined pair well depth :math:`\varepsilon_{ij}`:

.. math::


   \varepsilon_{ij}
   =
   \max\left(
   \varepsilon_{\mathrm{NN}},
   s_{\mathrm{IDR}}\varepsilon_{\mathrm{BT}}(i,j)
   \right).

It is the depth of the AH minimum,
:math:`U_{ij}(R_{ij})=-\varepsilon_{ij}`; it is not an additional model parameter.

The following values were calculated for actual residue pairs from
``bt_potential.csv``, using ``idr_scale = 0.10``,
:math:`\varepsilon_{\mathrm{EV}}=0.8368` kJ/mol, and :math:`T=300` K. The pair distance is
:math:`R_{ij}=R_{\min/2,i}+R_{\min/2,j}`, and the integral includes TOPO’s switching
region from 1.8 to 2.0 nm. These values describe the AH non-electrostatic term
only; Debye–Hückel electrostatics are not included.

+-----------------+--------------------------+----------------+---------------+-------------------+
| System or pair  | Pair well depth          | :math:`R_{ij}` | :math:`B_2`   | Interpretation    |
|                 | :math:`\varepsilon_{ij}` | (nm)           | (nm³)         |                   |
|                 | (kJ/mol)                 |                |               |                   |
+=================+==========================+================+===============+===================+
| Non-interacting | —                        | —              | 0             | No pair           |
| reference,      |                          |                |               | interaction       |
| :math:`U(r)=0`  |                          |                |               |                   |
+-----------------+--------------------------+----------------+---------------+-------------------+
| **CYS–CYS**     | **0.81170**              | 0.60613        | **−0.01890**  | Slight net        |
|                 |                          |                |               | attraction        |
+-----------------+--------------------------+----------------+---------------+-------------------+
| PHE–MET         | 0.62342                  | 0.69593        | +0.09878      | Repulsion still   |
|                 |                          |                |               | dominates         |
+-----------------+--------------------------+----------------+---------------+-------------------+
| TRP–MET         | 0.64434                  | 0.72399        | +0.09907      | Repulsion still   |
|                 |                          |                |               | dominates         |
+-----------------+--------------------------+----------------+---------------+-------------------+
| LEU–ASP         | 0.00837                  | 0.65664        | +0.38327      | Very weak         |
|                 |                          |                |               | attraction;       |
|                 |                          |                |               | excluded volume   |
|                 |                          |                |               | dominates         |
+-----------------+--------------------------+----------------+---------------+-------------------+
| ARG–LYS         | 0.04184                  | 0.72399        | **+0.49337**  | Largest positive  |
|                 |                          |                |               | non-electrostatic |
|                 |                          |                |               | :math:`B_2`       |
+-----------------+--------------------------+----------------+---------------+-------------------+

CYS–CYS is the strongest pair because its raw BT value is −1.34 kcal/mol:

.. math::


   \varepsilon_{\mathrm{BT}}(\mathrm{CYS,CYS})
   =4.184\left|-1.34-0.6\right|
   =8.11696\ \mathrm{kJ/mol},

and therefore

.. math::


   \varepsilon_{\mathrm{CYS,CYS}}
   =0.10\times8.11696
   =0.811696\ \mathrm{kJ/mol}.

.. note::
   **Pairwise behavior at the current default parameters.**

   With
   :math:`\varepsilon_{ij}=\max(\varepsilon_{\mathrm{NN}},
   0.10\,\varepsilon_{\mathrm{BT}}(i,j))` and
   :math:`\varepsilon_{\mathrm{EV}}=0.8368` kJ/mol, 209 of the 210 unique amino-acid
   pair types have positive second virial coefficients. Their values range from
   :math:`+0.09878` nm³ for PHE–MET to :math:`+0.49337` nm³ for ARG–LYS. CYS–CYS is the only
   negative pair, with :math:`B_2=-0.01890` nm³.

Thus, the default does not make every pair attractive. Most pairs remain
excluded-volume dominated, while the strongest pair is only slightly net
attractive. The overall behavior of an IDP still depends on the frequencies and
sequence arrangement of all pair types, electrostatics, bonded terms, and
collective conformational sampling.

As the AH attraction increases, :math:`B_2` decreases. The interaction therefore moves
in the expected thermodynamic direction: from excluded-volume-dominated behavior,
through the theta condition, and then toward net attraction.
The simulated chain dimensions follow the same trend: :math:`\nu` decreases from 0.637
at ``idr_scale: 0`` to 0.276 at ``idr_scale: 0.30``, passing through experimentally
relevant IDP dimensions before reaching the collapsed regime.

Sequence-dependent well depth
-----------------------------

Every nonlocal pair touching a disordered residue receives

.. math::


   \varepsilon_{ij}^{\mathrm{IDR}}
   =
   \max\left(
   \varepsilon_{\mathrm{NN}},
   s_{\mathrm{IDR}}\varepsilon_{\mathrm{BT}}(i,j)
   \right),

where:

- :math:`s_{\mathrm{IDR}}` is ``idr_scale``;
- :math:`\varepsilon_{\mathrm{NN}}` is the very small non-native floor;
- :math:`\varepsilon_{\mathrm{BT}}(i,j)` is the sidechain–sidechain BT energy for
  residue types :math:`i` and :math:`j`.

The BT energy is

.. math::


   \varepsilon_{\mathrm{BT}}(i,j)
   =
   4.184\left|\operatorname{raw}(i,j)-0.6\right|
   \quad\text{[kJ/mol]}.

This is the value returned by
``topo.utils.nonbonded.get_ss_interaction_energy``: the raw
``bt_potential.csv`` value is shifted by the 0.6 kcal/mol reference, converted to
an absolute magnitude, and converted from kcal/mol to kJ/mol.

The interaction is:

- **nonspecific in coverage**, because it acts on every eligible nonlocal pair
  rather than only native contacts;
- **chemically heterogeneous in depth**, because the BT energy depends on the
  two residue types.

Thus, the model can favor hydrophobic contacts more strongly without assigning
the IDR a predetermined fold.

Folded–IDR cross interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same depth rule is used for IDR–IDR and folded–IDR pairs. A disordered tail
can therefore form transient, sequence-weighted contacts with its folded domain
and may adsorb onto its surface.

This is a modeling choice rather than a result established by the calibration.
The default ``idr_scale = 0.10`` was fitted using fully disordered proteins, for
which every relevant pair was IDR–IDR. Folded–IDR pairs were not represented in
that benchmark. Excessive adsorption of an IDR onto its folded domain should
therefore be treated as a model-sensitivity question.

No ``nscale`` factor is applied
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The IDR well depth uses the unscaled BT energy matrix. It does not inherit a
folded domain’s ``nscale`` because ``nscale`` is a folding-stability parameter and an
IDR has no native fold to stabilize.

At ``idr_scale: 1.0``, the IDR attraction equals the unscaled sidechain–sidechain
BT energy for the same residue-type pair. The mean BT energy is approximately
2.36 kJ/mol, which is well beyond the AH theta point. The calibrated default
``idr_scale: 0.10`` is therefore approximately 10% of an unscaled native
sidechain-contact energy.

Pair distance and radius convention
-----------------------------------

For an IDR-involving pair,

.. math::


   R_{ij}=R_{\min/2,i}+R_{\min/2,j}.

An IDR residue uses the transferable amino-acid-specific :math:`R_{\min}/2` from the
parameter table. A folded residue keeps its structure-derived
Karanicolas–Brooks :math:`R_{\min}/2`, including when it interacts with an IDR
residue.

The radius is an :math:`R_{\min}/2` value, not a :math:`\sigma/2` value. TOPO performs the
conversion inside the AH force:

.. math::


   \sigma_{ij}=2^{-1/6}R_{ij}.

No additional :math:`2^{1/6}` conversion should be applied when populating the
per-residue radius array.

Overriding the per-residue radius, rather than only the pair matrix, also keeps
the nascent-chain–nascent-chain and nascent-chain–ribosome excluded-volume
channels consistent during continuous synthesis.

Configuration reference
-----------------------

All three top-level sections—``intra_domains``, ``inter_domains``, and
``disordered``—are optional. Only ``n_residues`` is required.

+--------------------------+-----------------+-----------------+--------------------+
| Key                      | Required?       | Type and        | Meaning            |
|                          |                 | default         |                    |
+==========================+=================+=================+====================+
| ``disordered``           | No              | Mapping; absent | Enables the IDR    |
|                          |                 |                 | treatment when     |
|                          |                 |                 | present.           |
+--------------------------+-----------------+-----------------+--------------------+
| ``disordered.residues``  | Yes, if         | List            | Residues to treat  |
|                          | ``disordered``  |                 | as disordered.     |
|                          | is present      |                 | Native contacts    |
|                          |                 |                 | involving these    |
|                          |                 |                 | residues are       |
|                          |                 |                 | removed.           |
+--------------------------+-----------------+-----------------+--------------------+
| ``disordered.idr_scale`` | No              | Float; ``0.10`` | Scales the         |
|                          |                 |                 | sequence-dependent |
|                          |                 |                 | BT attraction for  |
|                          |                 |                 | IDR–IDR and        |
|                          |                 |                 | folded–IDR pairs.  |
|                          |                 |                 | Increase to favor  |
|                          |                 |                 | compaction;        |
|                          |                 |                 | decrease to favor  |
|                          |                 |                 | expansion.         |
+--------------------------+-----------------+-----------------+--------------------+
| ``disordered.eps_ev_kj`` | No              | Float;          | Sets the AH        |
|                          |                 | ``0.8368``      | repulsive-core     |
|                          |                 |                 | energy in kJ/mol,  |
|                          |                 |                 | independently of   |
|                          |                 |                 | the pair well      |
|                          |                 |                 | depth.             |
+--------------------------+-----------------+-----------------+--------------------+

Residue numbering is one-based and must match the input PDB. Lists may contain
inclusive ranges, individual integers, or both:

.. code:: yaml

   disordered:
     residues: [1, 2, 5-10, 150-165]

TOPO applies exactly the residues supplied by the user. It does not predict
disorder. MobiDB annotations, experimental information, and low AlphaFold pLDDT
may help identify candidate regions, but the final residue definition remains a
modeling decision.

Low AlphaFold confidence should not be treated as proof of disorder. When used
as a practical screening rule, pLDDT below 70 has been used as an IDR proxy at
proteome scale [1].

Overlap with folded domains: disorder wins
------------------------------------------

A ``disordered:`` range may overlap an ``intra_domains`` range. The disorder
transformation runs after the folded nonbonded build. For every pair touching a
disordered residue, the domain scaling computed earlier in the build is discarded and the
IDR rule is applied.

.. code:: yaml

   n_residues: 100
   intra_domains:
     A: {residues: [1-100], nscale: 1.6871}
   disordered:
     residues: [40-50]

Here, residues 40–50 are disordered and their ``nscale`` has no effect. The
remaining portions of domain A still belong to one domain. The reader prints an
informational message listing overlapping residues so accidental overlap remains
visible.

Choosing and tuning ``idr_scale``
---------------------------------

The calibrated default is

.. code:: yaml

   idr_scale: 0.10

Use it unless system-specific experimental data support another value.

- Increasing ``idr_scale`` strengthens the sequence-dependent attractive well and
  generally lowers :math:`\nu`, favoring compaction.
- Decreasing ``idr_scale`` favors expansion.
- ``idr_scale: 0`` removes the BT-scaled attraction but retains the tiny
  non-native floor and the excluded-volume core. It is an approximately
  self-avoiding reference, not a noninteracting ghost chain.
- ``eps_ev_kj`` is primarily a repulsive-core parameter, not the preferred solvent-
  quality dial.

The benchmark theta point occurs near ``idr_scale ≈ 0.32``. Values much above 0.3
therefore enter the collapsed regime for that benchmark and should be justified
for the system being studied.

A self-avoiding reference may be appropriate for:

- an entropic linker whose reach is more important than its compaction;
- a strongly charged, highly expanded sequence;
- a deliberately minimal reference ensemble.

If SAXS, smFRET, or other ensemble measurements are available for the system of
interest, treat ``idr_scale`` as a parameter to recalibrate rather than assuming
that the global default is optimal.

Electrostatics must be interpreted separately. TOPO’s Debye–Hückel term can
favor expansion or compaction depending on charge composition and sequence
patterning; it is not universally repulsive.

.. _idr-validation:

Validation against SAXS radii of gyration
-----------------------------------------

The default was evaluated on 24 proteins of 24–273 residues with published SAXS
radii of gyration. Eighteen IDPs were used for calibration, and six foldable
globular proteins were reported separately as controls.

The excluded control proteins were CspTm, R15, R17, hCyp, Protein-L, and sNase.
Their published :math:`R_g` values describe folded states, so comparing those values
with fully disordered simulations does not test the IDR model. This
classification was made from protein identity rather than a label stored in the
benchmark dataset and should remain visible as a judgment call.

Each protein was simulated as a fully disordered chain for 90 ns of Langevin
dynamics at 300 K, starting from an expanded coil. The first 15 ns were
discarded. The reported value was the mass-weighted ensemble average

.. math::


   \sqrt{\left\langle R_g^2\right\rangle}

over the remaining 75 ns.

.. figure:: img/idr_validation.png
   :width: 100%
   :alt: TOPO and HPS-Urry radii of gyration versus SAXS measurements for 18 intrinsically disordered proteins

   **Left:** TOPO's Cα IDR model at the calibrated defaults
   ``idr_scale = 0.10`` and ``eps_ev_kj = 0.8368`` kJ/mol. **Right:** the
   HPS-Urry force field evaluated for the same 18 proteins as an external
   reference. The dashed line is :math:`y=x`, the green line is the
   ordinary-least-squares fit, and point color represents fractional deviation.
   TOPO follows :math:`y=x` more closely across the measured range; its fitted
   slope is 0.81, compared with 0.68 for HPS-Urry.

For the 18-IDP calibration set:

+------------+-------------+-------------+------------+-----------+-----------+
| Model      | :math:`\nu` | :math:`R_0` | RMS        | Pearson   | OLS slope |
|            |             |             | fractional | :math:`r` |           |
|            |             |             | error      |           |           |
+============+=============+=============+============+===========+===========+
| **TOPO Cα  | **0.566**   | **0.223**   | **12.0%**  | **0.89**  | **0.81**  |
| IDR**      |             |             |            |           |           |
+------------+-------------+-------------+------------+-----------+-----------+
| HPS-Urry   | 0.490       | 0.301       | 19.7%      | 0.70      | 0.68      |
| reference  |             |             |            |           |           |
+------------+-------------+-------------+------------+-----------+-----------+
| Experiment | 0.551       | 0.244       | —          | —         | —         |
+------------+-------------+-------------+------------+-----------+-----------+

The scaling relation was

.. math::


   R_g=R_0N^\nu.

The RMS value is the root-mean-square fractional deviation:

.. math::


   \sqrt{
   \frac{1}{N}
   \sum_i
   \left(
   \frac{R_{g,i}^{\mathrm{sim}}-R_{g,i}^{\mathrm{exp}}}
   {R_{g,i}^{\mathrm{exp}}}
   \right)^2
   }.

At ``idr_scale: 0.10``, TOPO gave :math:`\nu=0.566`. A three-seed check at
``idr_scale: 0.12`` gave :math:`\nu=0.559\pm0.009`.

Important limitation
~~~~~~~~~~~~~~~~~~~~

The 12.0% RMS deviation does not demonstrate that the model captures detailed
sequence-specific differences. A power law fitted directly to the same
experimental data already leaves a 9.5% residual. In addition, three chains of
length 185 span 36% in experimental :math:`R_g` but only 7% in the model, in the wrong
rank order.

The model is therefore supported for approximate ensemble dimensions across
different chain lengths. It should not currently be used to confidently rank
the dimensions of two IDPs of similar length. The remaining discrepancy may
involve bonded terms, charge patterning, or other omitted interactions rather
than only the contact potential.

Effect on native-contact analysis and stability optimization
------------------------------------------------------------

The IDR mask is applied consistently to the energy function and native-contact
definitions.

Native-contact Q analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

Any native contact touching a disordered residue is removed from:

- ``Q_protein``;
- every ``Q_domain``;
- interface-Q calculations.

Otherwise, contacts that cannot form would remain permanently in the
denominator and artificially lower Q.

Effective domain membership is therefore

.. math::


   \text{effective domain}=\text{declared domain}-\text{disordered residues}.

The ``nscale`` optimizer
~~~~~~~~~~~~~~~~~~~~~~~~

The optimizer evaluates folded domains and their interfaces. The IDR remains
present and physically active in every simulation round, but it is not a scoring
unit, does not enter the convergence check, and is never assigned an ``nscale``.

Always run the optimizer using the domain-definition file that already contains
the ``disordered:`` section. Declaring an IDR removes cross-boundary native
contacts and can therefore alter the stability of the remaining folded domain.
Optimizing the fully folded model and adding the IDR afterward would optimize a
different energy function.

Continuous synthesis
--------------------

CSP uses the same ``disordered:`` section without additional configuration. At
each nascent-chain length, the full disorder mask is restricted to the residues
that have emerged.

- Before any folded residue emerges, the folded 12–10–6 interaction group is
  empty.
- After folded residues emerge, the AH force evaluates IDR–IDR and folded–IDR
  pairs, while the Gō force evaluates folded–folded pairs.
- The per-residue radius array is shared with the nascent-chain–ribosome
  excluded-volume construction.

The separate nascent-chain–ribosome interaction remains the existing 12–10–6
form at ``RIBO_NC_EPS_KJ`` for every nascent bead, including disordered beads. This
preserves the O’Brien ribosome-interaction parameterization. It also means that
an IDR bead uses different effective interactions toward the ribosome and toward
the protein chain; users should keep this distinction in mind when interpreting
nascent-chain behavior.

Fully disordered proteins
-------------------------

For a fully disordered protein, list every residue and omit ``intra_domains``:

.. code:: yaml

   n_residues: 92
   disordered:
     residues: [1-92]

This is the configuration used to calibrate the defaults. All native contacts
are removed, and all eligible nonlocal pairs use the AH interaction.

The Q analysis returns an empty contact list and Q is ``NaN``. The ``nscale``
optimizer detects that no foldable native contacts remain and exits with a
“nothing to optimize” message.

Do not set both ``idr_scale`` and ``eps_ev_kj`` to zero as a way to create a
self-avoiding chain. Setting ``eps_ev_kj: 0`` removes the intended AH repulsive
core. To obtain the approximately self-avoiding reference while retaining
physical bead size, use:

.. code:: yaml

   disordered:
     residues: [1-92]
     idr_scale: 0
     eps_ev_kj: 0.8368

Starting structures and equilibration
-------------------------------------

The starting coordinates do not determine the equilibrium ensemble, provided
the simulation samples the equilibrium distribution adequately. An IDR that
starts from a compact or folded-looking structure can relax after its native
contacts are removed, so an extended-chain PDB is not strictly required.

Starting coordinates can nevertheless affect equilibration time. Discard the
initial relaxation period, check convergence, and consider simulations from
both compact and expanded starting structures when slow collapse, adsorption,
or barrier crossing is plausible.

Common pitfalls
---------------

- ``idr_scale: 0`` **does not remove excluded volume.** It retains the core set by
  ``eps_ev_kj`` and also retains TOPO’s minute non-native attraction floor.
- ``eps_ev_kj: 0`` **removes the intended repulsive core.** Do not use it to define
  a self-avoiding reference.
- **Numbering must match the input structure.** Residue numbering is one-based;
  an incorrect number disorders the wrong residue.
- **Domain overlap is legal.** If a residue appears in both a domain and the
  disorder mask, disorder wins and TOPO prints an informational message.
- **The current disorder mask is system-wide.** It is a flat residue set and has
  no per-chain qualifier.
- **Declaring an IDR weakens the original Gō network.** Every native contact with
  at least one disordered endpoint is deleted. Re-optimize folded-domain
  stability with the IDR declaration already present.
- **Use spaces in YAML.** Do not use tab indentation, and include a space after
  each colon.

Reference
---------

1. Tesei G. *et al.* Conformational ensembles of the human intrinsically
   disordered proteome. *Nature*. 2024;626:897–904.
   https://doi.org/10.1038/s41586-023-07004-5
