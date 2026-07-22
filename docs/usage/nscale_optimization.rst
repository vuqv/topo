.. _nscale-optimization:

Why optimize the contact nscale?
================================

The single most important *adjustable* quantity in a TOPO model is
:math:`n_\mathrm{scale}` — the ``nscale`` field in ``domain.yaml``. To see what
it controls, recall that each native contact between residues :math:`i` and
:math:`j` is a 12-10-6 well whose **depth** :math:`\varepsilon_{ij}` is the sum
of three physical contributions (:ref:`theory-contacts`):

.. math::

   \varepsilon_{ij}
   = \underbrace{E_\mathrm{HB}}_{\text{hydrogen bonds}}
   + \underbrace{E_\mathrm{BS}}_{\text{backbone–sidechain}}
   + \underbrace{n_\mathrm{scale}\; E_\mathrm{SS}}_{\text{scaled sidechain–sidechain}}\, .

:math:`E_\mathrm{HB}`, :math:`E_\mathrm{BS}`, and :math:`E_\mathrm{SS}` are fixed
per-contact energies read from the input structure; :math:`n_\mathrm{scale}` is
the **one free multiplier**, and it acts on the **sidechain–sidechain (SS) term
only**. Setting :math:`n_\mathrm{scale}` therefore tunes how deep the native
contacts are — i.e. how strongly the fold is held together — while leaving the
backbone hydrogen-bond and backbone–sidechain energies (which encode secondary
structure) untouched.

This page explains **why a raw, unscaled model is rarely usable as-is** and why
the right value is chosen **per domain and per interface** rather than globally.
For the file format see :doc:`domain_definition`; for the options that drive the
automatic search see :doc:`optimization_control`; for a worked run see the
tutorial :doc:`../tutorials/A5_opt_nscal`.


The raw model is under-stabilized
---------------------------------

Coarse-graining a protein to one Cα bead per residue discards the atomic detail
that, in an all-atom force field, sets the exact free-energy balance between the
folded and unfolded states. The SS contact energies that replace it are drawn
from a **transferable, knowledge-based parameterization** — one generic set of
well depths applied to every residue pair in every protein. That set cannot
simultaneously reproduce the folding stability of proteins that differ in size,
topology, and contact density.

In practice the transferable set is **too weak**: at :math:`n_\mathrm{scale} = 1`
(:math:`E_\mathrm{SS}` unscaled) most proteins sit only marginally folded, and at
a realistic simulation temperature (``ref_t`` ≈ 310 K) they partially or fully
unfold within a trajectory. Raising :math:`n_\mathrm{scale}` deepens the SS
contribution and restores the correct stability without re-parameterizing the
whole force field — it is the one knob calibrated per system.


Why per-domain and per-interface, not one global value
------------------------------------------------------

A multidomain protein is not uniformly stable. Its domains differ in

* **secondary-structure content** — an all-:math:`\alpha` domain, an all-\
  :math:`\beta` domain, and a mixed :math:`\alpha/\beta` domain reach the same
  marginal stability at *different* contact strengths (this is why the calibrated
  ladder is indexed by structural class; see :ref:`Table 1
  <nscale-ladder-table>`);
* **size and contact density** — a small, sparsely packed domain needs a deeper
  well to stay folded than a large, densely packed one;
* **interface character** — the contacts *between* two domains form their own
  structural unit, often much weaker than either domain interior, and must be
  scaled separately (``inter_domains`` in :doc:`domain_definition`).

A single global scale therefore cannot win: the value that finally folds the
**weakest** domain **over-stabilizes** the strongest, rigidifying it and
suppressing exactly the native fluctuations you want to observe. Assigning each
domain and each interface its own :math:`n_\mathrm{scale}` is what makes a
multidomain model faithful everywhere at once.


What the wrong value costs you
------------------------------

:math:`n_\mathrm{scale}` sets where each unit sits on the folding landscape, so a
mis-set value biases **every** downstream observable:

* **Too low** — the domain or interface unfolds spuriously during the run. The
  fraction of native contacts *Q* (:doc:`native_contacts`) drifts down, the
  ensemble is contaminated by non-native states, and any thermodynamic or kinetic
  quantity measured from it is meaningless. Under-stabilized chains are also more
  prone to mirror-image and other kinetic-trap artifacts (:doc:`mirror_detection`).
* **Too high** — the fold is over-stabilized. Native fluctuations, hinge motions,
  and conformational sampling are damped, the effective melting temperature is
  pushed unrealistically high, and the model reports a stiffness the real protein
  does not have.

The objective is therefore the **smallest** :math:`n_\mathrm{scale}` at which each
unit stays folded — *just stable enough*. This keeps native dynamics intact while
preventing spurious unfolding, and it is exactly the target the optimizer searches
for.


How TOPO chooses it
-------------------

Rather than tuning by hand, TOPO restricts :math:`n_\mathrm{scale}` to a small,
**pre-calibrated discrete ladder** — five levels per structural class plus a
median fallback, calibrated on a training set of 19 small single-domain proteins
(`Leininger et al., PNAS 116, 5523–5532, 2019
<https://www.pnas.org/doi/full/10.1073/pnas.1813003116>`_). A domain climbs the
ladder for its structural class (:math:`\alpha`, :math:`\beta`, or
:math:`\alpha/\beta`); every domain–domain interface uses the *Interface* ladder.

.. _nscale-ladder-table:

.. list-table:: Table 1. :math:`n_\mathrm{scale}` levels per structural class.
   :header-rows: 1
   :widths: 22 13 13 13 13 13 18

   * - Structural class
     - Level 1
     - Level 2
     - Level 3
     - Level 4
     - Level 5
     - Fallback (median)
   * - :math:`\alpha`
     - 1.1954
     - 1.4704
     - 1.7453
     - 2.0322
     - 2.5044
     - 1.7453
   * - :math:`\beta`
     - 1.4732
     - 1.8120
     - 2.1508
     - 2.5044
     - 2.5044
     - 2.1508
   * - :math:`\alpha/\beta`
     - 1.1556
     - 1.4213
     - 1.6871
     - 1.9644
     - 2.5044
     - 1.6871
   * - Interface
     - 1.2747
     - 1.5679
     - 1.8611
     - 2.1670
     - 2.5044
     - 1.8611

The optimizer starts every domain and interface at **level 1** and repeats one
**round** until every unit is stable:

#. **Build** a model with the current per-unit :math:`n_\mathrm{scale}` values.
#. **Simulate** ``ntraj`` independent trajectories at the reference temperature
   ``ref_t`` (produced as a single multi-copy run).
#. **Score** the fraction of native contacts *Q* (:doc:`native_contacts`) for
   every domain and interface in every trajectory.
#. **Decide.** A unit is *stable* when **all** ``ntraj`` trajectories keep its
   *Q* above ``q_threshold`` (default 0.6688) for at least ``frame_fraction``
   (default 98 %) of the frames. If every unit is stable the search is **done**;
   otherwise each **unstable** unit climbs one level of its ladder while the
   already-stable units stay frozen, and the next round runs.

A unit still unstable at level 5 is pinned at its class **median fallback** (the
level-3 value in Table 1) for the final model. Because domains and their shared
interfaces are coupled through the contact map, raising one unit's
:math:`n_\mathrm{scale}` can destabilize a unit that had already held, so
strongly coupled multidomain systems may need more rounds than the ladder is tall
— raise ``max_rounds`` for them (see :doc:`optimization_control`). The result is
the smallest per-unit :math:`n_\mathrm{scale}` that keeps the whole structure
folded across many independent trajectories.


Where to go next
----------------

* :doc:`domain_definition` — the ``domain.yaml`` format: per-domain and
  per-interface scaling, discontiguous domains, and decoupling.
* :doc:`../tutorials/A5_opt_nscal` — the hands-on optimizer tutorial (the ladder,
  the per-round search, and the recommended ``ntraj`` / ``min_contacts`` settings).
* :doc:`native_contacts` — the *Q* score that measures how folded each unit is.
* :doc:`model_theory` — the full potential energy function that
  :math:`n_\mathrm{scale}` scales.
