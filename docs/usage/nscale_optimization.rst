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
tutorial :doc:`../tutorials/05_opt_nscal`.


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
  ladder is indexed by structural class; see Table 1 in
  :doc:`../tutorials/05_opt_nscal`);
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
<https://www.pnas.org/doi/full/10.1073/pnas.1813003116>`_). The
:doc:`optimizer <../tutorials/05_opt_nscal>` then starts every domain and
interface at the bottom of its ladder and, round by round, raises only the units
that fail a stability test (all trajectories keeping *Q* above threshold for
nearly all frames), freezing each one as soon as it holds. The result is the
smallest per-unit :math:`n_\mathrm{scale}` that keeps the whole structure folded
across many independent trajectories.


Where to go next
----------------

* :doc:`domain_definition` — the ``domain.yaml`` format: per-domain and
  per-interface scaling, discontiguous domains, and decoupling.
* :doc:`../tutorials/05_opt_nscal` — the hands-on optimizer tutorial (the ladder,
  the per-round search, and the recommended ``ntraj`` / ``min_contacts`` settings).
* :doc:`native_contacts` — the *Q* score that measures how folded each unit is.
* :doc:`model_theory` — the full potential energy function that
  :math:`n_\mathrm{scale}` scales.
