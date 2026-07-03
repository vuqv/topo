What TOPO does
==============

TOPO builds a single **one-bead-per-residue, structure-based (Gō-like) model**
from a folded-protein structure and runs Langevin dynamics on
`OpenMM <https://openmm.org/>`_. That one model powers **two complementary
workflows**, each with its own set of tutorials:

* :ref:`A. Coarse-grained simulation of folded proteins <overview-simulation>`
  — take a fully-formed structure and study how it moves, unfolds, and comes
  apart.
* :ref:`B. Protein synthesis <overview-synthesis>` — grow the chain residue by
  residue on the ribosome and watch it fold *as it is made*.

Part B builds directly on the Part A model, so start with A if you are new here.


.. _overview-simulation:

A. Coarse-grained simulation of folded proteins
-----------------------------------------------

Start from a complete structure and run structure-based MD: folding and
unfolding, thermal and mechanical stability, and multidomain motions. Contact
energies can be scaled per domain and per interface, so different parts of a
protein can be made more or less stable.

.. rubric:: Tutorials

* :doc:`1 · Single-domain quickstart <tutorials/01_single_domain>` — build and
  run your first simulation of a small single-domain protein, then read its
  outputs.
* :doc:`2 · Multidomain & per-domain contact scaling <tutorials/02_multidomain>`
  — control contact-energy stability **within** each domain and **across**
  interfaces with a ``domain.yaml``.
* :doc:`3 · Restart & outputs <tutorials/03_restart>` — checkpoint and resume a
  run, and understand every file it writes.
* :doc:`4 · Many independent copies <tutorials/04_multicopy>` — replicate one
  protein into many non-interacting copies in a single job.
* :doc:`5 · Optimizing the contact scale (nscale) <tutorials/05_opt_nscal>` —
  automatically calibrate the per-domain and per-interface ``nscale`` so every
  unit stays folded.
* :doc:`6 · Anneal & quench <tutorials/06_anneal>` — drive temperature ramps to
  melt, quench, and observe (un)folding.

.. rubric:: Reference

* :doc:`The TOPO model: theory & force field <usage/model_theory>`
* :doc:`Simulation control options (md.ini) <usage/simulation_control>`
* :doc:`Domain definition file (domain.yaml) <usage/domain_definition>`
* :doc:`Native-contact (Q) analysis <usage/native_contacts>`
* :doc:`Output files & log format <usage/outputs>`


.. _overview-synthesis:

B. Protein synthesis (co-translational)
---------------------------------------

Grow the nascent chain **N→C, one residue at a time**, so the protein can fold
co-translationally as it emerges from the exit tunnel. Both tutorials use
codon-resolved kinetics with ``topo-csp`` on top of the Part A model; they differ
in how the ribosome exit tunnel is represented.

.. rubric:: Tutorials

* :doc:`7 · Synthesis through an analytic tunnel <tutorials/07_translation_cylinder>`
  — the exit tunnel is a cylindrical bore through a wall (no explicit ribosome
  beads); fast, never jams, and the chain folds co-translationally on egress.
* :doc:`8 · Synthesis on a coarse-grained ribosome <tutorials/08_ribosome_synthesis>`
  — the ribosome-based counterpart: grow the chain through TOPO's own truncated
  CG ribosome, then eject the completed protein (worked on 4c5c and P0CX28).

.. rubric:: Reference

* :doc:`Continuous synthesis protocol — CSP (topo-csp) <usage/continuous_synthesis>`
* :doc:`Using TOPO from Python <usage/python_api>`
