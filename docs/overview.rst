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

* :doc:`A.1 · Single-domain quickstart <tutorials/A1_single_domain>` — build and
  run your first simulation of a small single-domain protein, then read its
  outputs.
* :doc:`A.2 · Multidomain & per-domain contact scaling <tutorials/A2_multidomain>`
  — control contact-energy stability **within** each domain and **across**
  interfaces with a ``domain.yaml``.
* :doc:`A.3 · Restart & outputs <tutorials/A3_restart>` — checkpoint and resume a
  run, and understand every file it writes.
* :doc:`A.4 · Many independent copies <tutorials/A4_multicopy>` — replicate one
  protein into many non-interacting copies in a single job.
* :doc:`A.5 · Optimizing the contact scale (nscale) <tutorials/A5_opt_nscal>` —
  automatically calibrate the per-domain and per-interface ``nscale`` so every
  unit stays folded.
* :doc:`A.6 · Anneal & quench <tutorials/A6_anneal>` — drive temperature ramps to
  melt, quench, and observe (un)folding.
* :doc:`A.7 · Disordered / IDR regions <tutorials/A7_idr_mixed>` — mark part of a
  chain intrinsically disordered with a ``disordered:`` section; the folded core
  keeps its shape while the tails and loops stay flexible.

.. rubric:: Reference

* :doc:`The TOPO model: theory & force field <usage/model_theory>`
* :doc:`Simulation control options (md.ini) <usage/simulation_control>`
* :doc:`Domain definition file (domain.yaml) <usage/domain_definition>`
* :doc:`Native-contact (Q) analysis <usage/native_contacts>`
* :doc:`Output files & log format <usage/outputs>`


.. _overview-synthesis:

B. Protein synthesis
--------------------

Grow the nascent chain **N→C, one residue at a time**, so the protein can fold
co-translationally as it emerges from the exit tunnel. Both tutorials layer
codon-resolved kinetics on the Part A model but differ in how the ribosome exit
tunnel is represented — and so in which runner they use: Tutorial B.1 uses
``topo-cylinder`` (an analytic tunnel), Tutorial B.2 uses ``topo-csp`` (an explicit
coarse-grained ribosome).

.. rubric:: Tutorials

* :doc:`B.1 · Synthesis through an analytic tunnel <tutorials/B1_translation_cylinder>`
  — the exit tunnel is a cylindrical bore through a wall (no explicit ribosome
  beads); fast, never jams, and the chain folds co-translationally on egress.
* :doc:`B.2 · Synthesis on a coarse-grained ribosome <tutorials/B2_ribosome_synthesis>`
  — the ribosome-based counterpart: grow the chain through TOPO's own truncated
  CG ribosome, then eject the completed protein (worked on 4c5c and P0CX28).

.. rubric:: Reference

* :doc:`Synthesis in coarse-grained ribosome model <usage/continuous_synthesis>`
* :doc:`Synthesis in cylinder ribosome model <usage/cylinder_synthesis>`
* :doc:`Synthesis control options (csp.ini) <usage/synthesis_control>`
* :doc:`Using TOPO from Python <usage/python_api>`
