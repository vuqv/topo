.. TOPO documentation master file.

Welcome to the TOPO documentation
=================================

**TOPO** (*TOPOlogy-based coarse-grained model for folded prOteins*) is a Python
library and command-line toolkit for coarse-grained molecular dynamics of
**globular (folded) proteins**, built on `OpenMM <https://openmm.org/>`_. From a
single folded-protein structure it builds a one-bead-per-residue, structure-based
(Gō-like) model and runs Langevin dynamics — ideal for studying folding,
unfolding, thermal/mechanical stability, and multidomain motions.

**New here?** Read :doc:`overview` to see the two things TOPO does and jump to
the right tutorials, :doc:`modules/introduction` for installation and package
layout, and :doc:`usage/model_theory` for what the model is.

**Lineage.** TOPO implements the models developed in the **O'Brien lab** (E. P.
O'Brien and co-workers): the coarse-grained structure-based model of O'Brien
*et al.* (2012) and the per-codon co-translational-synthesis (CSP) protocol of
Jiang *et al.* (2023). Its physics and elongation kinetics reproduce that work.
If you publish with TOPO, please cite those references alongside the software —
see :doc:`citation` for the full list with DOIs.

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   overview
   diagrams
   modules/introduction
   usage/external_dependencies
   citation

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 1
   :caption: The TOPO model

   usage/model_theory

.. toctree::
   :maxdepth: 1
   :caption: Isolated-protein simulations

   usage/simulation_control
   usage/domain_definition
   usage/native_contacts
   usage/mirror_detection
   usage/outputs

.. toctree::
   :maxdepth: 1
   :caption: Protein synthesis

   usage/synthesis_overview
   usage/ribosome_preparation
   usage/codon_dwell_times
   usage/continuous_synthesis
   usage/cylinder_synthesis
   usage/synthesis_control

.. toctree::
   :maxdepth: 1
   :caption: Python & API reference

   usage/python_api
   modules/parameters
   modules/system
   modules/models
   modules/topo.reporter

.. toctree::
   :maxdepth: 2
   :caption: Full module index

   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
