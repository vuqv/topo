.. TOPO documentation master file.

Welcome to the TOPO documentation
=================================

**TOPO** (*TOPOlogy-based coarse-grained model for folded prOteins*) is a Python
library and command-line toolkit for coarse-grained molecular dynamics of
**globular (folded) proteins**, built on `OpenMM <https://openmm.org/>`_. From a
single folded-protein structure it builds a one-bead-per-residue, structure-based
(Gō-like) model and runs Langevin dynamics — ideal for studying folding,
unfolding, thermal/mechanical stability, and multidomain motions.

**New here?** Read :doc:`modules/introduction` for the feature overview, then
:doc:`usage/model_theory` for what the model is, and work through the
:doc:`tutorials/index`.

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   modules/introduction

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 1
   :caption: User guide

   usage/model_theory
   usage/simulation_control
   usage/domain_definition
   usage/native_contacts
   usage/outputs
   usage/python_api

.. toctree::
   :maxdepth: 1
   :caption: API reference

   modules/parameters
   modules/system
   modules/models
   modules/topo.reporter


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
