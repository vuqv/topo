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
   :caption: Background & theory

   usage/model_theory

.. toctree::
   :maxdepth: 1
   :caption: Running simulations

   usage/simulation_control
   usage/domain_definition

.. toctree::
   :maxdepth: 1
   :caption: Analysis & output

   usage/native_contacts
   usage/outputs

.. toctree::
   :maxdepth: 1
   :caption: Python & advanced features

   usage/python_api
   usage/continuous_synthesis

.. toctree::
   :maxdepth: 1
   :caption: API reference

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
