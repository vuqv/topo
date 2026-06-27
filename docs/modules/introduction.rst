Introduction 
=========================================================

The TOPO model is a Python library that offers flexibility to set up coarse-grained simulation of globular (folded) proteins using the MD framework of OpenMM toolkit.
The codebase was inherited from COSMO (for IDP simulation) and adapted for globular proteins.
It automates the creation of :code:`openmm.system` classes that contain the necessary force field parameters to run molecular dynamics simulations using a protein structure as the only necessary inputs.

TOPO is divided in three main classes:

1. :code:`system`
2. :code:`models`
3. :code:`geometry`

------------------

:code:`system`, is the main class that holds all the methods to define,
modify and create CG system to be simulated with OpenMM.
Class inheritance from :code:`openmm.system` with some more attributes for TOPO.

------------------

:code:`models`, class contains set of models, each model contains a collection of sequence of commands
to build model, allows to easily set up CG models.

------------------

:code:`geometry`, contains methods to calculate the geometrical parameters from the input structures.
It's not useful in current need of simulation method.

------------------

The library is open-source and offers flexibility to simulate globular proteins.

------------------

Installation
------------

From the repository root::

    pip install -e .

This installs TOPO as an editable package (your source edits take effect
immediately) and registers the ``topo-mdrun`` console command. STRIDE must be on
your ``PATH`` for the contact potential.

Running a simulation
--------------------

Every run is driven by a plain-text control file (``md.ini``) and launched with
any of the following equivalent forms::

    topo-mdrun -f md.ini              # installed console command
    python -m topo.mdrun -f md.ini    # module form
    python run_simulation.py -f md.ini  # the thin shim shipped in each tutorial

See :doc:`../usage/simulation_control` for every ``md.ini`` option, and
:doc:`../tutorials/index` for hands-on walkthroughs.

Beyond a single run
-------------------

* **Many copies in one run** — set ``n_copies`` to pack independent copies of a
  chain into one GPU run (:func:`topo.make_noninteracting_copies`), then split the
  combined trajectory into per-copy DCDs with :func:`topo.split_chains`
  (Tutorial 4).
* **Native-contact analysis (Q)** — :mod:`topo.analysis.native_contacts` computes
  the per-domain and per-interface fraction of native contacts from a trajectory.
* **Strength optimization** — ``optimization.py`` (Tutorial 5) automatically
  searches the per-domain/interface ``strength`` (*n*\ :sub:`scale`) that keeps
  each domain folded, instead of hard-coding it.

