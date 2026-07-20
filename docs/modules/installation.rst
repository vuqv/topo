How to install
==============

TOPO is a Python package built on `OpenMM <https://openmm.org/>`_. It depends on
OpenMM plus NumPy, ParmEd, MDAnalysis, mdtraj, pandas, and PyYAML, which are best
installed from ``conda-forge``. It also calls one external command-line program,
**STRIDE**, that is installed separately (plus **PULCHRA**, an opt-in extra only
needed for all-atom backmapping). The steps below target Linux.

Requirements
------------

Python dependencies
~~~~~~~~~~~~~~~~~~~~~

* **Python** >= 3.9
* **OpenMM** >= 7.7 (the MD engine; install from the ``conda-forge`` channel, and
  choose a CUDA toolkit compatible with your NVIDIA driver). ``getStepCount()`` is
  unreliable before 7.7 and is required to restart simulations; OpenMM 8.2 is
  recommended for better performance.
* **ParmEd** >= 3.4
* **NumPy** >= 1.22, **pandas** >= 1.4, **PyYAML** >= 6.0
* **MDAnalysis** >= 2.2, **mdtraj** >= 1.9.7 (trajectory/structure I/O)

These are the exact runtime dependencies of ``import topo``, declared in
``pyproject.toml`` and mirrored in ``requirements.txt``. Floors are the oldest
versions known to work; there are no upper caps (the package runs on current
releases, e.g. NumPy 2.x / OpenMM 8.x). The standalone tools under ``scripts/``
need a few extra packages (``scipy``, ``matplotlib``, ``numba``).

External program (STRIDE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TOPO calls one third-party command-line binary. It is a compiled C program, so it
is **not** installed by ``pip`` and **not** bundled in the wheel — you install it
once and TOPO locates it at run time.

* **STRIDE** (required) — secondary-structure / backbone H-bond assignment for the
  contact map. Only invoked when TOPO has to *build* the contact map; if you supply
  a precomputed STRIDE file (``stride_output_file=...``) it need not be installed
  for that run.

Install it with the bundled helper::

    scripts/install_deps.sh              # STRIDE, into $HOME/.local/bin

A second program, **PULCHRA**, is *optional* and installed only on request
(``scripts/install_deps.sh pulchra``). It backmaps a coarse-grained (Cα) structure
to all-atom coordinates; if you never do that, you do not need it.

TOPO resolves each program in this order: ``$TOPO_STRIDE`` / ``$TOPO_PULCHRA`` (an
explicit path) → the program on ``PATH`` → a copy vendored at ``topo/bin/``. See
:doc:`../usage/external_dependencies` for manual installs, the optional
backmapping tools (PULCHRA and the cg2all alternative), and details.

Steps
-----

1. **Create and activate a fresh conda/mamba environment** with the binary
   dependencies (``mamba`` is recommended for faster, more reliable solves)::

       mamba create -n topo -c conda-forge python">=3.9" openmm parmed \
           mdanalysis mdtraj numpy pandas pyyaml
       mamba activate topo

2. **Get the source code.** TOPO is not on PyPI, so clone (or download) the
   repository from GitHub to a target path, for example ``/path/to/topo``::

       git clone https://github.com/vuqv/topo.git /path/to/topo
       cd /path/to/topo

   Without ``git``, download the ZIP archive from
   `github.com/vuqv/topo <https://github.com/vuqv/topo>`_ ("Code" → "Download
   ZIP") and unpack it; the rest of the steps are the same.

3. **Install STRIDE**, e.g. with the bundled helper, run from the repository
   root (add ``pulchra`` as an argument only if you need backmapping)::

       scripts/install_deps.sh

4. **Install TOPO.** There are two main ways.

Two ways to install
-------------------

**(I) Add to** ``PYTHONPATH`` **(no install)**
    If you only need ``import topo`` and the module-form entry points (no console
    commands), add the repo root (the parent of ``topo/``) to ``PYTHONPATH``
    (persist it in ``.bashrc``)::

        export PYTHONPATH=$PYTHONPATH:/path/to/topo

    Invoke the tools as modules, e.g. ``python -m topo.mdrun -f md.ini``. You must
    still install the dependencies above.

**(II) Install with pip**
    From the repository root (the directory with ``pyproject.toml``). This
    additionally registers the ``topo-mdrun``, ``topo-optimize`` and other console
    commands on your CLI:

    * editable (recommended for development), reflects source edits immediately::

          pip install -e .

    * regular install, copies the package into ``site-packages`` (source edits
      require a reinstall)::

          pip install .

Console commands
----------------

``pip install`` registers these entry points (each also has a module form,
``python -m <module>``):

.. list-table::
   :header-rows: 1
   :widths: 22 26 52

   * - Command
     - Module
     - Purpose
   * - ``topo-mdrun``
     - ``topo.mdrun``
     - Run a folded-protein simulation from an ``md.ini`` control file.
   * - ``topo-optimize``
     - ``topo.optimize``
     - Calibrate per-domain / per-interface contact ``nscale``.
   * - ``topo-csp``
     - ``topo.csp.protocol``
     - Continuous synthesis on an explicit coarse-grained ribosome.
   * - ``topo-cylinder``
     - ``topo.csp.cylinder``
     - Continuous synthesis through an analytic (cylindrical) exit tunnel.
   * - ``topo-csp-movie``
     - ``topo.csp.movie``
     - Stitch per-residue/-stage synthesis trajectories into one VMD movie.
   * - ``topo-make-mrna``
     - ``topo.csp.synth_mrna``
     - Pre-generate a fastest/slowest synonymous-codon mRNA for a protein.

Verify
------

::

    topo-mdrun -h      # prints help if the console command is installed
    topo-csp -h        # prints help
    python -c "import topo; print(topo.__version__)"
