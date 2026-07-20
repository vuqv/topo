External dependencies
=====================

Beyond its Python dependencies, TOPO works alongside three **third-party
command-line programs**. Only **STRIDE** is part of the normal install; the other
two are backmapping tools you install only if you reconstruct all-atom
coordinates from a coarse-grained trajectory — pick **either** PULCHRA **or**
cg2all, you do not need both.

.. list-table::
   :header-rows: 1
   :widths: 14 18 68

   * - Program
     - Needed
     - Used for
   * - **STRIDE**
     - Required
     - Secondary-structure and backbone hydrogen-bond assignment. TOPO runs
       ``stride -h`` on your structure to build the native-contact / H-bond list
       (:doc:`native_contacts`) and to find helix/strand segments for
       mirror-image detection (:doc:`mirror_detection`).
   * - **PULCHRA**
     - Optional (opt-in)
     - Backmapping a coarse-grained (Cα) structure to an all-atom model. Only
       needed if you reconstruct atomistic coordinates from a TOPO trajectory;
       nothing in a standard simulation, optimization or synthesis run uses it.
   * - **cg2all**
     - Optional alternative
     - The same backmapping job done by a neural network instead of PULCHRA's
       rule-based rebuild, with native trajectory and CUDA support. Untested
       against TOPO output — see :ref:`cg2all-backmapping` below.

STRIDE and PULCHRA are compiled C binaries, so they are *not* installed by
``pip install topo`` and are *not* bundled in the wheel (STRIDE's redistribution
terms are restrictive and prebuilt binaries are platform-specific); you install
them once and TOPO locates them at run time. cg2all is a Python package with its
own heavy dependencies, installed in a **separate environment** and run by hand —
TOPO does not call it.

.. note::

   STRIDE is only invoked when TOPO has to *build* the contact map. If you pass a
   precomputed STRIDE file (``stride_output_file=...`` in the API, or ``-s`` on
   the mirror CLI), STRIDE need not be installed for that run.


How TOPO finds them
-------------------

For STRIDE and PULCHRA, :func:`topo.utils.external.find_executable` searches in
this order and uses the first hit:

#. **Environment override** — ``$TOPO_STRIDE`` / ``$TOPO_PULCHRA`` set to the
   full path of the executable. Use this to point at a specific build.
#. **On ``PATH``** — the program is found via ``shutil.which`` (the usual case).
#. **Vendored in the package** — ``topo/bin/stride`` or ``topo/bin/pulchra``, if
   you deliberately installed a copy there (off by default).

If none resolve, TOPO raises a ``RuntimeError`` naming the missing program and
how to supply it.


Installing STRIDE
-----------------

The convenience script installs STRIDE by default; PULCHRA is installed only if
you ask for it by name. For each tool it tries three sources in order — build
from upstream, build from a GitHub source mirror, then fall back to the vendored
binary — so a run succeeds even when the upstream sites are down. From the
repository root::

    scripts/install_deps.sh                 # STRIDE only, into $HOME/.local/bin
    scripts/install_deps.sh pulchra         # PULCHRA only (opt-in)
    scripts/install_deps.sh stride pulchra  # both
    PREFIX=~/bin scripts/install_deps.sh    # choose the install location

The STRIDE preference order (both in the script and for manual installs) is:

#. **Build from the upstream source** — the canonical STRIDE tarball
   (``$STRIDE_URL``, default ``https://webclu.bio.wzw.tum.de/stride/stride.tar.gz``):
   ``make`` and copy the resulting ``stride`` onto your ``PATH``. This site is
   often unreachable; if so, use the mirror below.
#. **Compile from the GitHub mirror** — the maintained
   `MDAnalysis/stride <https://github.com/MDAnalysis/stride>`_ source
   (``$STRIDE_GIT``): ``git clone`` it, ``make``, and copy ``stride`` onto your
   ``PATH``.
#. **Use the shipped binary** — if neither source build works, the validated
   STRIDE binary is vendored in the repository at ``assets/stride/`` (see its
   ``README.md`` for provenance). Put that directory on your ``PATH``, or set
   ``TOPO_STRIDE`` to ``assets/stride/stride``. This is a Linux x86-64 build.

   .. note::

      STRIDE is distributed for **academic use with restrictive redistribution
      terms**. The vendored binary is included for academic reproducibility only;
      it is kept out of the installable wheel. If you redistribute TOPO, review
      STRIDE's terms, and cite Frishman & Argos (*Proteins* **23**, 566–579,
      1995) in any work using its assignments.

Installing PULCHRA (optional)
-----------------------------

Skip this section unless you backmap coarse-grained trajectories to all-atom
models — nothing else in TOPO calls PULCHRA. When you do need it,
``scripts/install_deps.sh pulchra`` installs it, walking the same three sources
as STRIDE. To install it by hand, use the same preference order:

#. **Build from the upstream source** — ``pulchra_306.tgz``
   (``$PULCHRA_URL``, default ``http://www.pirx.com/downloads/pulchra_306.tgz``).
   Download, extract, and compile (PULCHRA builds from two C files and needs its
   ``*.h`` data files in the same directory, so compile from inside the extracted
   tree)::

       curl -fSLO http://www.pirx.com/downloads/pulchra_306.tgz
       tar -xzf pulchra_306.tgz
       cd pulchra_306 && cc -O3 pulchra.c pulchra_data.c -lm -o pulchra

   This site is often unreachable; if so, use the mirror below.

#. **Compile from the GitHub mirror** — the
   `euplotes/pulchra <https://github.com/euplotes/pulchra>`_ source
   (``$PULCHRA_GIT``)::

       git clone https://github.com/euplotes/pulchra.git
       cd pulchra && cc -O3 pulchra.c pulchra_data.c -lm -o pulchra

#. **Use the shipped binary** — if neither source build works, the PULCHRA binary
   is vendored in the repository at ``assets/pulchra/`` (MIT-licensed, so it is
   shipped freely; see its ``README.md``). Put that directory on your ``PATH``, or
   set ``TOPO_PULCHRA`` to ``assets/pulchra/pulchra``. This is a Linux x86-64
   build.

Then make TOPO see the binary: put it on your ``PATH`` or set ``TOPO_PULCHRA``
(see :ref:`below <making-topo-see-them>`).

.. note::

   Both **vendored binaries** (``assets/stride/`` and ``assets/pulchra/``) are
   **Linux x86-64** builds. On any other OS or architecture you must build the
   tool from source and point TOPO at your own binary.


.. _making-topo-see-them:

Making TOPO see them
--------------------

Put the install location on your ``PATH``::

    export PATH="$HOME/.local/bin:$PATH"

or point TOPO directly at each binary::

    export TOPO_STRIDE="$HOME/.local/bin/stride"
    export TOPO_PULCHRA="$HOME/.local/bin/pulchra"

Confirm resolution from Python::

    >>> from topo.utils.external import find_executable
    >>> find_executable("stride")
    '/home/you/.local/bin/stride'

.. _cg2all-backmapping:

cg2all — a deep-learning alternative to PULCHRA (optional)
-----------------------------------------------------------

`cg2all <https://github.com/huhlim/cg2all>`_ (Heo & Feig) reconstructs all-atom
coordinates from a coarse-grained model with an equivariant graph neural network
instead of PULCHRA's rule-based/statistical rebuild. It is a plausible drop-in
replacement for the backmapping step. Features relevant to TOPO output
(per its documentation):

* **Reads a Cα trace directly** — ``--cg CalphaBasedModel``, exactly what a TOPO
  trajectory is. Other CG representations (residue/sidechain centre of mass,
  backbone-only, Martini, Martini3, PRIMO) are supported by the same command.
* **Backmaps a whole trajectory in one call** — pass the CG topology with ``-p``
  and a DCD with ``-d``; the output is a DCD (``-opdb`` additionally writes the
  last frame as a PDB). ``--batch`` sets frames per forward pass and ``--proc``
  the input-loading threads. This replaces the per-frame loop in
  ``scripts/pulchra_backmap_no_optimize.py``.
* **Runs on CUDA** — ``--device cuda`` (or ``cpu``); with no flag it picks one
  automatically. Note the authors' own caveat: **CPU is usually the faster
  choice**, because loading the checkpoint onto the GPU dominates unless the
  system is large or the DCD has many frames. Benchmark before assuming the GPU
  wins.
* **Can pin the input coordinates** — ``--fix`` keeps the input Cα positions
  exactly, so the backmapped model stays on the TOPO trace instead of drifting
  off it.
* Output atom names are **CHARMM** by default; add ``--standard-name`` for IUPAC
  names if a downstream tool expects them.

.. warning::

   cg2all is **not integrated into TOPO and has not been validated against TOPO
   trajectories by us** — it is recorded here as a promising option. PULCHRA
   remains the path the ``scripts/`` helpers use. Check the reconstruction
   yourself (bond geometry, clashes, chirality) before trusting it downstream.

Install it into its own environment — it pulls in PyTorch/DGL, which you do not
want in the TOPO environment::

    # CPU only
    pip install git+https://github.com/huhlim/cg2all

    # GPU
    conda create -n cg2all pip cudatoolkit=11.3 dgl=1.0 -c dglteam/label/cu113
    conda activate cg2all
    pip install git+https://github.com/huhlim/cg2all

Backmap a TOPO structure or trajectory with its ``convert_cg2all`` command::

    # single structure
    convert_cg2all -p cg_model.pdb -o all_atom.pdb --cg CalphaBasedModel

    # trajectory (CG topology + DCD in, all-atom DCD out), on the GPU
    convert_cg2all -p cg_model.pdb -d traj.dcd -o all_atom.dcd \
        --cg CalphaBasedModel --device cuda --batch 8 --fix

cg2all is Apache-2.0 licensed. If you use it, cite L. Heo & M. Feig, *One bead
per residue can describe all-atom protein structures*,
`Structure 32, 97–111.e6 (2024) <https://doi.org/10.1016/j.str.2023.10.013>`_.


.. seealso::

   Python-dependency and package installation is covered in
   :doc:`../modules/introduction`.
