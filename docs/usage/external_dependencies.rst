External dependencies (STRIDE & PULCHRA)
========================================

Beyond its Python dependencies, TOPO calls two **third-party command-line
programs**. They are compiled C binaries, so they are *not* installed by
``pip install topo`` and are *not* bundled in the wheel (STRIDE's redistribution
terms are restrictive and prebuilt binaries are platform-specific). You install
them once and TOPO locates them at run time.

.. list-table::
   :header-rows: 1
   :widths: 16 16 68

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
     - Optional
     - Backmapping a coarse-grained (Cα) structure to an all-atom model. Only
       needed if you reconstruct atomistic coordinates from a TOPO trajectory.

.. note::

   STRIDE is only invoked when TOPO has to *build* the contact map. If you pass a
   precomputed STRIDE file (``stride_output_file=...`` in the API, or ``-s`` on
   the mirror CLI), STRIDE need not be installed for that run.


How TOPO finds them
-------------------

For each program, :func:`topo.utils.external.find_executable` searches in this
order and uses the first hit:

#. **Environment override** — ``$TOPO_STRIDE`` / ``$TOPO_PULCHRA`` set to the
   full path of the executable. Use this to point at a specific build.
#. **On ``PATH``** — the program is found via ``shutil.which`` (the usual case).
#. **Vendored in the package** — ``topo/bin/stride`` or ``topo/bin/pulchra``, if
   you deliberately installed a copy there (off by default).

If none resolve, TOPO raises a ``RuntimeError`` naming the missing program and
how to supply it.


Installing them
---------------

The convenience script installs both binaries. For each tool it tries three
sources in order — build from upstream, build from a GitHub source mirror, then
fall back to the vendored binary — so a run succeeds even when the upstream sites
are down. From the repository root::

    scripts/install_deps.sh              # both, into $HOME/.local/bin
    scripts/install_deps.sh stride       # just STRIDE
    PREFIX=~/bin scripts/install_deps.sh # choose the install location

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

Installing PULCHRA
------------------

PULCHRA is optional — you only need it to backmap a coarse-grained trajectory to
an all-atom model. ``scripts/install_deps.sh`` installs it automatically, walking
the same three sources as STRIDE. To install it independently, use the same
preference order:

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

.. seealso::

   Python-dependency and package installation is covered in
   :doc:`../modules/introduction`.
