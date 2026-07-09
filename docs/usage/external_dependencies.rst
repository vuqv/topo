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

The convenience script builds/installs both from their upstream source (and
prefers `bioconda <https://bioconda.github.io/>`_ for STRIDE when ``conda`` is
available). From the repository root::

    scripts/install_deps.sh              # both, into $HOME/.local/bin
    scripts/install_deps.sh stride       # just STRIDE
    PREFIX=~/bin scripts/install_deps.sh # choose the install location

.. note::

   Verify the upstream source URLs in ``scripts/install_deps.sh`` before running
   it; they are the canonical distribution points at time of writing but may
   move. Override them with the ``STRIDE_URL`` / ``PULCHRA_URL`` environment
   variables if needed.

Equivalent manual installs:

* **STRIDE via bioconda** — ``conda install -c bioconda stride``
* **STRIDE from source** — download the STRIDE tarball, ``make``, and copy the
  resulting ``stride`` binary onto your ``PATH``.
* **PULCHRA from source** — download ``pulchra_306.tgz`` and compile the single
  C file: ``cc -O3 -o pulchra pulchra*.c -lm``.


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
