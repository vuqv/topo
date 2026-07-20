Introduction
============

**TOPO** (*TOPOlogy-based coarse-grained model for folded prOteins*) is a Python
library and command-line toolkit for **coarse-grained molecular dynamics of
globular (folded) proteins**, built on the `OpenMM <https://openmm.org/>`_ engine.

Given only a folded-protein structure (a PDB or CIF file), TOPO automatically
builds a **one-bead-per-residue (alpha-carbon) structure-based model** — bonds,
angles, sequence-dependent torsions, screened electrostatics, and a Gō-like
native-contact potential — and runs Langevin dynamics. Because the native state
is the energy minimum by construction, TOPO is well suited to **protein folding
and unfolding, thermal/mechanical stability, and multidomain motions**.

If you are new here, read :doc:`../usage/model_theory` for what the model *is*,
then work through the :doc:`../tutorials/index`.

.. note::

   TOPO simulates **folded** proteins. For intrinsically disordered proteins,
   the model that TOPO's code base was originally adapted from (COSMO) is the
   more appropriate tool.


.. seealso::

   For a feature-by-feature map of what TOPO does and which tutorial teaches
   each, see :doc:`../overview`.


The model in one paragraph
--------------------------

TOPO keeps only the alpha-carbon of each residue and builds a potential with
**bonded** terms (rigid or harmonic Cα–Cα bonds, a bimodal Gaussian backbone
angle, and four-periodicity sequence-dependent torsions) and **non-bonded**
terms (Debye–Hückel screened electrostatics between charged residues, and a
12-10-6 structure-based contact potential). Native contacts — pairs in contact
in your input structure (via backbone hydrogen bonds detected by STRIDE, and
backbone–sidechain / sidechain–sidechain heavy-atom proximity) — get attractive
wells at their native Cα–Cα distances; all other pairs get a soft excluded-volume
repulsion. The full functional forms, constants, and parameter sources are in
:doc:`../usage/model_theory`.


Package layout
--------------

The codebase is organized into focused subpackages:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Module
     - Responsibility
   * - :mod:`topo.core`
     - The model: ``system`` (builds forces), ``models`` (the ``buildCoarseGrainModel`` entry point), ``geometry`` (distance helpers).
   * - :mod:`topo.parameters`
     - Force-field constants: per-residue mass/radii/charge, bond/angle constants, and the sequence-dependent dihedral table.
   * - :mod:`topo.mdrun`
     - The simulation runner (``topo-mdrun``): build → set up → run the temperature protocol → finalize.
   * - :mod:`topo.optimize`
     - The contact-nscale optimizer (``topo-optimize``).
   * - :mod:`topo.analysis`
     - Native-contact (*Q*) analysis.
   * - :mod:`topo.reporter`
     - The fixed-width state log writer and its parser.
   * - :mod:`topo.utils`
     - Config parsing, the non-bonded contact builder, multi-copy replication, and run provenance.


Installation
------------

TOPO depends on OpenMM (and ParmEd, MDAnalysis, mdtraj, NumPy, pandas, PyYAML),
which are best installed from ``conda-forge``, plus the external **STRIDE** (and
optional **PULCHRA**) binaries. In brief::

    mamba create -n topo -c conda-forge python">=3.9" openmm parmed \
        mdanalysis mdtraj numpy pandas pyyaml
    mamba activate topo
    pip install -e .

See :doc:`installation` for the full step-by-step guide (both the editable/regular
``pip`` install and the no-install ``PYTHONPATH`` alternative), and
:doc:`../usage/external_dependencies` for STRIDE and the optional backmapping
tools (PULCHRA, cg2all).


Running a simulation
--------------------

Every run is driven by a plain-text control file (``md.ini``) and launched with
any of the following equivalent forms::

    topo-mdrun -f md.ini                 # installed console command
    python -m topo.mdrun -f md.ini       # module form
    python run_simulation.py -f md.ini   # the thin shim shipped in each tutorial

See :doc:`../usage/simulation_control` for every ``md.ini`` option,
:doc:`../usage/outputs` for the files a run produces, and
:doc:`../tutorials/index` for hands-on walkthroughs.
