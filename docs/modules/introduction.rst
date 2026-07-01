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


What TOPO can do (feature overview)
-----------------------------------

A newcomer's map of the full feature set, with where each is documented:

.. list-table::
   :header-rows: 1
   :widths: 30 44 26

   * - Feature
     - What it gives you
     - Learn it in
   * - **Structure-based CA model**
     - One bead/residue; native fold = energy minimum; full force field built automatically from a PDB.
     - :doc:`../usage/model_theory`
   * - **Single-domain simulation**
     - The minimal workflow: one config file, one structure, run MD, read outputs.
     - :doc:`../tutorials/01_single_domain`
   * - **Per-domain / interface contact scaling**
     - Tune the stability of each domain and each domain–domain interface independently (incl. discontiguous domains).
     - :doc:`../tutorials/02_multidomain`, :doc:`../usage/domain_definition`
   * - **Restart from checkpoint**
     - Continue long runs across wall-clock limits; logs/trajectories append seamlessly.
     - :doc:`../tutorials/03_restart`
   * - **Many copies in one run**
     - Pack N non-interacting chains into one (GPU-filling) simulation → N independent trajectories; split afterwards.
     - :doc:`../tutorials/04_multicopy`
   * - **Nscale optimization**
     - Automatically search the contact nscale that keeps each domain/interface folded, instead of guessing.
     - :doc:`../tutorials/05_opt_nscal`
   * - **Temperature annealing / quenching**
     - Hold hot to unfold, then T-jump or slow-cool to study refolding; quench and production write separate trajectories.
     - :doc:`../tutorials/06_anneal`
   * - **Native-contact (Q) analysis**
     - Measure how folded the protein (and each domain/interface) is, frame by frame.
     - :doc:`../usage/native_contacts`
   * - **Flexible vs. rigid bonds, PBC, pressure coupling, custom I/O**
     - Constraint mode, periodic box + barostat, output frequency and log formatting.
     - :doc:`../usage/simulation_control`
   * - **Python API**
     - Build a model, inspect its forces, replicate/split trajectories from your own scripts.
     - :doc:`../usage/python_api`


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
which are best installed from ``conda-forge``. From the repository root::

    mamba create -n topo -c conda-forge python">=3.9" openmm parmed \
        mdanalysis mdtraj numpy pandas pyyaml
    mamba activate topo
    pip install -e .

This installs TOPO as an editable package (your source edits take effect
immediately) and registers the ``topo-mdrun`` and ``topo-optimize`` console
commands. **STRIDE** must be on your ``PATH`` for the contact potential (or
supply a precomputed STRIDE file via ``stride_output_file``). See the project
``README`` for a no-install (``PYTHONPATH``) alternative.


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
