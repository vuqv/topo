Output files and the run log
============================

Every run writes all of its files into **one self-contained folder**, set by
``output_dir`` (default ``traj/``) and named ``<outname>`` (default ``traj``).
This page lists every file TOPO produces, explains the state-log format and how
to parse it, and describes the run-provenance record.


The files a run produces
------------------------

After ``topo-mdrun -f md.ini`` finishes, the run folder contains:

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - File
     - Format
     - Contents / use
   * - ``<outname>.dcd``
     - binary
     - **Trajectory** — coordinates written every ``nstxout`` steps. Open with VMD, MDAnalysis, or mdtraj (load alongside the ``.psf``).
   * - ``<outname>.log``
     - text
     - **State log** — step, time, energies, temperature, speed, ETA, every ``nstlog`` steps. Fixed-width and machine-parsable (see below).
   * - ``<outname>.psf``
     - text
     - **Topology** of the CA model (atoms, masses, charges, bonds). Needed to load the ``.dcd`` in analysis tools.
   * - ``<outname>.chk``
     - binary
     - **Checkpoint** — full dynamical state (positions **and** velocities), every ``nstchk`` steps. Used to restart (:doc:`../tutorials/03_restart`).
   * - ``<outname>_final.pdb``
     - text
     - **Last conformation** (CA PDB). Reuse as ``init_position`` to seed a follow-up run.
   * - ``<outname>_runinfo.log``
     - text (INI)
     - **Run provenance** — package versions, hardware, GPU, timing, coordinate/velocity sources. See below.

**Conditional files:**

* ``<outname>_multi.psf`` — written when ``n_copies > 1``: the combined
  multi-chain topology matching the multi-copy ``.dcd`` (the plain ``.psf``
  remains the single-chain topology). See :doc:`../tutorials/04_multicopy`.
* ``<outname>_quench.dcd`` / ``<outname>_quench.log`` — written when
  ``anneal = yes``: the quench-phase trajectory and log, kept separate from the
  production ``.dcd`` / ``.log`` so the hot phase never contaminates your
  production ensemble. See :doc:`../tutorials/06_anneal`.
* ``<pdb_prefix>_stride.dat`` — cached STRIDE hydrogen-bond output, written **next
  to the input PDB** (not in the run folder) the first time the model is built.
  Reused automatically on later runs; delete to force regeneration.

.. tip::

   Each run **overwrites** the ``traj/`` folder. To keep several runs side by
   side, point each at its own folder (``output_dir = runs/P0CX28_T300``) or
   change ``outname``.


The state log format
--------------------

The ``.log`` is written by :class:`topo.reporter.topo_reporter.topoReporter`, a
subclass of OpenMM's ``StateDataReporter`` tuned for readability **and** easy
parsing. Each float column is rounded to ``log_precision`` decimals (default 4)
and right-justified to ``log_width`` characters (default 14), columns are
separated by two spaces, and the header is a ``#`` comment line:

.. code-block:: text

   #         Step       Time (ps)   Potential Energy (kJ/mole)   Kinetic Energy (kJ/mole)   ...   Temperature (K)   Speed (ns/day)   Time Remaining
              1000           15.0                   -1234.5678                    312.4567   ...          301.2345         123.4567       0:01:23
              2000           30.0                   -1241.0021                    305.1180   ...          299.8765         124.0012       0:01:10

Default columns: step, time (ps), potential energy, kinetic energy, total energy,
temperature, speed (ns/day), and estimated time remaining. The energies are in
**kJ/mol**.

* **``log_precision``** — decimals for float columns; ``None`` gives OpenMM's full
  ``repr`` precision.
* **``log_width``** — minimum column width for alignment; ``None`` disables
  fixed-width padding.

On a **restart**, the reporter *appends* to the existing ``.log`` and ``.dcd``
(no repeated header), so a multi-stage run stays one continuous record. In an
**annealing** run, the quench writes ``_quench.log`` and the production log's
step counter is reset to 0 (see :doc:`../tutorials/06_anneal`).

Parsing the log in Python
~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the companion reader, which returns a ``{column_name: [values...]}`` dict
(numeric columns as floats, others as strings) and respects the two-space
separator so multi-word headers like ``Potential Energy (kJ/mole)`` stay intact:

.. code-block:: python

   import numpy as np
   from topo.reporter.topo_reporter import readOpenMMReporterFile

   data = readOpenMMReporterFile("traj/traj.log")
   step = np.array(data["Step"])
   pe   = np.array(data["Potential Energy (kJ/mole)"])
   temp = np.array(data["Temperature (K)"])
   print("mean T:", round(temp.mean(), 1), "K")

This is the recommended way to plot energy or temperature vs. time, or to check
that a run is healthy (stable temperature near ``ref_t``, non-exploding
potential energy).


The run-provenance record
-------------------------

``<outname>_runinfo.log`` is an INI-format side channel (written by
:mod:`topo.utils.runinfo`) that does **not** affect the simulation but records
everything you need to reproduce or debug it:

* ``[run_start]`` — start timestamp; the control file and checkpoint paths;
  whether this was a restart; steps planned; **initial coordinate and velocity
  sources** (checkpoint vs. ``init_position`` vs. ``pdb_file``, and Boltzmann
  velocities at ``ref_t``); Python / NumPy / ParmEd / OpenMM versions; hostname,
  OS, CPU count, requested threads; and the selected OpenMM platform.
* ``[cuda_metadata]`` — on GPU runs: CUDA device name, driver/runtime versions,
  precision, and an ``nvidia-smi`` snapshot.
* ``[run_end]`` — end timestamp, wall-clock elapsed time, final step and
  simulated time, and the path to the final structure.

Keep this file with your trajectory: it answers "what exactly did I run, on what
hardware, from which coordinates?" long after the job is gone.


Loading a trajectory for analysis
---------------------------------

.. code-block:: python

   import MDAnalysis as mda
   u = mda.Universe("traj/traj.psf", "traj/traj.dcd")     # single chain
   print(u.atoms.n_atoms, "CA beads,", len(u.trajectory), "frames")

   # multi-copy run: use the combined topology
   u_multi = mda.Universe("traj/traj_multi.psf", "traj/traj.dcd")

From there, compute RMSD, radius of gyration, or the native-contact score *Q*
(:doc:`native_contacts`). For multi-copy runs, split the combined DCD into
per-copy trajectories first with :func:`topo.split_chains`
(:doc:`../tutorials/04_multicopy`).
