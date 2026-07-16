Using TOPO from Python
======================

Most users drive TOPO entirely through ``md.ini`` and the ``topo-mdrun`` command
(:doc:`simulation_control`). But every step is also a documented Python function,
so you can build a model, inspect its forces, script a custom protocol, or post-
process trajectories from your own code.

The most useful names are re-exported at the top level of the ``topo`` package:

.. code-block:: python

   import topo
   # model building
   topo.models.buildCoarseGrainModel
   topo.system
   # config + reporter
   topo.read_simulation_config, topo.SimulationConfig, topo.topoReporter
   # multi-copy + trajectory splitting
   topo.make_noninteracting_copies, topo.split_chains
   # native-contact (Q) analysis
   topo.build_native_contacts, topo.fraction_native_contacts, topo.load_domains


Build a model from a structure
------------------------------

:meth:`topo.models.buildCoarseGrainModel` is the single entry point. It reads a
structure, keeps the Cα atoms, and assembles the full force field (bonds, angles,
torsions, Yukawa electrostatics, and the structure-based contacts). It returns a
:class:`topo.core.system.system` object holding the OpenMM ``System``,
``Topology``, and ``positions``.

.. code-block:: python

   import topo

   cg = topo.models.buildCoarseGrainModel(
       "P0CX28_clean.pdb",
       domain_def="domain.yaml",       # optional: per-domain contact scaling
       stride_output_file=None,        # None -> run STRIDE automatically, cache it
       constraints="AllBonds",         # rigid bonds (default); None -> flexible
       box_dimension=None,             # None -> no PBC; float or [x,y,z] -> PBC
       minimize=False,                 # native PDB is already the minimum
   )

   print(cg.n_chains, "chains,", cg.n_atoms, "CA beads")
   print(type(cg.system))              # openmm.System with all forces added

Key arguments (full reference in :class:`topo.core.models.models`):

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - Argument
     - Default
     - Meaning
   * - ``structure_file``
     - *required*
     - Input PDB/CIF; defines topology, force field, and native contacts.
   * - ``domain_def``
     - ``None``
     - ``domain.yaml`` for per-domain/interface sidechain-contact scaling (:doc:`domain_definition`).
   * - ``stride_output_file``
     - ``None``
     - Precomputed STRIDE file; ``None`` runs STRIDE and caches ``<prefix>_stride.dat``.
   * - ``constraints``
     - ``'AllBonds'``
     - ``'AllBonds'`` = rigid (constraints, 15 fs step); ``None`` = flexible harmonic bonds.
   * - ``box_dimension``
     - ``None``
     - ``None`` = no PBC; a float = cubic box (nm); ``[x,y,z]`` = rectangular box.
   * - ``minimize``
     - ``False``
     - Energy-minimize the input geometry if large forces are found.
   * - ``check_forces``
     - ``True``
     - Run the build-time energy/large-force check (set ``False`` on restart).

Useful attributes and methods of the returned object:

* ``cg.system`` / ``cg.topology`` / ``cg.positions`` — the OpenMM objects to feed
  into a ``Simulation``.
* ``cg.n_atoms`` / ``cg.n_chains`` / ``cg.n_bonds`` / ``cg.n_angles`` /
  ``cg.n_torsions`` — geometry counts.
* ``cg.rmin_matrix`` / ``cg.energy_matrix`` — the contact :math:`R_{ij}` and
  :math:`\varepsilon_{ij}` matrices (nm, kJ/mol).
* ``cg.forceGroups`` — ordered map of force name → force object (one per log
  column).
* ``cg.dumpTopology("model.psf")`` — write the CA topology as PSF.
* ``cg.dumpStructure("model.pdb")`` — write the current CA coordinates as PDB.
* ``cg.dumpForceFieldData("ff.txt")`` — dump all force-field parameters to text.
* ``cg.reportEnergy(simulation, header=...)`` — print the total and per-force-group
  energy of a simulation's current state.


Run dynamics yourself
---------------------

If you want full control over the integrator and protocol, drive OpenMM directly.
This reproduces the core of what :mod:`topo.mdrun` does:

.. code-block:: python

   import openmm as mm
   from openmm import unit, app
   import topo

   cg = topo.models.buildCoarseGrainModel("P0CX28_clean.pdb", domain_def="domain.yaml")

   integrator = mm.LangevinIntegrator(300*unit.kelvin, 0.01/unit.picosecond,
                                      0.015*unit.picoseconds)
   sim = app.Simulation(cg.topology, cg.system, integrator)
   sim.context.setPositions(cg.positions)
   sim.context.setVelocitiesToTemperature(300*unit.kelvin)

   # TOPO's fixed-width log with one energy column per force group:
   sim.reporters.append(topo.topoReporter("traj.log", 1000, sbmObject=cg,
                                          step=True, time=True,
                                          potentialEnergy=True, temperature=True))
   sim.reporters.append(app.DCDReporter("traj.dcd", 1000))
   sim.step(100000)

Passing ``sbmObject=cg`` to :class:`~topo.reporter.topo_reporter.topoReporter`
adds one energy column per force group, so you can watch the contact energy,
angle energy, etc., separately.

In practice you usually don't need this — :func:`topo.read_simulation_config`
plus the :mod:`topo.engine` / :mod:`topo.mdrun` helpers already wrap build →
set-up → protocol → finalize. Reach for the manual route only for a non-standard
protocol the config file can't express.


Start from a hackable script
----------------------------

If you want a custom protocol but not a from-scratch OpenMM loop, start from the
ready-to-edit script in the repository's ``examples/`` directory rather than
writing one from nothing. ``examples/custom_md.py`` is the *open* version of the
``topo-mdrun`` command: it reproduces the runner's exact flow using the
:mod:`topo.engine` building blocks, with the customization points marked
``# === EDIT: ... ===``.

.. code-block:: bash

   cp examples/custom_md.py my_run.py     # copy, then edit the EDIT sections
   python my_run.py -f md.ini             # runs like topo-mdrun, reads the same md.ini

The engine layer it composes is the same one :mod:`topo.mdrun` uses, so an edited
copy behaves exactly like the console command except where you change it:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Helper
     - Does
   * - :func:`topo.engine.build_system`
     - Build the (optionally multi-copy) CG model; returns a ``BuiltSystem`` with ``.system`` / ``.topology`` / ``.positions``. Add your own forces to ``built.system`` here, before setup.
   * - :func:`topo.engine.setup_simulation`
     - Create the OpenMM ``Simulation`` (integrator, platform, restart, starting coordinates); returns a ``RunContext``.
   * - :func:`topo.engine.attach_reporters`
     - Attach the DCD / TOPO-log / checkpoint reporters.
   * - :func:`topo.mdrun.protocol.run_protocol`
     - Step the simulation through a schedule -- a list of ``(temperature, n_steps)`` stages you define.
   * - :func:`topo.engine.finalize_simulation`
     - Save the final checkpoint and structure and close out run metadata.

The three worked customizations in the script are: adding an OpenMM force to the
built system (e.g. a restraint), defining your own temperature schedule instead
of constant-\ *T* production, and running in segments with an analysis callback
between them. See ``examples/README.md`` for the copy-and-edit workflow.


Read a control file
-------------------

:func:`topo.read_simulation_config` parses an ``md.ini`` into a
:class:`~topo.utils.config.SimulationConfig` dataclass with OpenMM units already
applied — handy for scripting parameter sweeps:

.. code-block:: python

   import topo
   cfg = topo.read_simulation_config("md.ini")
   print(cfg.md_steps, cfg.ref_t, cfg.dt)
   kwargs = cfg.build_kwargs()          # ready to pass to buildCoarseGrainModel
   cg = topo.models.buildCoarseGrainModel(cfg.pdb_file, **kwargs)

``SimulationConfig`` also resolves output paths (``cfg.output_path('.dcd')``),
the checkpoint path (``cfg.checkpoint_path()``), the OpenMM platform
(``cfg.make_platform()``), and the annealing step counts (``cfg.quench_steps()``,
``cfg.total_steps()``).


Multi-copy replication and splitting
------------------------------------

Pack independent copies of a chain into one system, then split the combined
trajectory afterwards (:doc:`../tutorials/04_multicopy`):

.. code-block:: python

   import topo
   from openmm import unit

   cg = topo.models.buildCoarseGrainModel("P0CX28_clean.pdb", domain_def="domain.yaml")

   system, topology, positions = topo.make_noninteracting_copies(
       cg.system, cg.topology, cg.positions,
       n_copies=10, shift=5.0*unit.nanometer)

   # ... run, producing a combined traj.dcd ...

   # split back into one DCD per copy (memory-bounded streaming):
   topo.split_chains("traj.dcd",
                     [f"traj_{k}.dcd" for k in range(10)],
                     center=True)

The copies are guaranteed non-interacting: bonded terms are duplicated per copy
and each ``CustomNonbondedForce`` is restricted to intra-copy interaction groups,
so the total energy is exactly ``n_copies ×`` the single-chain energy.


Score native contacts (*Q*)
---------------------------

The :doc:`native_contacts` analysis is fully scriptable; see that page for a
complete example using :func:`topo.build_native_contacts` and
:func:`topo.fraction_native_contacts`.


Parse a log
-----------

:func:`topo.reporter.topo_reporter.readOpenMMReporterFile` reads a TOPO ``.log``
into a ``{column: [values]}`` dict (see :doc:`outputs`):

.. code-block:: python

   from topo.reporter.topo_reporter import readOpenMMReporterFile
   data = readOpenMMReporterFile("traj/traj.log")
   pe = data["Potential Energy (kJ/mole)"]


See also
--------

* :doc:`model_theory` — the physics behind every force the builder adds.
* :doc:`simulation_control` — the ``md.ini`` options each argument mirrors.
* The API reference pages: :doc:`../modules/models`, :doc:`../modules/system`,
  :doc:`../modules/parameters`.
