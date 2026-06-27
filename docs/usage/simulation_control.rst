Simulation control options
==========================

Simulation parameters are read from an ``.ini`` file (e.g. ``md.ini``) by
:func:`topo.read_simulation_config`, which returns a
:class:`~topo.utils.config.SimulationConfig`. The section title ``[OPTIONS]`` is
required.

* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Every option below has a default **except** the ones marked *required*; you
  only need to set the options you want to change.

Running a simulation
--------------------

The runner lives in the package as :mod:`topo.mdrun`. Once TOPO is installed
(``pip install -e .`` from the repo root) any of these are equivalent:

.. code-block:: bash

    topo-mdrun -f md.ini              # installed console command
    python -m topo.mdrun -f md.ini    # module form
    python run_simulation.py -f md.ini  # tutorial shim (calls topo.mdrun)

Example ``md.ini``:

.. code-block::

        [OPTIONS]
        md_steps = 500_000   ; number of steps (underscores allowed)
        dt = 0.01            ; time step in ps
        nstxout = 1000       ; steps between trajectory (DCD) writes
        nstlog = 1000        ; steps between log writes
        nstchk = 5000        ; steps between checkpoint writes (default: nstxout)
        nstcomm = 100        ; center-of-mass motion removal; omit for multi-chain runs
        log_precision = 4    ; decimal places for float columns in the .log
        log_width = 14       ; min column width for aligned, fixed-width .log
        model = topo         ; TOPO model (only option currently)
        constraints = AllBonds  ; AllBonds (rigid) or None (flexible bonds)

        ; temperature coupling
        tcoupl = yes
        ref_t = 310          ; Kelvin (also the low/refold temperature when annealing)
        tau_t = 0.01         ; ps^-1

        ; temperature protocol (annealing / quenching) -- off by default.
        ; Two phases: a quench (-> <outname>_quench.dcd/.log) then production
        ; (md_steps at ref_t -> <outname>.dcd/.log). anneal_steps is SEPARATE
        ; from md_steps; grand total = anneal_steps (+ ramp) + md_steps.
        ; anneal = yes       ; hold at t_high then quench/cool down to ref_t
        ; t_high = 600       ; Kelvin (high/unfolding temperature)
        ; anneal_steps = 1_000_000   ; quench-phase steps held at t_high
        ; anneal_ramp = jump ; jump (instant drop) or linear (gradual cool)
        ; anneal_ramp_steps = 500_000      ; linear only: quench-phase ramp t_high -> ref_t
        ; anneal_ramp_increments = 20      ; linear only: discrete T steps in the ramp

        ; pressure coupling
        pcoupl = no
        ref_p = 1
        frequency_p = 25

        ; periodic boundary condition
        pbc = yes
        box_dimension = 30   ; cubic 30 nm; or [30, 30, 60] for a box

        ; input
        pdb_file = 2ww4.pdb
        ; init_position = traj/traj_final.pdb   ; optional starting coordinates
        domain_def = domain.yaml      ; optional
        stride_output_file = stride.dat ; optional
        ; output  (all files -> <output_dir>/<outname>.*, default traj/traj.*)
        output_dir = traj
        outname = traj
        ; hardware
        device = GPU
        ppn = 4
        ; restart
        restart = no
        minimize = no


Parameter summary
+++++++++++++++++

"Required = yes" means the run cannot proceed without it. Options with a default
may be omitted. ``—`` in the *Default* column means there is no default (the
option is either required, or only meaningful in a specific mode noted in the
description).

.. list-table::
   :header-rows: 1
   :widths: 20 14 10 14 42

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``md_steps``
     - int
     - no
     - ``1000``
     - Total number of integration steps. Underscores are allowed (``500_000``).
   * - ``dt``
     - float [ps]
     - no
     - ``0.01``
     - Integration time step.
   * - ``nstxout``
     - int
     - no
     - ``10``
     - Steps between writing the trajectory (DCD).
   * - ``nstlog``
     - int
     - no
     - ``10``
     - Steps between writing the energy/temperature log.
   * - ``nstchk``
     - int
     - no
     - ``nstxout``
     - Steps between writing the checkpoint. Defaults to ``nstxout`` when unset, so checkpoint and trajectory frequency can be decoupled (e.g. frequent checkpoints, sparser frames).
   * - ``nstcomm``
     - int
     - no
     - ``—`` (off)
     - Steps between center-of-mass motion removals. **Unset by default → no COM removal.** COM removal suits a single chain but couples the drift of independent chains, so leave it off for multi-copy runs (``n_copies > 1``).
   * - ``log_precision``
     - int
     - no
     - ``4``
     - Decimal places for floating-point columns (energies, time, temperature, ...) in the ``.log``. Set to ``None`` for OpenMM's full ``repr`` precision.
   * - ``log_width``
     - int
     - no
     - ``14``
     - Minimum width (characters) of each ``.log`` column, right-justified, so columns line up. Each column uses ``max(header_length, log_width)``. Set to ``None`` to disable fixed-width formatting.
   * - ``model``
     - str
     - no
     - ``topo``
     - Force-field model. Only ``topo`` is currently supported.
   * - ``constraints``
     - str
     - no
     - ``AllBonds``
     - Bond treatment: ``AllBonds`` (rigid bonds via constraints) or ``None`` (flexible harmonic bonds). Mutually exclusive.
   * - ``constraint_tolerance``
     - float
     - no
     - ``1e-5``
     - Integrator relative constraint tolerance. Only meaningful with ``constraints = AllBonds``.
   * - ``tcoupl``
     - bool
     - no
     - ``yes``
     - Langevin thermostat on/off. (NVE is not used.)
   * - ``ref_t``
     - float [K]
     - no
     - ``300``
     - Reference temperature. Used when ``tcoupl = yes``.
   * - ``tau_t``
     - float [ps⁻¹]
     - no
     - ``0.01``
     - Friction coefficient coupling the system to the heat bath. Used when ``tcoupl = yes``.
   * - ``anneal``
     - bool
     - no
     - ``no``
     - Temperature protocol. ``no`` → constant-temperature equilibrium at ``ref_t`` (one phase). ``yes`` → a **quench** phase (hold at ``t_high``, then quench/cool to ``ref_t``; written to ``<outname>_quench.dcd``/``.log``) followed by a **production** phase (``md_steps`` at ``ref_t``; written to ``<outname>.dcd``/``.log``). Requires ``tcoupl = yes``. See :doc:`../tutorials/06_anneal`.
   * - ``t_high``
     - float [K]
     - if ``anneal = yes``
     - ``—``
     - High (unfolding) temperature held during the quench phase. ``ref_t`` is reused as the low/refold temperature (there is no ``t_low``).
   * - ``anneal_steps``
     - int (``> 0``)
     - no
     - ``0``
     - Quench-phase steps held at ``t_high``. **Separate from** ``md_steps`` (grand total = ``anneal_steps`` (+ ``anneal_ramp_steps``) + ``md_steps``). Must be many thermal relaxation times (``~1/tau_t``) and long enough to unfold. Used when ``anneal = yes``.
   * - ``anneal_ramp``
     - str
     - no
     - ``jump``
     - ``jump`` = instantaneous drop ``t_high → ref_t`` at the phase boundary (delta T-jump; best for folding kinetics). ``linear`` = gradual cool-down inside the quench phase (best for refolding yield). Used when ``anneal = yes``.
   * - ``anneal_ramp_steps``
     - int
     - no
     - ``0``
     - Additional quench-phase steps spent ramping ``t_high → ref_t`` (on top of ``anneal_steps``). Used only when ``anneal_ramp = linear``.
   * - ``anneal_ramp_increments``
     - int
     - no
     - ``20``
     - Number of discrete temperature steps in the ramp; the last lands exactly on ``ref_t``. Used only when ``anneal_ramp = linear``.
   * - ``pcoupl``
     - bool
     - no
     - ``no``
     - Monte Carlo barostat on/off. Requires ``pbc = yes``.
   * - ``ref_p``
     - float [bar]
     - no
     - ``1``
     - Reference pressure. Used when ``pcoupl = yes``.
   * - ``frequency_p``
     - int [steps]
     - no
     - ``25``
     - Barostat move attempt frequency. Used when ``pcoupl = yes``.
   * - ``pbc``
     - bool
     - no
     - ``no``
     - Periodic boundary conditions on/off.
   * - ``box_dimension``
     - float or [x,y,z] [nm]
     - if ``pbc = yes``
     - ``—``
     - Box size: a scalar ``L`` gives a cubic ``L×L×L`` box; a list ``[x, y, z]`` a rectangular box.
   * - ``pdb_file``
     - str
     - **yes**
     - ``—``
     - Input structure (``.pdb`` / ``.cif``) used to build the model (topology, force field) and, by default, the initial coordinates.
   * - ``init_position``
     - str
     - no
     - ``—``
     - Optional PDB of starting coordinates for a fresh run (atom count must match the system). If unset, the coordinates from ``pdb_file`` are used. Ignored on a successful restart (coordinates come from the checkpoint).
   * - ``output_dir``
     - str
     - no
     - ``traj``
     - Folder for all generated files; created if missing. One run = one self-contained folder.
   * - ``outname``
     - str
     - no
     - ``traj``
     - Basename for generated files: ``<output_dir>/<outname>.dcd``, ``.log``, ``.psf``, ``.chk``, ``_runinfo.log`` (and ``_multi.psf`` for multi-copy).
   * - ``checkpoint``
     - str
     - no
     - ``<output_dir>/<outname>.chk``
     - Explicit checkpoint path override. Normally leave unset so it lands in the run folder.
   * - ``domain_def``
     - str
     - no
     - ``—``
     - Path to a domain YAML for per-domain sidechain-contact scaling. If omitted, all SS contacts use scale 1.0. See :doc:`domain_definition`.
   * - ``stride_output_file``
     - str
     - no
     - ``—``
     - Path to a precomputed STRIDE output. If omitted, STRIDE is run automatically on the structure (and cached to ``{prefix}_stride.dat``).
   * - ``device``
     - str
     - no
     - ``CPU``
     - Compute platform: ``CPU`` or ``GPU`` (CUDA).
   * - ``ppn``
     - int
     - no
     - ``1``
     - Number of CPU threads. Used when ``device = CPU``.
   * - ``n_copies``
     - int
     - no
     - ``1``
     - Number of independent, **non-interacting** copies of the input chain to pack into one simulation (better GPU utilization; ``n_copies`` trajectories per run). ``1`` disables replication. See :doc:`../tutorials/04_multicopy`.
   * - ``copy_shift``
     - float [nm]
     - no
     - ``5.0``
     - Initial x-translation between successive copies. Only used when ``n_copies > 1``; since copies never interact, the exact value affects only the starting layout, not the physics.
   * - ``restart``
     - bool
     - no
     - ``no``
     - Restart from ``checkpoint`` instead of the PDB coordinates. Forces ``minimize = no``.
   * - ``minimize``
     - bool
     - no
     - ``yes``
     - Energy-minimize the input structure before dynamics. Forced ``no`` when ``restart = yes``.

.. note::

   Boolean options accept ``yes``/``no``, ``true``/``false``, ``1``/``0``
   (parsed by :func:`distutils.util.strtobool`).


Notes on individual options
+++++++++++++++++++++++++++

Output layout (``output_dir`` / ``outname``)
    Every generated file is written to ``<output_dir>/<outname><suffix>``, so a run
    is one self-contained folder (default ``traj/``): ``traj.dcd`` (trajectory),
    ``traj.log`` (state log), ``traj.psf`` (single-chain topology), ``traj.chk``
    (checkpoint), ``traj_final.pdb`` (last conformation — reusable as
    ``init_position`` for a follow-up run), ``traj_runinfo.log`` (provenance), and
    ``traj_multi.psf`` for a multi-copy run. ``output_dir`` is created
    automatically if missing. To keep
    several runs side by side, point each at its own folder (e.g.
    ``output_dir = runs/P0CX28_T300``) or change ``outname``.

Output frequency and log formatting (``nstxout`` / ``nstchk`` / ``nstcomm`` / ``log_precision`` / ``log_width``)
    ``nstxout`` controls trajectory (DCD) frames; ``nstchk`` controls checkpoint
    writes and defaults to ``nstxout`` when unset, so the two can be decoupled.
    ``nstcomm`` is **off unless set** — center-of-mass motion removal is only
    appropriate for a single chain (it couples the drift of independent chains),
    so leave it unset for multi-copy runs. The log reporter is
    :class:`topo.topoReporter`, which writes a fixed-width, aligned ``.log``:
    ``log_precision`` sets the decimals for float columns (default ``4``;
    ``None`` for full precision) and ``log_width`` sets the minimum column width
    (default ``14``; ``None`` to disable padding). Columns are separated by two
    spaces and the header is a ``#`` comment, so the log stays both human-readable
    and machine-parsable (see :func:`topo.reporter.readOpenMMReporterFile`).

Bond treatment (``constraints``)
    ``AllBonds`` (default) makes every bond a rigid distance constraint and adds
    no harmonic bond force — appropriate for the standard CA model and required
    for the usual 15 fs / 0.015 ps time step. ``None`` (also ``none`` / empty)
    instead adds a harmonic bond force and no constraints (flexible bonds). The
    two modes are mutually exclusive.

Temperature / pressure coupling
    ``ref_t`` and ``tau_t`` are only consumed when ``tcoupl = yes``; ``ref_p``
    and ``frequency_p`` only when ``pcoupl = yes``. Pressure coupling additionally
    requires ``pbc = yes`` (asserted at parse time).

Temperature protocol (``anneal`` / ``t_high`` / ``anneal_*``)
    By default a run is **equilibrium**: the Langevin thermostat is held at
    ``ref_t`` for all ``md_steps``, written to ``<outname>.dcd`` / ``.log``.
    Setting ``anneal = yes`` splits the run into **two phases**, each writing its
    own trajectory and log:

    * **Quench** — hold at ``t_high`` for ``anneal_steps`` and, for
      ``anneal_ramp = linear``, cool to ``ref_t`` over ``anneal_ramp_steps`` (in
      ``anneal_ramp_increments`` discrete temperature steps). Written to
      ``<outname>_quench.dcd`` / ``.log``.
    * **Production** — ``md_steps`` at ``ref_t``. Written to the usual
      ``<outname>.dcd`` / ``.log``, so the production trajectory is never
      contaminated by the hot hold.

    ``anneal_steps`` is therefore **separate** from ``md_steps``: the grand total
    is ``anneal_steps`` (+ ``anneal_ramp_steps`` for ``linear``) + ``md_steps``.
    For ``anneal_ramp = jump`` the drop to ``ref_t`` is instantaneous and lands
    exactly on the boundary between the two files. ``ref_t`` is always the low /
    refold temperature — there is no separate ``t_low`` key. The OpenMM step
    counter is continuous across phases (the production ``Step`` column starts at
    the quench length), and a single checkpoint ``<outname>.chk`` covers the whole
    run. The runner prints both phases at startup.

    Choose ``jump`` for clean folding kinetics (folding happens at a single
    temperature) and ``linear`` for maximum refolding yield. Crucially, a Langevin
    thermostat relaxes toward a new setpoint over roughly ``1/tau_t``, so
    ``anneal_steps`` must be **many** relaxation times (and long enough to actually
    unfold the protein) — with a production-typical ``tau_t = 0.01`` ps⁻¹ (≈100 ps
    relaxation) that means a hold of many nanoseconds. Annealing requires
    ``tcoupl = yes`` and composes with restarts: a restart resumes whichever phase
    it stopped in (appending to that phase's files). See
    :doc:`../tutorials/06_anneal` for a full walkthrough.

Periodic boundary conditions
    Turning ``pbc`` on affects the non-bonded forces and how coordinates are
    written to the PDB/DCD (handled internally). ``box_dimension`` must be given
    when ``pbc = yes``; an invalid/empty value silently disables PBC.

``domain_def`` and ``stride_output_file``
    Both are optional inputs to the TOPO structure-based non-bonded potential.
    Omit ``domain_def`` for no domain scaling (every SS contact scaled by 1.0).
    Omit ``stride_output_file`` to let the builder run STRIDE for you (STRIDE must
    be on ``PATH``); supply a path to reuse a precomputed file.

``restart`` and ``minimize``
    ``restart = yes`` loads positions **and** velocities from ``checkpoint`` and
    continues; reporters append to the existing log/trajectory, and ``minimize``
    is forced off. With ``restart = no`` you may choose ``minimize``. Note that a
    native input structure is already the energy minimum of the structure-based
    model, so minimization is usually unnecessary.

Initial coordinates and velocities (``init_position``)
    The starting state is resolved as follows (and recorded in
    ``<outname>_runinfo.log`` under ``initial_coordinates`` / ``initial_velocities``):

    * ``restart = yes`` **with** the checkpoint present → positions *and* velocities
      come from the checkpoint.
    * Otherwise (a fresh start, **including** ``restart = yes`` when the checkpoint
      is missing — which prints a warning) → coordinates come from ``init_position``
      if set, else from ``pdb_file``; velocities are drawn from the Boltzmann
      distribution at ``ref_t``.

    Every run also writes ``<outname>_final.pdb`` (the last conformation), which
    you can pass as ``init_position`` to seed a follow-up run.

Hardware (``device`` / ``ppn``)
    ``device = GPU`` runs on CUDA (mixed precision). ``device = CPU`` uses ``ppn``
    threads; ``ppn`` is ignored on GPU.

Multi-copy runs (``n_copies`` / ``copy_shift``)
    A single coarse-grained chain (a few hundred beads) badly underuses a GPU.
    Setting ``n_copies > 1`` packs that many independent copies of the input chain
    into one ``System`` and yields one trajectory per copy. The copies are
    guaranteed non-interacting: bonded terms and constraints are duplicated per
    copy, and every ``CustomNonbondedForce`` (Yukawa electrostatics and the
    structure-based contacts) is restricted to intra-copy interaction groups, so
    the total potential energy is exactly ``n_copies ×`` the single-chain energy.
    ``copy_shift`` sets the initial x-offset between copies (layout only). The
    runner (:mod:`topo.mdrun`) replicates automatically via
    :func:`topo.make_noninteracting_copies` whenever ``n_copies > 1``; afterwards
    split the combined trajectory into per-copy DCDs with
    :func:`topo.split_chains` (memory-bounded streaming; CLI
    ``python -m topo.utils.multichain``). See :doc:`../tutorials/04_multicopy` for
    a complete walkthrough.
