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

Example ``md.ini``:

.. code-block::

        [OPTIONS]
        md_steps = 500_000   ; number of steps (underscores allowed)
        dt = 0.01            ; time step in ps
        nstxout = 1000       ; steps between trajectory/checkpoint writes
        nstlog = 1000        ; steps between log writes
        nstcomm = 100        ; steps between center-of-mass motion removal
        model = topo         ; TOPO model (only option currently)
        constraints = AllBonds  ; AllBonds (rigid) or None (flexible bonds)

        ; temperature coupling
        tcoupl = yes
        ref_t = 310          ; Kelvin
        tau_t = 0.01         ; ps^-1

        ; pressure coupling
        pcoupl = no
        ref_p = 1
        frequency_p = 25

        ; periodic boundary condition
        pbc = yes
        box_dimension = 30   ; cubic 30 nm; or [30, 30, 60] for a box

        ; input
        protein_code = 2ww4
        pdb_file = 2ww4.pdb
        domain_def = domain.yaml      ; optional
        stride_output_file = stride.dat ; optional
        ; output
        checkpoint = 2ww4.chk
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
     - Steps between writing the trajectory (DCD) and checkpoint.
   * - ``nstlog``
     - int
     - no
     - ``10``
     - Steps between writing the energy/temperature log.
   * - ``nstcomm``
     - int
     - no
     - ``100``
     - Steps between center-of-mass motion removals.
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
     - Input structure (``.pdb`` / ``.cif``) for topology and initial coordinates.
   * - ``protein_code``
     - str
     - **yes**
     - ``—``
     - Output filename prefix, e.g. ``{protein_code}.dcd``, ``{protein_code}.log``.
   * - ``checkpoint``
     - str
     - **yes**
     - ``—``
     - Checkpoint file (``.chk``); written during the run and read on restart.
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
    standard ``run_simulation.py`` replicates automatically via
    :func:`topo.make_noninteracting_copies` whenever ``n_copies > 1``; afterwards
    split the combined trajectory into per-chain DCDs for analysis. See
    :doc:`../tutorials/04_multicopy` for a complete walkthrough.
