.. _optimization-control:

Optimization control options
============================

Optimizer parameters are read from an ``.ini`` file (e.g. ``optimize.ini``) by
:func:`topo.optimize.read_optimize_config`. For *why* the nscale needs
optimizing at all, see :doc:`nscale_optimization`; for a worked run see the
tutorial :doc:`../tutorials/05_opt_nscal`.

* The file is a flat list of ``key = value`` lines.
* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Only ``pdb_file`` and ``domain_def`` are *required*; everything else has a
  default.

.. important::

   **How this file differs from** ``md.ini``. The optimizer is not a simulation
   engine — it *drives* :mod:`topo.mdrun` in a loop. So ``optimize.ini`` is one
   flat key list that is **split two ways** when it is read:

   * the five **optimizer controls** (``ntraj``, ``q_threshold``,
     ``frame_fraction``, ``max_rounds``, ``min_contacts``) are consumed by the
     optimizer and *removed* — they never reach ``md.ini``;
   * **every other key** is treated as a simulation parameter and passed through
     verbatim into each round's generated ``round_N/md.ini``.

   There is no list of "allowed" simulation keys: the pass-through is
   open-ended, so any option from :doc:`simulation_control` may appear here and
   will reach the MD run. The flip side is that a **misspelled control key is
   not an error** — it is silently forwarded to ``md.ini`` as an unknown
   simulation parameter, and the control keeps its default. Check the header of
   ``optimization.log``, which echoes the controls actually in force.

   Two further consequences are covered below: the pass-through defaults are
   **not** the ``md.ini`` defaults (:ref:`Options <opt-options>`), and a few keys
   are **overridden every round** (:ref:`opt-overridden`).

Running the optimizer
---------------------

The optimizer lives in the package as :mod:`topo.optimize`. Once TOPO is
installed (``pip install -e .`` from the repo root) either of these works:

.. code-block:: bash

    topo-optimize -f optimize.ini -o opt_out        # installed console command
    python -m topo.optimize -f optimize.ini -o opt_out   # module form

A bare ``topo-optimize`` with no arguments prints help and exits.

Example ``optimize.ini``:

.. code-block::

        ; --- inputs (paths are resolved relative to THIS file) ---
        pdb_file   = P0A6E6.pdb    ; all-atom reference (native contacts + geometry)
        domain_def = domain.yaml   ; initial domains + class; nscales get optimized

        ; --- per-trajectory production (passed through to each round's md.ini) ---
        md_steps   = 50_000        ; steps per trajectory — SET THIS (see note below)
        nstxout    = 100           ; trajectory write frequency (these frames feed Q)
        nstlog     = 100
        ref_t      = 300           ; K; stability protocol temperature

        ; Override an implicit default by adding it here, e.g. to run on CPU:
        ; device = CPU

        ; Optional: skip STRIDE entirely by pointing at a precomputed file. If
        ; omitted, STRIDE runs once in round 1 and the result is reused.
        ; stride_output_file = P0A6E6_stride.dat

        ; --- optimizer controls (consumed here; never reach md.ini) ---
        ntraj          = 10        ; independent trajectories per round (= n_copies)
        q_threshold    = 0.6688    ; a frame is "folded" for a unit if Q > this
        frame_fraction = 0.98      ; a traj is "stable" if >= this fraction folded
        max_rounds     = 6         ; 5 ladder levels + median fallback
        min_contacts   = 0         ; 0 disables the too-few-contacts freeze

.. note::

   **Set** ``md_steps`` **explicitly.** Its implicit default is only **10 000
   steps** (150 ps at ``dt = 0.015``) — a smoke-test length, far too short to
   distinguish a folded domain from an unfolded one, and the optimizer will
   happily report "converged" on meaningless trajectories. The published
   protocol uses ~0.5 µs per trajectory (≈ ``3.3e7`` steps at ``dt = 0.015``).
   Lower ``nstxout`` alongside it: at its ``5000``-step default a 10 000-step
   round yields two frames per trajectory to compute Q from.

.. _opt-options:

Options
-------

Everything worth deciding about, in one table. The keys are of three kinds:

* **Inputs** (``pdb_file``, ``domain_def``, ``stride_output_file``) — the
  structure and domains the model is built from. The first two are the only
  required keys in the file.
* **Controls** (``ntraj``, ``q_threshold``, ``frame_fraction``, ``max_rounds``,
  ``min_contacts``) — consumed by the optimizer itself, and the only keys that do
  **not** reach the per-round ``md.ini``.
* **Simulation parameters** (``md_steps``, ``nstxout``, ``device``, ``ref_t``,
  ``minimize``) — passed through to each round's ``md.ini``. All are filled in
  from the optimizer's own protocol defaults, but ``md_steps`` and ``nstxout``
  are the two you should set yourself, since their defaults only describe a
  smoke test. ``device`` and ``minimize`` deliberately differ from what the same
  key would default to in a plain :doc:`simulation_control` run.

.. list-table::
   :header-rows: 1
   :widths: 18 8 10 12 52

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``pdb_file``
     - str
     - **yes**
     - ``—``
     - All-atom reference structure. Defines both the native contacts used for Q scoring and the geometry the model is built from. Resolved relative to the ``.ini`` file.
   * - ``domain_def``
     - str
     - **yes**
     - ``—``
     - Seed domain YAML. Its ``residues`` and ``class`` fields are read and kept; its ``nscale`` values are **ignored and overwritten** each round (they are what is being optimized). Resolved relative to the ``.ini`` file. See :doc:`domain_definition`.
   * - ``stride_output_file``
     - str
     - no
     - ``—``
     - Precomputed STRIDE output. If omitted, STRIDE runs during round 1 and the cached file is reused for every later round (STRIDE depends only on the fixed structure, so it never changes). Resolved to an absolute path, since each round's ``md.ini`` lives in a different directory.
   * - ``md_steps``
     - int
     - no
     - ``10_000``
     - **Set this explicitly** (see the note above). Steps per trajectory, passed through to ``md.ini``. The default is a smoke-test length — far too short to tell a folded domain from an unfolded one.
   * - ``nstxout``
     - int
     - no
     - ``5000``
     - Steps between trajectory (DCD) writes, passed through to ``md.ini``. These frames are what Q is computed from, so this sets the sampling resolution of the stability decision — lower it alongside ``md_steps``. ``nstlog`` and ``nstchk`` default to ``5000`` the same way.
   * - ``ntraj``
     - int
     - no
     - ``10``
     - Independent trajectories per round. Produced as a single **multi-copy** MD run (``n_copies = ntraj``) that is then split into per-copy DCDs, so this scales the cost of every round roughly linearly.
   * - ``q_threshold``
     - float
     - no
     - ``0.6688``
     - A frame counts as *folded* for a unit if its fraction of native contacts Q > this value. See :doc:`native_contacts`.
   * - ``frame_fraction``
     - float
     - no
     - ``0.98``
     - A trajectory is *stable* for a unit if at least this fraction of its frames are folded. A unit is stable only when **all** ``ntraj`` trajectories pass — one unstable trajectory is enough to promote the unit.
   * - ``max_rounds``
     - int
     - no
     - ``6``
     - Hard cap on rounds. The default covers the normal case exactly: 5 ladder levels + the median fallback, climbing one level per round. More rounds are needed only if raising one unit's nscale destabilizes an already-stable unit (rare).
   * - ``min_contacts``
     - int
     - no
     - ``0``
     - A domain or interface with fewer than this many native contacts is considered too weakly structured to fold: it is **pinned at the first ladder level, frozen, and never optimized** (it cannot otherwise stabilize, so it would drag every round to the fallback). ``0`` disables the check.
   * - ``device``
     - str
     - no
     - ``GPU``
     - Passed through to ``md.ini``, where the default would be ``CPU``. Optimization is a production-scale task (``ntraj`` × ``max_rounds`` trajectories).
   * - ``ref_t``
     - float
     - no
     - ``300`` K
     - Passed through to ``md.ini``. The stability protocol temperature; pinned here so each round's generated file records it explicitly.
   * - ``minimize``
     - bool
     - no
     - ``no``
     - Passed through to ``md.ini``, where the default would be ``yes``. The native structure is already the model's energy minimum.

The table stops at the keys worth a decision. The remaining implicit defaults
match the model's standard parameterization and are rarely touched:
``dt = 0.015`` (ps), ``model = topo``, ``tcoupl = yes``, ``tau_t = 0.05``,
``pcoupl = no``, ``pbc = no``, ``restart = no`` (every round starts fresh),
``ppn = 4``, ``nstlog = 5000``, ``nstchk = 5000``. Setting any of these in ``optimize.ini`` overrides the implicit
value. The full set is :data:`topo.optimize.IMPLICIT_DEFAULTS`. Beyond these,
any other key from :doc:`simulation_control` may appear and is passed through
verbatim.

.. _opt-overridden:

Keys the optimizer overrides every round
----------------------------------------

The optimizer owns the round layout and the multi-copy fan-out, so it rewrites
these keys when it expands ``optimize.ini`` into ``round_N/md.ini``. **Setting
them in** ``optimize.ini`` **has no effect** — your value is silently replaced:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Key
     - Value forced by the optimizer
   * - ``n_copies``
     - ``ntraj`` — the round's trajectories are one multi-copy run.
   * - ``output_dir``
     - ``<outdir>/round_N/traj``
   * - ``outname``
     - ``traj``
   * - ``pdb_file``, ``domain_def``
     - Rewritten to absolute paths; ``domain_def`` points at the round's freshly written ``round_N/domain.yaml``, not your seed file.

``device`` and ``md_steps`` are also overridden, but only when you pass the
corresponding command-line flag (below); otherwise your ``.ini`` value stands.

Command-line flags
------------------

.. list-table::
   :header-rows: 1
   :widths: 22 14 20 44

   * - Flag
     - Type
     - Default
     - Description
   * - ``-f``, ``--config``
     - str
     - ``—`` (**required**)
     - Path to ``optimize.ini``.
   * - ``-o``, ``--outdir``
     - str
     - ``opt_out``
     - Optimization root directory. Round subdirectories and the final model are written here.
   * - ``--device``
     - str
     - ``—``
     - Override ``device`` (``CPU``/``GPU``) for every round, without editing the ``.ini``.
   * - ``--md-steps``
     - int
     - ``—``
     - Override ``md_steps`` for every round. Handy for a quick smoke test before committing to a production run.
   * - ``--python``
     - str
     - current interpreter
     - Python interpreter used to launch :mod:`topo.mdrun` subprocesses.

What the optimizer writes
-------------------------

Everything lands under ``--outdir`` (default ``opt_out/``):

.. code-block:: text

    opt_out/
      optimization.log          # full report: per-round nscales, Q verdicts, convergence
      domain_optimized.yaml     # <- the calibrated model; this is the deliverable
      round_1/
        md.ini                  # generated: optimize.ini + implicit defaults + overrides
        domain.yaml             # this round's nscales
        traj/
          traj.dcd, traj.psf    # the multi-copy run
          traj_0.dcd ...        # split per-copy trajectories
          Q_0.csv ...           # per-frame Q per unit, paired with traj_<k>.dcd
      round_2/ ...

``optimization.log`` is line-buffered and flushed per line, so ``tail -f`` works
during a long run. The per-round ``Q_<k>.csv`` files hold the Q time series
behind every stability decision — read them when a unit refuses to stabilize.

.. warning::

   **A non-converged run still writes** ``domain_optimized.yaml``. If a protein
   does not stabilize within ``max_rounds``, the optimizer leaves the unstable
   units at their highest level (or the median fallback), writes the file
   anyway, and flags the run with a boxed ``WARNING`` naming the unstable units.
   The file's presence is therefore **not** evidence of a calibrated model —
   check the log for ``CONVERGED`` before using it.

Limitations
-----------

The optimization is **not resumable**: each invocation starts from ladder level
1 and re-runs every round. A long run that is interrupted must start over.
