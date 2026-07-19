.. _optimization-control:

Optimization control options
============================

Optimizer parameters are read from an ``.ini`` file (e.g. ``optimize.ini``) by
:func:`topo.optimize.read_optimize_config`. For *why* the nscale needs
optimizing at all, see :doc:`nscale_optimization`; for a worked run see the
tutorial :doc:`../tutorials/A5_opt_nscal`.

* The file is a flat list of ``key = value`` lines.
* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Only ``pdb_file`` and ``domain_def`` are *required*; everything else has a
  default.

.. important::

   **How this file differs from** ``md.ini``. The optimizer is not a simulation
   engine — it *drives* :mod:`topo.mdrun` in a loop. So ``optimize.ini`` is one
   flat key list that is **split two ways** when it is read:

   * the six **optimizer controls** (``ntraj``, ``q_threshold``,
     ``frame_fraction``, ``max_rounds``, ``min_contacts``, ``outdir``) are
     consumed by the optimizer and *removed* — they never reach ``md.ini``;
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

    topo-optimize -f optimize.ini        # installed console command
    python -m topo.optimize -f optimize.ini   # module form

``-f``/``--config`` is the only required argument: it names the ``optimize.ini``
below, and is the one setting that cannot live *inside* that file. Every other
setting is a key in :ref:`Options <opt-options>`. Three of those keys also have a
command-line flag for overriding the file without editing it — ``-o``/``--outdir``,
``--device`` and ``--md-steps`` — each documented on its own row in that table.

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
        ; omitted, STRIDE runs once up front (written to the optimization root)
        ; and every round reuses that file.
        ; stride_output_file = P0A6E6_stride.dat

        ; --- optimizer controls (consumed here; never reach md.ini) ---
        ntraj          = 10        ; independent trajectories per round (= n_copies)
        q_threshold    = 0.6688    ; a frame is "folded" for a unit if Q > this
        frame_fraction = 0.98      ; a traj is "stable" if >= this fraction folded
        max_rounds     = 6         ; 5 ladder levels + median fallback
        min_contacts   = 0         ; 0 disables the too-few-contacts freeze
        ; outdir       = opt_out   ; optimization root; relative to THIS file.
                                   ; -o on the command line overrides it.

.. note::

   **Set** ``md_steps`` **explicitly.** Its implicit default is only **10 000
   steps** (150 ps at ``dt = 0.015``) — a smoke-test length, far too short to
   distinguish a folded domain from an unfolded one, and the optimizer will
   happily report "converged" on meaningless trajectories. The published
   protocol uses ~1 µs per trajectory (≈ ``6.7e7`` steps at ``dt = 0.015``).
   Lower ``nstxout`` alongside it: at its ``5000``-step default a 10 000-step
   round yields two frames per trajectory to compute Q from.

.. _opt-options:

Options
-------

The keys fall into four categories:

* **System input & force field** (``pdb_file``, ``domain_def``,
  ``stride_output_file``) — the structure and domains the model is built from.
  The first two are the only required keys in the file.
* **Per-trajectory production** (``md_steps``, ``nstxout``) — passed through to
  each round's ``md.ini``. These are the two you should set yourself, since
  their defaults only describe a smoke test.
* **Optimizer controls** (``ntraj``, ``q_threshold``, ``frame_fraction``,
  ``max_rounds``, ``min_contacts``, ``outdir``) — consumed by the optimizer
  itself, and the only keys that do **not** reach the per-round ``md.ini``.
* **Protocol & hardware** (``device``, ``ref_t``, ``minimize``) — also passed
  through, but filled in from the optimizer's own protocol defaults. ``device``
  and ``minimize`` deliberately differ from what the same key would default to
  in a plain :doc:`simulation_control` run.

.. list-table::
   :header-rows: 1
   :widths: 18 8 10 12 52

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - **System input & force field**
     -
     -
     -
     -
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
     - Precomputed STRIDE output. If omitted, STRIDE runs once up front, is written to the optimization root, and every round reuses that file (STRIDE depends only on the fixed structure, so it never changes). Resolved to an absolute path, since each round's ``md.ini`` lives in a different directory.
   * - **Per-trajectory production** — passed through to each round's ``md.ini``
     -
     -
     -
     -
   * - ``md_steps``
     - int
     - no
     - ``10_000``
     - **Set this explicitly** (see the note above). Steps per trajectory, passed through to ``md.ini``. **Recommended: ~1 µs per trajectory** (≈ ``6.7e7`` steps at ``dt = 0.015``), the published protocol length. The ``10_000`` default is a smoke-test length — far too short to tell a folded domain from an unfolded one. ``--md-steps`` on the command line overrides it for every round, handy for a smoke test before committing to a production run.
   * - ``nstxout``
     - int
     - no
     - ``5000``
     - Steps between trajectory (DCD) writes, passed through to ``md.ini``. These frames are what Q is computed from, so this sets the sampling resolution of the stability decision — lower it alongside ``md_steps``. ``nstlog`` and ``nstchk`` default to ``5000`` the same way.
   * - **Optimizer controls** — consumed here; never reach ``md.ini``
     -
     -
     -
     -
   * - ``ntraj``
     - int
     - no
     - ``10``
     - Independent trajectories per round. **Recommended: 10 trajectories**, which is the default. Produced as a single **multi-copy** MD run (``n_copies = ntraj``) that is then split into per-copy DCDs, so this scales the cost of every round roughly linearly.
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
   * - ``outdir``
     - str
     - no
     - ``opt_out``
     - Optimization root directory — the ``round_N/`` subdirectories and the final model are written here. A relative path is resolved **against this** ``.ini`` **file** (like ``pdb_file`` and ``domain_def``), so the file gives the same layout from any working directory. ``-o``/``--outdir`` on the command line overrides it; when neither is set, ``opt_out`` is created in the *current* directory. Not to be confused with ``output_dir``, a per-round ``md.ini`` key the optimizer :ref:`overrides every round <opt-overridden>`.
   * - **Protocol & hardware**
     -
     -
     -
     -
   * - ``device``
     - str
     - no
     - ``GPU``
     - Passed through to ``md.ini``, where the default would be ``CPU``. Optimization is a production-scale task (``ntraj`` × ``max_rounds`` trajectories). ``--device`` on the command line overrides it for every round, without editing the ``.ini``.
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

The table covers the keys you are likely to set. The remaining implicit defaults
fix the stability protocol's ensemble and the model's standard parameterization,
and are rarely touched.

**The ensemble.** Every round samples constant-temperature dynamics: a Langevin
integrator at ``ref_t = 300`` K with friction ``tau_t = 0.05`` ps⁻¹
(``tcoupl = yes``), no pressure coupling (``pcoupl = no``), and no periodic box
(``pbc = no``). Nscale is therefore calibrated at fixed particle number and
temperature — this is canonical, constant-*T* sampling rather than a strict NVT
box run, since with ``pbc = no`` the chain is unbounded and there is no
simulation volume being held fixed. The thermostat is what makes the stability
verdict meaningful: Q is scored on thermal fluctuations at 300 K, so a unit
counts as stable only if it resists unfolding at the protocol temperature.

**Parameterization and bookkeeping.** ``dt = 0.015`` (ps), ``model = topo``,
``restart = no`` (every round starts fresh), ``ppn = 4``, ``nstlog = 5000``,
``nstchk = 5000``. Setting any of these — or any ensemble key above — in
``optimize.ini`` overrides the implicit value. The full set is
:data:`topo.optimize.IMPLICIT_DEFAULTS`. Beyond these, any other key from
:doc:`simulation_control` may appear and is passed through verbatim.

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
corresponding command-line flag; otherwise your ``.ini`` value stands.

What the optimizer writes
-------------------------

Everything lands under the optimization root — ``--outdir``, else ``outdir`` in
the ``.ini``, else ``opt_out/``:

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
