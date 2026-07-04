Synthesis control options
=========================

Continuous-synthesis (co-translational folding) runs are configured from a
separate ``.ini`` control file (e.g. ``csp.ini``), read by
:func:`topo.csp.protocol.read_csp_config`, which returns a
:class:`~topo.csp.protocol.CSPConfig` (inputs, length schedule, and a populated
:class:`~topo.csp.core.RunParams`). This is the synthesis analogue of the
single-protein :doc:`simulation_control` page: the same ``[OPTIONS]`` /
``key = value`` grammar, a different set of options. The section title
``[OPTIONS]`` is required.

* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Every option below has a default **except** the ones marked *required*; you
  only need to set the options you want to change.
* Units are OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm² — and **dwell times
  are in seconds**. Integers may use ``_`` digit separators (``200_000``).

For the *physics* behind these knobs (the three-stage elongation cycle, the
codon-resolved kinetics, the rigid ribosome and tunnel wall, the stability
guard), see :doc:`continuous_synthesis`. This page is the parameter reference.

Running a synthesis
-------------------

The runner lives in the package as :mod:`topo.csp`. Once TOPO is installed
(``pip install -e .`` from the repo root) either of these works:

.. code-block:: bash

    topo-csp -f csp.ini              # installed console command
    python -m topo.csp -f csp.ini    # module form

    # stitch the per-stage trajectories into one VMD movie afterwards
    topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb

A GPU is recommended: the truncated ribosome adds several thousand rigid beads.
See :doc:`../tutorials/08_ribosome_synthesis` for a runnable, validated example.

Example ``csp.ini``:

.. code-block::

        [OPTIONS]
        ; --- inputs (pdb_file, ribosome, domain_def are required) ---
        pdb_file           = 4c5c_model_clean.pdb   ; full native PDB; the CG model is built from it
        ribosome           = ribosome_trunc.pdb     ; truncated CG ribosome (P-/A-anchors + rigid scenery)
        domain_def         = domain.yaml            ; per-domain contact-nscale definition
        stride_output_file = 4c5c_model_clean_stride.dat  ; optional; else STRIDE is run for you

        ; --- length schedule ---
        L0    = 1            ; first nascent length (default 1)
        L_max = 306          ; final nascent length (default: full residue count)

        ; --- codon-resolved kinetics ---
        mrna         = 4c5c_mrna.txt   ; one codon per residue (required for per-codon timing)
        ; codon_times = trans_times.txt  ; table path = per-codon; a positive number of s = uniform;
                                          ; omit -> bundled E. coli 310 K table
        scale_factor = 216564650       ; in-vivo seconds -> in-silico ns compression (larger = faster)
        time_stage_1 = 0.000340        ; mean peptidyl-transfer dwell (s)
        time_stage_2 = 0.004201        ; mean translocation dwell (s)
        random_seed  = 20240629        ; reproducible first-passage-time schedule

        ; --- test clamps (leave UNSET in production; see the warning below) ---
        max_steps_per_stage = 2000
        min_steps_per_stage = 100

        ; --- integrator / output ---
        dt     = 0.015       ; timestep (ps)
        ref_t  = 310         ; temperature (K)
        tau_t  = 0.05        ; Langevin friction (1/ps)
        nstout = 100         ; trajectory/log output interval (steps)

        ; --- ribosome / PTC mechanics ---
        optimize_ptc_geometry = yes       ; seed at the equilibrium peptide-bond geometry ...
        constraints           = AllBonds  ; ... and pair with rigid bonds (stable 15 fs path)
        restraint_k = 83680  ; C-terminus harmonic restraint (kJ/mol/nm^2)
        buffer      = 0.4    ; A-site seed offset into the tunnel (nm)
        minimize    = yes    ; energy-minimize each seeded structure before its MD
        ; tunnel_wall = yes  ; one-sided tunnel floor (default on; plane auto-derived)

        ; --- post-synthesis phases (steps; 0 = skip) ---
        ejection_steps     = 20000   ; release the tether; chain diffuses out (+x)
        dissociation_steps = 0       ; free run away from the ribosome

        ; --- hardware / output ---
        device = GPU
        ppn    = 4
        outdir = synth_out


Parameter summary
+++++++++++++++++

"Required = yes" means the run cannot proceed without it. Options with a default
may be omitted. ``—`` in the *Default* column means there is no default (the
option is required, or only meaningful in a specific mode noted in the
description). The options are grouped by role; they all live in the same single
``[OPTIONS]`` section.

Inputs & length schedule
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 14 12 16 38

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``pdb_file``
     - str
     - **yes**
     - ``—``
     - Full native (all-atom) PDB of the target protein. TOPO builds the one-bead-per-residue CG model, force field, and native-contact map from it. The nascent chain at length ``L`` uses residues ``1..L`` of this structure.
   * - ``ribosome``
     - str
     - **yes**
     - ``—``
     - Truncated CG ribosome PDB. Source of the P-/A-anchors and the rigid (mass-0) scenery. Supplying it *is* the signal to load it as rigid — there is no ``rigid_ribosome`` key. Must carry the tRNA beads under the expected names (segids ``PtR``/``AtR``, resid 76).
   * - ``domain_def``
     - str
     - **yes**
     - ``—``
     - Domain YAML defining the protein's per-domain / per-interface contact-nscale (Gō well-depth scaling) — the structure-based analogue of O'Brien's ``nscal``. See :doc:`domain_definition`.
   * - ``stride_output_file``
     - str
     - no
     - ``—``
     - Precomputed STRIDE output for the contact build. If omitted, STRIDE is run automatically on the structure (STRIDE must be on ``PATH``).
   * - ``mrna``
     - str
     - for per-codon timing
     - ``—``
     - mRNA sequence file (raw nucleotides; one codon per residue plus one stop). Drives the per-codon timing, so it is required unless ``codon_times`` is set to a number (uniform timing).
   * - ``codon_times``
     - str or float
     - no
     - bundled E. coli 310 K
     - Codon-timing key. A **path** to a codon→mean-time table (``CODON  seconds``) selects **per-codon** timing; a **positive number of seconds** selects **uniform** timing (that mean dwell for every codon, no ``mrna`` needed). Omit it to use the bundled Fluitt *E. coli* 310 K table (organism-universal). A table filename must **not** be a bare number (a numeric value is always read as a uniform time).
   * - ``L0``
     - int
     - no
     - ``1``
     - First nascent length to synthesize (start from a single residue).
   * - ``L_max``
     - int
     - no
     - full length
     - Final nascent length. Omit/blank to synthesize the whole chain. Must satisfy ``1 <= L0 <= L_max <= N_full``.
   * - ``outdir``
     - str
     - no
     - ``synth_out``
     - Output root; each residue writes ``L_<L>/stage_<1,2,3>/`` beneath it.

Kinetics & schedule length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 12 12 14 42

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``scale_factor``
     - float
     - no
     - ``4331293``
     - In-vivo-seconds → in-silico-ns compression factor (``t_sim_ns = t_s · 1e9 / scale_factor``). Larger ⇒ fewer MD steps per residue ⇒ faster run, while preserving the *relative* timing of fast vs. slow codons.
   * - ``time_stage_1``
     - float [s]
     - no
     - ``0.00034``
     - Mean peptidyl-transfer (peptide-bond) dwell. The actual per-residue dwell is an exponential draw with this mean.
   * - ``time_stage_2``
     - float [s]
     - no
     - ``0.004201``
     - Mean translocation dwell. Stage 3's mean is the remainder: ``τ(next codon) − time_stage_1 − time_stage_2``.
   * - ``random_seed``
     - int
     - no
     - ``—`` (nondet.)
     - Seed for the first-passage-time sampler, making the whole dwell schedule reproducible. Accepts ``_`` separators.
   * - ``max_steps_per_stage``
     - int
     - no
     - ``—`` (uncapped)
     - **Testing only.** Upper clamp on each stage's MD step count. Breaks the physical timescale mapping — leave **unset** in production. See the warning below.
   * - ``min_steps_per_stage``
     - int
     - no
     - ``1``
     - **Testing only.** Lower clamp on each stage's MD step count. See the warning below.
   * - ``ejection_steps``
     - int
     - no
     - ``0``
     - Post-synthesis ejection phase length (steps); ``0`` = skip. Releases the C-terminus tether so the finished chain diffuses out of the tunnel (+x).
   * - ``dissociation_steps``
     - int
     - no
     - ``0``
     - Post-synthesis dissociation phase length (steps); ``0`` = skip. A further free run that lets the released protein drift off the ribosome.

.. warning::

   **``max_steps_per_stage`` / ``min_steps_per_stage`` are testing-only knobs.**
   They clamp the MD step count per stage so tutorials finish quickly, which
   **breaks the physical timescale mapping** — the integrated MD per stage no
   longer matches the sampled dwell time. In a production run leave them
   **unset** so step counts are driven entirely by the kinetics
   (``scale_factor``, the codon times, and ``dt``). The sampled dwell **times in
   seconds** are always written to ``dwell_times.dat`` regardless of any clamp.

Integrator, ribosome & PTC mechanics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These configure the shared per-length engine (:class:`~topo.csp.core.RunParams`,
consumed by :func:`~topo.csp.core.run_length`).

.. list-table::
   :header-rows: 1
   :widths: 22 12 16 50

   * - Option
     - Type
     - Default
     - Description
   * - ``dt``
     - float [ps]
     - ``0.015``
     - Integration timestep. A stage that diverges is transparently re-run with ``dt`` halved and the step count doubled (dwell time unchanged) — see the stability guard in :doc:`continuous_synthesis`.
   * - ``ref_t``
     - float [K]
     - ``310``
     - Langevin temperature. Defaults to **310 K** — the temperature of the bundled codon-time table — so the kinetics and the thermostat are consistent by default.
   * - ``tau_t``
     - float [ps⁻¹]
     - ``0.05``
     - Langevin friction coefficient.
   * - ``nstout``
     - int
     - ``5000``
     - Trajectory (DCD) and log output interval, in steps, for every stage.
   * - ``constraints``
     - str
     - ``None``
     - Bond treatment. ``None`` (flexible harmonic bonds) is required for the default far-seed path. Set ``AllBonds`` (rigid) **only together with** ``optimize_ptc_geometry = yes``, which seeds at the equilibrium bond length so rigid constraints stay stable.
   * - ``restraint_k``
     - float [kJ/mol/nm²]
     - ``83680``
     - Stiffness of the C-terminus harmonic restraint to the A/P target point (= 200 kcal/mol/Å²). Switching the target A→P is how translocation is reproduced.
   * - ``buffer``
     - float [nm]
     - ``0.4``
     - Offset into the tunnel (+x) at which a newly added bead is seeded from the A-anchor. Ignored when ``optimize_ptc_geometry = yes`` (which seeds at the optimal A-site target instead).
   * - ``minimize``
     - bool
     - ``yes``
     - Energy-minimize each seeded structure before running that stage's MD.
   * - ``tunnel_wall``
     - bool
     - ``yes``
     - Apply the one-sided half-harmonic tunnel wall (a floor below the synthesis point that keeps the chain extruding forward). The plane ``x₀`` is **auto-derived** from the ribosome structure and the stiffness is a fixed model constant — neither is an INI key; only this on/off toggle is exposed.
   * - ``ptc_offset``
     - float [nm]
     - auto ``0.476``
     - Offset into the tunnel (+x) from the P-anchor bead at which the C-terminus is held, so it clears the P-tRNA bead. Defaults to the tRNA-tether bond length.
   * - ``optimize_ptc_geometry``
     - bool
     - ``no``
     - Seed each new residue at the optimal A-site target — one equilibrium peptide bond (0.381 nm) from the previous C-terminus — and place the A/P restraint targets and tunnel-wall plane there too. This lets a rigid ``constraints = AllBonds`` build run cleanly at 15 fs without the dt-halving guard firing. Pair it with ``constraints = AllBonds``.
   * - ``nascent_ev_radii``
     - str
     - ``kb``
     - Nascent per-residue excluded-volume radius (Rmin/2) for the nascent-chain ↔ ribosome force: ``kb`` = per-residue Karanicolas–Brooks radii from the native structure (O'Brien's actual values); ``per_aa`` = per-amino-acid sidechain radii (fallback). Must be one of ``kb`` / ``per_aa``.
   * - ``trna_tether``
     - bool
     - ``no``
     - C-terminus restraint mode. ``no`` (default via ``csp.ini``) uses the validated position restraint to the A/P target point (a→a→p migration). ``yes`` switches to O'Brien's full tRNA tether (bond + two orienting angles + improper; A-site in stages 1–2, P-site in stage 3).
   * - ``device``
     - str
     - ``CPU``
     - Compute platform: ``CPU`` or ``GPU`` (CUDA). A GPU is strongly recommended given the rigid-ribosome bead count.
   * - ``ppn``
     - int
     - ``1``
     - Number of CPU threads. Used when ``device = CPU``.

.. note::

   Boolean options accept ``yes``/``no``, ``true``/``false``, ``1``/``0``
   (parsed by :func:`topo.utils.config.strtobool`).


Notes on individual options
+++++++++++++++++++++++++++

Required inputs (``pdb_file`` / ``ribosome`` / ``domain_def``)
    A synthesis run cannot start without all three. ``pdb_file`` supplies the
    target fold (topology, force field, and the full native contact map, which is
    computed **once** and sliced to the top-left ``L×L`` block per length —
    STRIDE and the contact analysis are never re-run per length). ``ribosome`` is
    the truncated CG ribosome that provides the P-/A-anchors and the rigid
    scenery; it must carry tRNA beads under the fixed names (segids
    ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``) or anchor lookup fails.
    ``domain_def`` sets the per-domain contact-nscale (the ``nscal`` analogue);
    see :doc:`domain_definition`.

Per-codon vs. uniform timing (``mrna`` / ``codon_times``)
    A single key, ``codon_times``, selects the timing mode by its **value type**:

    * **A path** (e.g. ``codon_times = trans_times.txt``) → **per-codon** timing:
      the ``mrna`` file is required and each codon's mean translation time is read
      from that ``CODON  seconds`` table.
    * **A positive number of seconds** (e.g. ``codon_times = 0.05``) → **uniform**
      timing: every codon gets that mean dwell and no ``mrna`` is needed — a
      control/testing mode.
    * **Omitted** → per-codon timing with the bundled *E. coli* 310 K table (the
      table is organism-universal, so only the protein-specific ``mrna`` is
      mandatory).

    Because a numeric value is always read as a uniform time, a codon-time **table
    filename must not be a bare number** (give it a name like ``trans_times.txt``).

Kinetics and step counts (``scale_factor`` / ``time_stage_1`` / ``time_stage_2``)
    Each residue is added over three sub-stages whose dwell times are exponential
    draws with means ``time_stage_1`` (peptidyl transfer), ``time_stage_2``
    (translocation), and ``τ(next codon) − time_stage_1 − time_stage_2`` (the
    decoding wait). Those seconds are mapped to MD steps through ``scale_factor``
    and ``dt`` (``t_sim_ns = t_s · 1e9 / scale_factor``; ``n_steps = t_sim_ns /
    dt_ns``). A larger ``scale_factor`` shortens the run while keeping the
    *relative* fast/slow-codon timing that matters for co-translational folding.
    ``random_seed`` fixes the whole schedule. Full derivation in
    :doc:`continuous_synthesis`.

Bond treatment and PTC geometry (``constraints`` / ``optimize_ptc_geometry`` / ``buffer``)
    Two consistent paths exist:

    * **Default (flexible bonds).** ``constraints = None`` with
      ``optimize_ptc_geometry = no``. The new bead is seeded far (``buffer`` nm
      from the A-anchor) and the harmonic bond absorbs the stretch, but stiff
      native contacts can force the dt-halving stability guard to fire.
    * **Optimized (rigid bonds).** ``optimize_ptc_geometry = yes`` together with
      ``constraints = AllBonds`` seeds each residue at the equilibrium peptide-bond
      length, so rigid constraints stay stable at 15 fs and the guard never fires.
      This is the path the ribosome-synthesis tutorial uses.

    Do **not** set ``constraints = AllBonds`` without ``optimize_ptc_geometry =
    yes`` — a rigid constraint cannot represent the far A-site seed.

C-terminus restraint (``trna_tether`` / ``restraint_k`` / ``ptc_offset``)
    The nascent C-terminus is held at the PTC and translocated A→P each residue.
    The default ``trna_tether = no`` uses a simple harmonic position restraint
    (stiffness ``restraint_k``) to the A/P target point — the validated path.
    ``trna_tether = yes`` replaces it with O'Brien's full tRNA tether (bond + two
    orienting angles + improper), which additionally aims the chain down the
    tunnel. ``ptc_offset`` sets how far past the anchor bead the hold point sits
    (so the C-terminus clears the tRNA bead); it defaults to the tether bond
    length (``0.476`` nm). (Note: the underlying ``RunParams`` dataclass field
    defaults ``trna_tether`` to ``True``, but ``read_csp_config`` makes the
    effective ``csp.ini`` default ``no``.)

Tunnel wall (``tunnel_wall``)
    A one-sided half-harmonic wall supplies the "floor" that the truncated 50S
    body no longer provides, so the chain can only extrude forward (+x) and cannot
    slip back into the void below the synthesis point. The wall **plane** is
    auto-placed at the lower C-terminus hold plane (recomputed for whatever
    ribosome PDB you supply, so it can never go stale) and its stiffness is a
    fixed model constant — hence only the ``tunnel_wall`` on/off toggle is a knob.
    Leave it on for any truncated ribosome.

Post-synthesis phases (``ejection_steps`` / ``dissociation_steps``)
    After the last residue, ``ejection_steps > 0`` runs a phase with the
    C-terminus tether **released** (rigid ribosome and tunnel wall still present),
    so the finished chain diffuses out of the tunnel — the model's analogue of
    termination. ``dissociation_steps > 0`` then continues the free protein away
    from the ribosome. Both write their own ``ejection/`` and ``dissociation/``
    output folders; ``0`` skips the phase.

Output layout and progress log
    Every stage writes a standalone, **nascent-only** trajectory to
    ``<outdir>/L_<L>/stage_<s>/`` (``traj.dcd``, ``traj_final.pdb``, ``traj.log``,
    …); different lengths have different bead counts, so they are separate files.
    A per-residue ``<outdir>/dwell_times.dat`` records the codon and the three
    sampled dwell times (seconds, ns, and MD steps) — the physical schedule,
    independent of any step clamp. Stitch the per-stage trajectories into one
    VMD movie with ``topo-csp-movie -o <outdir> --ribosome <ribosome>.pdb``. Set
    ``TOPO_CSP_VERBOSE=1`` to restore the full per-stage banners. See
    :doc:`continuous_synthesis` for the output tree in full.
