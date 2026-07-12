Synthesis control options
=========================

Both synthesis runners are configured from a ``.ini`` control file with a single
``[OPTIONS]`` section:

* **``topo-csp``** — the coarse-grained-ribosome runner (``csp.ini``), read by
  :func:`topo.csp.protocol.read_csp_config`;
* **``topo-cylinder``** — the analytic-tunnel runner (``cylinder.ini``), read by
  :func:`topo.csp.cylinder.read_cylinder_config`.

Both return a populated :class:`~topo.csp.core.RunParams`, so **most keys are shared**.
**This page is the single reference for every control key** — grouped into *shared*,
*coarse-grained-ribosome-only*, and *cylinder-only* sets below; the
:doc:`continuous_synthesis` and :doc:`cylinder_synthesis` pages cover the *physics* and
point back here for the keys.

* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Every option has a default **except** the ones marked *required*; set only what you
  want to change.
* Units are OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm² — and **dwell times are in
  seconds**. Integers may use ``_`` digit separators (``200_000``).

Running a synthesis
-------------------

Once TOPO is installed (``pip install -e .`` from the repo root):

.. code-block:: bash

    topo-csp      -f csp.ini         # coarse-grained-ribosome runner
    topo-cylinder -f cylinder.ini    # analytic-tunnel runner

    # stitch the per-stage trajectories into one VMD movie afterwards
    topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb

A GPU is recommended for ``topo-csp``: the truncated ribosome adds several thousand rigid
beads. See :doc:`../tutorials/08_ribosome_synthesis` for a runnable, validated example.

Example control files
---------------------

``csp.ini`` (coarse-grained ribosome):

.. code-block::

        [OPTIONS]
        ; --- inputs (pdb_file, ribosome, domain_def are required) ---
        pdb_file           = 4c5c_model_clean.pdb   ; all-atom native PDB; the CG model is built from it
        ribosome           = ribosome_trunc.pdb     ; truncated CG ribosome (P-/A-anchors + rigid scenery)
        domain_def         = domain.yaml            ; per-domain contact-nscale definition
        stride_output_file = 4c5c_model_clean_stride.dat  ; optional; else STRIDE is run for you

        ; --- length schedule ---
        L0    = 1            ; first nascent length (default 1)
        L_max = 306          ; final nascent length (default: full residue count)

        ; --- codon-resolved kinetics ---
        mrna         = 4c5c_mrna.txt   ; one codon per residue; or fastest/slowest/median
        codon_times  = trans_times.txt ; table path = per-codon; a number of s = uniform
        scale_factor = 216564650       ; in-vivo s -> in-silico ns compression (larger = faster)
        time_stage_1 = 0.000340        ; mean peptidyl-transfer dwell (s)  [CSP only]
        time_stage_2 = 0.004201        ; mean translocation dwell (s)      [CSP only]
        random_seed  = 20240629        ; reproducible first-passage-time schedule

        ; --- integrator / ribosome mechanics ---
        dt     = 0.015 ; ref_t = 310 ; tau_t = 0.05 ; nstout = 100
        constraints = AllBonds ; minimize = yes
        ; tunnel_wall = yes                              ; one-sided tunnel floor (default on)
        ; ribo_free_mask = L24 : 42 - 59                 ; free the L24 exit-tunnel loop (optional)
        ; ribo_free_pdb  = .../ecoli/L24_atomistic.pdb   ; all-atom chain for native contacts

        ; --- post-synthesis phases + resume + hardware ---
        ejection_steps = 20000 ; dissociation_steps = 0 ; resume = auto
        device = GPU ; ppn = 4 ; outdir = synth_out

``cylinder.ini`` (analytic tunnel) — the shared keys are identical; it drops ``ribosome``
/ the tRNA-tether / flexible-loop keys and adds the ``tunnel_*`` geometry:

.. code-block::

        [OPTIONS]
        pdb_file   = P0CX28_clean.pdb ; domain_def = domain.yaml   ; (no `ribosome` PDB)
        L0 = 1 ; L_max =              ; mrna = P0CX28_mrna.txt ; codon_times = trans_times.txt
        scale_factor = 216564650 ; random_seed = 20240629
        constraints = None ; minimize = yes  ; (cylinder allows None or AllBonds)
        dt = 0.015 ; ref_t = 300 ; tau_t = 0.05 ; nstout = 100

        ; --- analytic exit tunnel (cylinder only) ---
        tunnel_radius = 0.9 ; tunnel_length = 10.0 ; tunnel_x_lo = 0.0
        tunnel_center = 0.0, 0.0 ; tunnel_k = 8368

        ejection_steps = 300000 ; dissociation_steps = 0 ; resume = auto
        device = GPU ; ppn = 4 ; outdir = synth_out


Option reference
++++++++++++++++

"Required = yes" means the run cannot proceed without it; options with a default may be
omitted. ``—`` in the *Default* column means there is no default. The keys are grouped
into three sets — all in the same single ``[OPTIONS]`` section:

* **Shared** — accepted (with identical meaning) by *both* runners.
* **Coarse-grained ribosome runner only** (``topo-csp`` / ``csp.ini``).
* **Cylinder runner only** (``topo-cylinder`` / ``cylinder.ini``).

Shared options (both runners)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Inputs & length schedule.* Both runners build a length-``L`` structure-based (Gō) model
from ``pdb_file`` + ``domain_def`` (the native contact map is computed once and sliced
per length — STRIDE / contact analysis are never re-run).

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
     - Full native **all-atom** PDB of the target protein — TOPO's Gō model derives secondary structure (STRIDE) and the native-contact map from the heavy atoms, so a Cα-only CG PDB is *not* sufficient (the sibling cosmo IDP model does accept one). The nascent chain at length ``L`` uses residues ``1..L``.
   * - ``domain_def``
     - str
     - **yes**
     - ``—``
     - Domain YAML defining the per-domain / per-interface contact-nscale (Gō well-depth scaling) — the structure-based analogue of O'Brien's ``nscal``. See :doc:`domain_definition`.
   * - ``stride_output_file``
     - str
     - no
     - ``—``
     - Precomputed STRIDE output for the contact build. If omitted, STRIDE is run automatically (must be on ``PATH``).
   * - ``L0``
     - int
     - csp: no / cyl: **yes**
     - ``1``
     - First nascent length. Optional for ``topo-csp`` (default 1); **required** for ``topo-cylinder``.
   * - ``L_max``
     - int
     - no
     - full length
     - Final nascent length. Omit/blank to synthesize the whole chain. Must satisfy ``1 ≤ L0 ≤ L_max ≤ N_full``.
   * - ``outdir``
     - str
     - no
     - ``synth_out``
     - Output root; each residue writes one ``L_<L>/`` folder (CSP: shared ``traj.psf`` + per-stage ``traj_s<1,2,3>.dcd``; cylinder: a single ``traj.dcd``).

*Codon kinetics & schedule length.*

.. list-table::
   :header-rows: 1
   :widths: 20 12 12 14 42

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``mrna``
     - str
     - for per-codon timing
     - ``—``
     - mRNA sequence file (one codon per residue + one stop), **or** ``fastest``/``slowest``/``median`` to auto-build a synonymous-codon mRNA (written next to the PDB). Required unless ``codon_times`` is a number. A real filename must not be ``fastest``/``slowest``/``median``. See the notes below.
   * - ``codon_times``
     - str or float
     - for per-codon timing
     - ``—``
     - A **path** to a ``CODON  seconds`` table selects **per-codon** timing (no bundled default — pick one under ``assets/csp/codon_dwell_times/``); a **positive number of seconds** selects **uniform** timing (no ``mrna`` needed). A table filename must not be a bare number.
   * - ``scale_factor``
     - float
     - no
     - ``4331293``
     - In-vivo-seconds → in-silico-ns compression (``t_sim_ns = t_s · 1e9 / scale_factor``). Larger ⇒ fewer MD steps per residue ⇒ faster, preserving the *relative* fast/slow-codon timing.
   * - ``random_seed``
     - int
     - no
     - ``—`` (nondet.)
     - Seed for the first-passage-time sampler — makes the whole dwell schedule reproducible.
   * - ``max_steps_per_stage``
     - int
     - no
     - ``—`` (uncapped)
     - **Testing only.** Upper clamp on each stage's MD step count. See the warning below.
   * - ``min_steps_per_stage``
     - int
     - no
     - ``1``
     - **Testing only.** Lower clamp on each stage's MD step count.

.. warning::

   **``max_steps_per_stage`` / ``min_steps_per_stage`` are testing-only knobs.** They
   clamp the MD step count so tutorials finish quickly, which **breaks the physical
   timescale mapping** — the integrated MD no longer matches the sampled dwell time. In
   production leave them **unset** so step counts come entirely from the kinetics
   (``scale_factor``, the codon times, ``dt``). The sampled dwell **times in seconds** are
   always written to ``dwell_times.dat`` regardless of any clamp.

*Integrator & MD mechanics* (the shared per-length engine,
:class:`~topo.csp.core.RunParams`).

.. list-table::
   :header-rows: 1
   :widths: 22 14 16 48

   * - Option
     - Type
     - Default
     - Description
   * - ``dt``
     - float [ps]
     - ``0.015``
     - Integration timestep. A diverging stage is transparently re-run with ``dt`` halved and the step count doubled (dwell unchanged) — the stability guard in :doc:`continuous_synthesis`.
   * - ``ref_t``
     - float [K]
     - ``310``
     - Langevin temperature. Defaults to 310 K (the E. coli codon-time table temperature) so kinetics and thermostat are consistent; set it to match your ``codon_times`` table.
   * - ``tau_t``
     - float [ps⁻¹]
     - ``0.05``
     - Langevin friction coefficient.
   * - ``nstout``
     - int
     - ``5000``
     - Trajectory (DCD) and log output interval, in steps.
   * - ``constraints``
     - str
     - ``AllBonds``
     - Bond treatment. Rigid ``AllBonds`` is the default and is stable because each new residue is seeded at the equilibrium peptide-bond length. Set ``None`` for flexible harmonic bonds.
   * - ``restraint_k``
     - float [kJ/mol/nm²]
     - ``83680``
     - Stiffness of the C-terminus harmonic restraint to the PTC target (= 200 kcal/mol/Å²).
   * - ``minimize``
     - bool
     - ``yes``
     - Energy-minimize each seeded structure before running that stage's MD.
   * - ``device``
     - str
     - ``CPU``
     - Compute platform: ``CPU`` or ``GPU`` (CUDA). GPU recommended for ``topo-csp``.
   * - ``ppn``
     - int
     - ``1``
     - Number of CPU threads (``device = CPU``).

*Post-synthesis phases & resume.*

.. list-table::
   :header-rows: 1
   :widths: 22 10 12 56

   * - Option
     - Type
     - Default
     - Description
   * - ``ejection_steps``
     - int
     - ``0``
     - Post-synthesis ejection phase length (steps); ``0`` = skip. Releases the C-terminus restraint so the finished chain diffuses out of the tunnel (+x).
   * - ``dissociation_steps``
     - int
     - ``0``
     - Post-synthesis dissociation phase length (steps); ``0`` = skip. A further free run that lets the released protein drift off.
   * - ``resume``
     - str
     - ``auto``
     - Resume policy for interrupted runs: ``auto`` (resume iff a ``progress.log`` is present under ``outdir``, else fresh), ``yes`` (require a resumable run, else error), ``no`` (always fresh). CLI ``--fresh`` / ``--no-resume`` forces ``no``. See :doc:`synthesis_resume`.

Coarse-grained ribosome runner only (``topo-csp``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These keys are read only from ``csp.ini``. (``time_stage_1`` / ``time_stage_2`` are
*accepted* by the cylinder INI but have **no effect** there — with a single MD segment
per residue the whole codon dwell is one segment, not a three-way split.)

.. list-table::
   :header-rows: 1
   :widths: 22 12 14 52

   * - Option
     - Type
     - Default
     - Description
   * - ``ribosome``
     - str
     - **required**
     - Truncated CG ribosome PDB — the P-/A-anchors and the rigid (mass-0) scenery. Supplying it *is* the signal to load it as rigid (no ``rigid_ribosome`` key). Must carry tRNA beads under the fixed names (segids ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``).
   * - ``time_stage_1``
     - float [s]
     - ``0.00034``
     - Mean peptidyl-transfer (stage 1) dwell; the actual dwell is an exponential draw with this mean.
   * - ``time_stage_2``
     - float [s]
     - ``0.004201``
     - Mean translocation (stage 2) dwell. Stage 3's mean is the remainder ``τ(next codon) − time_stage_1 − time_stage_2``.
   * - ``tunnel_wall``
     - bool
     - ``yes``
     - Apply the one-sided half-harmonic tunnel wall (a floor below the synthesis point). The plane ``x₀`` is **auto-derived** from the ribosome structure and the stiffness is a fixed constant — neither is a key; only this on/off toggle is exposed.
   * - ``nascent_ev_radii``
     - str
     - ``kb``
     - Nascent per-residue excluded-volume radius (Rmin/2) for the nascent ↔ ribosome force: ``kb`` = per-residue Karanicolas–Brooks radii from the native structure; ``per_aa`` = per-amino-acid sidechain radii (fallback).
   * - ``trna_tether``
     - bool
     - ``no``
     - How the growing C-terminus is held at the PTC. ``no`` (default) pins the C-terminal bead to a **fixed point in space** (the A-/P-site target, harmonic, stiffness ``restraint_k``) — position only, no orientation. ``yes`` instead couples it to the actual tRNA beads with O'Brien's full peptidyl-tRNA linkage (bond + 2 orienting angles + improper + a backbone angle), which controls **position *and* orientation** so the chain threads out N-first instead of balling up at the PTC. See the note below.
   * - ``ribo_free_mask``
     - str
     - *(none)*
     - Free a portion of one ribosomal protein as a mobile Gō loop, ``SEG : lo - hi`` (E. coli ``L24``, eukaryotes ``L26``). Requires ``ribo_free_pdb``. Absent → whole ribosome frozen. See :ref:`the flexible-loop notes <flexible-loop>`.
   * - ``ribo_free_pdb``
     - path
     - *(none)*
     - All-atom PDB of the ``ribo_free_mask`` chain — the native-contact source the Cα-only ribosome cannot supply. Required whenever ``ribo_free_mask`` is set.

Cylinder runner only (``topo-cylinder``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The analytic exit tunnel — a bore of radius ``r`` drilled through an infinite wall — that
replaces the explicit ribosome beads. Read only from ``cylinder.ini``.

.. list-table::
   :header-rows: 1
   :widths: 24 14 62

   * - Option
     - Default
     - Description
   * - ``tunnel_radius``
     - ``0.9``
     - Bore radius ``r`` (nm); ~3 CG beads wide.
   * - ``tunnel_length``
     - ``10.0``
     - Bore length (nm). The exit face sits at ``x_exit = tunnel_x_lo + tunnel_length`` (derived, not a key).
   * - ``tunnel_x_lo``
     - ``0.0``
     - PTC / closed end of the bore (nm); the C-terminus is restrained on-axis here.
   * - ``tunnel_center``
     - ``0.0, 0.0``
     - Tunnel axis ``(y0, z0)`` (nm). The axis runs along +x.
   * - ``tunnel_k``
     - ``8368``
     - Wall stiffness (kJ/mol/nm² = 20 kcal/mol/Å²).
   * - ``tunnel_mouth_round``
     - ``0.2``
     - Mouth-corner fillet radius ``rho`` (nm); rounds the 90° inner corner so the potential is continuous.

.. note::

   Boolean options accept ``yes``/``no``, ``true``/``false``, ``1``/``0``
   (parsed by :func:`topo.utils.config.strtobool`).


Notes on individual options
+++++++++++++++++++++++++++

Required inputs (``pdb_file`` / ``ribosome`` / ``domain_def``)
    ``pdb_file`` supplies the target fold (topology, force field, and the full native
    contact map, computed **once** and sliced to the top-left ``L×L`` block per length).
    It **must be all-atom** — STRIDE (secondary structure) and the native-contact map are
    derived from the heavy atoms; a Cα-only CG PDB is not sufficient (the sibling cosmo
    IDP runner does accept one). ``ribosome`` (``topo-csp`` only) provides the P-/A-anchors
    and rigid scenery; it must carry tRNA beads under the fixed names (segids
    ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``) or anchor lookup fails.
    ``domain_def`` sets the per-domain contact-nscale (the ``nscal`` analogue); see
    :doc:`domain_definition`.

Per-codon vs. uniform timing (``mrna`` / ``codon_times``)
    A single key, ``codon_times``, selects the timing mode by its **value type**:

    * **A path** (e.g. ``codon_times = trans_times.txt``) → **per-codon** timing: ``mrna``
      is required and each codon's mean time is read from that ``CODON  seconds`` table.
    * **A positive number of seconds** (e.g. ``codon_times = 0.05``) → **uniform** timing:
      every codon gets that mean dwell and no ``mrna`` is needed (a control/testing mode).

    There is no bundled default table: per-codon timing needs an explicit path (pick one
    under ``assets/csp/codon_dwell_times/``). Because a numeric value is always read as a
    uniform time, a codon-time **table filename must not be a bare number**.

    **Fastest / slowest / median mRNA.** Setting ``mrna = fastest``, ``slowest`` or
    ``median`` (instead of a file) auto-builds a synonymous-codon mRNA that holds the
    protein sequence fixed and reassigns each residue's codon — a controlled comparison in
    which only the elongation timing changes. The sequence is read from ``pdb_file`` and
    each residue is encoded by its shortest-``τ`` (``fastest``), longest-``τ``
    (``slowest``) or middle-``τ`` (``median``) synonymous codon per the ``codon_times``
    table, written next to the PDB as ``mrna_<mode>.txt``. This needs a ``codon_times``
    **table** (a uniform numeric value is rejected). Pre-generate one standalone with
    ``topo-make-mrna --pdb P.pdb --codon-times T.txt --mode fastest``.

Kinetics and step counts (``scale_factor`` / ``time_stage_1`` / ``time_stage_2``)
    In ``topo-csp`` each residue is added over three sub-stages whose dwell times are
    exponential draws with means ``time_stage_1`` (peptidyl transfer), ``time_stage_2``
    (translocation), and ``τ(next codon) − time_stage_1 − time_stage_2`` (the decoding
    wait). Those seconds map to MD steps through ``scale_factor`` and ``dt``. ``random_seed``
    fixes the whole schedule. **``topo-cylinder`` uses a single segment per residue**, so
    the whole codon dwell ``τ`` is one MD segment and ``time_stage_1`` / ``time_stage_2``
    are ignored. Full derivation in :doc:`continuous_synthesis`.

Bond treatment and PTC geometry (``constraints``)
    The PTC-geometry optimization is **always on**: each new residue is seeded one
    equilibrium peptide bond (0.381 nm) from the previous C-terminus, clear of the
    excluded volume, so rigid ``constraints = AllBonds`` (the CSP default) seeds and
    minimizes cleanly without the dt-halving guard firing. The cylinder runner seeds each
    new residue one equilibrium bond deeper than its restraint point, which likewise keeps
    ``AllBonds`` stable (its example uses ``None`` historically; either works). Set
    ``None`` for flexible harmonic bonds.

C-terminus restraint (``trna_tether`` / ``restraint_k``)
    Each new residue is added at the PTC and translocated A→P over the three stages
    (A-site in stages 1–2, P-site in stage 3). ``trna_tether`` selects **how** the
    growing C-terminus is held there — the two modes differ in whether they control the
    chain's *orientation* as well as its *position*:

    * ``trna_tether = no`` **(default) — position restraint to a fixed point.** The
      C-terminal bead is pulled by a single harmonic spring (stiffness ``restraint_k``,
      default 83680 kJ/mol/nm² = 200 kcal/mol/Å²) toward a **fixed point in space** — the
      A-site or P-site target computed once per residue by ``optimal_ptc_targets``
      (the two points sit one equilibrium peptide bond, 0.381 nm, apart and clear of the
      ribosome). This constrains only *where* the C-terminus sits, not which way the
      residue faces; it is simpler, needs no tRNA orientation data, and is the right
      choice for routine runs.

    * ``trna_tether = yes`` **— full O'Brien peptidyl-tRNA linkage.** Instead of a fixed
      point, the C-terminus is coupled to the **actual tRNA beads** of the loaded
      ribosome, reproducing ``continuous_synthesis_v6.py``. Four terms act on the
      restrained bead ``N`` and the current-site tRNA beads (``R`` = acceptor, ``P`` =
      R−1, ``PU2`` = R+2):

      - a **bond** ``N–R`` (sets the tRNA–residue distance);
      - two **orienting angles** ``N–R–P`` and ``N–R–PU2`` (fix the residue's bearing in
        the tRNA frame);
      - an **improper** ``N–R–P–PU2`` (fixes the out-of-plane sense);
      - a **backbone orienting angle** ``prev–N–R`` (aims the chain down the tunnel).

      Together these control **position *and* orientation**, so the nascent chain
      extrudes N-terminus-first down the exit tunnel rather than curling up at the PTC —
      the behaviour to pick when faithfully reproducing O'Brien's directional mechanism.
      It requires the tRNA beads to be present in the ``ribosome`` PDB (segids
      ``PtR``/``AtR``, resid 76, beads ``R``/``P``/``BR2``). The always-present peptide
      bond ``N-1↔N`` stays in place either way (topo keeps the chain permanently bonded —
      a deliberate deviation, see :doc:`continuous_synthesis`).

    ``restraint_k`` sets the spring stiffness of the **position restraint only** (the
    ``no`` mode). The tether-mode ``N–R`` bond uses a fixed internal constant of the same
    magnitude (200 kcal/mol/Å²) and does **not** read ``restraint_k``, so changing that
    key has no effect under ``trna_tether = yes``. (``topo-cylinder`` has no tRNA tether —
    it always restrains the C-terminus to the on-axis PTC point.)

Tunnel wall vs. analytic tunnel (``tunnel_wall`` / ``tunnel_*``)
    ``topo-csp`` uses a one-sided half-harmonic **wall** (``tunnel_wall``) that supplies
    the "floor" the truncated 50S no longer provides; its plane is auto-placed and its
    stiffness fixed, so only the on/off toggle is exposed. ``topo-cylinder`` instead uses a
    fully **analytic tunnel** — a cylindrical bore in an infinite wall — whose geometry is
    set by the ``tunnel_radius`` / ``tunnel_length`` / ``tunnel_x_lo`` / ``tunnel_center``
    / ``tunnel_k`` / ``tunnel_mouth_round`` keys (see :doc:`cylinder_synthesis`).

.. _flexible-loop:

Flexible exit-tunnel loop (``ribo_free_mask`` / ``ribo_free_pdb``)
    (``topo-csp`` only.) By default every ribosome bead is frozen. ``ribo_free_mask`` frees
    a portion of the exit-tunnel protein (uL24 family: E. coli **L24**, eukaryotic **L26**)
    as a mobile, topo-style Gō loop — reproducing O'Brien's mobile-L24 setup, where the
    β-hairpin that lines the tunnel can breathe and contact the nascent chain while the
    rest of the ribosome stays rigid. Its bonded terms and native contacts are built from
    an **all-atom** copy of the chain (``ribo_free_pdb``) — the Cα-only ribosome cannot
    supply the STRIDE H-bonds and heavy-atom contacts. The freed residues get mass; the two
    boundary residues just outside the range stay frozen and anchor the loop. Requires
    ``minimize = yes``. Prepare the all-atom chain with
    ``assets/csp/prepare_ribosome/helpers/carve_flexible_protein.py``; ready-carved chains
    and per-organism mask ranges are in ``assets/csp/prepare_ribosome/MANIFEST.md``.
    Example::

        ribo_free_mask = L24 : 42 - 59
        ribo_free_pdb  = assets/csp/prepare_ribosome/structures/ecoli/L24_atomistic.pdb

Post-synthesis phases (``ejection_steps`` / ``dissociation_steps``)
    After the last residue, ``ejection_steps > 0`` runs a phase with the C-terminus
    restraint **released** (ribosome/tunnel still present), so the finished chain diffuses
    out — the analogue of termination. ``dissociation_steps > 0`` then continues the free
    protein away from the ribosome. Both write their own ``ejection/`` / ``dissociation/``
    folders; ``0`` skips the phase.

Output layout, resume and movies
    Each residue writes one ``<outdir>/L_<L>/`` folder (CSP: shared ``traj.psf`` +
    per-stage ``traj_s<1,2,3>.dcd``, a folded ``traj_runinfo.log`` and a single
    ``traj_final.pdb``; cylinder: a single ``traj.dcd`` + ``traj_final.pdb``), plus a
    per-residue ``dwell_times.dat`` (the schedule, with a ``#PTC`` header for CSP) and an
    append-only ``progress.log``. Re-invoking a runner on an interrupted ``outdir``
    **resumes** from the last completed residue (``resume = auto``); see
    :doc:`synthesis_resume`. Stitch a movie with ``topo-csp-movie -o <outdir>``; set
    ``TOPO_CSP_VERBOSE=1`` to restore the full per-stage banners.
