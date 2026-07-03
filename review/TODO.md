# Consolidated TODO — topo source (collected 2026-07-01)

Every open item found by sweeping the source tree: inline code markers
(`TODO`/`FIXME`/`XXX`/`HACK`), dedicated TODO/TASK files, and "open items / remaining /
not yet" sections in docs. Grouped by source, deduplicated where the same item appears in
several places (cross-refs noted). Branch `tut15-claude-fix`, HEAD `e9667c1`.

Legend: `[ ]` open · `[~]` partial/mapped · `[x]` done (kept for context).

---

## A. Inline code markers (`topo/` package)

The package is nearly marker-free — only these:

- [ ] [`topo/core/system.py:189`](../topo/core/system.py#L189) — check that residue indices
  are consecutive (validation not implemented).
- [ ] [`topo/core/system.py:309`](../topo/core/system.py#L309) — bond length is
  protein-only; handle Protein+DNA+RNA complexes.
- [ ] [`topo/core/system.py:812`](../topo/core/system.py#L812) — `threshold=0.5` is only
  "safe" for protein systems; revisit for other systems.
- Pointers (not tasks): [`topo/csp/kinetics.py:24`](../topo/csp/kinetics.py#L24) and
  [`topo/csp/protocol.py:119`](../topo/csp/protocol.py#L119) both refer to the CSP TODO
  items, now consolidated here in §B (the former `topo/csp/TODO.md` was merged in and
  deleted 2026-07-02).

### Default value fixes
- [x] **Change the default `tau_t` (Langevin friction) from `0.01` → `0.05` ps⁻¹** —
  *Done.* 0.05 ps⁻¹ is the production value; 0.01 was only ever a placeholder. Updated
  all three hard-coded defaults — [`topo/utils/config.py:308`](../topo/utils/config.py#L308)
  (runner fallback), [`topo/optimize/optimize.py:123`](../topo/optimize/optimize.py#L123)
  (nscale-optimizer per-round default), and [`topo/csp/core.py:827`](../topo/csp/core.py#L827)
  (`RunParams` dataclass default) — plus every tutorial INI / demo script pinned to `0.01`
  and all docs that cited `0.01` as the production default (relaxation time `1/tau_t`
  restated 100 ps → 20 ps). Friction is thermodynamically neutral (does not shift Tm /
  folded populations / the nscale calibration); it changes only kinetics and per-round
  equilibration speed.

---

## B. CSP — validation, production, extensions

*(Absorbed the former `topo/csp/TODO.md`, deleted 2026-07-02; revised against the current
codebase — the standalone `elongate.py` runner + Tutorial 7 are gone, `run_continuous_synthesis`
is the only synthesis path.)*

### Validation & production
- [ ] **Validate the C-terminus tether** *(marked highest priority)* — establish whether
  the tether (bond + CA–CA–tRNA angle) improves extrusion vs. a plain position restraint:
  ≥5–10 replicas per `trna_tether` on/off, production-relevant dwell (`n_steps` ≳ 50k,
  ideally 840k) with truncated `L_max`, robust per-residue radial/x metric averaged over
  frames *and* replicas. Keep tether default on only if it measurably helps.
- [ ] **Realistic production parameters** — define/validate a production INI: dwell
  `n_steps = 840_000` (12.6 ns @ 15 fs), full `L_max`, `device = GPU`, many independent
  trajectories (O'Brien runs 50/protein — add a replica driver or document launching),
  benchmark wall-clock per length for the ~4,600-particle v2 system.
- [ ] **Analysis layer** (DESIGN §6 phase 4 — none exists) — Q-vs-length co-translational
  folding curves per domain (reuse `topo.analysis.native_contacts`); ejection-time
  observable (steps from tether release to tunnel exit); wire a worked example into the
  tutorial.
- [ ] **Full-length / threading validation run** — one full `L0 → N_full` (P0CX28) at a
  realistic dwell confirming the chain threads the tunnel and the N-terminal domain folds
  outside. (Can combine with the two items above.)
- [ ] **Tunnel wall → 5.8 nm in ejection phase** (short note item in the file).

### Model extensions
- [ ] **tRNA presence / naming robustness** — the P-/A-anchors and `optimal_ptc_targets`
  (used when `optimize_ptc_geometry = yes`) assume the ribosome PDB carries tRNA beads
  under hardcoded names (segids `PtR`/`AtR`, resid 76, beads `R`/`P`/`BR2`); see
  [`topo/csp/protocol.py`](../topo/csp/protocol.py) (anchor block) and
  [`topo/csp/core.py:optimal_ptc_targets`](../topo/csp/core.py#L128). A ribosome with no
  tRNA — or with differently-named tRNA segments — fails with a generic "expected exactly
  one bead" `ValueError` from `anchor_coord`. Revise: (a) detect this up front and raise
  an actionable error naming the expected segids/resid/beads, and/or (b) make the tRNA
  segids/resid/bead names configurable (INI keys). `optimize_ptc_geometry` in particular
  depends on tRNA presence+naming, so it must be gated/validated alongside this.
- [x] **Variable per-codon elongation schedule** — *Obsolete.* The 3-stage CSP kinetics
  (`kinetics.stage_steps`, driven by `run_continuous_synthesis`) already samples a
  per-residue dwell from the codon decoding time; the single-`n_steps` `elongate.py` runner
  this item targeted was removed (2026-06-30).
- [ ] **Uniform translation** — support a uniform (constant per-residue dwell) elongation
  mode as an alternative to the codon-dependent variable schedule (the 3-stage kinetics
  above). Scope/spec TBD.
- [ ] **Restart / resume across lengths** (DESIGN §4) — skip lengths whose
  `L_<L>/traj_final.pdb` exists; resume an interrupted length from checkpoint. *(Dup of
  §D "restart=1".)*
- [ ] **`ejection_steps = auto`** — stop ejection once the C-terminus clears the ribosome
  (`x(Cterm) − max(x_ribosome) > cutoff`, default 2.0 nm), chunked-stepping loop + safety
  cap; add `ejection_cutoff_nm`, `ejection_check_every`; consider same for
  `dissociation_steps`. *(Overlaps §D "quantitative validation / ejection".)*

### Revision list
- [ ] **Remove per-stage step clamps for production** — `max_steps_per_stage` /
  `min_steps_per_stage` are testing-only and break the physical timescale mapping. Decide:
  drop the `RunParams` fields + `read_csp_config` parsing + `stage_steps` clamp args, or
  keep them but loudly warn when set. *(Also flagged inline in tutorials 12/13/14 `csp.ini`
  "delete both for production" — see §E.)*
- [ ] **Ribosome-traffic correction** (hidden since 2026-06-30) — decide: (a) finish it
  (vendor/reimplement the `ribosome_traffic` upstream-queue correction, validate the
  intrinsic→real stretch, re-expose in docs+config) or (b) drop it entirely (remove
  `RunParams` fields, `read_csp_config` parsing, `kinetics.ribosome_traffic_times`). Also
  unhide/remove the `ribosome_traffic=off` runtime banner in `protocol.py`. *(Dup of §D
  "external ribosome_traffic binary".)*
- [x] **Restore rigid `AllBonds`** — *Done for the surviving 3-stage runner.*
  `optimize_ptc_geometry = yes` places the new bead one CG bond (0.381 nm) from its
  restraint target (`optimal_ptc_targets`), so rigid `constraints = AllBonds` seeds at
  equilibrium (the exact fix this item proposed). The plain `elongate.py` runner it also
  named was removed. Whether to make it the *default* (vs. the validated flexible-bond +
  dt-halving Tutorials 12/13 path) is a design choice, not a TODO — see the step-2 doc
  TODO in §E for retiring the dt-halving guard, and DIFFERENCES §"Chain chemistry".

---

## C. `CHANGELOG.md §8` — standalone-ribosome migration open items

- [ ] **Commit the working tree** — nothing committed yet; suggested logical commits
  (a) C5′ + regenerated structures, (b) `model_parameters` Rmin_2, (c) K-B nascent radius
  wiring, (d) naming renames, (e) remove O'Brien loaders, (f) repoint tut15/P0CX28,
  (g) docstrings + docs.
- [ ] **Docs duplicate-object warnings (160)** — flat apidoc stubs (`docs/topo.*.rst`) and
  curated pages describe the same modules. Decide: drop the overlapping curated pages, or
  `:no-index:` one set, or `napoleon_use_ivar = True`. Build still succeeds.
- [ ] **PtR:76 A76 P-anchor** — kept 3.45 Å off O'Brien; splice O'Brien's A76 R coord into
  the CG PDB only if a bit-exact P-site is ever needed. *(Accepted deviation — DIFFERENCES
  item 6.)*
- [ ] **L24** — approximated per-AA instead of O'Brien's per-residue B-types (102 beads);
  mean protein Rmin/2 error 0.013 nm. *(Accepted — DIFFERENCES item 2.)*
- [ ] **Egress D5b** — the in-run 20k-step ejection is too short for a 306-mer; run
  `eject_demo.py` for the extended egress (apples-to-apples with tut15's PASS).
- [ ] *(note, not a task)* `interaction_details.md` documents the NC↔ribosome model.
- [ ] *(perf note, not a task)* CSP throughput is CPU-bound (per-stage OpenMM context
  rebuild + JIT), ~900 steps/s on a contended node — a node with free CPU cores is the win.

---

## D. `tutorials/10_csp_obrien/TASK.md` — "Remaining"

- [ ] **CHARMM PSF/TOP/PRM/COR ingestion** — run O'Brien's exact 4c5c + 50S systems.
- [ ] **Multi-trajectory multiprocessing** (`num_traj`/`tpn`/`ppn`, GPU device fan-out).
  *(Overlaps §F per-chain split + §B "many trajectories".)*
- [~] **Literal 3-stage mechanics** (peptide-bond toggling, explicit A/P tRNA bonded
  geometry) — currently mapped via the A→P restraint switch. *(This is DIFFERENCES §"Chain
  chemistry" — intentional KEEP; the `trna_tether` path now covers the bonded geometry.)*
- [ ] **`restart = 1`** resume of a partial trajectory. *(Dup of §B restart/resume.)*
- [ ] **Working external `ribosome_traffic` binary** (exits 127 here). *(Dup of §B ribosome
  traffic.)*
- [ ] **Quantitative validation vs. v6 reference outputs.**
- [ ] **Docs site page** for the tutorial.
- [ ] **Stage-1 blow-up at some residue additions** — 5/306 stages explode (PotE ~1e13)
  when a new bead is seeded on top of an existing bead; self-recovers next stage but
  corrupts those frames (see `OBSERVATIONS.md` #1). Fix: robustify seed placement /
  minimization (soft-core or large-force re-seed). *(Largely mitigated by the
  dt-halving guard + optimal PTC seeding; verify it's fully gone at full length.)*

---

## E. Tutorial 14/15 doc TODOs + INI debug markers

- [ ] [`tutorials/14_obrien_topo_consistency/step2_optimal_ptc_geometry.md:153`](../tutorials/14_obrien_topo_consistency/step2_optimal_ptc_geometry.md#L153)
  — run a **longer `L_max` baseline** to capture an actual dt-halving event for a full
  before/after, then **delete the stability guard** once equilibrium geometry is the
  default for AllBonds runs.
- [ ] **DEBUG-profile INIs must be switched to production before real runs** — `csp.ini`
  in tutorials 12/13/14 (and tut15 `csp.ini`/`csp_tether.ini`, `P0CX28/csp.ini`,
  `tutorials/16_csp_standalone/4c5c/csp_debug.ini`) use short `L_max`, inflated `scale_factor`
  (×50), and the step clamp with an explicit "delete both for production" comment. Not a
  code bug — a reminder that these are smoke profiles. *(Ties to §B "remove step clamps".)*

---

## F. Trajectory output & serialization (general runner, not CSP)

*(Merged in from the former root `todo.md`, 2026-07-01 — full detail preserved here.)*

- [ ] **Auto-split per-chain trajectories** — **multi-copy runs only** (`n_copies > 1`,
  non-interacting copies):
  - [ ] At end of run, optionally emit per-chain DCDs `traj_<k>.dcd` (+ appendix) so the
    manual `split_chains` step isn't needed. Make it a config flag (e.g.
    `split_chains = yes|no`), not unconditional.
  - [ ] Keep the combined `traj.dcd` as the canonical/source-of-truth output; per-chain
    files are *derived* (don't replace the combined one — it's needed for some analyses and
    is the restart/append target).
  - [ ] Per-chain centering applies to independent copies only (current
    `split_chains(center=True)`).
  - [ ] Do **NOT** auto-split or auto-center genuine *interacting* multichain
    systems/complexes (see "interacting chains" below) — the inter-chain arrangement is the
    physics; splitting + independent centering would destroy it.
  - [ ] Decide naming convention: code/optimizer currently use 0-based `traj_0..N-1`; note
    proposed 1-based `_1..N`. Pick one and document the mapping (e.g. `traj_1.dcd` = copy
    index 0).
  - [ ] Note: if split runs at finalize, a restart re-splits the whole grown DCD
    (idempotent but redundant) — acceptable, or split incrementally.
- [ ] **Interacting chains** (support, distinct from the split-copies case above).
- [ ] **OpenMM XML serialization** *(optional/later; discussed, revisit later)*:
  - [ ] **System XML export** (`XmlSerializer.serialize(system)`): optional flag (e.g.
    `write_system_xml = yes`) to dump the fully-built System (all forces + tabulated contact
    R/eps matrices + exclusions) to `<outname>_system.xml`. Value: reproducibility/provenance
    (pins every force-field number), decouples build from run, and lets reruns load the model
    with no STRIDE installed.
  - [ ] **State XML** as a *portable* restart alternative to the binary `.chk` (OpenMM
    checkpoints are not portable across GPU/OpenMM build/platform). Highest practical payoff
    for long resubmitted runs.
  - [ ] *(optional)* run-from-XML import path in the runner: `system_xml = path` -> skip
    build, load forces from XML (still build the cheap CA topology from PDB for the Topology
    + positions).
  - [ ] Caveats to handle: Topology not stored in System XML (pair with .psf/.pdb); n×n
    tabulated-function size for large complexes (few MB now, tens of MB for
    ribosome-nascent-chain); reconstruct topoReporter force-group *names* from group order
    after deserialize.
  - Skip ForceField-template XML (amber14.xml-style): residue templates don't fit a
    structure-based contact model.

---

## G. `DESIGN.md §5` — open questions (mostly resolved)

- [~] **Starting length `L0`** — numeric start value still a per-study choice (layout
  agreed). Not a blocker.
- [~] **Effective timescale mapping** — marked deferred in DESIGN; now largely addressed by
  the kinetics module (`scale_factor` × codon time ÷ dt). Superseded by §B "production
  parameters".

---

## Cross-cutting themes (the same work under several headings)

1. **Restart/resume** — §B, §D.
2. **Ribosome-traffic correction: finish or delete** — §B, §D.
3. **Production vs. debug run config** (remove step clamps, real dwell, GPU, many
   replicas) — §B, §E, §D.
4. **Rigid `AllBonds` as default + delete the dt-halving guard** once equilibrium seeding
   is proven at full length — §B, §E; interacts with DIFFERENCES §"Chain chemistry".
5. **Commit the working tree + fix docs duplicate-object warnings** — §C.

### Highest-leverage next actions (opinion)
- Commit the working tree (§C-1) — large, uncommitted, and blocks everything downstream.
- Decide the ribosome-traffic feature's fate (§B/§D) — dead-ish code carried in three
  places.
- Run the longer-`L_max` AllBonds baseline (§E) so the dt-halving guard can be retired and
  rigid bonds become the default — closes the biggest remaining ⚪ item in DIFFERENCES.
