# Consolidated TODO — topo source (collected 2026-07-01)

Every open item found by sweeping the source tree: inline code markers
(`TODO`/`FIXME`/`XXX`/`HACK`), dedicated TODO/TASK files, and "open items / remaining /
not yet" sections in docs. Grouped by source, deduplicated where the same item appears in
several places (cross-refs noted). Branch `tut15-claude-fix`, HEAD `e9667c1`.

Legend: `[ ]` open · `[~]` partial/mapped. (Done items are removed once completed.)

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

---

## B. CSP — validation, production, extensions

*(Absorbed the former `topo/csp/TODO.md`, deleted 2026-07-02; revised against the current
codebase — the standalone `elongate.py` runner + Tutorial 7 are gone, `run_continuous_synthesis`
is the only synthesis path.)*

### Validation & production

- [ ] **Analysis layer** (DESIGN §6 phase 4 — none exists) — Q-vs-length co-translational
  folding curves per domain (reuse `topo.analysis.native_contacts`); ejection-time
  observable (steps from tether release to tunnel exit); wire a worked example into the
  tutorial.


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
- [ ] **`ejection_steps = auto`** — stop ejection once the C-terminus clears the ribosome
  (`x(Cterm) − max(x_ribosome) > cutoff`, default 2.0 nm), chunked-stepping loop + safety
  cap; add `ejection_cutoff_nm`, `ejection_check_every`; consider same for
  `dissociation_steps`.

### Revision list
- [ ] **Ribosome-traffic correction** (hidden since 2026-06-30) — decide: (a) finish it
  (vendor/reimplement the `ribosome_traffic` upstream-queue correction, validate the
  intrinsic→real stretch, re-expose in docs+config) or (b) drop it entirely (remove
  `RunParams` fields, `read_csp_config` parsing, `kinetics.ribosome_traffic_times`). Also
  unhide/remove the `ribosome_traffic=off` runtime banner in `protocol.py`.
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

### Highest-leverage next actions (opinion)

- Decide the ribosome-traffic feature's fate (§B) — dead-ish code carried in three
  places.
