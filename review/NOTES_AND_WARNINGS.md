# Review notes & warnings — tick when addressed

Every observation, warning, and recommendation I raised while reviewing `topo/csp/` against
the O'Brien reference and sweeping the tree for TODOs (2026-07-01, branch `tut15-claude-fix`,
HEAD `e9667c1`). Tick `[x]` as each is handled. Sources: this review + companions
[`DIFFERENCES.md`](DIFFERENCES.md), [`REVIEW_NOTES.md`](REVIEW_NOTES.md), [`TODO.md`](TODO.md).

Severity: ⚠️ = correctness / could bite silently · 📌 = decision needed · 📝 = doc/accuracy ·
💡 = recommendation.

---

## Warnings (⚠️ — could bite silently)

- [ ] ⚠️ **`trna_tether` default mismatch (footgun).** The `RunParams` dataclass
  defaults `trna_tether = True` ([`core.py:848`](../topo/csp/core.py#L848)), but the CSP
  runner `read_csp_config` defaults it to **False**
  ([`protocol.py:654-655`](../topo/csp/protocol.py#L654-L655)). So `topo-csp` runs the
  position-restraint path, but **any direct caller** of `run_length` / `RunParams`
  that bypasses `read_csp_config` silently gets the *tether* path. Fix: unify the default
  (make the dataclass default `False`), or add a comment + assertion at the call sites.

- [ ] ⚠️ **DEBUG-profile INIs will silently produce non-physical timescales if run as-is.**
  `csp.ini` in tutorials 12/13/14, tut15 (`csp.ini`, `csp_tether.ini`, `P0CX28/csp.ini`),
  and `tutorials/16_csp_standalone/4c5c/csp_debug.ini` use short `L_max`, inflated `scale_factor` (×50),
  and the per-stage step clamp — all marked "delete for production". A production run must
  clear these. See [`TODO.md`](TODO.md) §E.

- [ ] ⚠️ **Stage-1 blow-up on seed-on-top-of-bead (tut10).** 5/306 stages historically
  exploded (PotE ~1e13) when a new bead seeded onto an existing one; self-recovers next
  stage but corrupts those frames (`tutorials/10_csp_obrien/OBSERVATIONS.md` #1). The
  dt-halving guard + optimal PTC seeding should mask it now — **verify it is fully gone at
  full length** before trusting per-frame data.

- [ ] ⚠️ **Working tree is large and uncommitted.** The entire standalone-ribosome
  migration lives only in the working tree (CHANGELOG §8-1). Anything downstream (further
  edits, runs, a lost node) risks losing it. Commit before continuing.

---

## Decisions needed (📌 — "finish or delete" / pick one)

- [ ] 📌 **Ribosome-traffic correction: finish or delete.** Off by default and hidden since
  2026-06-30, but the code is still carried in `RunParams`, `read_csp_config`,
  `kinetics.ribosome_traffic_times`, and a `ribosome_traffic=off` banner. It depends on an
  external binary that is not bundled (exits 127). Decide: (a) vendor/reimplement +
  validate + re-expose, or (b) remove the fields/parsing/function/banner. (TODO §B, §D.)

- [ ] 📌 **Per-stage step clamps: drop or guard.** `max_steps_per_stage` /
  `min_steps_per_stage` are testing-only and break the dwell→steps timescale mapping.
  Decide: remove the fields + parsing + `stage_steps` clamp args, or keep with a loud
  production warning. (TODO §B.)

- [ ] 📌 **Rigid `AllBonds` as the default + retire the dt-halving guard.** Equilibrium PTC
  seeding (`optimize_ptc_geometry`) makes rigid `AllBonds` seed cleanly, but flexible
  bonds + the guard are still the default. Decision blocked on the baseline below.

- [ ] 📌 **A76 P-anchor (3.45 Å off O'Brien).** Kept as topo's; PTC-seeding absorbs it.
  Splice O'Brien's single A76 R coord in only if a bit-exact P-site is ever required.
  (Accepted deviation — confirm the decision stands.)

- [ ] 📌 **L24 radii (per-AA vs O'Brien per-residue B-types).** Approximated per-AA (mean
  Rmin/2 error 0.013 nm). Accepted; confirm.

---

## Recommendations (💡 — highest leverage first)

- [ ] 💡 **Commit the working tree** in the logical chunks CHANGELOG §8-1 lists (C5′ +
  structures / `model_parameters` Rmin_2 / K-B nascent wiring / renames / loader removal /
  repoint / docs). Unblocks everything.
- [ ] 💡 **Run the longer-`L_max` AllBonds baseline** to capture a real dt-halving event for
  a before/after, then retire the guard and make rigid bonds the default — closes the
  biggest remaining ⚪ item in DIFFERENCES §"Chain chemistry".
  (`tutorials/14.../step2_optimal_ptc_geometry.md:153`.)
- [ ] 💡 **Fix the 160 docs duplicate-object warnings** (flat apidoc stubs vs curated pages):
  drop the overlaps, `:no-index:` one set, or set `napoleon_use_ivar = True`. (CHANGELOG §8-2.)
- [ ] 💡 **Add the C-terminal mobility window** (freeze all but the last N=15 nascent
  residues) — the last high-value 🟢 in DIFFERENCES still open. NB tut15 D9 flagged mass-0
  conflicts with `AllBonds` and with topo's diffusion-extrusion route, so this needs design
  thought, not a drop-in mask.

---

## Doc / accuracy notes (📝)

- [ ] 📝 **L24 free-loop residue range was wrong in the original DIFFERENCES.md.** It said
  42–**59**; the O'Brien reference is `ribo_free_mask = L24 : 42-56` (line 37). Corrected in
  the revised doc — propagate the fix if the number appears elsewhere.
- [ ] 📝 **"Placement geometry" is superseded, not matched.** topo's `optimize_ptc_geometry`
  solves for a least-buried seed rather than reproducing O'Brien's fixed 4.27 Å + 10° tilt.
  Same intent, different mechanism → listed as "still different" but not a regression to
  close. Only add the fixed tilt if exact O'Brien seed reproduction is required.
- [ ] 📝 **Tether item is closed only in the tether path.** "Restrain previous AA +
  orientation control" is done when `trna_tether = yes`, still absent in the default
  position-restraint path — the split verdict must be kept when quoting status.
- [ ] 📝 **Original DIFFERENCES.md predates the alignment commits.** Anyone reading the
  root/tutorial copy should use `review/DIFFERENCES.md` instead — the root copy is a
  verbatim 2026-06-30 baseline kept only for diffing.

---

## Scope caveats of this review (what I did *not* verify)

- [ ] 📝 Verdicts are from **static code + git-history reading**, not fresh runs. The "still
  different / done" tags reflect what the source does, not a re-executed validation. Re-run
  the tut13-style debug check if a tag is load-bearing for a decision.
- [ ] 📝 I swept tracked text/code files for `TODO/FIXME/XXX/HACK/BUG` + dedicated
  TODO/TASK/open-items docs. I did **not** inspect binary assets, `__pycache__`, `docs/_build`,
  or untracked run directories for embedded notes.

---

## Suggestions for improvement (🔧 — beyond the existing TODOs)

Forward-looking ideas I noticed that aren't tracked anywhere yet. Grouped; tick if adopted
or consciously declined.

### Architecture / code quality
- [ ] 🔧 **Single source of truth for defaults.** `read_csp_config` re-hardcodes defaults
  that also live on the `RunParams`/`RunParams` dataclasses (the `trna_tether`
  mismatch is the symptom). Have the INI reader fall back to the dataclass default when a
  key is absent, so a default is defined once.
- [ ] 🔧 **Warn on unknown INI keys.** `opt()` returns `None` for absent keys, so a typo'd
  key (`trna_teather`, `L_maX`) silently falls back to a default with no error. Validate the
  parsed INI against the known key set and warn/raise on unknowns — cheap, prevents silent
  mis-runs.
- [ ] 🔧 **Quarantine experimental/off-by-default code.** Ribosome-traffic and the step
  clamps are carried live across several modules while disabled. If not deleted, move them
  behind an `experimental`/feature-flag boundary so the main path reads clean.
- [ ] 🔧 **`optimal_ptc_targets` uses a constant `aa_rmin_2_nm = 0.5`** (the conservative
  max K-B radius) for the seed EV, while the simulation uses the *actual* per-residue K-B
  radius. Consider passing the seeded residue's real radius for full seed↔sim consistency
  (the seeding-consistency fix already aligned the *form*; this aligns the *value*).

### Testing / reproducibility
- [ ] 🔧 **Add numeric regression tests.** The CHANGELOG verifies "numerically identical
  across renames" by hand (native contacts 819, K-B Rmin/2 0.251–0.531 nm,
  `rmin_matrix[i,j] = rmin_2[i]+rmin_2[j]`, ribosome radii/charges vs O'Brien). Freeze these
  as pytest assertions (there's already an untracked `test/`) so a future rename can't
  silently drift the physics.
- [ ] 🔧 **Automate the D5/D6 validation.** Turn the manual threading/dwell/Rg checks into a
  script that emits the metrics vs the reference bundle and PASS/FAIL — makes "quantitative
  validation vs v6" (tut10 remaining) repeatable rather than a one-off.
- [ ] 🔧 **Log and surface the RNG seed** for FPT sampling + velocity re-draw, per run and
  per replica. Essential once replica studies start (the tether-validation item needs
  reproducible independent seeds).

### Performance
- [ ] 🔧 **Reuse the OpenMM Context across the 3 stages of one length.** CHANGELOG §8-7
  notes CSP is CPU-bound on per-stage context rebuild + JIT. The topology is identical for
  a given `L` across stages 1→2→3 — only the restraint target point and step count change.
  Rebuilding once per length and just updating the restraint parameter (and skipping
  minimize, which you already do for stage 2) could cut the rebuild/JIT cost ~3×. Biggest
  throughput lever available.
- [ ] 🔧 **Narrow the dt-halving retry.** On divergence the guard re-runs the *whole* stage
  at halved dt. Consider retrying only from the last good checkpoint / a sub-window, so a
  late blow-up doesn't repay the whole stage.

### Docs
- [ ] 🔧 **Consolidate the overlapping doc set.** `DESIGN.md`, `FILES.md`, `PROMPT.md`,
  `TODO.md`, `README.md`, plus per-tutorial `TASK.md`/`NOTES.md` and the top-level
  `CHANGELOG.md` carry overlapping state that drifts (the 42-56 vs 42-59 error is a symptom).
  Pick one canonical status doc per concern and cross-link rather than duplicate.

### Science
- [ ] 🔧 **Benchmark whole-chain vs C-terminal-window mobility as a physics question, not
  just fidelity.** Moving the entire chain (topo) vs only the C-terminal 15 (O'Brien)
  changes co-translational relaxation, not only cost. A controlled A/B (folding order,
  Q-vs-length, Rg) would tell you whether the deviation actually matters for results —
  turning a 🟢 "should match" into an evidence-based keep-or-fix.
- [ ] 🔧 **Report a confidence interval on dwell/Rg from replicas.** Current D6 ratios
  (1.01×, 1.06×) are single-run point estimates; O'Brien runs ~50 trajectories/protein.
  Once the replica driver exists, quote mean ± CI so "matches reference" is statistically
  grounded.
