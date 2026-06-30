# AGENTS.md — Tutorial 15 (`claude_fix`) orientation + goal for AI agents

> **Read me first.** This file tells an agent **what this tutorial is for**, **what code
> change it validates**, **where the reference code and the topo port live**, and **how to
> work the task**. It merges the *autonomous-directive* style of
> [`../12_auto/Goal.md`](../12_auto/Goal.md) with the *orientation map* style of
> [`../13_validate_claude_fix12/AGENTS.md`](../13_validate_claude_fix12/AGENTS.md) and
> [`../14_obrien_topo_consistency/AGENTS.md`](../14_obrien_topo_consistency/AGENTS.md).
>
> **Agent directive (autonomous).** This is an autonomous task. Before any edits, **make a
> git commit** to checkpoint the repo (see §0). Then work the loop in
> [§6](#6-execution-loop) until the **final goal** below is achieved (or a
> [§8 stop condition](#8-stop-conditions) fires). Don't end a turn just because one step
> passed — re-read §5, find the first unmet criterion, and continue.
>
> **Maintain [`TASK.md`](TASK.md)** in this folder: a checklist of every concrete sub-task,
> **ticked `[x]` as it is finished** (with a one-line result), `[ ]` while open. Create it on
> the first turn from §5 + §1b, keep it the single source of truth for progress, and update it
> every milestone. Also log decisions/numbers in `NOTES.md`.
>
> ### FINAL GOAL (the loop runs until this is true)
> A **stable** continuous-synthesis run in which:
> 1. the protein is **synthesized successfully** to full length (4c5c: 1→306; then P0CX28: 1→106),
>    with finite per-stage energy throughout (no blow-ups);
> 2. the nascent chain **does not overlap / clash with any ribosome bead** at any frame, and
>    **never collapses back into the ribosome** — i.e. no nascent bead reaches **`x < 0`** along
>    the aligned exit-tunnel (+x) axis;
> 3. on tether release the chain **extrudes out through the exit tunnel** with **axial egress
>    (+x out the exit port)** and fully clears the ribosome;
> 4. it does so **without passing through the tunnel wall** of the **truncated CG ribosome**.
>    The wall is **radial (y–z)**: the chain must not leak sideways — within the tunnel x-window,
>    no nascent bead's radial distance `r = √(y² + z²)` may cross outward past the ribosome wall
>    beads. It exits via the real tunnel lumen, respecting the ribosome's excluded volume as a
>    **biological fact**, not by leaking through the truncation's open shell or back-face.
>
> This goal supersedes "all D-boxes ticked": the D-checks (§5) are how you *verify* it.

---

## 0. Before you start — checkpoint & repo scope

- **Commit first (mandatory).** On your **first turn**, make a git checkpoint commit so every
  change afterward is reversible:
  ```bash
  cd /storage/work/qzv5006/src/topo
  git add -A && git commit -m "tut15 claude_fix: checkpoint before O'Brien-consistency work"
  ```
  Commit again at each milestone (after a feature lands + validates) so progress is bisectable.
- **You may touch ANYTHING in this repo** to reach the final goal — `topo/csp/*`, other
  `topo/*`, `pyproject.toml`, tests, tutorials — not just this folder. Tutorial 15 is
  **independent of tutorials 12–14** (§9): you need not preserve their behavior or the legacy
  path. The only off-limits targets are the **raw inputs** and the **read-only reference**
  (§9). When a change is non-trivial, commit it on its own so it can be reverted in isolation.

## 1. What this tutorial is

Tutorial 15 validates **"the claude fix"** — the **equilibrium-bond PTC geometry + rigid
`AllBonds`** elongation path that became the **`topo-csp` default** in Tutorial 14
(steps 2–4 of [`../14_obrien_topo_consistency/step2_optimal_ptc_geometry.md`](../14_obrien_topo_consistency/step2_optimal_ptc_geometry.md)).

It is the **next link in the `validate_claude_fix` chain**:

| Tutorial | Fix it validated |
|----------|------------------|
| 12_auto | First O'Brien reproduction on 4c5c + the **dt-halving stability guard**. |
| 13_validate_claude_fix12 | Full-length (L=1→306) stress test of that **dt-halving guard**. |
| **15_claude_fix (this)** | The **equilibrium-PTC-geometry + `AllBonds` default** that makes the guard unnecessary. |

**The science question.** Tutorials 12 & 13 reproduced O'Brien using the *legacy* path
(far-seed at `AtR:76@R + 0.4 nm·x̂`, **flexible** harmonic bonds, propped up by the
dt-halving guard). The claude fix instead seeds every new residue **exactly 3.81 Å — one
peptide bond — from the current C-terminus** at the optimal PTC target points, so a **rigid
`AllBonds`** build is stable at 15 fs with **no dt-halving**. This is **more physically
faithful to O'Brien** (his bonds are rigid constraints, not springs).

> **Goal.** Show the new default **reproduces O'Brien's continuous-synthesis results on
> 4c5c** — physically sane synthesis + clean tunnel ejection (D-checks below) and
> quantitative agreement (dwell times / geometry) with the reference run — **and** that it
> does so on the *unpinned, current default* path (no `equil_peptide_geometry = no`, no
> `constraints = None`). **If validation exposes a defect in the default path, fix it in
> `topo/csp/*`** (this is a `claude_fix` tutorial — see §9 for the guardrails), then
> re-validate. **If 4c5c succeeds, repeat for P0CX28** in the [`P0CX28/`](P0CX28/) subfolder.

**What "the fix" is, precisely** (synthesized from Tutorial 14's notes — read those for the
derivation):

- **Optimal PTC target points.** `core.optimal_ptc_targets(ribosome)` solves for an A-site
  target `a_target` and P-site target `p_target` that are a **hard 0.381 nm apart**
  (= the `AllBonds` peptide-bond constraint length), with the two tRNA distances
  (0.427 / 0.476 nm to `AtR:76@R` / `PtR:76@R`) as **soft** harmonic penalties at O'Brien's
  `k = 83680 kJ/mol/nm²`, an **exit-side** inequality (`x ≥` the tRNA R-bead x, so the
  residue sits between tRNA and exit port), and ribosome excluded-volume clash minimized
  over a full-sphere Fibonacci multistart. Solved **once per run** (~17 s), in nm/kJ/rad.
- **Seed + migration.** Residue **L** is seeded at `a_target`, restrained to `a_target` in
  stages 1–2, then to `p_target` in stage 3 (the A→P hand-off). Because the previous
  residue rests at `p_target`, the new bond is born at exactly 3.81 Å — never stretched.
- **Rigid bonds, no guard.** Build with `constraints = "AllBonds"`; the dt-halving guard in
  `run_length` never fires (verified in Tutorial 14: max |PotE| ≈ 30–85 kJ/mol vs the 1e9
  divergence limit; old far-seed path peaked ~500 kJ/mol from a 1.9×-stretched seed bond).
- **Now the default.** `protocol.csp_default_elong()` sets `equil_peptide_geometry=True,
  constraints="AllBonds"`; `CSPParams.elong` / `read_csp_config` use it when the INI is
  silent. Legacy path stays reachable via `equil_peptide_geometry = no` + `constraints =
  None` (how 12 & 13 are pinned). **This tutorial deliberately leaves the default ON.**

---

## 1b. Consistency targets for this tutorial (DECIDED with the user)

Beyond validating the current default, Tutorial 15 will **close these specific O'Brien gaps**
(chosen from [`DIFFERENCES.md`](DIFFERENCES.md) — read it for the full side-by-side; **do not
edit it**). Already-closed items (rigid `AllBonds`, equilibrium-PTC seeding) are *not* relisted.

**Approach = VALIDATE-FIRST, then decide.** Do **not** implement everything up front:

1. **Baseline first.** Run the *current* AllBonds default on 4c5c and compare dwell times /
   geometry to the reference (`../12_auto/continuous_synthesis/output/`). Record the gaps the
   comparison *actually* shows (D6).
2. **Then implement** the selected features below — **one at a time**, validating dwell/geometry
   vs the reference after each (debug-then-scale) — **prioritizing whatever the baseline
   comparison shows matters most.** A feature that doesn't move the numbers can be deferred.

**Selected features (✅ in scope):**

| # | Feature | What to change | Code location |
|---|---------|----------------|---------------|
| ✅1 | **C-terminal mobility window** | Freeze all but the last **N=15** nascent residues each stage (mass-0); the rest of the chain currently moves. | `run_length` ([`core.py`](../../topo/csp/core.py)) — add an optional freeze mask. |
| ✅2 | **tRNA tether + orientation** | Replace the isotropic position restraint with O'Brien's **bond + 2 angles + improper** to tRNA beads (A-site stages 1–2, P-site stage 3). | `add_trna_tether` ([`ribosome.py`](../../topo/csp/ribosome.py), currently a bond + 1 angle to `PtR:76@R`, **unused by CSP**); `protocol.py:191` forces `trna_tether=False`. |
| ✅3 | **Restrain previous AA (L−1)** | Also restrain residue **L−1** (P-site) and use the **`P`/`PU2`** sub-beads for the angle/improper, not just the `R` bead. Pairs with ✅2. | `add_cterm_restraint` ([`core.py`](../../topo/csp/core.py)) — currently restrains only residue L. |
| ✅4 | **Ribosome L24 free loop** | Let ribosome residues **42–59** (L24 tunnel loop) move instead of freezing the whole ribosome. | rigid-scenery build in `run_length`/[`ribosome.py`](../../topo/csp/ribosome.py) — add a free mask. |
| ✅5 | **Placement 10° off-axis tilt** | Add O'Brien's 10° xy tilt to the seed direction. Minor (largely superseded by optimal-PTC targeting). | `seed_positions` ([`core.py`](../../topo/csp/core.py)). |

**NOT in scope (❌ — user excluded):** ribosome-traffic correction (kinetics only; leave off),
tunnel-wall-post-only (keep topo's persistent anti-leak wall).

> Because tut 15 is **independent of tut 12–14** (§9), these may change shared CSP defaults
> directly — no need to keep the legacy path reachable.

## 2. Where the reference code is (read-only — semantics only, never run/modify as deliverable)

| Path | What it is |
|------|-----------|
| `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` | O'Brien lab CG model + protocol suite (Jiang, Nissley, O'Brien). |
| `…/Continuous_synthesis_protocol/continuous_synthesis_v6.py` | **The protocol being reproduced.** Key fns: `run_elongation`, `elongation`, `A_site_tRNA_binding`, `peptide_bond_formation`, `translocation_AtR`, `create_elongation_system`. |
| `…/Continuous_synthesis_protocol/ribosome_traffic` | Upstream-queue delay binary. |
| `…/CG_protein_parameterization/`, `…/CG_ribosome_parameterization/` | How the reference CG `.psf/.top/.prm` were built. |

**Reference *run* to validate against** (inputs + outputs) lives in Tutorial 12, not here:
[`../12_auto/continuous_synthesis/`](../12_auto/continuous_synthesis/)
(`input/cont_synth_ecoli.cntrl`, `input/setup/`, `output/info.log`, `output/output/1.out`).
The reference covers **4c5c only**; P0CX28 has no O'Brien reference run (see §4).

## 3. The topo port (reuse — do **not** re-implement)

| Path | Role |
|------|------|
| [`../../topo/csp/core.py`](../../topo/csp/core.py) | MD building blocks: `optimal_ptc_targets` (the fix's solver), `seed_positions` (`seed_point` arg), `add_cterm_restraint`, `run_length` (+ dt-halving guard, now dormant on the default path), `ElongationParams`. |
| [`../../topo/csp/protocol.py`](../../topo/csp/protocol.py) | 3-stage outer loop; `csp_default_elong()`, `CSPParams`, `read_csp_config`; forces `trna_tether=False` (position-restraint path). |
| [`../../topo/csp/ribosome.py`](../../topo/csp/ribosome.py) | Rigid ribosome scenery, tunnel wall, EV constants. |
| [`../../topo/csp/kinetics.py`](../../topo/csp/kinetics.py) | FPT/dwell-time sampling + ribosome-traffic correction. |
| [`../../topo/csp/DESIGN.md`](../../topo/csp/DESIGN.md) | Design rationale + invariants. |
| `topo-csp` (console script) | The runner: `topo-csp -f <ini>`. |

Background on the fix: Tutorial 14's
[`step1_allbonds_no_dt_halving.md`](../14_obrien_topo_consistency/step1_allbonds_no_dt_halving.md),
[`step2_optimal_ptc_geometry.md`](../14_obrien_topo_consistency/step2_optimal_ptc_geometry.md),
[`DIFFERENCES.md`](../14_obrien_topo_consistency/DIFFERENCES.md).

## 4. This folder's structure

```
15_claude_fix/
├── AGENTS.md                 # this file
├── TASK.md                   # (you create + maintain) ticked checklist of every sub-task — source of truth for progress
├── DIFFERENCES.md            # O'Brien-vs-topo gap menu — READ ONLY, do not edit
├── NOTES.md                  # (you create) decisions, deviations, validation tables
├── README.md                 # (you create) exact reproduce commands + result tables
├── 4c5c_model_clean.pdb      # raw input — 4c5c all-atom (306 aa)   [copy from tut 14]
├── 4c5c_model_clean_stride.dat
├── 4c5c_mrna.txt             # one codon / residue
├── trans_times.txt           # Fluitt E. coli codon-time table
├── ribosome_trunc.pdb        # truncated CG 50S + tRNAs (X-aligned); P-/A-anchors
├── domain.yaml               # 4c5c 3-domain map + Go-scale strengths
├── analyze_validation.py     # D5 energy scan + D5b ejection check [copy from tut 12/13]
├── csp.ini                   # DEBUG profile (short L, large scale_factor) -> synth_out_debug/
├── csp_val.ini               # FULL-LENGTH profile (L=1->306)            -> synth_out/
├── synth_out_debug/  synth_out/        # OUTPUTS ONLY — safe to regenerate
└── P0CX28/                   # the second target, run ONLY after 4c5c passes
    ├── P0CX28_clean.pdb  P0CX28_clean_stride.dat  P0CX28_mrna.txt
    ├── trans_times.txt  ribosome_trunc.pdb  domain.yaml  (single-domain, strength 2.5044)
    ├── csp.ini  csp_val.ini        (L_max = 106)
    └── synth_out*/                 # OUTPUTS ONLY
```

**Inputs are copies** of existing validated inputs — 4c5c from
[`../14_obrien_topo_consistency/`](../14_obrien_topo_consistency/) (or `../13_*`), P0CX28
from [`../14b_P0CX28/`](../14b_P0CX28/). **Never overwrite** any raw input. The fixed-point
INI shape to mirror: tut 14's `csp_step2_allbonds.ini` for the **debug** profile and tut
14c / 14b's `csp_*.ini` for the **full** profile (note: those leave the default ON by simply
*not* pinning `equil_peptide_geometry`/`constraints`).

**P0CX28 caveat.** There is **no O'Brien reference run** for P0CX28, so D6 (quantitative
match) does **not** apply there — the P0CX28 deliverable is D3/D4/D5/D5b (runs, sane
energies, clean ejection) + internal consistency (e.g. folds to a sensible R_g), **not** a
dwell-time match. State this in `P0CX28/NOTES.md`.

## 5. Definition of Done

Verify each — don't assume it. 4c5c first; the P0CX28 subset only after 4c5c's D1–D7 pass.

- [ ] **D0 — Scaffold.** Inputs copied into this folder (and `P0CX28/`); raw inputs never
      overwritten. `analyze_validation.py` present.
- [ ] **D1 — Configs.** `csp.ini` (debug) + `csp_val.ini` (full) exist and **use the current
      default fix path** — i.e. they do **not** pin `equil_peptide_geometry = no` and do
      **not** set `constraints = None`. Confirm the run banner reports the
      equilibrium-geometry + `AllBonds` path (not "flexible (harmonic)"). Kinetics
      (`time_stage_1/2`, `scale_factor`, `mrna`, `trans_times`) match the reference mapping.
- [ ] **D2 — Fix is active.** A debug run prints `optimal_ptc_targets` output and shows the
      **seed peptide bond ≈ 3.79–3.81 Å** (not ~7.4 Å) and the **dt-halving guard never
      fires**. Record the target points + max |PotE|.
- [ ] **D3 — Run completes.** `topo-csp -f csp_val.ini` runs end-to-end, exit 0, full
      length (4c5c: 1→306). A documented short demo length is acceptable *only* with a
      stated reason.
- [ ] **D4 — Outputs.** Trajectory + per-residue `dwell_times.dat` under `synth_out/`
      (mirroring the reference layout).
- [ ] **D5 — Physically sane.** No stage |PotE| ≳ 1e9; chain threads the tunnel (monotonic-ish
      +x egress), never collapses back into the ribosome. **Collapse is defined geometrically:
      NO nascent bead may reach `x < 0`** along the aligned exit-tunnel axis (the tunnel central
      line is the +x axis; `x < 0` means a bead has gone *behind the PTC, back into the ribosome
      body*). Check: `min nascent x` over the whole run **≥ 0**.
- [ ] **D5b — Clean ejection (axial egress, no radial wall penetration).** On tether release the
      chain must leave **axially** — CoM and beads advance in **+x out the exit port** (x ≈
      `x_exit`) — and **never radially through the tunnel wall**. **The tunnel wall is defined in
      the y–z plane**, not just by an x-plane: the wall is the ribosome beads forming the lumen
      surface around the x-axis, so a bead *penetrates the wall* if — **while still inside the
      tunnel x-window (`0 ≤ x ≲ x_exit`)** — its **radial distance `r = √(y² + z²)` from the
      tunnel axis crosses outward past the ribosome wall beads / lumen radius** (i.e. it leaks
      sideways through the truncated ribosome's open shell instead of threading out the port).
      Verify **all three**: (a) `min nascent x ≥ 0` (no back-collapse, D5); (b) **radial check**
      — within the tunnel x-window, no nascent bead's `r` exceeds the local ribosome-wall radius
      (no sideways leak); (c) `min nascent–ribosome distance` ≥ bead-contact σ (no steric
      overlap) throughout, energy finite. Record min nascent x, the radial-leak check, min
      nascent–ribosome distance, and CoM-x vs frame.
      > ⚠️ **`analyze_validation.py` is currently incomplete for D5b:** its wall test is
      > **x-plane only** (`min_x ≥ x0_A`, ~line 164) and the min-distance test can *miss* a
      > radial leak (a bead escaping past the truncated shell is *far* from all beads → looks
      > clean). **Extend it** with the (b) radial-leak check above (per-frame max nascent `r`
      > within the tunnel x-window vs the local ribosome-wall radius), keeping the x-plane and
      > min-distance checks. Commit that as its own change.
- [ ] **D6 — Quantitative match (4c5c only).** vs `../12_auto/continuous_synthesis/output/`:
      (a) total length matches; (b) per-codon dwell / total synthesis time agree within a
      stated tolerance (~2× on summed in-vivo time given stochastic FPT); (c) final geometry
      (R_g) in range. **Also compare the fix path vs the legacy 12/13 path** — the fix should
      reproduce the reference *at least as well* while being more physically faithful (rigid
      bonds, no guard). Numbers in `NOTES.md`.
- [ ] **D7 — Documented.** `README.md` gives exact reproduce commands; `NOTES.md` logs
      decisions + the validation table. Add the tutorial to `tutorials/README.md`.
- [ ] **D8 — P0CX28 extension.** *Only after D0–D7 pass for 4c5c:* repeat D0–D5/D5b/D7 in
      [`P0CX28/`](P0CX28/) (L=1→106). D6 is N/A (no reference) — report internal consistency
      instead.
- [ ] **D9 — O'Brien-consistency features (§1b).** *Validate-first:* after D6 establishes the
      baseline gap, implement the selected features (✅1–✅5) **only as the comparison shows
      they matter**, one at a time, re-validating dwell/geometry vs the reference after each.
      For each: record in `NOTES.md` what it changed (before/after numbers) and whether it
      moved the run closer to O'Brien. Excluded features (traffic, post-only wall) stay off.

## 6. Execution loop

Repeat until §5 is checked or §8 fires:

1. **Orient.** Re-read §5; pick the first unmet criterion.
2. **Plan** the smallest action that advances it.
3. **Act.** Copy/adapt inputs, write/adjust the INI (leave the fix default ON), or run `topo-csp`.
4. **Verify.** Inspect the banner, energies, seed-bond length, trajectory; compare to the
   reference. A finished run is a hypothesis to check, not a success.
5. **Record** the milestone + findings in `NOTES.md`.
6. **Iterate.** On failure, diagnose (traceback / energy log), fix, retry.

**Debug-then-scale.** Get a short run green (D1/D2/D5 at small `L_max`, large `scale_factor`)
before the full-length validation (D3/D6).

> **Fast-test `scale_factor`.** `steps = dwell_s · 1e9 / scale_factor / dt`, so a *larger*
> `scale_factor` → *fewer* MD steps = quicker runs. For debug use
> `scale_factor = 4331293 × 50 = 216564650` (50× the production `4331293`). Restore `4331293`
> for the D6 validation run; note which value produced each result in `NOTES.md`.

### Quick-test INI options (technical)

For the **debug/short** profile (`csp.ini`), pin the per-stage step bounds and output
stride so each stage is cheap but still long enough to expose a blow-up:

```ini
max_steps_per_stage = 2000          # documented cap (see header) — hard ceiling per stage
min_steps_per_stage = 400           # floor: enough steps to surface an instability before the cap
nstout              = 200           # write every 200 steps -> ~5-10 frames/stage for the D5/D5b checks
```

Notes:
- These bound only the **MD step count per stage**; they do **not** change the sampled
  in-vivo dwell times (`dwell = n_steps × dt` is reported separately). A stage's real step
  count is `clamp(dwell·1e9/scale_factor/dt, min_steps_per_stage, max_steps_per_stage)`.
- `min_steps_per_stage = 400` (vs the `100` used in tut 14b/14c) makes the quick run a
  **more honest stability probe** — too small a floor can finish a stage before a stiff Go
  well has time to diverge, hiding a problem the full run would hit.
- `nstout = 200` keeps enough frames per stage for `analyze_validation.py` (D5 energy scan,
  D5b ejection) without bloating `synth_out_debug/`.
- Combine with a small `L_max` (e.g. 5–10) and the 50× `scale_factor` above for the smoke
  test; the **full** `csp_val.ini` restores production `scale_factor`, full `L_max`, and
  larger step caps.

## 7. How to run

From this folder (4c5c), then from `P0CX28/` after 4c5c passes:

```bash
topo-csp -f csp.ini        # debug smoke test -> synth_out_debug/   (watch the geometry banner + seed-bond length)
topo-csp -f csp_val.ini    # full-length      -> synth_out/         (watch for "[stability]" lines — should be NONE)
python analyze_validation.py   # D5 energy scan + D5b ejection check
```

Headline checks: **no stage |PotE| exceeds 1e9 kJ/mol**, **no `[stability]` dt-halving
line**, **seed peptide bond ≈ 3.81 Å**, **clean +x ejection**.

## 8. Stop conditions (pause and ask the user)

- A required input is missing/corrupt and cannot be derived from this repo or the reference dirs.
- An external dependency is unavailable (GPU/OpenMM platform, `scipy` for the solver, the
  `ribosome_traffic` binary) **and** no documented fallback suffices.
- The same run fails the same way after **3** distinct fix attempts — report the traceback,
  what you tried, and your best hypothesis.
- A validation criterion (D5/D6) can't be met and the gap looks like a genuine **science**
  question (model mismatch), not a bug — surface it with the numbers (e.g. if the rigid
  `AllBonds` default reproduces the reference *worse* than the legacy path, that's a finding).
- Any step would overwrite a raw input or the reference data under
  `../12_auto/continuous_synthesis/`. **Write only under `synth_out*/`.**

When stopping, state: which milestone you reached, what's blocking, and the one decision or
input you need.

## 9. Safety rules

- The reference source in `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` is
  **read-only** — consult it; never modify or run it as the deliverable.
- **Never overwrite** raw inputs (`*_clean.pdb`, `*_stride.dat`, `*mrna*`, `trans_times.txt`,
  `ribosome_trunc.pdb`, `domain.yaml`) or the reference data under
  `../12_auto/continuous_synthesis/`. Write run outputs only under `synth_out*/`.
- **Editing anything in the repo is allowed and expected** (§0) — `topo/csp/*`, other
  `topo/*`, `pyproject.toml`, tests, tutorials. This is a `claude_fix` tutorial whose goal is
  an **O'Brien-consistent** synthesis path; change whatever it takes. Reuse/understand first,
  but a genuine gap or bug gets a real code fix here. **Checkpoint with a commit first** and
  commit each non-trivial change on its own so it stays revertible.
- **Tutorial 15 is INDEPENDENT of Tutorials 12–14.** You are **not** required to keep the
  legacy far-seed / flexible-bond / dt-halving path reachable, and you are **not** required
  to keep 12/13's pinned runs passing. Make whatever change brings the path closer to
  O'Brien, including changing shared defaults — do not constrain the work for back-compat.
- Record every edit + its before/after evidence (energies, seed-bond length, dwell/R_g) in
  `NOTES.md`, and which O'Brien feature (§ DIFFERENCES.md) it closes.
