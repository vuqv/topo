# Tutorial 15 (`claude_fix`) — consolidated changelog

A single narrative of everything changed while building/validating Tutorial 15, tying together
the incremental commits (see `git log` for the per-step history). Final state: the O'Brien
continuous-synthesis path is reproduced for **4c5c (L=1→306)** and **P0CX28 (L=1→106)** with the
**FINAL GOAL met on both** — full-length synthesis, finite energy, **no nascent–ribosome clash**,
clean +x egress, and no leak through the truncated ribosome.

## 1. Validate the "claude fix" (equilibrium-PTC + rigid `AllBonds`)
- Wrote debug/full INIs on the unpinned fix path; validated D0–D6 for 4c5c: seed peptide bond
  3.810 Å, **no dt-halving**, dwell **1.01×** ref, R_g in range. Extended `eject_demo.py` shows
  clean +x egress. D8: P0CX28 L=1→106 (D6 N/A — no O'Brien reference).
- Files: `csp.ini`, `csp_val.ini`, `eject_demo.py`, `analyze_validation.py` (REF → `../12_auto/`),
  `README.md`, `NOTES.md`, `TASK.md`, `P0CX28/*`.

## 2. §1b O'Brien-consistency features (validate-first)
- **✅2/✅3 implemented** (behind `trna_tether = yes`): full O'Brien tRNA tether — bond + 2 angles
  + improper, A-site stages 1–2 / P-site stage 3, previous-AA tether in stage 1
  (`ribosome.py:add_trna_tether`, `core.py:run_length`, `protocol.py`). Stable; kinetics identical;
  does not by itself fix the clash → kept as an option, position restraint stays default.
- **✅1/✅4 deferred** (incompatible with topo's diffusion-extrusion / scenery-only ribosome),
  **✅5 superseded** by equil-PTC seeding. Documented in `NOTES.md`.

## 3. Force-field audit — `TOPO_OBrien_params_compare.md`
- Side-by-side of every energy term. **Angles, dihedrals, electrostatics, 12-10-6 Gō match O'Brien
  exactly** (dihedrals verified ratio 1.000; the ×0.756 reconstructs the prm). Flagged two real gaps
  → both since fixed (§5, §6).

## 4. Ribosome CG-coordinate fix (user-directed)
- topo's home-built ribosome placed the **ribose (R) bead off O'Brien's C5′** (P-anchor off 3.6 Å)
  and used a **uniform 0.71 nm** radius vs O'Brien's per-type Rmin/2. Added
  `load_obrien_ribosome(.cor/.psf/.prm)` + `load_ribosome_auto` (dispatch on `.cor`); wired CSP,
  `eject_demo`, and analyzers to O'Brien's authentic truncated ribosome (`ribosome_obrien.*`).
  Improved fidelity (R_g 0.96×), but revealed the clash was **not** ribosome-geometry alone.

## 5. NC↔ribosome nonbonded interaction — root-cause fix (user-directed)
- Deep comparison (`TOPO_OBrien_NCribosome_nonbonded_compare.md`): topo's NC↔ribosome excluded
  volume was **~1000× too soft** — pure `ε(σ/r)¹²` + **average** rule vs O'Brien's **12-10-6**
  (`ε[13(R/r)¹²−18(R/r)¹⁰+4(R/r)⁶]`) + **sum** rule `R_ij = Rmin/2_i + Rmin/2_j`. Confirmed the
  user's "Rmin/2 = σ/2" point from code.
- **Fix (`ribosome.py:append_ribosome`)**: 12-10-6 form; sum rule; matched ε; **nascent per-AA
  `SA..SY` Rmin/2** (`OBRIEN_SC_RMIN2_NM`, user "Option B"); ribosome per-type Rmin/2 from the prm.
- **Stability**: the stiff (correct) wall violently ejects previously-embedded beads → NaN, so the
  **dt-halving guard was made NaN-robust** (`core.py:run_length` catches blow-ups + minimize
  failures, retries at halved dt, preserving dwell time).
- **Result — clash RESOLVED on both proteins:**
  - 4c5c: hard clashes **36→0**, min NC dist 1.98→**3.32 Å**, extrusion x→**179 Å**, **D5b all PASS**
    (96 dt-halving recoveries; max|PotE| 1.9e6 ≪ 1e9).
  - P0CX28: hard clashes **18→0**, threading corr −0.64→**−0.969**, extrusion x→**98 Å**,
    R_g **2.31 nm ≈ native 2.11**, **D5b all PASS** (22 recoveries).
  - No y-z tunnel wall needed — the NC↔ribosome EV was the root cause.

## 6. Bond force-constant fix (repo-wide)
- `model_parameters['topo']['bond_force_constant']` **20920 → 41840** (added the missing
  CHARMM→OpenMM ×2). Now E = 209.2 kJ/mol at a 1 Å stretch = O'Brien/CHARMM (50 kcal/mol/Å²). Moot
  under `AllBonds` (bonds are constraints; tut15 unaffected), corrects flexible-bond runs.

## 7. Nascent NC↔ribosome Rmin/2: Option A vs B, seeding fix, speedup (2026-06-30, after `0c8f7d0`)

The initial NC↔ribosome fix (§5) used **Option B** for the *nascent* side (per-AA `SA..SY`). Verified
from O'Brien's source (`rnc.prm`/`rnc.xml`: nascent = per-position `A1..An`, NBFIX only A–A + B–B,
nascent↔ribosome = LB of self-Rmin) that O'Brien actually uses **per-position Karanicolas–Brooks**
radii — **Option A**. Switched to A (default), kept B optional, then removed A's cost with a seeding
fix and a stage-2 speedup. Commit-by-commit:

| commit | change |
|--------|--------|
| `1dc64d6` | Record Option B baseline metrics (for A/B compare) → `OPTION_A_vs_B.md`. Commit `0c8f7d0` **is** the Option B code state. |
| `83dfd6a` | **Option A implemented**: nascent Rmin/2 = per-residue K–B (structure-derived). Exposed via `build_nonbonded_interaction(return_rmin2=True)` → `precompute_contacts` (now 3-tuple) → `run_length(nascent_rmin2=)` → `append_ribosome(nascent_rmin2=)`. Reproduces O'Brien `.prm` A-values (76/306 exact, mean 0.06 nm). tut09 `cylinder.py` updated for the 3-tuple. |
| `efcf44d` | 4c5c A-vs-B recorded: A better (excl 4.12 vs 3.32 Å, energy 1.94e4 vs 1.9e6, egress +161 vs +22 Å); ~2.7× more dt-halving. |
| `a6c50e0` | P0CX28 A-vs-B + overall verdict: A ≈ B on P0CX28; adopt A. |
| `40f54e3` | **Optional key** `nascent_ev_radii = kb` (Option A, DEFAULT) \| `per_aa` (Option B). Both kept (`ElongationParams.nascent_ev_radii`; read in `read_csp_config`). |
| `2c01dec` | **Seeding-consistency fix**: `optimal_ptc_targets` now minimizes the *same* 12-10-6 sum-rule EV (`R = aa_rmin2_nm + Rmin/2_ribo`, default `aa_rmin2_nm=0.5`) as the simulation, so the A/P seed sits at the least-buried point of the real wall (a_target clearance 2.37 → 4.27 Å). Param renamed `aa_radius_nm` → `aa_rmin2_nm`. |
| `7bda1c6` | Seeding-fix **result** (4c5c): dt-halving **263/534 → 0/0**, worst \|PotE\| **1.94e4 → 1.02e3**, min NC-dist 4.12 → 4.79 Å, still 0 clash / no leak. Recorded in `OPTION_A_vs_B.md`. |
| `18af082` | Run **all profiles on GPU** (debug INIs `csp.ini`/`csp_tether.ini`/`P0CX28/csp.ini` switched CPU→GPU; full + eject already GPU). |
| `a2b8519` | **Low-risk speedup**: skip the redundant **stage-2 minimize** (continues from stage-1's relaxed final at the same A-site target). `run_length(minimize_override=)`; protocol passes `False` for stage 2. |

**Net (4c5c, Option A + seeding fix):** clash-free (0 hard clashes, min NC-dist 4.79 Å), no leak,
finite energy (worst 1.02e3 kJ/mol), **0 dt-halving**, D5b all PASS, dwell 1.01×. Measured cost split:
the seed solve (`optimal_ptc_targets`) runs **once per run** (not per length); per-stage cost is
minimize+MD (rebuild/context ≈ 0.5 s ≈ 80% of a fast 0.6 s stage). The "single build per L" refactor
was measured to save only ~40% (not 3×) and was declined in favor of the stage-2-minimize skip.

**Investigation records (read for the reasoning):** `TOPO_OBrien_NCribosome_nonbonded_compare.md`
(NC↔ribosome force comparison + root cause), `OPTION_A_vs_B.md` (A/B + seeding-fix tables),
`TOPO_OBrien_params_compare.md` (full force-field audit).

**Note on outputs (2026-07-01):** superseded run dirs (`synth_out_oldribo/`, `synth_out_softEV/`,
`synth_out_optB/`, `synth_out_optA_noseed/` and their P0CX28 twins) were deleted to stay under the
disk quota — all their metrics are preserved in the `.md` records above and they are regenerable via
the INIs.

## Code touched (outside this folder)
- `topo/csp/ribosome.py` — `load_obrien_ribosome`, `load_ribosome_auto`, `anchor_coord`,
  `OBRIEN_SC_RMIN2_NM`, full O'Brien tRNA tether, O'Brien 12-10-6 NC↔ribosome EV.
- `topo/csp/core.py` — NaN-robust dt-halving guard; `tether_segid`/`tether_prev_segid` in `run_length`.
- `topo/csp/protocol.py` — `.cor` ribosome dispatch, per-stage tether wiring, `trna_tether` INI key.
- `topo/parameters/model_parameters.py` — bond_force_constant fix.

## Reproduce
See `README.md`. Outputs (`synth_out*/`) are git-ignored; soft-EV baselines kept as
`synth_out_softEV/` for before/after visual comparison. Static ribosome overlay for VMD/PyMOL:
`ribosome_obrien_view.pdb` (committed convenience copy; the DCDs are nascent-only, so load the
ribosome separately). Tunnel/exit axis = +x.
