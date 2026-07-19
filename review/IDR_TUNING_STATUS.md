# IDR Rg-tuning — status & next steps (2026-07-19)

Working log for the effort to make topo's Cα-only IDR model reproduce experimental
SAXS Rg across 24 fully-disordered proteins. Sandbox (gitignored):
`sandbox/Tune_idr_scale/` — `scan/` (idr_scale α-scan), `scan2d/` (eps_gen scan),
`SOP-IDP/` (paper+SI), `exp_data.csv`, `hps_urry.csv`.

## Objective
Find IDR energy parameters so simulated Rg matches SAXS Rg for 24 IDPs. topo's IDR
term (`apply_disorder`, `topo/utils/nonbonded.py`): every disordered pair gets a
12-10-6 well at the sum-rule collision distance with depth
`max(NON_NATIVE_KJ, eps_gen_kj + idr_scale·|w−0.6|·KCAL_TO_KJ)`.

## What we did & found
1. **α-scan** (idr_scale ∈ {0,0.03,0.06,0.15,0.30,0.50,1,2,3}, folded start, 90 ns).
   - Rg *increases* with α up to ~0.3–0.5, then a **collapse branch** sets in (~1 kBT
     well). No single global α fits: at α=2 hCyp is +32% while IN is −22% (turnovers
     cross). Best full-set point is α=0 (RMS 33%). See `scan/REPORT.md`, `scan/figs/`.
   - α=3 (and most α=2) folded runs *hung* — traced to a **minimizer bug** (below).
2. **Reference/model comparisons** (memory: `idr-scale-expansion-mechanism`,
   `sop-idp-vs-topo-energy-scale`):
   - **O'Brien Cα-SCM** (generic-bt, JPCB 2019 9b02575): topo's SS term is a faithful
     port — same `0.03·|w−0.6|`, treat-BT-as-kcal, ~identical eps/kT (310 vs 300 K).
   - **SOP-IDP** (Baul/Thirumalai JPCB 2019): 2-bead, 3 channels BB 0.2 / BS 0.4 /
     SS 0.3 kBT, SS shift **0.7**. topo (SS-only, ~0.034 kBT at α=0.03) is ~7× weaker
     on SS and lacks the generic BB+BS cohesion.
   - **HPS-Urry** (`hps_urry.csv`): RMS 29%, slope 0.54 (compresses range). topo's
     error is a *systematic offset* (slope→0.95 at α=1, intercept ~+1.2 nm) — the
     "fixable" kind. `scan/figs/model_compare_*.png`.
3. **eps_gen scan** (idr_scale=1 fixed, eps_gen_kj ∈ {0,2.51,3.77,5.02} kJ, extended
   start, 40 ns; `scan2d/`). RMS stayed 43–48% — no compaction. **BUT the runs are NOT
   equilibrated** (see below), so this is INCONCLUSIVE, not a model verdict. The
   implementation is verified correct: eps_gen deepens every well by exactly +eps_gen_kj
   (~1→3 kBT). Deep wells *should* compact at equilibrium.
4. **THE blocking problem — equilibration** (`scan2d/figs/rg2d_rgt.png`): canonical MD
   does not equilibrate Rg for **N ≳ 100**. Rg(t) drifts / hops between compact and
   expanded basins (hCyp eps_gen=5: compact ~2.7 nm until ~18 ns, then jumps to ~4 nm),
   and even the 90 ns α=0-baseline swings 4↔7 nm for K19/osteopontin. Only the small
   chains (Histatin N=24) converge. **This confounds BOTH scans** — Rg values for the
   large chains are unreliable in all canonical-MD runs so far.
5. **Power/equilibration analysis** (`scan/figs/power_analysis.png`): Rg autocorrelation
   is short *for the expanding* α-scan (folded→coil is fast; 15 ns equil + ~30 ns prod
   is plenty). It is the *collapsing* regime (extended→globule) that traps.

## Committed model changes (this session)
- **`eps_gen_kj`** additive generic-cohesion term in `apply_disorder` + config plumbing
  (`domain.yaml: eps_gen_kj`, kJ/mol, default 0.0 = byte-identical to old build) +
  `tests/test_idr.py::test_eps_gen_additive`. Full IDR suite 11/11 green.
- **Minimizer fix** (`topo/core/system.py checkLargeForces`): the force-tolerance loop
  decremented `10→…→0.1` by float subtraction and never matched the `== 0.1` exit, so
  it **spun forever** whenever max force could not reach ≤10 kJ/mol/nm (any stiff CG
  IDP; deep wells especially). Now it **breaks at the tolerance floor and proceeds**
  (residual force > threshold with converged energy is benign — thermostat relaxes it),
  and **prints `[EM]` progress** each iteration (force_tol, max force, PE; flushed).
  This was the true cause of the earlier "hangs" (not the folded-start overlap).
- Not committed: `topo/parameters/model_parameters.py` mass corrections (pre-existing,
  unrelated — left for a separate decision).

## Open decisions / next steps (prioritized)
> **Enhanced sampling (replica exchange / simulated tempering) is OUT OF SCOPE — user
> decision 2026-07-19.** The equilibration limitation (§4 above) is therefore an
> ACCEPTED constraint, not a blocker: large-chain (N ≳ 100) Rg from canonical MD carries
> a non-equilibrium caveat and is used only qualitatively. Analysis leans on the regimes
> that DO converge — the small chains, and the folded→coil expansion (fast).

1. **Combining-rule experiment (main next step; UNDECIDED on adoption).** topo places the
   IDR well at `Rmin_2_i + Rmin_2_j`; O'Brien's generic-bt uses
   `rvdw_i + rvdw_j = (Rmin_2_i+Rmin_2_j)/2^(1/6)` — ~11% tighter contacts (verified:
   mean 0.639→0.569 nm for hCyp), which systematically compacts chains and could cut the
   ~+1.2 nm over-expansion offset. Change is a one-liner in `apply_disorder` step 2:
   divide the `rmin_matrix[involves_dis]` assignment by `RMIN_SCALE_FACTOR` (update
   `tests/test_idr.py::test_idr_folded_excluded_only` + docstring). **Run in a NEW folder
   (`scan2d_rvdw/`); do NOT clobber `scan2d/`.** Read the effect primarily on the SMALL
   chains that equilibrate cleanly (Histatin, IN, ACTR, Protein-L) — those isolate the
   geometric effect from the sampling caveat.
2. **Model decoupling (open question).** Whether a Cα-only model ultimately needs true
   decoupling (fixed purely-repulsive core + separately-scaled attraction, à la SOP-IDP's
   separate beads) vs whether tighter contacts + a generic term on the current well
   suffice — judged on the equilibrating cases only, given sampling is capped.

## Key artifacts
- `scan/REPORT.md`, `scan/figs/{rg_curves,rg_correlation_panels,power_analysis,model_compare_*}.png`
- `scan2d/REPORT2D.md`, `scan2d/rg2d_results.csv`,
  `scan2d/figs/{rg2d_correlation_panels,rg2d_curves,rg2d_rgt}.png`
- `review/IDR_EPS_GEN_PLAN.md` (eps_gen design + retracted-conclusion note)
- Memory: `idr-scale-expansion-mechanism`, `sop-idp-vs-topo-energy-scale`, `idr-design-finalized`
