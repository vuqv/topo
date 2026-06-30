# NOTES.md — Tutorial 15 (`claude_fix`) decisions & validation log

Chronological log of decisions, deviations, and validation numbers. See `TASK.md` for the
ticked checklist and AGENTS.md for the directive.

## 2026-06-30 — Orientation & setup

**Environment.** `bioenv` conda env: OpenMM 8.2, scipy, MDAnalysis, numpy all present;
`topo-csp` on PATH; GPU = Tesla T4. No blockers (AGENTS.md §8 deps satisfied).

**Code-state finding (important).** AGENTS.md §1 step 4 says the fix is "now the default" via a
`protocol.csp_default_elong()` that sets `equil_peptide_geometry=True, constraints="AllBonds"`.
**That function does not exist.** The actual code default is `ElongationParams.equil_peptide_geometry
= False` (core.py:737) and `constraints` is read from the INI (default `None` if silent). So the
fix path is **opt-in via the INI**, exactly as tut14c/14b enable it. Decision: the tut15 INIs set
`equil_peptide_geometry = yes` + `constraints = AllBonds` **explicitly** (matching the validated
tut14c `csp_4c5c.ini`). This satisfies D1 literally (no `=no`, no `=None`) and guarantees the fix
path regardless of code defaults. Flipping the code default to ON is a possible later cleanup
(noted, not required for the goal).

**Reference run scope.** `12_auto/continuous_synthesis/output/output/1.out` synthesized only
**L=1→10** ("Elongation finished at length 10"). So D6 dwell comparison is over residues 1–10;
there is no full-length 306 O'Brien reference. Mean in-vivo per-codon dwell times are deterministic
from the Fluitt codon table (RNG only affects the *sampled* values), so the robust D6 comparison is
mean in-vivo dwell per codon + summed total (~2× tolerance) + final R_g range.

**Inputs present.** 4c5c: `4c5c_model_clean.pdb` (306 aa, 2305 atoms), `_stride.dat`, `4c5c_mrna.txt`
(wrapped nt, ~307 codons), `trans_times.txt`, `ribosome_trunc.pdb`, `domain.yaml`,
`analyze_validation.py`. P0CX28 subfolder has its inputs (no `analyze_validation.py` yet; needed for D8).

**Plan.** Validate-first: (D1) write debug + full INIs → (D2) debug smoke test confirms fix active →
(D5/D6) baseline AllBonds default vs reference → then (D9) §1b features one at a time by measured need.

## 2026-06-30 — D1/D2 debug smoke test (4c5c)

`csp.ini` (debug): L=1→8, scale_factor=216564650 (50×), max/min steps 2000/400, nstout 200, CPU.
`topo-csp -f csp.ini` → exit 0. Banner confirms the fix path:
- `[equil_peptide_geometry] optimal PTC restraint targets (|A-P| = 0.3810 nm)`;
  A-target `[0.871,-0.165,-0.236]`, P-target `[1.007,0.115,-0.016]` nm.
- All 24 stages built `rigid (AllBonds)` (0 "flexible (harmonic)").
- Tunnel wall plane auto `x ≥ 0.8711 nm`. Anchors: P `[0.5705,0.298,-0.059]`, A `[0.828,-0.584,-0.137]`.

**D2 PASS:** seed peptide bond (L−1↔L) = **3.810 Å** at L=2,4,8 (exact equilibrium; AllBonds holds it,
not the ~7.4 Å far-seed). max |PotE| over the whole run = **42.78 kJ/mol** (limit 1e9 → 2.3e7× margin).
**Zero `[stability]` dt-halving lines.** Confirms the claude fix makes the guard unnecessary.

## 2026-06-30 — D3/D4/D5/D6 full-length baseline (4c5c, L=1→306)

`csp_val.ini`: L=1→306, scale_factor=216564650 (50×), max/min steps 2000/100, nstout 100, GPU,
ejection_steps=20000. `topo-csp -f csp_val.ini` → "Done. Synthesized 1 → 306", ejection ran.

- **D3 PASS** — exit clean, full length 1→306. **Zero `[stability]` dt-halving lines** across the
  whole run (919 stages). The claude fix makes the guard unnecessary at full length too.
- **D4 PASS** — per-stage trajectories + `synth_out/dwell_times.dat` (306 rows) + ejection/.
- **D5 PASS** — 919 stage logs scanned; worst max|PotE| = **1.48e3 kJ/mol** (at ejection) ≪ 1e12.
- **D6 (4c5c) PASS on the meaningful metrics** vs ref (L=1→10):
  - in-vivo total dwell ratio **1.01×** (per-codon nearly exact, e.g. L2 0.0386 vs 0.0386 s); ≪2× tol.
  - final R_g (L=10) **0.799 vs 0.750 nm = 1.06×** — in range.
  - in-silico sampled ns ratio 0.02× — expected: our 50× scale_factor → ~50× fewer in-silico ns.
    Dwell comparison is scale-independent (compared in seconds), so this is not a discrepancy.
  - vs legacy 12/13 path: the fix reproduces the reference at least as well (1.01× dwell, 1.06× Rg)
    while being more physically faithful (rigid AllBonds, **no dt-halving guard ever firing**).

**Chain topology (end of synthesis, L=306):** corr(residue index, x) = **−0.926** — the chain threads
the tunnel correctly: N-terminus extruded far out (x≈109 Å), C-terminus at the PTC (x≈10 Å), monotonic.
**Not balled up.** No nascent bead is deep in the ribosome (none at x<0); only 1 bead is below the
8.71 Å wall plane (grazing). So FINAL-GOAL #4 (no leak through the truncated tunnel wall) holds.

**Analyzer fix:** `analyze_validation.py` computed the tunnel wall at the *legacy* x0 (min(anchor)+
tether = 10.46 Å); the equil-PTC run actually used x0 = min(a_target,p_target) = **8.71 Å**. Fixed the
analyzer to recompute the wall via `optimal_ptc_targets` exactly as the runner does → **D5b wall check
now PASS** (min nascent x 8.37 Å vs wall 8.71 Å, within slack).

**D5b remaining gaps (the real work):**
1. *Clash:* min nascent–ribosome distance dips to **2.41 Å** during the free ejection (freed C-terminus
   relaxes into the PtR/PTC beads). At end-of-synthesis, 6 of 306 residues (218–276, all threading the
   23S tunnel) sit at 2.87–2.96 Å from rRNA beads (radius 7.1 Å) — soft-EV interpenetration, finite
   energy. **Calibration:** at the only comparable length L=10, ref min dist = 4.57 Å, ours = 3.36 Å
   (both clash-free). The sub-3 Å contacts are a *full-length tunnel-packing* effect, ours slightly
   tighter than O'Brien — the signal to add the §1b orientation/mobility features.
2. *Egress not demonstrated:* in-run ejection is only 20000 steps (~0.3 ns) — far too short for a folded
   306-mer to diffuse off; CoM-x barely moves. Need a proper extended-ejection demo for FINAL-GOAL #3.

### Validation table (fill as runs complete)
| Run | path | L range | scale_factor | max\|PotE\| (kJ/mol) | seed bond (Å) | dt-halving? | min NC dist (Å) | dwell ratio | Rg ratio | notes |
|-----|------|---------|--------------|----------------------|---------------|-------------|-----------------|-------------|----------|-------|
| debug | equil-PTC + AllBonds | 1→8 | 216564650 | 42.78 | 3.810 | none | — | — | — | D1/D2 PASS (CPU) |
| baseline full | equil-PTC + AllBonds | 1→306 | 216564650 | 1.48e3 | 3.810 | none | 2.41 (eject) / 2.87 (synth) | 1.01× | 1.06× | D3/D4/D5/D6 PASS; D5b clash+egress open |
