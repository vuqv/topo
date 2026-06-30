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

### Validation table (fill as runs complete)
| Run | path | L range | scale_factor | max\|PotE\| (kJ/mol) | seed bond (Å) | dt-halving? | notes |
|-----|------|---------|--------------|----------------------|---------------|-------------|-------|
| debug | equil-PTC + AllBonds | 1→8 | 216564650 | 42.78 | 3.810 | none | D1/D2 PASS (CPU) |
