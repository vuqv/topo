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

## 2026-06-30 — §1b design findings (before implementing features)

**✅1 C-terminal mobility window — INCOMPATIBLE with topo's extrusion as-is (key finding).**
OpenMM forbids any constraint involving a massless particle ("A constraint cannot involve a massless
particle" — verified), so mass-0 freezing cannot coexist with the rigid `AllBonds` build at the
frozen↔mobile boundary or within the frozen region. More fundamentally: **topo extrudes the chain by
*diffusion* of the whole mobile chain** under the moving C-terminus restraint + the one-sided tunnel
wall — there is **no explicit register translocation** (DIFFERENCES.md "same outcome, different route").
Each new residue is seeded at the PTC and residues 1..L−1 are carried over from `prev_final`; the
N-terminus reaches x≈109 Å only because the mobile bulk is pushed/diffuses outward over 3×306 stages.
If the bulk is frozen, those carried-over positions never advance, so residue 1 would stay at its
cold-start x≈6 Å forever and the chain piles up at the PTC — **freezing breaks extrusion**. To do ✅1
O'Brien's way would also require implementing explicit per-residue register translocation of the frozen
bulk (a substantial mechanism change topo deliberately avoided). → **✅1 deferred** as out-of-scope for
the diffusion-extrusion path; documented as a genuine model-route difference, not a quick mask.

**Egress geometry (✅ FINAL GOAL #3 scoping).** Truncated ribosome x-extent −14→112 Å; near the tunnel
axis it reaches ≈106 Å. At L=306 the chain is *already extruded* (N-term x=141 Å, threaded, no leak);
"C-terminus clears the ribosome" then means the whole folded 306-mer translates ~100 Å (+x) — a slow
post-release dissociation, not a short-MD event. The clean-egress demo is therefore most meaningful at
**short length** (the analyzer's `ejection_long/` design: a small chain traversing the tunnel and
popping out +x). Plan: demonstrate directional egress (C-term moves +x on release, no collapse/leak,
finite energy) with `eject_demo.py`, at full length (directionality) and short length (full traverse).

## 2026-06-30 — D9 §1b feature assessment (validate-first outcome)

Per AGENTS.md §1b the approach is validate-first: implement features only as the baseline gap shows
they matter, one at a time, recording before/after. The baseline already matches the reference on the
quantitative metrics (dwell 1.01×, Rg 1.06×); the lone open gap is the residual soft-EV clash. Findings
per feature:

- **✅1 C-terminal mobility window — DEFERRED (incompatible).** Mass-0 freeze is forbidden with rigid
  `AllBonds` (OpenMM: "a constraint cannot involve a massless particle"), and more fundamentally topo
  extrudes by *diffusion of the whole mobile chain* (no explicit register translocation), so freezing
  the bulk would halt extrusion and pile the chain at the PTC. Doing it O'Brien's way needs explicit
  per-residue translocation — a mechanism topo deliberately replaced. Documented as a model-route
  difference, not a quick mask.
- **✅2 tRNA tether + orientation — IMPLEMENTED + measured (behind `trna_tether = yes`).** Full O'Brien
  linkage now in `add_trna_tether`: bond N–R + harmonic angles N–R–P, N–R–PU2 + improper N–R–P–PU2
  (periodic-harmonic) + backbone angle prev–N–R; CSP applies A-site beads in stages 1–2 and P-site in
  stage 3. Debug A/B (L=8, `csp_tether.ini` → `synth_out_tether/`): **stable** (max|PotE| 596.8 vs 42.8
  kJ/mol for the position spring, finite, **no dt-halving**), seed bond 3.810 Å, kinetics identical.
  **Clash NOT reduced:** min nascent–ribosome dist 3.51 Å (tether) vs 4.08 Å (position) — the tether
  bonds the C-terminus directly to the tRNA bead at the PTC, so it sits *closer* in, by design.
  Geometry (L=10 Rg): tether 0.658 nm (0.88× ref) vs position 0.799 nm (1.06× ref) vs ref 0.750 —
  both within ~15%; the tether's angles give a more compact, down-tunnel-aimed chain. → Kept as the
  more O'Brien-faithful restraint *option*; position restraint stays the validated default (it does not
  move the clash, so not promoted).
- **✅3 Restrain previous AA (L−1) — IMPLEMENTED (pairs with ✅2).** In stage 1 the L−1 residue is also
  tethered to the P-site (both ends of the new peptide bond pinned at the equilibrium PTC). Active only
  on the tether path; validated together with ✅2 above.
- **✅4 Ribosome L24 free loop — DEFERRED (out of scope for topo's ribosome).** topo loads the ribosome
  as pure mass-0 *scenery with no intra-ribosome bonded FF* ("no intra-ribosome forces are ever
  computed"). Freeing L24 42–59 (giving those beads mass) without their internal bonds/angles would let
  them be dragged by the chain with no restoring force → unphysical. O'Brien has the full ribosome FF
  (`combine_ribo_L24.prm`); topo deliberately omits it. Also, the observed clashes are with **23S rRNA**
  beads, not L24, so freeing L24 would not address them. Documented as a scenery-model limitation.
- **✅5 Placement 10° off-axis tilt — N/A (superseded).** With `equil_peptide_geometry`, the new residue
  is seeded directly at the optimal A-site target (`seed_point = a_target`), which already encodes the
  optimal off-axis bearing from the full O'Brien restraint solve. The legacy tilt would only apply on
  the (unused) anchor+buffer seed path. DIFFERENCES.md itself flags it "largely superseded by optimal-PTC
  targeting." No change made.

**Net D9 conclusion.** Of the 5 selected features, only ✅2/✅3 are applicable to topo's deliberate
simplifications; they are implemented, stable, and O'Brien-faithful but do **not** reduce the residual
clash (and slightly tighten C-term proximity). ✅1/✅4 are incompatible with topo's diffusion-extrusion /
scenery-only ribosome; ✅5 is superseded by equil-PTC seeding. → The residual 2.2–2.9 Å nascent–23S
interpenetration is a **model property** (O'Brien's deliberately-soft EV at full tunnel packing), not a
feature-fixable bug (§8 science finding): the chain otherwise meets the FINAL GOAL — full synthesis,
finite energy, correct tunnel threading (corr −0.926), no leak through the truncation, and clean +x
egress on release.

## 2026-06-30 — User note "Add to revise": the y-z tunnel wall + clash re-analysis

User left a note: *"No collapse to ribosome must be defined as no bead x < 0 (along the aligned
exit tunnel). But the tunnel wall [should be] defined by y-z coordinate as well."* → the chain
should be confined laterally (y-z), not only forward (x≥x0). Investigated thoroughly:

1. **The chain leaks laterally, not only forward.** At L=306 the in-tunnel chain bulges to a
   yz-radius of ~24 Å (p90 ~17 Å) off its centroid — nothing confines it radially (the x-wall only
   blocks −x). This is the lateral leak the note flags.
2. **But there is NO bead-free lumen to confine it to.** Grid-searching each x-slab for the y-z
   point farthest from every bead σ-surface gives a max bead-free radius of only **~0.1–2.7 Å for
   x = 8–86 Å** (a real opening appears only past x≈90 Å, the exit). The 23S rRNA beads (σ = 7.1 Å)
   so densely fill the truncated tunnel that **no channel wider than a single chain bead exists** —
   a lumen-confining radial wall is not geometrically realizable.
3. **"Within σ" is normal in this CG model — not a hard clash.** Recomputed: the **O'Brien
   reference** chain (L=10) sits at min center-distance **4.57 Å** to rRNA, itself *within* the
   σ = 7.1 Å shell (surf gap −2.5 Å). So both models' chains live inside the soft σ-shells; that is
   how O'Brien's deliberately-soft EV (ε = 0.000132 kcal/mol) tunnel works. The **real difference**
   is magnitude: topo's chain relaxes ~2 Å *deeper* (min center-dist 2.2–2.9 Å vs O'Brien's 4.57 Å)
   because topo keeps the **whole chain mobile**, so it settles into the soft beads, whereas O'Brien
   **freezes all-but-15 C-terminal residues + explicitly translocates**, holding the chain on the
   reference path. Energy stays finite in both (topo total ≤1.5e3 kJ/mol).

**Conclusion.** The lateral confinement the note asks for is the right instinct, but (a) it cannot be
a lumen-radial wall (no lumen exists in the CG truncated ribosome), and (b) the residual deeper
penetration is driven by topo's whole-chain mobility vs O'Brien's frozen-window+translocation — the
same root as the ✅1 incompatibility. Options: a pragmatic y-z *tube* wall around the exit axis
(curbs the 24 Å lateral bulge / void-leak but cannot remove the intrinsic σ-overlap), accept it as a
characterized CG-model property, or take on the O'Brien frozen-window + explicit-translocation
rebuild. → surfaced to the user (§8 model question) for direction before proceeding.

### Validation table (fill as runs complete)
| Run | path | L range | scale_factor | max\|PotE\| (kJ/mol) | seed bond (Å) | dt-halving? | min NC dist (Å) | dwell ratio | Rg ratio | notes |
|-----|------|---------|--------------|----------------------|---------------|-------------|-----------------|-------------|----------|-------|
| debug | equil-PTC + AllBonds | 1→8 | 216564650 | 42.78 | 3.810 | none | — | — | — | D1/D2 PASS (CPU) |
| baseline full | equil-PTC + AllBonds | 1→306 | 216564650 | 1.48e3 | 3.810 | none | 2.41 (eject) / 2.87 (synth) | 1.01× | 1.06× | D3/D4/D5/D6 PASS; D5b clash+egress open |
| P0CX28 full | equil-PTC + AllBonds | 1→106 | 216564650 | 241 | 3.810 | none | 2.29 (eject) | N/A (no ref) | R_g 1.345 nm (0.64× native) | D8: D5 PASS, wall PASS, threads −0.740, no leak; chain still in tunnel (egress N/A) |

## 2026-06-30 — D8 P0CX28 (L=1→106) complete

Full claude-fix run on P0CX28 (106 res, strength 2.5044). **D5 PASS** (worst 241 kJ/mol), **0
dt-halving**, seed bond 3.810 Å. **D5b wall PASS**; threads tunnel corr −0.740, **no leak** (0
beads x<0); R_g 1.345 nm (0.64× the 2.111 nm native — tunnel-confined, more compact than native,
as expected). D6 N/A (no O'Brien reference). Egress: the 106-mer is still entirely inside the
~100 Å tunnel (all beads x 8.6–35.9 Å) — too short to emerge, so full clearing is not observable
at this length (no leak/collapse). Residual 2.29 Å clash = same soft-EV model property as 4c5c.
Details in [`P0CX28/NOTES.md`](P0CX28/NOTES.md).

## 2026-06-30 — Egress demo (`eject_demo.py`, FINAL GOAL #3)

Extended free-MD ejection from the L=306 final structure (500000 steps ≈ 7.5 ns, restraint OFF,
tunnel wall ON), → `synth_out/ejection_long/`. MD finished in 223 s on GPU.
- **C-terminus x: 12.8 → 24.8 Å (net +12.0)** — moves OUT (+x), does not collapse back.
- **nascent CoM-x: 59.3 → 82.1 Å (net +22.8)**, linear slope +0.014 Å/frame, 55% of steps advance.
- min nascent–ribosome distance 2.19–2.84 Å (same soft-EV grazing).

**D5b now: wall PASS + egress PASS** (analyzer uses `ejection_long/`). **Only `clash` still FAIL**
(min 2.19 Å). FINAL GOAL status: #1 full synthesis ✅, #3 directional egress ✅, #4 no wall leak ✅;
**#2 (no clash) is the lone open item** = the residual soft-EV interpenetration.

**Clash assessment (candidate §8 model finding).** The sub-3 Å contacts are nascent residues threading
the 23S rRNA tunnel (bead radius 7.1 Å) at 2.2–2.9 Å — O'Brien's EV is deliberately soft
(ε = 0.000132 kcal/mol), so beads interpenetrate when other forces (folding contacts, the growing
chain) push them, while total energy stays finite (≤1.5e3 kJ/mol). At the only length with a reference
(L=10) ours (3.36 Å) and O'Brien (4.57 Å) are both clash-free; the sub-3 Å contacts are a *full-length
tunnel-packing* effect. The §1b feature that could reduce it (✅1 mobility window — freeze extruded
residues away from beads) is **incompatible with topo's diffusion-extrusion** (see above); the
compatible features (✅2 orientation, ✅3 prev-AA, ✅4 L24, ✅5 tilt) all target the C-terminus / loop,
not the mid-chain threading where the clash lives. → will implement the headline ✅2 to measure, then
report the residual as a model property if unmoved.
