# TASK.md — Tutorial 15 (`claude_fix`) progress checklist

Single source of truth for progress. Built from AGENTS.md §5 (Definition of Done) + §1b
(consistency features). Tick `[x]` with a one-line result as each item finishes; `[ ]` while open.
Full decisions/numbers live in `NOTES.md`.

## FINAL GOAL
Stable continuous synthesis: (1) full length (4c5c 1→306, then P0CX28 1→106) with finite energy;
(2) no clash with any ribosome bead; (3) clean +x tunnel ejection clearing the ribosome;
(4) without passing through the truncated-ribosome tunnel wall.

---

## Phase 0 — setup
- [x] **P0a** Read AGENTS.md, DIFFERENCES.md, core.py, protocol.py, ribosome.py. — done
- [x] **P0b** Verify env (bioenv: openmm 8.2, scipy, MDAnalysis, Tesla T4 GPU). — done
- [x] **P0c** Checkpoint commit before any code edits. — done (see git log)
- [x] **P0d** Create TASK.md + NOTES.md. — this file + NOTES.md

## Phase 1 — Definition of Done (4c5c)
- [x] **D0** Scaffold: 4c5c + P0CX28 inputs present; `analyze_validation.py` present; raw inputs
      untouched. *(P0CX28 needs its own analyze_validation.py copy for D8 — deferred to Phase 2)*
- [x] **D1** Configs: `csp.ini` (debug, L=1→8) + `csp_val.ini` (full, L=1→306) written; both set
      `equil_peptide_geometry=yes`+`constraints=AllBonds` (no `=no`/`=None`); banner confirms fix path;
      kinetics (time_stage_1/2, scale_factor, mrna, trans_times) match O'Brien reference.
- [x] **D2** Fix active (debug): optimal_ptc_targets printed (|A−P|=0.3810 nm); seed peptide bond
      = **3.810 Å**; max|PotE| = **42.78 kJ/mol** (2.3e7× under limit); **0 dt-halving lines**.
- [x] **D3** Run completes: `topo-csp -f csp_val.ini` → "Done. Synthesized 1 → 306", exit clean,
      zero dt-halving lines across 919 stages.
- [x] **D4** Outputs: per-stage trajectories + `synth_out/dwell_times.dat` (306 rows) + ejection/.
- [x] **D5** Physically sane: worst max|PotE| = 1.48e3 kJ/mol ≪ 1e12; chain threads tunnel
      (corr(residue,x)=−0.926, monotonic +x egress); no collapse into ribosome (no bead x<0).
- [~] **D5b** Ejection: **wall PASS** (analyzer wall fixed 8.71 Å) + **egress PASS** (`eject_demo.py`
      ejection_long: C-term +12.0 Å, CoM +22.8 Å along +x over 7.5 ns). **clash OPEN** — min 2.19 Å
      (soft-EV grazing of 23S tunnel beads; candidate §8 model finding, see NOTES).
- [x] **D6** Quantitative match (4c5c): in-vivo dwell ratio **1.01×**, final R_g **1.06×** (both well
      within tolerance); fix path reproduces ref ≥ as well as legacy 12/13 while being more faithful
      (rigid AllBonds, no guard). in-silico ns ratio 0.02× = 50× scale_factor (documented, scale-indep).
- [x] **D7** Documented: `README.md` (reproduce + results + §1b + clash finding); `NOTES.md` log;
      `tutorials/README.md` rows 14 + 15 added.

## Phase 2 — P0CX28 (only after 4c5c D0–D7)
- [x] **D8** P0CX28 L=1→106 on the claude-fix path: D2 fix active (seed 3.810 Å), D3 completes
      (0 dt-halving), D5 PASS (worst 241 kJ/mol), D5b wall PASS, threads corr −0.740, no leak,
      R_g 1.345 nm (0.64× native, tunnel-confined). D6 N/A (no reference). Chain still inside the
      ~100 Å tunnel (egress not observable at 106 res). Docs: `P0CX28/NOTES.md`.

## Phase 3 — §1b O'Brien-consistency features (VALIDATE-FIRST, one at a time)
- [x] **D9** Validate-first §1b assessment done; before/after recorded in NOTES. Baseline already
      matches ref (dwell 1.01×, Rg 1.06×); the lone gap (clash) is a soft-EV model property.
  - [x] ✅1 mobility window — DEFERRED: mass-0 forbidden with AllBonds + breaks topo's diffusion-extrusion
        (needs explicit translocation). Documented model-route difference.
  - [x] ✅2 tRNA tether + orientation — IMPLEMENTED behind `trna_tether=yes` (bond+2 angles+improper,
        A-site st1–2/P-site st3). Stable, kinetics identical; does NOT reduce clash (tether sits closer
        to PTC by design). Kept as O'Brien-faithful option; position restraint stays default.
  - [x] ✅3 Restrain previous AA (L−1) at P-site in stage 1 — IMPLEMENTED (pairs with ✅2).
  - [x] ✅4 L24 free loop — DEFERRED: topo's scenery-only ribosome has no intra-ribosome FF; clashes are
        23S not L24. Documented scenery-model limitation.
  - [x] ✅5 10° tilt — N/A: superseded by equil-PTC seeding (seed = optimal A-target). No change.
  - NOT in scope: ribosome-traffic correction, post-only tunnel wall.

## Phase 4 — Ribosome CG-coordinate fix (user direction) + force-field audit
- [x] **R1** Confirmed inconsistency: topo R bead ≠ O'Brien C5′ (P-anchor off 3.6 Å); uniform 0.71 nm
      radii vs O'Brien per-type. Added `load_obrien_ribosome(.cor/.psf/.prm)` + `load_ribosome_auto`.
- [x] **R2** Re-validated 4c5c + P0CX28 on O'Brien's ribosome: D5/D6 PASS, R_g fidelity ↑ (0.96×),
      surface penetration shallower, no leak; clash persists (weaker EV) → y-z wall confirmed needed.
- [x] **R3** Force-field comparison `TOPO_OBrien_params_compare.md` (angles/dihedrals/electrostatics/
      12-10-6 Gō match exactly; bond k 2× soft = latent unit bug; ribosome fixed).
- [x] **R4→superseded** y-z lumen wall NOT needed — the root cause was the NC↔ribosome nonbonded
      interaction (see R5).
- [x] **R5 NC↔ribosome EV fix (user-directed, all 4 points)** — `append_ribosome` now O'Brien-
      consistent: 12-10-6 form, sum rule R=Rmin/2_i+Rmin/2_j, matched ε, nascent per-AA `SA..SY`
      Rmin/2 (Option B); dt-halving guard made NaN-robust. **RESOLVES the clash on both proteins:**
      4c5c hard clashes 36→0 (D5b all PASS, extrusion x→179); P0CX28 18→0 (D5b all PASS, threads
      −0.969, R_g 2.306≈native, x→98). FINAL GOAL met. Docs: `TOPO_OBrien_NCribosome_nonbonded_compare.md`.

## Stop conditions (pause + ask): see AGENTS.md §8.
