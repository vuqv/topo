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
- [~] **D5b** Ejection: **wall not penetrated PASS** (analyzer wall bug fixed: 8.71 Å not 10.46 Å);
      **clash + net-egress OPEN** — min NC dist 2.41 Å (eject) / 2.87 Å (synth, 6/306 res); 20k-step
      ejection too short to show a 306-mer clearing. → addressed by §1b features + extended ejection demo.
- [x] **D6** Quantitative match (4c5c): in-vivo dwell ratio **1.01×**, final R_g **1.06×** (both well
      within tolerance); fix path reproduces ref ≥ as well as legacy 12/13 while being more faithful
      (rigid AllBonds, no guard). in-silico ns ratio 0.02× = 50× scale_factor (documented, scale-indep).
- [ ] **D7** Documented: `README.md` reproduce commands; `NOTES.md` table; add to `tutorials/README.md`.

## Phase 2 — P0CX28 (only after 4c5c D0–D7)
- [ ] **D8** Repeat D0–D5/D5b/D7 in `P0CX28/` (L=1→106). D6 N/A (no reference) → internal consistency.

## Phase 3 — §1b O'Brien-consistency features (VALIDATE-FIRST, one at a time)
- [ ] **D9** After D6 baseline gap established, implement selected features only as comparison shows
      they matter; re-validate dwell/geometry after each; record before/after in NOTES.md.
  - [ ] ✅1 C-terminal mobility window (freeze all but last N=15 nascent residues; mass-0) — `run_length`
  - [ ] ✅2 tRNA tether + orientation (bond + 2 angles + improper; A-site st1–2, P-site st3) — `add_trna_tether` + protocol
  - [ ] ✅3 Restrain previous AA (L−1) at P-site; use P/PU2 sub-beads — `add_cterm_restraint`
  - [ ] ✅4 Ribosome L24 free loop (residues 42–59 movable) — rigid-scenery build
  - [ ] ✅5 Placement 10° off-axis tilt in seed direction — `seed_positions`
  - NOT in scope: ribosome-traffic correction, post-only tunnel wall.

## Stop conditions (pause + ask): see AGENTS.md §8.
