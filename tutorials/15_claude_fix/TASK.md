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
- [ ] **D0** Scaffold: inputs copied (4c5c + P0CX28); `analyze_validation.py` present;
      raw inputs untouched. *(inputs present; need csp.ini + csp_val.ini)*
- [ ] **D1** Configs: `csp.ini` (debug) + `csp_val.ini` (full) exist, use the fix path
      (equil PTC geometry + AllBonds, no `=no`/`=None` pinning); banner confirms; kinetics match ref.
- [ ] **D2** Fix active: debug run prints `optimal_ptc_targets`; seed peptide bond ≈ 3.79–3.81 Å;
      dt-halving guard never fires. Record targets + max|PotE|.
- [ ] **D3** Run completes: `topo-csp -f csp_val.ini` exit 0, full length 1→306.
- [ ] **D4** Outputs: trajectory + per-residue `dwell_times.dat` under `synth_out/`.
- [ ] **D5** Physically sane: no stage |PotE| ≳ 1e9; monotonic-ish +x egress; no collapse.
- [ ] **D5b** Clean ejection: chain diffuses +x, clears ribosome, no wall penetration / bead overlap.
      Record min nascent–ribosome distance + CoM-x vs frame.
- [ ] **D6** Quantitative match (4c5c only) vs `../12_auto/.../output/` (ref = L 1→10):
      length, per-codon dwell / total time (~2×), final R_g. Also compare fix vs legacy 12/13 path.
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
