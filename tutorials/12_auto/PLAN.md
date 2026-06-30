# PLAN — reproduce O'Brien CSP on 4c5c with `topo-csp`

The design and status for this tutorial. It is the *reproduction* counterpart of
`tutorials/10_csp_obrien` (the fast machinery demo). Full decisions and the validation
table are in [`NOTES.md`](NOTES.md); the elongation mechanics are in
[`STAGES.md`](STAGES.md); the Tutorial-10 post-mortem is in
[`WHY_10_FAILS.md`](WHY_10_FAILS.md).

## Approach

1. **Map** `continuous_synthesis/input/cont_synth_ecoli.cntrl` → `csp.ini` /
   `csp_val.ini` (kinetics, ribosome, traffic, initiation) — see `NOTES.md` §1.
2. **Prepare inputs** as a topo-native CG build (no CHARMM ingest): `4c5c` CG model +
   `domain.yaml` (mapped from `protein_cg_model/`) + topo-ready truncated ribosome +
   the reference's exact mRNA + Fluitt codon-time table — see `NOTES.md` §2.
3. **Debug→scale**: get a fast end-to-end run green (`csp.ini`, L=1→6, clamped), then
   run the production validation (`csp_val.ini`, L=1→10).
4. **Validate** energies (D5), ejection (D5b) and dwell/length/Rg vs the reference
   (D6) with `analyze_validation.py` + `eject_demo.py`.

## Reused (not re-implemented)

- `topo.csp` (`csp.py`, `kinetics.py`) + `topo-csp` console script — the runner.
- `topo.csp.core.run_length` and the rigid-ribosome / tunnel-wall /
  build-once-subset machinery — the per-stage MD.

## Done

- [x] Config mapping (`csp.ini` debug, `csp_val.ini` validation) — D1.
- [x] Inputs prepared (topo CG 4c5c + ribosome + identical mRNA + Fluitt table) — D2.
- [x] End-to-end debug run + production validation (L=1→10) — D3/D4.
- [x] **Per-stage stability guard** in `run_length` (dt-halving on divergence) — fixes
      the blow-up that breaks D5; verified worst stage 524 kJ/mol — D5.
- [x] `dwell_times.dat` per-residue dwell log (topo analog of reference `1.out`) — D4.
- [x] Ejection verified: clean +x egress, no wall penetration / collapse — D5b.
- [x] Quantitative validation vs reference (length, dwell 1.01×, Rg 1.06×) — D6.
- [x] Docs (`README.md`, `NOTES.md`, `WHY_10_FAILS.md`) + tutorials index — D7/D8.

## Remaining (out of scope here; see `NOTES.md` §7)

- [ ] CHARMM `.psf/.top/.prm/.cor` ingestion to run O'Brien's exact `go_bt` force field.
- [ ] Integrate with rigid `AllBonds` (remove only the new bond's constraint) as the
      clean alternative to the dt-halving guard.
- [ ] Full-length (306 aa) multi-trajectory ensemble; working `ribosome_traffic` binary.
