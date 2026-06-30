# SUMMARY — Tutorial 12_auto (autonomous task)

**Goal** ([`Goal.md`](Goal.md)): reproduce O'Brien's `continuous_synthesis_v6.py`
protocol *in topo-style* (an `csp.ini` + a `topo-csp` run) on the **4c5c** inputs, get
a physically sane, quantitatively consistent result vs the bundled reference run, and
explain why the earlier `tutorials/10_csp_obrien` does not work.

**Outcome: achieved — all Definition-of-Done criteria (D0–D8) pass.**

---

## What I did

1. **Oriented** in the existing scaffold (a prior session had produced `csp.ini`,
   `csp_val.ini`, inputs, and an L=1→10 run) and verified it against the reference.
   The reference `cont_synth_ecoli.cntrl` itself synthesizes only **1 → 10**
   (`total_nascent_chain_length = 10`), so that is the validation length.

2. **Found the blocking bug.** The prior validation run had a real **energy blow-up**:
   `L_010/stage_2` reached PotE ≈ 2×10¹³ kJ/mol (fails D5; the same class of bug as
   `10_csp_obrien/OBSERVATIONS.md` #1).

3. **Root-caused it by direct measurement** (not the guess in OBSERVATIONS #1):
   - The seed minimizes **cleanly** (≈243 kJ/mol, no steric overlap) — so it is *not*
     a seed-on-bead near-singularity.
   - It is a **deterministic integrator instability**: at `dt = 0.015 ps` with flexible
     bonds, L=10 diverges on **all** velocity seeds (max PotE 6–8×10¹²), L=9 is stable
     (516), and L=10 at **half the timestep** is stable on all seeds (≈260). A new
     native (Go) contact at L=10 forms a stiff well too fast for 15 fs. The reference
     avoids this with rigid `AllBonds` constraints; topo needs flexible bonds to seed
     the far-placed A-site residue.

4. **Fixed it at the source** (`topo/csp/elongate.py`, `run_length`): a
   **per-stage stability guard** that runs each stage in chunks tracking the *maximum*
   |PotE|, and if a stage diverges (> 10⁹ kJ/mol) re-runs it with a **halved timestep
   and double the steps** — which keeps the physical dwell time `n_steps·dt` identical
   while stabilising the integration (up to 6 halvings). The common case is unchanged
   (one pass at `dt = 0.015 ps`). Judging from the *max* (not final) energy matters:
   a diverged stage can cool back under threshold by its last frame yet still have
   corrupted frames — which is exactly how Tutorial 10's blow-ups "self-recover".

5. **Added the missing per-residue dwell-time log** (`topo/csp/protocol.py`):
   `synth_out/dwell_times.dat` (codon, sampled stage dwell times in s, in-silico ns,
   step counts) — the topo analog of the reference `output/1.out` (D4).

6. **Re-ran the production validation** (`csp_val.ini`, L=1→10, production
   `scale_factor = 4331293`). Exactly two L=10 stages auto-stabilised at the halved
   timestep; **no stage exceeds 524 kJ/mol** (D5 pass).

7. **Verified ejection (D5b)** with a dedicated long run (`eject_demo.py`,
   2×10⁶ steps): the released chain diffuses out **+x** (CoM 11 nm → 122 nm), **never
   penetrates the tunnel wall**, **never clashes with / collapses into** the ribosome
   (min distance grows from 11.9 Å), and fully clears it.

8. **Validated quantitatively (D6)** with `analyze_validation.py`: length 10 = 10;
   total in-vivo dwell **0.256 s vs 0.253 s (1.01×)**; final Rg **0.80 nm vs 0.75 nm
   (1.06×)**.

9. **Documented everything**: `README.md` (repro commands), `NOTES.md` (mapping,
   deviations, validation table), `STAGES.md` + `PLAN.md` (scaffold), `WHY_10_FAILS.md`
   (D8 post-mortem), and added a Tutorial-12 row to `tutorials/README.md`.

## Files changed / added

- **Package code** (the substantive fixes):
  `topo/csp/elongate.py` — per-stage stability guard (dt-halving on
  divergence) + constants. `topo/csp/protocol.py` — write `dwell_times.dat`.
- **Tutorial** (`tutorials/12_auto/`): `README.md`, `NOTES.md`, `WHY_10_FAILS.md`,
  `PLAN.md`, `STAGES.md`, `SUMMARY.md`; `analyze_validation.py`, `eject_demo.py`;
  outputs `synth_out/` (+ `ejection_long/`), `synth_out_debug/`, run logs.
  `tutorials/README.md` index updated. No raw input was overwritten.

## Key decisions (rationale in `NOTES.md`)

- **topo CG build, not CHARMM ingest** (CHARMM ingestion is the open gap in
  `10_csp_obrien/TASK.md`); the mRNA used is **byte-identical** to the reference, so
  the kinetics are exactly the reference's.
- **`ribosome_traffic = no`** (the binary is unavailable; effect ≲1 % at L≤10).
- **`max_steps_per_stage = 200000`** caps MD steps only, not the dwell times in
  seconds, so the D6 comparison is unaffected.
- **Validation at L=10** (matches the reference), not the full 306 residues.

## Open questions / follow-ups

- **CHARMM `.psf/.top/.prm/.cor` ingestion** would let `topo.csp` run O'Brien's exact
  `go_bt` force field instead of a topo CG rebuild (last gap from Tutorial 10).
- **Rigid `AllBonds`** (removing only the new bond's constraint during A-site delivery)
  is the cleaner long-term alternative to the dt-halving guard.
- **Full-length (306 aa) ensemble** with many `random_seed`s and a working
  `ribosome_traffic` binary for a production study.
- The dt-halving fix also protects **Tutorial 10** (shared `run_length`); a confirming
  full-length re-run of Tutorial 10 would close the loop on its 5/306 blow-ups.
