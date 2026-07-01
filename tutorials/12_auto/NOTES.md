# NOTES — decisions, deviations, and the validation table (Tutorial 12_auto)

Companion to [`Goal.md`](Goal.md). Records how the reference
`cont_synth_ecoli.cntrl` was mapped to `csp.ini`, the one code change needed to make
the run physically sane, the deliberate deviations, and the quantitative validation
against `continuous_synthesis/output/`.

Run on: NVIDIA GPU (CUDA), OpenMM 8.2, `topo` editable install. Profiles:
`csp.ini` (debug) and `csp_val.ini` (validation, the numbers below).

---

## 1. Config mapping: `cont_synth_ecoli.cntrl` → `csp_val.ini`

| reference key (`cont_synth_ecoli.cntrl`) | value | topo key (`csp_val.ini`) | value | note |
|---|---|---|---|---|
| `start_nascent_chain_length` | 1 | `L0` | 1 | |
| `total_nascent_chain_length` | 10 | `L_max` | 10 | reference only goes to 10 |
| `temp_prod` | 310 | `ref_t` | 310 | |
| `Time step` | 0.015 ps | `dt` | 0.015 | |
| `scale_factor` | 4331293.0 | `scale_factor` | 4331293 | production value |
| `time_stage_1` | 0.00034 s | `time_stage_1` | 0.000340 | peptidyl transfer mean |
| `time_stage_2` | 0.004201 s | `time_stage_2` | 0.004201 | translocation mean |
| `mrna_seq` | `setup/4c5c_mrna_sequence_fast.txt` | `mrna` | `4c5c_mrna.txt` | **byte-identical** |
| `trans_times` | `setup/Fluitt_ecoli_...` | `trans_times` | `trans_times.txt` | Fluitt *E. coli* table |
| `ribosome_traffic` | 1 | `ribosome_traffic` | **no** | binary unavailable — see §3 |
| `initiation_rate` | 0.083333 | `initiation_rate` | 0.083333 | (traffic only) |
| `x_eject` | 60 Å | (analysis) | — | egress target, see D5b |
| `File save steps` | 5000 | `nstout` | 2000 | finer output |
| `prot_psf/top/param`, `ribo_*`, `starting_strucs` | CHARMM `setup/*` | `pdb_file` + `ribosome` + `domain_def` | topo CG build | see §2 |

`random_seed = 20240629` (topo) is our own reproducible FPT seed; the reference draws
a fresh per-length seed (e.g. 371754059 for L=1), so the *sampled* dwell times differ
between the two runs while the per-codon *means* are identical (same mRNA + table).

---

## 2. Inputs (D2): topo CG build vs CHARMM ingestion

The reference system is CHARMM `.psf/.top/.prm/.cor`. `topo.csp` does **not** ingest
CHARMM (the known gap in `10_csp_obrien/TASK.md`), so — as Goal §2 permits — we use a
topo-native CG build from the all-atom PDB plus a topo-ready truncated ribosome:

- **Protein:** `4c5c_model_clean.pdb` (306 CA) built by topo's contact builder, with
  `domain.yaml` carrying the 3-domain map + Go-scale nscales mapped from
  `protein_cg_model/domain_def.dat` (4 cg_simtk segments → 3 structural domains; the
  collapse is information-preserving — see the header of `domain.yaml`).
- **Ribosome:** `ribosome_trunc.pdb` (4,576 CG beads, X-aligned) — the topo-ready
  equivalent of `setup/50S_tRNA_cg_truncated.*`, supplying the P-/A-anchors
  (PtR/AtR residue-76 R beads) and the rigid scenery.
- **mRNA:** `4c5c_mrna.txt` is **byte-identical** to the reference fast mRNA, so the
  codon sequence (and therefore the per-residue kinetics) is exactly the reference's.

This is the principal deviation from the reference: a topo-rebuilt CG model rather
than O'Brien's exact CHARMM `bt`/`go_bt` force field. Energies are therefore in
topo's CG units, not directly comparable bead-for-bead; the validation (D6) compares
*derived* quantities (dwell times, length, Rg) that do not depend on the absolute
force-field scale.

---

## 3. Deliberate deviations (and why)

1. **`ribosome_traffic = no`** (reference `= 1`). O'Brien's compiled `ribosome_traffic`
   binary is not on the PATH here; `topo.csp.kinetics.ribosome_traffic_times` degrades
   gracefully to `real == intrinsic` (no upstream-queue correction). Impact on D6: the
   traffic correction only adds to stage-2 when `real − intrinsic > 0`; at this short
   length its effect on the summed in-vivo time is ≲1 % (our 0.256 s vs reference-mean
   0.253 s, §5), so the comparison is unaffected.

2. **`max_steps_per_stage = 200000`** (reference uncapped; its longest stage is
   1,124,463 steps). The cap keeps the rare ~10⁶-step stage tractable on one GPU. It
   truncates **integration steps only**, never the sampled dwell **times in seconds**
   (`dwell_times.dat` records the untruncated seconds), so the D6 dwell comparison in
   seconds is exact; only the *in-silico ns actually integrated* is shortened.

3. **Flexible bonds (`constraints = None`)** vs the reference's rigid `AllBonds`. topo
   must seed the new residue ~1 nm from its bond partner (A-site delivery), which a
   rigid distance constraint cannot represent; a harmonic bond absorbs the stretch and
   the minimizer relaxes it. The cost of flexible bonds at 15 fs is the stability bug
   in §4.

4. **3-stage mechanics are mapped, not literal** (inherited from the port): A→P
   translocation is the restraint-target switch; the peptide bond is present from
   stage 1; explicit A/P tRNA bonded geometry is not modeled. The *timing* is faithful.

---

## 4. Code change: per-stage stability guard (D5)

**Symptom.** In the first validation run, stage `L_010/stage_2` diverged to
PotE ≈ 2×10¹³ kJ/mol and stayed there — fatal for D5 (and exactly the self-recovering
stage-1/2 blow-up logged in `10_csp_obrien/OBSERVATIONS.md` #1).

**Root cause (measured, deterministic — not a random clash).** Building the L=10
stage in isolation, the seed minimizes cleanly to ~243 kJ/mol (new bond relaxes to
0.374 nm, no steric overlap). The divergence is purely an **integrator instability**:
at `dt = 0.015 ps` with flexible bonds, the L=10 configuration diverged on **all 3**
velocity seeds (max PotE 7.8e12 / 7.4e9 / 6.5e12), while L=9 was stable (516) and L=10
at **half the timestep** (`dt = 0.0075 ps`) was stable on all seeds (257 / 264 / 263).
A new native (Go) contact at L=10 forms a stiff well whose vibrational period is too
short for 15 fs. The reference avoids this with rigid `AllBonds` constraints.

**Fix** (`topo/csp/elongate.py`, `run_length`): a per-stage stability guard.
The stage is run in chunks while tracking the **maximum** |PotE| (a mid-run blow-up
can cool back under the limit by the final frame yet still corrupt those frames, so
the final-state energy alone is *not* a sufficient test). If a stage's max |PotE|
exceeds `STABILITY_POTE_LIMIT_KJ = 1e9` (sane CG stages are 10¹–10⁴), it is re-run
with **a halved timestep and double the steps** — which keeps the physical dwell time
`n_steps · dt` identical while stabilising the integration — up to
`STABILITY_MAX_ATTEMPTS = 6` halvings. The fresh `setup_simulation` + `attach_reporters`
truncate the per-stage output, so a successful retry cleanly overwrites the aborted one.
The common case runs once at the configured `dt = 0.015 ps` (matching the reference);
only the rare unstable stage deviates, and it logs a `[stability] ...` line.

**Result in the validation run.** Exactly two stages auto-stabilised — `L_010/stage_1`
(952 → 1904 steps @ 0.0075 ps) and `L_010/stage_2` (113936 → 227872 steps @ 0.0075 ps).
Their `dwell_times.dat` rows still record the *physical* dwell (seconds) and the
nominal step counts; the doubled steps are an internal integration detail.

A secondary fix (D4): `run_continuous_synthesis` now writes
`synth_out/dwell_times.dat` (per-residue codon, sampled stage dwell times in s, the
in-silico ns, and step counts) — the topo analog of the reference `output/1.out`.

---

## 5. Validation table (D3–D6) — `csp_val.ini`, L = 1 → 10

Reproduce with: `python analyze_validation.py` (D5/D6) and `python eject_demo.py 2000000` (D5b).

### D3 — run completes
`topo-csp -f csp_val.ini` → exit 0; synthesized 1 → 10 (10 residues × 3 stages + ejection).

### D4 — outputs
Per-residue/-stage trajectories under `synth_out/L_<L>/stage_<s>/`, plus
`synth_out/dwell_times.dat` and `synth_out/ejection/`. Layout mirrors the reference
`output/traj/1/`.

### D5 — finite energies (no ≳10¹² blow-up)
Scanned **31** stage logs. Worst max|PotE| = **524 kJ/mol** (`L_009/stage_2`). No stage
exceeds 10⁹. **PASS.** Two stages reached the guard and were re-run at half timestep.

### D5b — clean ejection through the tunnel
In-run ejection (50 k steps) + extended demo (`ejection_long`, 2×10⁶ steps ≈ 30 ns):

| quantity | in-run (50 k) | extended (2 M) |
|---|---|---|
| nascent CoM-x | 24.5 → 23.4 Å | **110.9 → 1227.5 Å** (max 1522 Å), net **+1117 Å** |
| min nascent x (wall x0 = 10.5 Å) | 13.9 Å | 12.2 Å — wall **never penetrated** |
| min nascent–ribosome distance | 3.0 Å | 11.9 Å (frame 0), grows as it leaves — **no clash / no collapse** |
| Rg | — | 0.83 → 0.54 nm (stays compact, not exploded) |

The released chain diffuses out along **+x** (CoM 11 nm → 122 nm), **fully clears the
ribosome** (~10 nm extent), and at no frame penetrates the tunnel wall or overlaps a
ribosome bead. **PASS.** (50 k steps is too short for a 10-bead chain to traverse the
~4 nm tunnel; the extended run shows the egress.)

### D6 — quantitative consistency vs `continuous_synthesis/output/`

(a) **Length:** ours max L = 10, reference max L = 10. ✓

(b) **Per-residue in-vivo dwell (s)** — ours sampled `t_total` vs reference per-residue
mean (`mean_pt + mean_tl + mean_tr` from `output/1.out`):

| L | codon | ours t_total (s) | ref mean (s) |
|---|---|---|---|
| 1 | AUG | 2.757e-2 | 2.129e-2 |
| 2 | ACU | 3.859e-2 | 3.863e-2 |
| 3 | GAU | 3.809e-2 | 3.714e-2 |
| 4 | AAA | 4.862e-2 | 4.738e-2 |
| 5 | AUU | 1.955e-2 | 1.801e-2 |
| 6 | GCU | 1.303e-2 | 1.366e-2 |
| 7 | GUU | 1.754e-2 | 2.036e-2 |
| 8 | CUG | 1.754e-2 | 1.913e-2 |
| 9 | CUG | 1.754e-2 | 1.878e-2 |
| 10 | GGU | 1.754e-2 | 1.882e-2 |

- **Total in-vivo dwell:** ours (sampled) **0.2556 s** vs reference (mean) **0.2532 s**
  → **ratio 1.01×** (tolerance was ~2×). ✓
- **Total in-silico dwell:** ours **66.8 ns** vs reference **76.4 ns** → **0.87×** (both
  sample exponentials about the same means with different seeds). ✓

(c) **Final nascent radius of gyration (L=10):** ours **0.798 nm** vs reference
**0.750 nm** → **ratio 1.06×**. ✓

> Note: `dwell_times.dat` `t_total` = `intrinsic[L]` is the *next* codon's mean time —
> this is O'Brien's stage-3 ("wait for the next aminoacyl-tRNA") indexing, faithfully
> reproduced (`topo/csp/kinetics.py` docstring). CUG and GGU both have time 0.017542 s
> in the Fluitt table, so the repeated values for L=7–10 are correct, not a bug.

---

## 6. Status vs Goal §3

D0 scaffold ✓ · D1 config ✓ · D2 inputs ✓ · D3 run ✓ · D4 outputs ✓ · D5 finite ✓ ·
D5b ejection ✓ · D6 quantitative ✓ · D7 docs ✓ (this file + `README.md` + tutorials
index) · D8 post-mortem ✓ (`WHY_10_FAILS.md`).

## 7. Open questions / follow-ups

- **CHARMM ingestion** of the exact `setup/*` system would let `topo.csp` run O'Brien's
  precise `go_bt` force field (closes the last gap in `10_csp_obrien/TASK.md`); the
  topo CG rebuild is a faithful but not bit-identical stand-in.
- **Rigid bonds.** The stability guard halves the timestep where needed; a cleaner
  long-term fix is to integrate with `AllBonds` constraints (as the reference does),
  removing only the newly added bond's constraint during A-site delivery.
- **Full length (306 aa) + ensemble.** This run matches the reference's 10-residue
  length; a production study needs the full chain and many `random_seed`s.
- **Ribosome traffic binary.** Building/locating `ribosome_traffic` would let us match
  the reference's traffic-on kinetics exactly (small effect at L≤10).
