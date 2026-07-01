# Tutorial 12 (auto) — O'Brien Continuous Synthesis reproduced in topo-style (4c5c)

This tutorial **reproduces O'Brien's `continuous_synthesis_v6.py` protocol** on the
**4c5c** protein using the topo runner `topo-csp` driven by an `csp.ini` — *not* by
editing the legacy script. It maps the reference control file
`continuous_synthesis/input/cont_synth_ecoli.cntrl` onto topo names, runs the
per-codon, three-stage elongation cycle on the 4c5c system, and validates the result
against the bundled reference run in `continuous_synthesis/output/`.

It is the autonomous worked example asked for in [`Goal.md`](Goal.md); the decisions,
deviations and the full validation table are in [`NOTES.md`](NOTES.md), the elongation
mechanics are in [`STAGES.md`](STAGES.md), and a post-mortem on why the earlier
`tutorials/10_csp_obrien` demo does not reproduce O'Brien is in
[`WHY_10_FAILS.md`](WHY_10_FAILS.md).

> **What was fixed to make this work.** The CSP port had a latent numerical bug: a
> handful of stages diverge (PotE → ~10¹³ kJ/mol) at the 15 fs timestep with flexible
> bonds when a new native contact forms a stiff Go well (see
> `10_csp_obrien/OBSERVATIONS.md` #1). `topo.csp.core.run_length` now has
> a **per-stage stability guard**: it detects a diverging stage and re-runs it with a
> halved timestep and proportionally more steps — preserving the physical dwell time
> (`dwell = n_steps · dt`) while stabilising the integration. See
> [`WHY_10_FAILS.md`](WHY_10_FAILS.md).

---

## Files in this folder

### Raw inputs (never overwritten)
| File | Role |
|------|------|
| `4c5c_model_clean.pdb` | All-atom 4c5c structure (306 residues) = the nascent chain at full length. |
| `protein_cg_model/cg.cntrl`, `protein_cg_model/domain_def.dat` | O'Brien CG model control + domain/Go-scale definitions (source for `domain.yaml`). |
| `continuous_synthesis/` | The **reference** run: `input/cont_synth_ecoli.cntrl` + `input/setup/*` (CHARMM CG system, mRNA, Fluitt codon times) and `output/` (reference trajectory + `info.log` + `output/1.out`). |

### Prepared inputs for `topo-csp`
| File | Role | Provenance |
|------|------|-----------|
| `domain.yaml` | 3-domain map + contact nscale for the topo CG build. | Mapped from `protein_cg_model/domain_def.dat` (see header in the file). |
| `4c5c_mrna.txt` | mRNA, one codon/residue (+stop). | **Identical** to `continuous_synthesis/input/setup/4c5c_mrna_sequence_fast.txt`. |
| `trans_times.txt` | Codon → mean in-vivo time (s). | Fluitt *E. coli* table (= `setup/Fluitt_ecoli_trans_time_310K_avg_16.5.txt`). |
| `4c5c_model_clean_stride.dat` | Precomputed STRIDE (skips running STRIDE). | From the all-atom PDB. |
| `ribosome_trunc.pdb` | Truncated CG 50S + tRNAs (4,576 beads), X-aligned. | topo-ready equivalent of `setup/50S_tRNA_cg_truncated.*`. |
| `csp.ini` | **Debug** profile (fast end-to-end check). | maps `cont_synth_ecoli.cntrl`. |
| `csp_val.ini` | **Validation** profile (production kinetics, L=1→10). | maps `cont_synth_ecoli.cntrl`. |

### Outputs (written here; safe to delete/regenerate)
`synth_out/` (validation run), `synth_out_debug/` (debug run), `synth_out/ejection_long/`
(extended ejection, D5b), `val_run.log`, `debug_run.log`, `eject_demo.log`,
`analyze_validation.py`, `eject_demo.py`.

---

## Reproduce it

All commands run **from this folder** (paths in the `.ini` files are relative):

```bash
cd tutorials/12_auto

# 1. Fast debug run (~1 min on GPU): L=1→6, scale_factor 50× the production value,
#    step counts clamped. Sanity-checks the whole pipeline end to end.
topo-csp -f csp.ini                       # -> synth_out_debug/

# 2. Validation run (a few minutes on GPU): L=1→10 with PRODUCTION kinetics,
#    matching the reference cont_synth_ecoli.cntrl. Two L=10 stages auto-stabilise
#    at a halved timestep (printed as "[stability] ..." lines).
topo-csp -f csp_val.ini                   # -> synth_out/

# 3. Validate against the bundled reference output (D5/D5b/D6 report).
python analyze_validation.py

# 4. Extended ejection (D5b): release the tether on the L=10 final structure and
#    let the chain diffuse out the tunnel (+x), longer than the in-run ejection.
python eject_demo.py 2000000              # -> synth_out/ejection_long/
```

`topo-csp` echoes the config, loads the ribosome, precomputes contacts **once**, then
prints a per-residue kinetic schedule and grows the chain stage by stage. The
per-residue dwell-time table is written to `synth_out/dwell_times.dat`.

### Output layout (mirrors the reference `output/traj/1/` per-length layout)

```
synth_out/
├── dwell_times.dat          per-residue sampled dwell times (s) + steps  (= ref 1.out)
├── L_001/ … L_010/
│   ├── stage_1/   peptidyl transfer (A-site delivery)   traj.dcd/.log/.psf/.chk, traj_final.pdb
│   ├── stage_2/   translocation begins
│   └── stage_3/   tRNA binding / wait  ← traj_final.pdb seeds the next residue
├── ejection/                in-run ejection (50 k steps)
└── ejection_long/           extended ejection demo (D5b)
```

---

## What this run shows (validation summary)

The validation run (`csp_val.ini`, L = 1 → 10, production `scale_factor = 4331293`)
was checked against `continuous_synthesis/output/`. Headline numbers (full table in
[`NOTES.md`](NOTES.md)):

| Criterion | Result |
|-----------|--------|
| **D3** run completes | exit 0, synthesized 1 → 10 |
| **D5** finite energies | worst stage PotE **524 kJ/mol** (no ≳10¹² blow-ups) |
| **D5b** no collapse / no wall penetration | min nascent–ribosome distance **3.0 Å**, all beads stay `x ≥` tunnel wall |
| **D6a** length | 10 = 10 ✓ |
| **D6b** total in-vivo dwell | ours 0.256 s vs ref 0.253 s (**1.01×**) |
| **D6c** final radius of gyration | ours 0.80 nm vs ref 0.75 nm (**1.06×**) |

---

## Production / full-length notes

- `csp_val.ini` caps stages at `max_steps_per_stage = 200000` so the rare ~10⁶-step
  stage stays tractable on one GPU. The cap truncates **MD steps only**, not the
  sampled **dwell times in seconds**, so the D6 dwell comparison is unaffected.
- The reference itself synthesizes only **1 → 10** (`total_nascent_chain_length = 10`
  in `cont_synth_ecoli.cntrl`), so this tutorial matches that length. For the full
  306-residue chain leave `L_max` blank and remove the cap; expect a long multi-GPU
  run and use many `random_seed`s for an ensemble.
- `ribosome_traffic = no`: O'Brien's external `ribosome_traffic` binary is not on the
  PATH here, so the upstream-queue correction is off (the reference used it; its
  effect on the summed dwell time is ≲1 % at this length — see `NOTES.md`).

See the module reference `topo/csp/README.md` for every option and the kinetics API.
