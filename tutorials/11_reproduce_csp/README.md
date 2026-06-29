# Tutorial 11 — O'Brien CSP reference (legacy / CHARMM)

> **Status: ❌ does not synthesize in this environment.** The reference script runs
> to completion without error (`slurm-53953059.out` → "All Done", exit 0) but **only
> sets up nascent length 1 and never elongates** — `traj/1/` contains just
> `rnc_l1.psf` and `output/1.out` logs a single "Elongation at length 1". This is the
> motivation for the topo port; for a *working* continuous-synthesis run use
> **[Tutorial 12](../12_auto/)** (`topo-csp`, validated) and
> **[Tutorial 13](../13_validate_claude_fix12/)** (full 306-residue chain).

This folder is **not a topo tutorial** in the usual sense — it is the **original
O'Brien Continuous Synthesis Protocol**, kept here as the reference that the topo
port (`topo.csp`, Tutorials 10/12/13) reproduces and is validated against. It is the
upstream `continuous_synthesis_v6.py` plus the exact CHARMM-style coarse-grained
system it needs.

## What it is

`continuous_synthesis_v6.py` builds a ribosome–nascent-chain complex from a **CHARMM
CG force field** (`.psf` topology + `.top`/`.prm` parameters, loaded via `parmed` into
OpenMM) and grows the nascent chain with codon-resolved, three-stage elongation
kinetics (peptidyl transfer → translocation → tRNA binding). It is the protocol whose
*timing core* and *3-stage mechanics* `topo.csp` re-implements in topo style; the
difference is the force field (CHARMM `go_bt` here vs a topo structure-based CG build
in 10/12/13) and the runner.

## Files

| File / dir | Role |
|------------|------|
| `continuous_synthesis_v6.py` | The O'Brien protocol (the program). |
| `setup/` | The CHARMM CG system: protein `4c5c_model_clean_ca.{psf,top,cor}` + `..._nscal1_fnn1_go_bt.prm`; truncated ribosome `50S_tRNA_cg_truncated.{psf,cor}` + `50S_tRNA_cg.top` + `combine_ribo_L24_Yang.prm`; mRNA `4c5c_mrna_sequence_fast.txt`; Fluitt codon-time table; domain/secondary-structure defs; `rnc.{prm,xml}`. |
| `cont_synth_test.cntrl` | Control file — quick test (`total_nascent_chain_length = 3`, `scale_factor = 433129300`, no traffic). |
| `cont_synth_ecoli.cntrl` | Control file — the reference run (`= 10`, `scale_factor = 4331293`, `ribosome_traffic = 1`, `restart = 1`). Matches `12_auto/continuous_synthesis/input/cont_synth_ecoli.cntrl`. |
| `cont_synth_full.cntrl` | Control file — full chain (`= 306`, `scale_factor = 43312930`, no traffic). |
| `ribosome_traffic` | O'Brien's compiled upstream-queue-delay binary (the traffic correction). |
| `parse_cg_prm.py` | Helper to parse the CG `.prm`. |
| `submit_csp.slurm` | SLURM job: picks `cont_synth_full.cntrl`, auto-sets `restart`, runs the script on a GPU. |
| `render_*_RNC.tcl`, `visualize_cont_synth.py` | VMD/Python visualization of the ribosome–nascent complex. |
| `output/`, `traj/`, `info.log` | Output of the **local** (incomplete — length-1-only) run here. |
| `example_info.log` | An `info.log` from a complete upstream reference run, for comparison. |
| `test_run_bak/` | A backup of an earlier local run. |

> The **validated reference trajectory** (O'Brien's original 2021 `cont_synth_ecoli`
> run, lengths 1→10) is **not** here — it is bundled in
> [`../12_auto/continuous_synthesis/output/`](../12_auto/continuous_synthesis/) and is
> what Tutorial 12's `analyze_validation.py` compares against. This folder's `setup/`
> is the source of that reference system.

## How to run (and what currently happens)

```bash
cd tutorials/11_reproduce_csp
conda activate bioenv                       # environment with OpenMM + parmed
export PATH="$PWD:$PATH"                     # so `ribosome_traffic` is found

# Pick a control file and run the reference protocol directly:
python continuous_synthesis_v6.py -f cont_synth_test.cntrl     # total length 3
# or cont_synth_ecoli.cntrl (10) / cont_synth_full.cntrl (306)

# Or submit the bundled SLURM job (runs cont_synth_full.cntrl, auto restart):
sbatch submit_csp.slurm
```

**Observed result here (A100, OpenMM 8.2, `cont_synth_full.cntrl`):** the job finishes
with exit 0 and prints `All Done`, but produces **only** `traj/1/rnc_l1.psf` and a
single length-1 record in `output/1.out` — no per-stage MD trajectories, no elongation
past length 1. The run sets up the length-1 system, reaches the "minimizing" status
(`info.log`), then the worker exits without synthesizing. So the reference script does
**not** reproduce a synthesis run in this environment as-is.

### To come back later / debug

- The script dispatches trajectories through a multiprocessing pool ("Setup process
  pool containing 1 processors") — check whether the pool worker is silently
  swallowing an exception (e.g. a `parmed`/OpenMM template clash; note the
  `CGG`/`CGA`, `CGC`/`CGU` residue-collision warnings in `slurm-*.out`) and returning
  before running the per-length MD. Run with a single trajectory inline (no pool) to
  surface the real error.
- Confirm the `ribosome_traffic` binary runs (`./ribosome_traffic` with no args) when
  a `.cntrl` sets `ribosome_traffic = 1`.
- Compare `info.log` against `example_info.log` (a complete upstream run) to see where
  this run stops short.

## Relationship to the topo port

| | Tutorial 11 (here) | Tutorials 10 / 12 / 13 (topo) |
|--|--|--|
| program | O'Brien `continuous_synthesis_v6.py` | `topo-csp` (`topo.csp`) |
| force field | CHARMM CG (`go_bt`, via parmed) | topo structure-based CG build from the PDB |
| status here | ❌ length-1 only | ✅ 10 (demo, validated) / ✅ 306 (full) |
| use it for | the protocol definition + the CHARMM `setup/` system + the reference data | actually running and validating continuous synthesis |

If CHARMM `.psf/.top/.prm/.cor` ingestion is added to `topo.csp` (the open item in
`10_csp_obrien/TASK.md`), this folder's `setup/` is the system to ingest so topo can
run O'Brien's exact force field rather than a topo rebuild.
