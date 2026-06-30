# AGENTS.md — Tutorial 13 orientation for AI agents

> **Read me first.** This file tells an agent **where the reference code lives**, **what
> this tutorial is**, and **how its files are organized**, so you don't re-discover any of
> it. It is the Tutorial-13 analogue of [`../12_auto/Goal.md`](../12_auto/Goal.md) (which
> was an *autonomous task* directive); Tutorial 12's task is **done**, so this file is an
> **orientation/context map**, not a to-do loop.

---

## 1. What this tutorial is

Tutorial 13 is a **full-length stress test** of the per-stage **stability guard** that
Tutorial 12 added to `topo.csp`. Tutorial 12 reproduced O'Brien's continuous-synthesis
protocol on **4c5c** at the reference length (L = 1 → 10); Tutorial 13 runs the **entire
306-residue chain** (L = 1 → 306) to confirm the dt-halving guard in
[`topo/csp/core.py`](../../topo/csp/core.py) `run_length` eliminates the energy blow-ups
across the full chain. See [`README.md`](README.md) for the result table and
[`THEORY.md`](THEORY.md) for how the three MD sub-stages map onto the biological
elongation cycle.

**This is a validation/robustness run — not a quantitative dwell-time match.** For the
dwell-time / R_g comparison against the reference, use **Tutorial 12** (it bundles the
reference run; this folder does not).

## 2. Where the reference code is (read-only — semantics only, never run as deliverable)

The O'Brien lab CG synthesis source is **outside this repo**, on the filesystem:

| Path | What it is |
|------|-----------|
| `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` | O'Brien lab CG model + protocol suite (Jiang, Nissley, O'Brien). |
| `…/Continuous_synthesis_protocol/continuous_synthesis_v6.py` | **The protocol being reproduced.** OpenMM driver: per-codon, 3-stage elongation (A-site tRNA binding → peptide-bond formation → translocation). Key fns: `run_elongation`, `elongation`. **New-residue placement lives here** — read `elongation()` (≈ line 142+: it slices `prot_psf[':1-L']`, combines with the ribosome psf, and seeds coordinates for residue L). |
| `…/Continuous_synthesis_protocol/ribosome_traffic` | Upstream-queue delay binary (ribosome-traffic model). |
| `…/Continuous_synthesis_protocol/visualize_cont_synth.py`, `render_*_RNC.tcl` | Reference VMD visualization scripts. |
| `…/CG_protein_parameterization/`, `…/CG_ribosome_parameterization/` | How the reference CG `.psf/.top/.prm` were built (model-building reference). |

The **reference *run*** (inputs + outputs to validate against) is **not** bundled here —
it is in **Tutorial 12**:

| Path | What it is |
|------|-----------|
| [`../12_auto/continuous_synthesis/input/cont_synth_ecoli.cntrl`](../12_auto/continuous_synthesis/input/cont_synth_ecoli.cntrl) | Reference protocol parameters (kinetics, ribosome, traffic, `scale_factor`). |
| [`../12_auto/continuous_synthesis/input/setup/`](../12_auto/continuous_synthesis/input/setup/) | Reference CG `.psf/.top/.prm/.cor`, mRNA, codon trans-times (CHARMM format). |
| [`../12_auto/continuous_synthesis/output/`](../12_auto/continuous_synthesis/output/) | Reference trajectory + `info.log` / `1.out` to compare numbers against. |

## 3. The topo reimplementation (reuse — do **not** re-implement)

The protocol is already ported into the package. Use it; don't rewrite it.

| Path | Role |
|------|------|
| [`../../topo/csp/`](../../topo/csp/) | The port. `core.py` (MD building blocks + `run_length` stability guard), `protocol.py` (`.ini` → params), `kinetics.py` (FPT sampling), `ribosome.py`, `cg_ribosome.py`. |
| `topo-csp` (console script) | The runner. `topo-csp -f <ini>`. |
| [`../../topo/csp/DESIGN.md`](../../topo/csp/DESIGN.md), [`PROMPT.md`](../../topo/csp/PROMPT.md), [`README.md`](../../topo/csp/README.md) | Design rationale, invariants, config reference. |
| [`../12_auto/WHY_10_FAILS.md`](../12_auto/WHY_10_FAILS.md) | Post-mortem on the original blow-up that the guard fixes. |

## 4. This folder's file structure

| File / dir | Role |
|------------|------|
| `4c5c_model_clean.pdb`, `4c5c_model_clean_stride.dat` | 4c5c all-atom structure (306 aa) + precomputed STRIDE. **Raw input — do not overwrite.** |
| `domain.yaml` | 3-domain map + per-domain Go-scale strengths. |
| `4c5c_mrna.txt`, `trans_times.txt` | mRNA (one codon/residue) + Fluitt *E. coli* codon-time table. **Raw input.** |
| `ribosome_trunc.pdb` | Truncated CG 50S + tRNAs (X-aligned); P-/A-anchors + rigid scenery. **Raw input.** |
| `csp.ini` | **Debug** profile (short L) → `synth_out_debug/`. |
| `csp_val.ini` | **Full-length** profile (L = 1 → 306) → `synth_out/`. |
| `THEORY.md`, `THEORY.vi.md` | Protocol theory (EN / VI). Read for the biology↔MD mapping. |
| `analyze_validation.py` | D5 energy scan + D5b ejection check (copied from Tutorial 12). |
| `eject_demo.py`, `make_movie.py` | Ejection demo + VMD movie stitcher (wraps `topo-csp-movie`). |
| `synth_out/`, `synth_out_debug/` | **Outputs only.** Per-residue `L_<L>/stage_<s>/`, `ejection/`, `dwell_times.dat`, `movie.*`. Safe to regenerate. |

## 5. How to run / verify

See [`README.md`](README.md) for full commands. In short, from this folder:

```bash
topo-csp -f csp.ini        # debug smoke test  -> synth_out_debug/
topo-csp -f csp_val.ini    # full-length test  -> synth_out/  (watch for "[stability]" lines)
```

The headline check: **no stage's |PotE| exceeds 1e9 kJ/mol** (the README has the scan
one-liner; `analyze_validation.py` prints the D5/D5b checks).

## 6. Safety rules

- **Never overwrite raw inputs** (`4c5c_model_clean.pdb`, `domain.yaml`, `*mrna*`,
  `trans_times.txt`, `ribosome_trunc.pdb`) or the reference data under
  `../12_auto/continuous_synthesis/`. Write only under `synth_out*/`.
- The reference source in `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` is
  **read-only** — consult it for semantics; do not modify or run it as the deliverable.
