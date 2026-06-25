# Tutorial 4 — Many copies in one run (better GPU utilization)

**Goal:** run **many non-interacting copies** of a single chain in one simulation,
so a GPU (which is wasted on one ~100-bead chain) stays busy and you collect
**N independent trajectories per run** for better sampling. Then split the
multi-chain trajectory back into per-chain DCDs for normal analysis.

**Time:** the 10-copy demo finishes in a few seconds on a CPU; the real win is on
a GPU.

**Prerequisite:** [Tutorial 1](https://vuqv.github.io/topo/tutorials/01_single_domain.html).

---

## Why do this?

A coarse-grained protein is tiny — `P0CX28` is 106 CA beads. A modern GPU can
integrate tens of thousands of particles per step at almost the same wall-clock
cost as a few hundred, so simulating **one** small chain leaves the device almost
idle. Packing `n_copies` independent chains into a single `System` fills the GPU
and turns one run into `n_copies` trajectories — exactly what you want when a
study needs many trajectories for statistics (folding rates, free-energy
estimates, etc.).

The copies must be **truly independent** (no chain should feel another). This is
guaranteed two ways:

- bonded terms (bonds/angles/torsions) and constraints are duplicated **within**
  each copy with offset atom indices, and
- every `CustomNonbondedForce` (the Yukawa electrostatics and the structure-based
  contacts) is restricted to **intra-copy** interactions via OpenMM
  `addInteractionGroup(group_k, group_k)`.

A direct check confirms it: the potential energy of `N` copies equals exactly
`N ×` the single-chain energy.

## Files in this folder

| File | Role |
|------|------|
| `P0CX28_clean.pdb` | Single-chain input structure (106 residues). |
| `domain.yaml` | Calibrated single-domain strength (2.5044), as in Tutorial 1. |
| `md.ini` | Config; note the **`n_copies`** and `copy_shift` options. |
| `run_simulation.py` | The standard runner — replicates automatically when `n_copies > 1`. |
| `split_chains.py` | Post-process: split the multi-chain DCD into per-chain DCDs. |

## How it works in the run script

Multi-copy support is controlled entirely by `n_copies` (default `1`). The
standard `run_simulation.py` keeps its normal structure and just adds one
conditional right after building the single-chain model:

```python
cgModel = topo.models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())

if cfg.n_copies > 1:
    system, topology, positions = topo.make_noninteracting_copies(
        cgModel.system, cgModel.topology, cgModel.positions,
        n_copies=cfg.n_copies, shift=cfg.copy_shift * unit.nanometer)
else:
    system, topology, positions = cgModel.system, cgModel.topology, cgModel.positions
# ... everything below runs on (system, topology, positions) unchanged ...
```

So Tutorials 1–3 and this one use the **same** `run_simulation.py`; only
`n_copies` in `md.ini` differs.

`make_noninteracting_copies` is the package primitive that takes the single-chain
model and returns the replicated `(system, topology, positions)`. It wraps
`replicate_system_intra_only`, `replicate_topology` and `replicate_positions`
(all on `topo`). Copy *k* is the contiguous atom block `[k*n : (k+1)*n]`.

## Step-by-step

### 1. Set `n_copies` in `md.ini`
```ini
n_copies = 10          ; independent chains in one run
copy_shift = 5.0       ; nm, x-offset between copies at the start
protein_code = P0CX28_x10
device = CPU           ; <-- set to GPU on a CUDA machine; that's the point
```

### 2. Run
```bash
python run_simulation.py -f md.ini
```
The **same** standard runner is used as in Tutorials 1–3; because `n_copies = 10`
the run script's `if cfg.n_copies > 1` branch replicates the model via
`topo.make_noninteracting_copies`. Output (DCD + PSF only, no PDB):
`P0CX28_x10.dcd` (all chains), `P0CX28_x10.log`, `P0CX28_x10.chk`,
`P0CX28_x10.psf` (single-chain topology) and `P0CX28_x10_multi.psf`
(multi-chain topology, matches the combined DCD).

### 3. Split into per-chain trajectories
```bash
python split_chains.py -f md.ini
```
This writes `P0CX28_x10_chain{0..9}.dcd` — one ordinary single-chain trajectory
per copy. Load each with the single-chain PSF:
```python
import MDAnalysis as mda
u = mda.Universe("P0CX28_x10.psf", "P0CX28_x10_chain0.dcd")
# RMSD, Q, Rg, ... per independent trajectory
```

## Notes & tips

- **Pick `n_copies` for your GPU.** Start with 10–50 small chains and increase
  until the per-step time stops being flat (you've saturated the device). Memory
  is the eventual limit.
- **Independent sampling.** Each copy gets different random Langevin forces, so
  the `n_copies` trajectories are independent samples of the same ensemble —
  ideal for averaging observables and estimating error bars.
- **`copy_shift`** only sets the starting separation; since the chains never
  interact, its exact value does not affect the physics (it keeps the initial
  structure tidy and easy to view).
- **No PBC here.** With `pbc = no` the chains simply occupy different regions of
  space. If you enable PBC, make the box large enough to hold all shifted copies.

## Try next
- Switch `device = GPU` and compare ns/day for `n_copies = 1` vs `50` — the
  per-copy cost drops sharply once the GPU is filled.
- Apply the same pattern to the multidomain system from Tutorial 2.
