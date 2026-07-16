# Parallel tempering (REMD) — from scratch

A **dependency-free** temperature replica-exchange driver:
[`remd.py`](remd.py) implements parallel tempering directly on TOPO's public
building blocks (`topo.engine` + OpenMM), with no external sampling library. It is
the readable reference for *how REMD actually works* — copy it when you want full
control of the exchange loop, or when you can't install
[openmmtools](../parallel_tempering_openmmtools) (the managed alternative).

This folder is **self-contained** — inputs, [`md.ini`](md.ini), and the script all
live here. Run from this directory.

## How it works

- Builds the coarse-grained system **once** and shares it across N OpenMM
  `Context`s, one per replica.
- Lays out a **geometric** temperature ladder between `--tmin` and `--tmax`.
- Runs **serially on one device**: each cycle steps every replica for
  `--exchange-interval` steps, then attempts nearest-neighbour swaps
  (alternating even/odd pairs).
- Accepts swaps with the standard Sugita–Okamoto Metropolis criterion,
  `Δ = (β_i − β_j)(U_i − U_j)`. On acceptance it swaps **coordinates** and
  rescales velocities by `√(T_new/T_old)`, so each `Context` stays at a fixed
  temperature — meaning `out/remd_rep00.dcd` is the coldest canonical ensemble,
  `_rep01` the next rung, and so on, **ready for analysis without demultiplexing**.

Restart/checkpointing and multi-GPU/MPI are intentionally omitted to keep the
loop legible; use the openmmtools version for those.

## Contents

| File | Purpose |
| --- | --- |
| [`remd.py`](remd.py) | the from-scratch REMD driver |
| [`md.ini`](md.ini) | base MD settings (system, dt, friction, device) |
| `P0CX28_clean.pdb`, `domain.yaml`, `P0CX28_clean_stride.dat` | the demo system |

No extra dependencies beyond OpenMM + TOPO.

## Run

```bash
python remd.py -f md.ini --tmin 300 --tmax 500 --nreplicas 8 \
    --exchange-interval 1000 --seed 1
```

- `--tmin/--tmax/--nreplicas` — geometric ladder (K); `ref_t` in the `md.ini` is ignored.
- `--exchange-interval` — MD steps between exchange attempts.
- `--cycles` — number of cycles (default: `md_steps // exchange-interval`).
- `--seed` — RNG seed for the Metropolis accept/reject stream (reproducibility).

## Outputs (under `out/`)

| File | Contents |
| --- | --- |
| `remd_repNN.dcd` / `.log` | one trajectory + log per temperature rung (fixed-T ensembles) |
| `remd_remd.log` | exchange log: per-cycle attempts + per-pair acceptance ratios |
| `remd_repNN_final.pdb` | last frame of each replica (seed for a follow-up run) |
| `remd.psf`, `P0CX28_clean_CA.pdb` | topology / CG reference |

The per-pair acceptance summary is printed at the end and appended to
`remd_remd.log`. Aim for **neighbour acceptance ≈ 0.2–0.4** across the ladder; if
a pair drops below ~0.2, add replicas there or narrow the range.

---

This script needs no OpenMM *extension* at all — it drives the raw OpenMM
`Context` API that `topo.engine` exposes. When you'd rather not hand-roll the
exchange loop, the same TOPO `openmm.System` drops straight into
[openmmtools](../parallel_tempering_openmmtools) or any other OpenMM-compatible
sampler.
