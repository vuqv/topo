# Simulated tempering

[`simulated_tempering.py`](simulated_tempering.py) runs OpenMM's built-in
`SimulatedTempering` on a TOPO coarse-grained model. It reuses the standard TOPO
build → setup → reporter → finalize flow (`topo.engine`) and wraps the resulting
`Simulation` in `openmm.app.SimulatedTempering`.

This folder is **self-contained** — inputs, [`md.ini`](md.ini), and the script all
live here. Run from this directory. No extra dependencies beyond OpenMM + TOPO.

## Method — NOT replica exchange

Simulated tempering is **not** parallel tempering. It runs **one** copy of the
system whose *temperature* is a dynamical variable: the single walker periodically
attempts to hop up/down a temperature ladder, accepting with a Metropolis rule
that uses per-temperature free-energy **weights** `g_k`:

```
p_acc = min(1, exp[ -(β_n − β_m) U + (g_n − g_m) ])
```

Those weights are unknown a priori. By default this script passes `weights=None`,
so OpenMM **learns them on the fly** (Wang–Landau style) — convenient, but it
means early sampling is not yet canonical while the weights are still moving. For
a clean production run, learn the weights first (a short run, or a parallel-
tempering pass) and then supply them fixed.

Trade-offs vs. replica exchange:

| | Simulated tempering (here) | Parallel tempering ([openmmtools](../parallel_tempering_openmmtools) / [scratch](../parallel_tempering_scratch)) |
| --- | --- | --- |
| Copies | 1 walker | N replicas |
| Weights needed | **yes** (learned/estimated) | no |
| Device | single, cheap | N replicas (scales to multi-GPU) |
| Output | 1 trajectory, temperature varies | N fixed-temperature ensembles |

Use simulated tempering when you are GPU/memory-limited or want a lightweight
single-trajectory run; otherwise prefer parallel tempering.

## Contents

| File | Purpose |
| --- | --- |
| [`simulated_tempering.py`](simulated_tempering.py) | the tempering driver |
| [`md.ini`](md.ini) | base MD settings (run length = `md_steps`, dt, friction, device) |
| `P0CX28_clean.pdb`, `domain.yaml`, `P0CX28_clean_stride.dat` | the demo system |

## Run

```bash
python simulated_tempering.py -f md.ini --tmin 300 --tmax 500 --nreplicas 8 \
    --temp-change-interval 100
```

- `--tmin/--tmax/--nreplicas` — geometric ladder (K); `ref_t` in the `md.ini` is ignored.
- `--temp-change-interval` — MD steps between temperature-jump attempts.
- `--report-interval` — steps between `_st.log` lines (default: `nstlog`).

Run length is the `md.ini` `md_steps` (one continuous trajectory).

## Outputs (under `out/`)

| File | Contents |
| --- | --- |
| `remd.dcd` / `.log` | the single trajectory + standard TOPO log |
| `remd_st.log` | temperature index, temperature, and learned weights over time |
| `remd.psf`, `P0CX28_clean_CA.pdb` | topology / CG reference |

At the end the script prints the visit histogram and learned weight per rung —
flat-ish occupancy means the weights have converged. To recover a canonical
ensemble at one temperature, filter the trajectory frames by the temperature
column in `remd_st.log`.

---

This example shows TOPO composing with an **OpenMM built-in** enhanced-sampling
class: `build_system` yields a standard `openmm.System`/`Simulation`, and
`openmm.app.SimulatedTempering` wraps it unchanged. The same interface opens the
door to the rest of the ecosystem (openmmtools, openmm-plumed metadynamics,
openmm-torch ML potentials, custom forces, …).
