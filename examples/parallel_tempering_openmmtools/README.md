# Parallel tempering (REMD) with openmmtools

The **recommended** way to run temperature replica-exchange MD on a TOPO
coarse-grained model. [`remd_openmmtools.py`](remd_openmmtools.py) builds the
system with TOPO's public API (`topo.engine.build_system`) and hands the plain
`openmm.System` to the [openmmtools](https://openmmtools.readthedocs.io/)
`ParallelTemperingSampler`, which manages the N replicas, the Gibbs all-pairs
exchange scheme, online free-energy estimation, and NetCDF storage /
checkpointing / restart.

This folder is **self-contained** — the input structure, domain definition,
STRIDE file, and [`md.ini`](md.ini) all live here. Run everything from this
directory.

## Contents

| File | Purpose |
| --- | --- |
| [`remd_openmmtools.py`](remd_openmmtools.py) | the REMD driver |
| [`md.ini`](md.ini) | base MD settings (system, dt, friction, device) |
| [`analyze.ipynb`](analyze.ipynb) | notebook to explore the output (mixing, walker diffusion, energy overlap, MBAR free energies, structure demux) |
| `P0CX28_clean.pdb`, `domain.yaml`, `P0CX28_clean_stride.dat` | the demo system (single domain, calibrated) |

## Requirements

```bash
conda install -c conda-forge openmmtools
```

(OpenMM and TOPO are already in your environment.)

## Run

```bash
python remd_openmmtools.py -f md.ini --tmin 300 --tmax 500 --nreplicas 8 \
    --exchange-interval 1000 --cycles 100
```

- `--tmin/--tmax/--nreplicas` — geometric temperature ladder (K). `ref_t` in the
  `md.ini` is ignored.
- `--exchange-interval` — MD steps per iteration between exchange attempts.
- `--cycles` — number of exchange iterations (default: `md_steps //
  exchange-interval`). Total steps per replica = `cycles × exchange-interval`.

Set `device = GPU` in the `md.ini` for large systems; for this 106-bead demo the
run is exchange-overhead-bound, so CPU is just as fast.

## Outputs (under `out/`)

| File | Contents |
| --- | --- |
| `remd_remd.nc` | all replicas, all iterations (energies, states, periodic coordinates) |
| `remd_remd_checkpoint.nc` | restart checkpoint |
| `remd.log` | openmmtools per-iteration log (progress, timing, ETA) |
| `remd.psf`, `P0CX28_clean_CA.pdb` | topology / CG reference from the build step |

## Analyze

Open [`analyze.ipynb`](analyze.ipynb) after a run — it reads `out/remd_remd.nc`
and walks through the standard REMD diagnostics. The essentials directly:

```python
from openmmtools.multistate import MultiStateReporter, ParallelTemperingAnalyzer
r = MultiStateReporter("out/remd_remd.nc", open_mode="r")
a = ParallelTemperingAnalyzer(r)
dF, dF_err = a.get_free_energy()                            # MBAR free energies (kT)
mixing, eigenvalues, tau = a.generate_mixing_statistics()  # exchange mixing
r.read_replica_thermodynamic_states()                      # replica -> temperature map
r.read_sampler_states(iteration=i)                         # demux one temperature
```

> **Note:** coordinates are stored only every `checkpoint_interval` iterations, so
> the structural ensemble is coarser than the energy/state records. The notebook
> handles this.

## Tuning the ladder

Aim for **neighbour exchange acceptance ≈ 0.2–0.4** across the whole ladder and
check that walkers make several full **round trips** (bottom → top → bottom). If
acceptance collapses somewhere (typically near a melting transition), add
replicas there or narrow the range. A short run is enough to read these
diagnostics before committing to a long one.

---

TOPO composes here without any change to TOPO itself: `build_system` returns a
standard `openmm.System`, which is exactly what openmmtools consumes. The same
handoff works for the rest of the OpenMM ecosystem — see the sibling
[`../simulated_tempering/`](../simulated_tempering) (OpenMM built-in) and
[`../parallel_tempering_scratch/`](../parallel_tempering_scratch) (hand-written)
examples.
