# Custom MD — the open, hackable folded-protein workflow

This folder is the **from-scratch starting point** for building and running a
TOPO simulation in Python. It is the *open* version of the `topo-mdrun` console
command: it spells out every step of the `build → setup → run → finalize` flow
using TOPO's public building blocks (`topo.engine`), so you can insert your own
logic between them — a custom force, a bespoke temperature schedule, an
on-the-fly analysis callback — instead of driving the whole thing through a
config file.

Everything here uses only TOPO's **public API** (`topo`, `topo.engine`); nothing
reaches into private internals. The more specialized examples
([`../parallel_tempering_openmmtools/`](../parallel_tempering_openmmtools),
[`../simulated_tempering/`](../simulated_tempering), …) build on exactly this
same flow.

This folder is **self-contained** — the input structure, domain definition,
STRIDE file, and [`md.ini`](md.ini) all live here. Run from this directory.

## Contents

| File | Purpose |
| --- | --- |
| [`custom_md.py`](custom_md.py) | the hackable script, with marked `# === EDIT: ... ===` points |
| [`custom_md.ipynb`](custom_md.ipynb) | the same flow as a notebook — run it cell by cell to inspect each intermediate object (`built`, `sim`, …) interactively |
| [`md.ini`](md.ini) | control file (system, dt, thermostat, output, device) |
| `P0CX28_clean.pdb`, `domain.yaml`, `P0CX28_clean_stride.dat` | the demo system (single domain, calibrated) |

The `.py` and `.ipynb` are two views of the identical workflow — the script for
batch runs and copy-and-edit, the notebook for interactive exploration.

## The flow

```
cfg   = read_simulation_config('md.ini')      # parse the control file
built = engine.build_system(cfg)              # coarse-grain model -> openmm.System
ctx   = engine.setup_simulation(cfg, built)   # integrator, platform, coords/vels, restart
engine.attach_reporters(cfg, ctx.simulation)  # DCD / log / checkpoint
run_protocol(sim, schedule, ...)              # step through a temperature schedule
engine.finalize_simulation(cfg, ctx, ...)     # checkpoint + final structure + metadata
```

The three edit points demonstrated in the script:

1. **Add your own OpenMM force** to `built.system` before the Context is created
   (e.g. a harmonic positional restraint or a pulling force).
2. **Define a custom temperature schedule** — a list of `(temperature, n_steps)`
   stages — instead of plain constant-`ref_t` production.
3. **Run in segments with a callback** for on-the-fly analysis or adaptive control.

## Run

Script:

```bash
python custom_md.py -f md.ini
```

Notebook:

```bash
jupyter lab custom_md.ipynb      # then run cells top to bottom
```

Both read the same [`md.ini`](md.ini), so all the usual options (`pdb_file`,
`ref_t`, `md_steps`, `dt`, `device`, `n_copies`, output naming, …) apply.

## Outputs (under `traj/`)

`<outname>.dcd` / `.log` / `.chk`, `<outname>_final.pdb`, `<outname>_runinfo.log`,
plus the topology (`.psf`) and the CG reference (`*_CA.pdb`) from the build step.

## Where to go next

For the fully-featured runner (anneal/quench, restart, multi-copy plumbing) that
this example is distilled from, read the canonical sources:
[`../../topo/mdrun/mdrun.py`](../../topo/mdrun/mdrun.py) and
[`../../topo/engine.py`](../../topo/engine.py).
