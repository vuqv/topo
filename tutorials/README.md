# TOPO Tutorials

Hands-on, ready-to-run tutorials for the **TOPO** package — a topology-based
(structure-based / Gō-like) coarse-grained model for folded proteins, built on
[OpenMM](https://openmm.org/).

Each subfolder is **one self-contained example**: it ships the input files you
need, and its `README.md` walks you through the run step by step and explains
the concepts involved. Work through them in order.

| # | Folder | What you learn |
|---|--------|----------------|
| 1 | [`01_single_domain_quickstart/`](01_single_domain_quickstart/) | The minimal workflow: a config file, one PDB, run an MD simulation, read the outputs. |
| 2 | [`02_multidomain_domain_scaling/`](02_multidomain_domain_scaling/) | Multidomain proteins: per-domain contact scaling via `domain.yaml`, including a discontiguous domain. |
| 3 | [`03_restart_and_outputs/`](03_restart_and_outputs/) | Continuing a run from a checkpoint, and a tour of every output file. |

---

## What is TOPO? (the 1-minute version)

TOPO turns a folded-protein structure into a **one-bead-per-residue** (alpha-carbon,
"CA") coarse-grained model and simulates it in OpenMM. It is a **structure-based
model**: the native crystal structure you give it *defines* the energy minimum.
The force field has these terms:

- **Bonds, angles, torsions** — geometry of the CA chain (bonds are rigid
  constraints by default).
- **Yukawa electrostatics** — Debye–Hückel screened Coulomb between charged
  residues.
- **Structure-based non-bonded (contacts)** — the heart of the model. Native
  contacts (pairs of residues close together in the folded structure) get an
  attractive well; everything else is repulsive. Contact strengths come from
  hydrogen bonds (via **STRIDE**), backbone–sidechain and sidechain–sidechain
  geometry, and optional **per-domain scaling** (Tutorial 2).

Because the native state is the energy minimum, TOPO is ideal for studying
folding/unfolding, domain motions, and mechanical/thermal stability.

## How you run it

Every simulation is driven by a plain-text config file (conventionally
`md.ini`) and launched with:

```bash
python run_simulation.py -f md.ini
```

`run_simulation.py` reads `md.ini`, builds the CA model from your PDB
(`topo.models.buildCoarseGrainModel`), and runs Langevin dynamics. The same
script is copied into each tutorial folder so every example is self-contained.

A full reference for `md.ini` options lives in
[`../docs/usage/simulation_control.rst`](../docs/usage/simulation_control.rst);
the `domain.yaml` format is documented in
[`../docs/usage/domain_definition.rst`](../docs/usage/domain_definition.rst).

---

## Prerequisites

1. **Python environment with TOPO + OpenMM installed.** From the repo root the
   `topo` package must be importable:
   ```bash
   python -c "import topo, openmm; print('OpenMM', openmm.__version__)"
   ```
2. **STRIDE** on your `PATH`. TOPO calls `stride` to detect backbone hydrogen
   bonds for the contact potential. Check:
   ```bash
   which stride        # should print a path
   ```
   If `stride_output_file` is not set in `md.ini`, TOPO runs STRIDE for you and
   caches the result to `<pdb_prefix>_stride.dat`. If STRIDE is *not* installed,
   either install it or pre-compute the file (`stride -h yourprotein.pdb > stride.dat`)
   and point `stride_output_file` at it.
3. **(Optional) A GPU.** The tutorials default to `device = CPU` so they run
   anywhere; switch to `device = GPU` in `md.ini` if you have a CUDA device.

## Conventions used in the tutorials

- Settings are deliberately **small and fast** (`md_steps = 5000`, `device = CPU`)
  so each example finishes in seconds. They are demos, not production runs — for
  real science you would increase `md_steps` to millions and likely use a GPU.
- `minimize = no`: the input PDBs are native structures (already the model's
  energy minimum), so no pre-minimization is needed.
- Output files use the `protein_code` prefix from `md.ini` and are **not**
  committed — you generate them by running the tutorial.
