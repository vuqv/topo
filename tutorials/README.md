# TOPO Tutorials

Hands-on, ready-to-run tutorials for the **TOPO** package — a topology-based
(structure-based / Gō-like) coarse-grained model for folded proteins, built on
[OpenMM](https://openmm.org/).

Each subfolder is **one self-contained example**: it ships the input files you
need, and its `README.md` walks you through the run step by step and explains
the concepts involved. The tutorials fall into two parts — **(A)** simulating a
folded protein with the coarse-grained model, and **(B)** co-translational
synthesis (growing the chain on the ribosome). Work through them in order; Part B
builds on Part A.

## Part A — Coarse-grained protein simulation

Model a folded protein as a one-bead-per-residue structure-based (Gō-like) model
and run / analyze its dynamics: the basic workflow, multidomain scaling, restarts
and outputs, many copies at once, contact-nscale optimization, and temperature
protocols.

| # | Tutorial | What you learn |
|---|----------|----------------|
| 1 | [Single-domain quickstart](01_single_domain.md) | The minimal workflow: a config file, one PDB, run an MD simulation, read the outputs. |
| 2 | [Multidomain & domain scaling](02_multidomain.md) | Multidomain proteins: per-domain contact scaling via `domain.yaml`, including a discontiguous domain. |
| 3 | [Restart & outputs](03_restart.md) | Continuing a run from a checkpoint, and a tour of every output file. |
| 4 | [Many copies in one run](04_multicopy.md) | Run N non-interacting chains at once to fill the GPU, then split into per-chain trajectories. |
| 5 | [Optimizing the contact nscale](05_opt_nscal.md) | Automatically search the per-domain/interface `nscale` (*n*<sub>scale</sub>) that keeps each domain folded, instead of hard-coding it. |
| 6 | [Temperature annealing & quenching](06_anneal.md) | Run a temperature protocol — hold the protein hot to unfold it, then T-jump (or slow-cool) back to `ref_t` to study refolding. |

## Part B — Translation (protein synthesis)

Grow the nascent chain **residue by residue**, so it can fold *as it is
synthesized*. Both tutorials build on the Part A model and use codon-resolved
kinetics; they differ only in how the ribosome exit tunnel is represented — an
analytic bore (7) vs an explicit truncated ribosome (8). See the
[protein synthesis overview](../usage/synthesis_overview.md) for the
biology and which model to use, and
[the ribosome structure](../usage/ribosome_preparation.md) for where tutorial 8's
`ribosome_trunc.pdb` comes from (four organisms ship ready-made).

| # | Tutorial | What it is |
|---|----------|------------|
| 7 | [Synthesis through an analytic tunnel](07_translation_cylinder.md) | The exit tunnel is modelled analytically as a cylindrical bore through an infinite wall (a "hole in a wall") — **no explicit ribosome beads**. The nascent chain is the only system, so it is fast and never jams, and the protein **folds co-translationally** as it extrudes and clears the bore. |
| 8 | [Synthesis on a coarse-grained ribosome](08_ribosome_synthesis.md) | The ribosome-based counterpart: grow the chain through TOPO's own truncated **coarse-grained ribosome** (`ribosome_trunc.pdb`, standalone — no CHARMM inputs), with codon-resolved kinetics and equilibrium-PTC seeding, then eject the completed protein. Worked on 4c5c (306 aa) and P0CX28 (106 aa). |

The **ready-to-run files** for each tutorial (PDB, `md.ini`, `run_simulation.py`,
…) live in the matching folder under
[`tutorials/`](https://github.com/vuqv/topo/tree/main/tutorials) in the repository.

**Reference docs** that complement the tutorials:

- [The TOPO model: theory & force field](../usage/model_theory.rst) — every energy term, its formula, constants, and parameter sources. Read this to understand *why* the model behaves as it does.
- [Simulation control options](../usage/simulation_control.rst) — every `md.ini` option.
- [Domain definition file](../usage/domain_definition.rst) — the `domain.yaml` format.
- [Native-contact (Q) analysis](../usage/native_contacts.rst) — measuring how folded the protein is.
- [Output files & log format](../usage/outputs.rst) — every file a run writes, and how to parse the log.
- [Using TOPO from Python](../usage/python_api.rst) — the scripting API.
- [Continuous synthesis protocol (CSP)](../usage/continuous_synthesis.md) — the `topo-csp` codon-resolved synthesis runner: `csp.ini` options, the rigid-ribosome force field, the tRNA tether, and the movie tool.

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
  attractive well; everything else is repulsive. Contact nscales come from
  hydrogen bonds (via **STRIDE**), backbone–sidechain and sidechain–sidechain
  geometry, and optional **per-domain scaling** (Tutorial 2).

Because the native state is the energy minimum, TOPO is ideal for studying
folding/unfolding, domain motions, and mechanical/thermal stability.

For the full functional forms, constants, and where each parameter comes from,
see [The TOPO model: theory & force field](../usage/model_theory.rst).

## How you run it

Every simulation is driven by a plain-text config file (conventionally
`md.ini`) and launched with either of:

```bash
python -m topo.mdrun -f md.ini          # the canonical package runner
python run_simulation.py -f md.ini    # same thing, from inside a tutorial folder
```

The runner reads `md.ini`, builds the CA model from your PDB
(`topo.models.buildCoarseGrainModel`), and runs Langevin dynamics. It lives in
the package as `topo.mdrun`; each tutorial folder keeps a tiny `run_simulation.py`
that just calls it (`from topo.mdrun import mdrun`), so every example stays
self-contained while the runner has a single canonical implementation.

A full reference for `md.ini` options lives in
[Simulation control options](../usage/simulation_control.rst);
the `domain.yaml` format is documented in
[Domain definition file](../usage/domain_definition.rst).

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
- Output files are written to the run folder `output_dir` (default `traj/`),
  named `<outname>.*` (default `traj.*`), and are **not** committed — you
  generate them by running the tutorial.
