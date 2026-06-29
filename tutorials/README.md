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

### Part A — Coarse-grained protein simulation

Model a folded protein as a one-bead-per-residue structure-based (Gō-like) model
and run / analyze its dynamics: the basic workflow, multidomain scaling, restarts
and outputs, many copies at once, contact-strength optimization, and temperature
protocols.

| # | Tutorial | What you learn |
|---|----------|----------------|
| 1 | [Single-domain quickstart](https://vuqv.github.io/topo/tutorials/01_single_domain.html) | The minimal workflow: a config file, one PDB, run an MD simulation, read the outputs. |
| 2 | [Multidomain & domain scaling](https://vuqv.github.io/topo/tutorials/02_multidomain.html) | Multidomain proteins: per-domain contact scaling via `domain.yaml`, including a discontiguous domain. |
| 3 | [Restart & outputs](https://vuqv.github.io/topo/tutorials/03_restart.html) | Continuing a run from a checkpoint, and a tour of every output file. |
| 4 | [Many copies in one run](https://vuqv.github.io/topo/tutorials/04_multicopy.html) | Run N non-interacting chains at once to fill the GPU, then split into per-chain trajectories. |
| 5 | [Optimizing the contact strength](https://vuqv.github.io/topo/tutorials/05_opt_nscal.html) | Automatically search the per-domain/interface `strength` (*n*<sub>scale</sub>) that keeps each domain folded, instead of hard-coding it. |
| 6 | [Temperature annealing & quenching](https://vuqv.github.io/topo/tutorials/06_anneal.html) | Run a temperature protocol — hold the protein hot to unfold it, then T-jump (or slow-cool) back to `ref_t` to study refolding. |

### Part B — Translation (co-translational synthesis)

Grow the nascent chain **residue by residue** on (or through) the ribosome, so it
can fold *as it is synthesized*. These build on the Part A model and add the
elongation drivers (`topo-elongate`, `topo-csp`): vectorial synthesis, an analytic
exit tunnel, codon-resolved kinetics, and a reproduction + validation of O'Brien's
continuous-synthesis protocol.

Status legend: ✅ works · ⚠️ works with a caveat · ❌ does not run to completion here.

| # | Tutorial | Status | What it is |
|---|----------|--------|------------|
| 7 | [Protein synthesis](https://vuqv.github.io/topo/tutorials/07_protein_synthesis.html) | ✅ works | **The foundation.** Synthesize a protein vectorially on a rigid coarse-grained ribosome with `topo-elongate`: grow the chain N→C one residue per step at a **fixed** `n_steps`, restraining the C-terminus to the ribosome P-site, then eject it — and make a movie of the chain emerging from the exit tunnel. |
| 9 | [Co-translational synthesis through an analytic tunnel](https://github.com/vuqv/topo/tree/main/tutorials/09_translation_cylinder) | ✅ works | A variant of Tutorial 7 with **no ribosome beads**: the exit tunnel is modelled analytically as a cylindrical bore drilled through an infinite wall (a "hole in a wall"). The nascent chain is the only system, so it is fast and never jams, and the folded protein **folds co-translationally** as it extrudes and clears the bore. A standalone `cylinder.py` reusing the unchanged `topo.translation` machinery. |
| 10 | [Continuous synthesis (O'Brien CSP)](https://github.com/vuqv/topo/tree/main/tutorials/10_csp_obrien) | ⚠️ demo only | **The kinetic upgrade of Tutorial 7.** Synthesize with **codon-resolved kinetics** using `topo-csp`: time each residue from its mRNA codon and add it through O'Brien's three sub-stages (peptidyl transfer → translocation → tRNA binding). The short clamped demo runs fine, but a **full-length run blows up** (≈5/306 stages → PotE ~10¹³; see its `OBSERVATIONS.md`). That bug is **fixed in Tutorials 12/13** (the dt-halving stability guard in `run_length`); re-run against the patched `topo` to make it work. |
| 11 | [O'Brien CSP reference (legacy/CHARMM)](https://github.com/vuqv/topo/tree/main/tutorials/11_reproduce_csp) | ❌ does not synthesize here | **The reference, not a topo runner.** The raw O'Brien `continuous_synthesis_v6.py` + its CHARMM `setup/` (`.psf/.top/.prm/.cor`), `.cntrl` files, the `ribosome_traffic` binary, SLURM script and outputs. In this environment the script exits cleanly ("All Done", exit 0) but **only sets up length 1 and never elongates** (`traj/1/` holds just `rnc_l1.psf`) — which is exactly why the topo port (10/12/13) exists. Kept as the protocol Tutorials 12/13 are validated against and the source of the CHARMM-ingest path. See its `README.md`. |
| 12 | [Reproducing O'Brien CSP on 4c5c](https://github.com/vuqv/topo/tree/main/tutorials/12_auto) | ✅ works (validated) | A **validated reproduction** of O'Brien's protocol on the 4c5c system with `topo-csp`: maps the reference `cont_synth_ecoli.cntrl` onto `csp.ini`, runs L=1→10 with production kinetics, and checks energies, ejection and per-codon dwell times against the bundled reference (Tutorial 11) run. Adds the per-stage **stability guard** that fixes the Tutorial-10 blow-ups (see its `WHY_10_FAILS.md`). |
| 13 | [Full-length validation of the Tutorial-12 fix](https://github.com/vuqv/topo/tree/main/tutorials/13_validate_claude_fix12) | ✅ works (validated) | **Validates Tutorial 12 at full scale** — synthesizes the *entire* 306-residue 4c5c chain, the regime where Tutorial 10 blew up. **Works as expected:** all 306 residues synthesize with **zero** stage blow-ups (919 stages scanned, worst max\|PotE\| ≈ 1755 kJ/mol), confirming the dt-halving stability guard holds across the whole chain. Includes a VMD movie of the full chain growing out of the exit tunnel. |

The **ready-to-run files** for each tutorial (PDB, `md.ini`, `run_simulation.py`,
…) live in the matching folder under
[`tutorials/`](https://github.com/vuqv/topo/tree/main/tutorials) in the repository.

**Reference docs** that complement the tutorials:

- [The TOPO model: theory & force field](https://vuqv.github.io/topo/usage/model_theory.html) — every energy term, its formula, constants, and parameter sources. Read this to understand *why* the model behaves as it does.
- [Simulation control options](https://vuqv.github.io/topo/usage/simulation_control.html) — every `md.ini` option.
- [Domain definition file](https://vuqv.github.io/topo/usage/domain_definition.html) — the `domain.yaml` format.
- [Native-contact (Q) analysis](https://vuqv.github.io/topo/usage/native_contacts.html) — measuring how folded the protein is.
- [Output files & log format](https://vuqv.github.io/topo/usage/outputs.html) — every file a run writes, and how to parse the log.
- [Using TOPO from Python](https://vuqv.github.io/topo/usage/python_api.html) — the scripting API.
- [Protein synthesis](https://vuqv.github.io/topo/usage/protein_synthesis.html) — the `topo.translation` elongation runner: `elongate.ini` options, the rigid-ribosome force field, the tRNA tether, and the movie tool.

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

For the full functional forms, constants, and where each parameter comes from,
see [The TOPO model: theory & force field](https://vuqv.github.io/topo/usage/model_theory.html).

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
[Simulation control options](https://vuqv.github.io/topo/usage/simulation_control.html);
the `domain.yaml` format is documented in
[Domain definition file](https://vuqv.github.io/topo/usage/domain_definition.html).

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
