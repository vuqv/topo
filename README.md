# TOPO: TOPOlogy-based coarse-grained model for folded prOteins

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21360706.svg)](https://doi.org/10.5281/zenodo.21360706)

A coarse-grained molecular-dynamics engine for **globular (folded) proteins**,
built on [OpenMM](https://openmm.org/). From a single folded-protein structure,
TOPO builds a **one-bead-per-residue, structure-based (Gō-like) model** — bonds,
angles, sequence-dependent periodic torsions, Debye–Hückel electrostatics, and a
native-contact potential — and runs Langevin dynamics.

That one model powers **two complementary workflows**:

- **A · Folded-protein simulation** — take a complete structure and study how it
  moves, unfolds, and comes apart: folding/unfolding, thermal and mechanical
  stability, and multidomain motions. Contact energies can be scaled per domain
  and per interface.
- **B · Protein synthesis** — grow the nascent chain N→C, one residue at a time,
  with codon-resolved kinetics, so the protein folds *co-translationally* as it
  emerges from the ribosome exit tunnel (analytic-tunnel or explicit CG-ribosome
  variants).

Part B builds directly on the Part A model, so start with A if you are new here.

## Requirements

### Python dependencies

- **Python** ≥ 3.9
- **OpenMM** ≥ 7.7 (the MD engine; install from the `conda-forge` channel)
- **ParmEd** ≥ 3.4
- **NumPy** ≥ 1.22, **pandas** ≥ 1.4, **PyYAML** ≥ 6.0
- **MDAnalysis** ≥ 2.2, **mdtraj** ≥ 1.9.7 (trajectory/structure I/O and analysis)

These are the exact runtime dependencies of `import topo`, declared in
[`pyproject.toml`](pyproject.toml) and mirrored in
[`requirements.txt`](requirements.txt). Floors are the oldest versions known to
work; there are no upper caps — the package runs on current releases (e.g.
NumPy 2.x, OpenMM 8.x). The standalone tools under `scripts/` need a few extra
packages (`scipy`, `matplotlib`, `numba`); see the commented section in
`requirements.txt`.

> OpenMM ships compiled CUDA/OpenCL kernels, so install it (and ideally the rest)
> with **conda/mamba** from `conda-forge` rather than pip. `mamba` is recommended
> for faster, more reliable solves.

### External program (STRIDE)

TOPO also calls a third-party command-line binary, **STRIDE**. It is a compiled C
program, so it is **not** installed by `pip` and **not** bundled in the wheel —
you install it once and TOPO locates it at run time.

| Program     | Needed            | Used for                                                              |
| ----------- | ----------------- | --------------------------------------------------------------------- |
| **STRIDE**  | Required          | Secondary-structure / backbone H-bond assignment for the contact map. |
| **PULCHRA** | Optional (opt-in) | Backmapping a coarse-grained (Cα) structure to all-atom coordinates.  |

STRIDE is only invoked when TOPO has to *build* the contact map; if you supply a
precomputed STRIDE file (`stride_output_file=...`), it need not be installed for
that run. PULCHRA is installed only if you ask for it — nothing but backmapping
uses it, and [cg2all](https://github.com/huhlim/cg2all) is a deep-learning
alternative for that step (untested here; see the docs). Install with the bundled
helper:

```bash
scripts/install_deps.sh              # STRIDE, into $HOME/.local/bin
scripts/install_deps.sh pulchra      # PULCHRA only (opt-in)
```

TOPO resolves each program in this order: `$TOPO_STRIDE` / `$TOPO_PULCHRA` (an
explicit path) → the program on `PATH` → a copy vendored at `topo/bin/`. See
[docs/usage/external_dependencies.rst](docs/usage/external_dependencies.rst) for
manual installs and details.

## Install

Two supported ways, depending on whether you want the console commands.

### 1. pip install (recommended)

Installs the package and registers the console commands (below). Use an editable
install (`-e`) so changes to the source take effect immediately.

```bash
# Create an environment with the binary dependencies first (OpenMM etc.):
mamba create -n topo -c conda-forge python">=3.9" openmm parmed mdanalysis \
    mdtraj numpy pandas pyyaml
mamba activate topo

# Then install TOPO itself from the repo root (the directory with pyproject.toml):
pip install -e .
```

Use the editable install (`-e`, above) for development. For a stable deployment,
a regular install copies the package into `site-packages` instead (source edits
then require a reinstall):

```bash
pip install .
```

Verify:

```bash
topo-mdrun        # prints help
topo-csp          # prints help
```

### 2. Add to PYTHONPATH (no install)

If you only need `import topo` and the module-form entry points (no console
commands), add the repo root (the parent of `topo/`) to `PYTHONPATH`. You must
still install the dependencies above (e.g. with conda/mamba).

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/topo      # add to ~/.bashrc to persist
```

With this method, invoke the tools as modules (the console scripts are created
only by `pip install`), e.g. `python -m topo.mdrun -f md.ini`.

## Console commands

`pip install` registers these entry points (each has a module-form equivalent,
`python -m <module>`):

| Command            | Module               | Purpose                                                              |
| ------------------ | -------------------- | ------------------------------------------------------------------- |
| `topo-mdrun`       | `topo.mdrun`         | Run a folded-protein simulation from an `md.ini` control file.       |
| `topo-optimize`    | `topo.optimize`      | Calibrate per-domain / per-interface contact `nscale`.              |
| `topo-csp`         | `topo.csp.protocol`  | Continuous synthesis on an explicit coarse-grained ribosome.         |
| `topo-cylinder`    | `topo.csp.cylinder`  | Continuous synthesis through an analytic (cylindrical) exit tunnel.  |
| `topo-csp-movie`   | `topo.csp.movie`     | Stitch per-residue/-stage synthesis trajectories into one VMD movie. |
| `topo-make-mrna`   | `topo.csp.synth_mrna`| Pre-generate a fastest/slowest synonymous-codon mRNA for a protein.  |

## Usage

### A · Folded-protein simulation

Run a simulation from a control file (`md.ini`):

```bash
topo-mdrun -f md.ini                   # installed console command
python -m topo.mdrun -f md.ini         # module form
```

An `md.ini` sets options such as `pdb_file`, `model` (use `topo`), `md_steps`,
`dt`, `device`, `n_copies`, and output naming. To calibrate contact scales so
every domain/interface stays folded:

```bash
topo-optimize -f optimize.ini -o opt_out
```

See [tutorials/A5_opt_nscal/](tutorials/A5_opt_nscal/) for the optimizer.

Part of a chain (or all of it) can be marked **intrinsically disordered** by adding
an optional `disordered:` section to the same domain-definition YAML — those
residues lose their native contacts and instead feel a weak, non-specific
attraction whose default strength is calibrated against SAXS radii of gyration:

```yaml
disordered:
  residues: [1-24, 150-165]   # 1-based, same syntax as a domain's residues
```

See [tutorials/A7_idr_mixed_protein/](tutorials/A7_idr_mixed_protein/) and the
[Disordered / IDR regions](https://vuqv.github.io/topo/usage/disordered_regions.html)
guide.

### B · Protein synthesis

Grow the chain co-translationally, either through an explicit CG ribosome or an
analytic tunnel:

```bash
topo-csp -f csp.ini                    # explicit coarse-grained ribosome
topo-cylinder -f cylinder.ini          # analytic cylindrical tunnel
```

See [tutorials/B2_ribosome_synthesis/](tutorials/B2_ribosome_synthesis/) and
[tutorials/B1_translation_cylinder/](tutorials/B1_translation_cylinder/).

## Tutorials

Ready-to-run, ordered examples live in [`tutorials/`](tutorials/):

| #  | Tutorial                         | Topic                                                    |
| -- | -------------------------------- | ------------------------------------------------------- |
| 1  | `A1_single_domain_quickstart`    | Build and run your first single-domain simulation.       |
| 2  | `A2_multidomain_domain_scaling`  | Per-domain and per-interface contact scaling.            |
| 3  | `A3_restart_and_outputs`         | Checkpoint/resume and the files a run writes.            |
| 4  | `A4_multicopy`                   | Many independent, non-interacting copies in one job.     |
| 5  | `A5_opt_nscal`                   | Automatically calibrate the contact `nscale`.            |
| 6  | `A6_anneal_quench`               | Temperature ramps to melt/quench and observe (un)folding.|
| 7  | `A7_idr_mixed_protein`           | Mark tails/loops intrinsically disordered (IDR regions).  |
| 8  | `B1_translation_cylinder`        | Co-translational synthesis through an analytic tunnel.    |
| 9  | `B2_ribosome_synthesis`          | Co-translational synthesis on a CG ribosome.             |

## Documentation

Full documentation (model theory, control-file references, Python API, and the
tutorials above) is a Sphinx site built from [`docs/`](docs/):

```bash
cd docs && ./build_docs.sh     # output in docs/_build/html/index.html
```

## Repository layout

```
topo/            The importable package
  core/          The model: system (forces), models (build entry point), geometry
  parameters/    Force-field constants (masses, radii, charges, dihedral table)
  utils/         Non-bonded/contact building, config parsing, external-binary lookup
  csp/           Continuous synthesis: protocol, ribosome, cylinder, mRNA kinetics
  analysis/      Post-processing (native-contact Q, mirror-image detection)
  mdrun/         Folded-protein simulation runner
  optimize/      Contact-scale (nscale) optimizer
  reporter/      Trajectory/energy reporters
assets/          Codon dwell-time tables and the ribosome-preparation pipeline
docs/            Sphinx documentation source
tutorials/       Ordered, ready-to-run examples (.ini-driven console commands)
examples/        Hackable Python workflow scripts (copy and edit for custom runs)
scripts/         Standalone analysis tools and install_deps.sh
tests/           Test suite
```

## Citation

If TOPO contributed to your work, please cite it (and give the repo a ⭐):

> Vu, Q. (2026). *TOPO: a topology-based coarse-grained model for folded
> proteins* (Version 2026.2) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.21360706

The DOI above is the **concept DOI** — it always resolves to the latest release.
GitHub's "Cite this repository" button reads
[`CITATION.cff`](CITATION.cff). Please also cite the underlying O'Brien-lab
models and other methods — see the
[citation guide](docs/citation.rst) for the full list.

## License

TOPO is released under the **GNU General Public License v3.0** — see
[`LICENSE`](LICENSE).

## Contact

Report issues to Quyen Vu (`vuqv.phys@gmail.com`) or open an issue on the
[GitHub repository](https://github.com/vuqv/topo).
