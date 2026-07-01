# TOPO: TOPOlogy-based coarse-grained model for folded prOteins

A coarse-grained simulation engine for **globular (folded) proteins** built on
OpenMM. TOPO builds topology/structure-based (Go-like) models with bonds, angles,
periodic torsions, electrostatics, and optional contact-based non-bonded
interactions.

## Requirements

- **Python** ≥ 3.9
- **OpenMM** ≥ 7.7 (the MD engine; install from the `conda-forge` channel)
- **ParmEd** ≥ 3.4
- **NumPy** ≥ 1.22, **pandas** ≥ 1.4, **PyYAML** ≥ 6.0
- **MDAnalysis** ≥ 2.2, **mdtraj** ≥ 1.9.7 (trajectory/structure I/O and analysis)

These are the exact runtime dependencies of `import topo`. They are declared in
[`pyproject.toml`](pyproject.toml) and mirrored in
[`requirements.txt`](requirements.txt). Floors are the oldest versions known to
work; there are no upper caps — the package runs on current releases (e.g.
NumPy 2.x, OpenMM 8.x). The standalone tools under `scripts/` need a few extra
packages (`scipy`, `matplotlib`, `numba`); see the commented section in
`requirements.txt`.

> OpenMM ships compiled CUDA/OpenCL kernels, so install it (and ideally the rest)
> with **conda/mamba** from `conda-forge` rather than pip. `mamba` is recommended
> for faster, more reliable solves.

## Install

Two supported ways, depending on whether you want the console commands.

### 1. pip install (recommended)

Installs the package and registers the `topo-mdrun` / `topo-optimize` console
commands. Use an editable install (`-e`) so changes to the source take effect
immediately.

```bash
# Create an environment with the binary dependencies first (OpenMM etc.):
mamba create -n topo -c conda-forge python">=3.9" openmm parmed mdanalysis \
    mdtraj numpy pandas pyyaml
mamba activate topo

# Then install TOPO itself from the repo root (the directory with pyproject.toml):
pip install -e .
```

Verify:

```bash
topo-mdrun        # prints help
topo-optimize     # prints help
```

### 2. Add to PYTHONPATH (no install)

If you only need `import topo` and the module-form entry points (no console
commands), add the repo root (the parent of `topo/`) to `PYTHONPATH`. You must
still install the dependencies above (e.g. with conda/mamba).

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/topo      # add to ~/.bashrc to persist
```

With this method, invoke the tools as modules (the `topo-mdrun`/`topo-optimize`
console scripts are created only by `pip install`):

```bash
python -m topo.mdrun -f md.ini
python -m topo.optimize -f optimize.ini
```

## Usage

**Run a simulation** from a control file (`md.ini`) — equivalent forms:

```bash
topo-mdrun -f md.ini                   # installed console command
python -m topo.mdrun -f md.ini         # module form
python run_simulation.py -f md.ini     # thin shim shipped in each tutorial
```

A `md.ini` sets options such as `pdb_file`, `model` (use `topo`), `md_steps`,
`dt`, `device`, `n_copies`, output naming, etc. See the
[tutorials](tutorials/) for ready-to-run templates.

**Optimize interaction nscales** (per-domain / per-interface `nscale`) from a
minimal `optimize.ini` — equivalent forms:

```bash
topo-optimize -f optimize.ini -o opt_out   # installed console command
python -m topo.optimize -f optimize.ini -o opt_out
```

See [tutorials/05_opt_nscal/](tutorials/05_opt_nscal/) for details.

## Contact

Report issues to Quyen Vu (`vuqv.phys@gmail.com`).
