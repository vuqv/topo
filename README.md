# TOPO: TOPOlogy-based coarse-grained model for folded prOteins

A coarse-grained simulation engine for **globular (folded) proteins** using OpenMM. TOPO builds topology/structure-based models with bonds, angles, periodic torsions, electrostatics, and optional contact-based non-bonded interactions.

## Requirements

- **OpenMM** ≥ 7.7
- **Parmed**
- See `requirements.txt` for full list.

## Install

```bash
# Clone or download the repo, then add to Python path (project root = parent of topo/)
export PYTHONPATH=$PYTHONPATH:/path/to/topo
```

Install dependencies (e.g. conda): `conda install -c conda-forge openmm parmed` (and others from `requirements.txt`).
`mamba` is recommemded for performance and speed to install

## Usage

**From an example directory:**

```bash
cd examples/standard_example
python run_simulation.py -f md.ini
```

**Or run the top-level script:**

```bash
python topo/topo-simulation.py -f md.ini
```

Use a control file (e.g. `md.ini`) with options such as: `md_steps`, `dt`, `pdb_file`, `model` (use `topo`), `box_dimension`, `device`, etc. See `examples/standard_example/md.ini` for a template.

## Contact

Report issues to Quyen Vu (`vuqv.phys@gmail.com`).
