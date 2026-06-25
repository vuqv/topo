# Tutorial 1 — Single-domain quickstart

**Goal:** run your first TOPO simulation end-to-end and understand the inputs and
outputs. We use a small single-domain protein (`P0CX28`, 106 residues) and the
simplest possible configuration.

**Time:** the simulation finishes in ~2 seconds on a CPU.

---

## Files in this folder

| File | Role |
|------|------|
| `P0CX28_clean.pdb` | Input structure (all-atom PDB; TOPO keeps only the CA atoms). |
| `domain.yaml` | Defines the single domain and its **calibrated** contact strength. |
| `md.ini` | Simulation configuration (steps, temperature, I/O, hardware). |
| `run_simulation.py` | The runner script (reads `md.ini`, builds the model, runs MD). |

## Background concepts

- **Coarse-graining.** TOPO reads your all-atom PDB but keeps **one bead per
  residue** at the alpha-carbon position. A 106-residue protein → 106 particles.
- **Structure-based force field.** The native contacts in `P0CX28_clean.pdb`
  define the attractive interactions.
- **Contact strength is a calibration parameter.** Even for a single domain,
  `domain.yaml` carries a `strength` — the global scale factor on the
  sidechain–sidechain contact energies. It is tuned so the model reproduces a
  target property (here, realistic stability at 300 K); for `P0CX28` the
  calibrated value is **2.5044**. The raw, unscaled value (1.0) would leave the
  protein under-stabilized and only marginally folded. Tutorial 2 generalizes
  this to *different* strengths per domain and across domain interfaces.
- **STRIDE.** On the first run TOPO calls `stride` on the PDB to find backbone
  hydrogen bonds and writes `P0CX28_clean_stride.dat`. You don't manage this file
  by hand.

## Step-by-step

### 1. Check your environment
From this folder:
```bash
python -c "import topo, openmm; print('OK', openmm.__version__)"
which stride
```
Both must succeed (see the [tutorials overview](https://vuqv.github.io/topo/tutorials/index.html) for setup).

### 2. Look at `md.ini`
Open `md.ini`. The important lines:
```ini
md_steps = 5000          # how long to run (short, for a demo)
ref_t = 300              # temperature in Kelvin
pbc = no                 # no periodic box (single protein, no solvent)
pdb_file = P0CX28_clean.pdb
protein_code = P0CX28    # prefix for ALL output files
domain_def = domain.yaml # single domain with the calibrated contact strength (2.5044)
device = CPU             # runs anywhere; switch to GPU if you have CUDA
minimize = no            # native structure is already the energy minimum
```
There is no `stride_output_file` line, so STRIDE is run automatically (next step).

### 3. Run it
```bash
python run_simulation.py -f md.ini
```
You'll see TOPO build the model (it prints the number of chains, adds each force
term, runs STRIDE, builds the contact matrices) and then step the dynamics. It
ends with `--- Finished in … seconds ---`.

### 4. Inspect the outputs
After a successful run you get (all prefixed with `protein_code = P0CX28`):

| File | What it is |
|------|------------|
| `P0CX28.log` | Tab-separated energy/temperature log (one line every `nstlog` steps). |
| `P0CX28.dcd` | Trajectory (coordinates every `nstxout` steps) — open with VMD/MDAnalysis. |
| `P0CX28.chk` | Binary checkpoint (positions + velocities) for restarting (Tutorial 3). |
| `P0CX28.psf` | Topology of the CA model (for loading the DCD in analysis tools). |
| `P0CX28_init.pdb` | The CA-only starting structure TOPO actually simulated. |
| `P0CX28_final.pdb` | The last frame. |
| `P0CX28_clean_stride.dat` | Cached STRIDE output (auto-generated). |

Peek at the log:
```bash
column -t P0CX28.log | head
```
The columns are step, time (ps), potential / kinetic / total energy (kJ/mol),
temperature (K), speed, and remaining time. A stable temperature near 300 K and
a non-exploding potential energy mean the run is healthy.

## Try next

- Bump `md_steps` to `50000` and watch the trajectory in VMD:
  `vmd P0CX28.psf P0CX28.dcd`.
- Raise `ref_t` toward the protein's unfolding temperature and look for the
  potential energy climbing as native contacts break.
- Move on to **Tutorial 2** to handle a multidomain protein.

> Tip: each run overwrites the `P0CX28.*` outputs. To keep a run, copy them
> aside or change `protein_code`.
