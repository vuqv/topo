# Changelog

All notable changes to TOPO are documented here. The format is loosely based on
[Keep a Changelog](https://keepachangelog.com/); releases correspond to git tags.

## [2026.1] — 2026-07-12

**Initial tagged release.** Releases use CalVer `YEAR.N`; this first 2026 release
(`2026.1`) keeps TOPO in lockstep with its
sibling project [COSMO](https://github.com/vuqv/cosmo) (whose `cosmo.csp` mirrors
`topo.csp`). This entry summarizes TOPO's feature set at first release.

TOPO is a topology-based (Gō-like) coarse-grained MD engine for **folded proteins**,
built on OpenMM. From one folded-protein structure it builds a one-bead-per-residue,
structure-based model and powers two complementary workflows.

### Part A — Folded-protein simulation
- **Structure-based (Gō) model**: bonds, angles, sequence-dependent periodic torsions,
  Debye–Hückel electrostatics, and a native-contact potential (STRIDE secondary
  structure; PULCHRA back-mapping).
- **Per-domain contact scaling** — `nscale` per domain and per interface via a
  `domain.yaml` definition file.
- **`topo-mdrun`** simulation runner (`md.ini`): PBC, Langevin thermostat + barostat,
  minimization, and **checkpoint/restart**.
- **Temperature protocols** — annealing / quenching (jump and linear-ramp).
- **Multi-copy runs** — pack `n_copies` independent, non-interacting chains into one
  system to saturate a GPU and collect N trajectories per run (`split_chains`).
- **`topo-optimize`** — calibrate per-domain `nscale` from stability criteria.
- **Analysis** — native-contact and mirror-image (chirality) detection utilities.

### Part B — Co-translational protein synthesis (`topo.csp`)
- Grows the nascent chain N→C with codon-resolved kinetics, in two confinement models:
  - **Cylinder model** (`topo-cylinder`) — an analytic cylindrical exit tunnel (no
    ribosome beads; fast, never jams).
  - **CG-ribosome model** (`topo-csp`) — an explicit coarse-grained ribosome as rigid
    scenery, with the 3-stage elongation cycle and the A→P C-terminus switch.
- **Standalone CG ribosome** — no external CHARMM `.cor/.psf/.prm` at runtime; the
  ribosome loads from a truncated CG PDB. Structures are provided for **four organisms
  — E. coli, yeast, N. crassa, and human** — swapped via the `ribosome` key.
- **Per-codon kinetics** with dwell-time tables (E. coli, human, yeast, N. crassa);
  **mRNA generator** (`topo-make-mrna`) in fastest / slowest / median synonymous-codon
  modes.
- **tRNA tether**, **rigid-bond (`AllBonds`) support**, **resume** (continue from the
  last completed residue), and post-synthesis **ejection / dissociation** phases.
- **Movie tool** (`topo-csp-movie`) — stitch per-residue/-stage trajectories into a VMD
  movie, with `--tunnel` (analytic tunnel) or `--ribosome` (CG ribosome) scenery.

### Tooling & docs
- **Console commands**: `topo-mdrun`, `topo-optimize`, `topo-csp`, `topo-cylinder`,
  `topo-csp-movie`, `topo-make-mrna`.
- **Eight hands-on tutorials** (single-domain, multidomain, restart, multi-copy,
  nscale optimization, anneal/quench, cylinder synthesis, ribosome synthesis) and full
  Sphinx documentation (usage guides + API reference), including a
  "Visualizing the synthesis process" page.
- The CG ribosome pipeline is standalone with model parameters consolidated into a
  single source of truth (`model_parameters`), and 100% docstring coverage.

### Credit
The `topo.csp` continuous-synthesis subsystem follows the published continuous-synthesis
protocol (Jiang, Nissley & O'Brien), adapted to TOPO's structure-based Gō model.
