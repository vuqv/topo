# Changelog

All notable changes to TOPO are documented here. The format is loosely based on
[Keep a Changelog](https://keepachangelog.com/); releases correspond to git tags.

## [2026.2] — 2026-07-18

Headline: **intrinsically disordered regions (IDRs)** in the Cα model, together with
a per-part tutorial reorganization and a round of optimizer / config refinements.

### Added
- **Disordered / IDR regions.** An optional `disordered:` section in the
  `domain_def` YAML marks residues whose native contacts are removed and replaced by
  a weak, non-specific attraction. One knob — `idr_scale` (default `0.03`; `0` = a
  pure self-avoiding chain). The mask reaches the energy path (`apply_disorder`), the
  native-contact (Q) analysis, and the nscale optimizer (IDR-involving contacts are
  excluded from Q), and continuous synthesis gets it for free — all from the single
  section. Folded-only runs (no `disordered:` section) are byte-for-byte identical.
  `intra_domains` is now optional, so a fully-disordered chain needs only
  `disordered:` + `n_residues`.
- **Tutorial A.7 — a protein with IDRs** (mixed folded + disordered, P28566) and a
  **Disordered / IDR regions** user guide (how the model treats each pair class, how
  to write the `disordered:` section, and sourcing IDR ranges via MobiDB or the
  AlphaFold pLDDT < 70 proxy).
- IDR build-log reporting: a compact, range-collapsed domain/IDR overlap notice, and
  an "IDR masking" line reporting how many input-structure native contacts are
  excluded vs kept.

### Changed
- **Tutorials renumbered to a per-part scheme:** Part A `A.1–A.7`, Part B `B.1–B.2`
  (folders and doc stubs renamed accordingly, cross-references updated).
- **Optimizer / config:** control files no longer need an `[OPTIONS]` header; each
  round's `md.ini` is headerless; `outdir` is a control key; STRIDE is computed once
  into the optimization root and reused every round; the `--python` flag was dropped.
- **Default protocol:** `md_steps` 10 000, output frequency 5000, optimizer
  `ref_t` 300 K.
- Internal refactors: a single `run_stride` in `topo.utils.external`; the non-bonded
  energy constants flattened into named UPPER_CASE module constants.

### Fixed
- The build-time large-force check prints a concise notice instead of a multi-line
  `RuntimeWarning` for a benign, unminimized input.

### Docs
- README Citation section with Zenodo DOI badges; the enhanced-sampling and
  custom-MD examples split into self-contained folders; an optimizer control-options
  page plus option-table cleanup; stale references to an earlier tutorial layout
  removed; meta-commentary trimmed across the doc set.

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
  last completed residue), and a post-synthesis **ejection** phase.
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
