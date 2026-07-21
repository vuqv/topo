# Changelog

All notable changes to TOPO are documented here. The format is loosely based on
[Keep a Changelog](https://keepachangelog.com/); releases correspond to git tags.

## [Unreleased]

## [2026.3] — 2026-07-21

Headline: **SAXS-calibrated IDR defaults** — a generic cohesion term and transferable
van der Waals radii bring disordered-chain dimensions in line with experiment, and TOPO
is now described as a unified model for globular *and* disordered proteins.

### Added
- **Generic cohesion for disordered regions.** The IDR–IDR well depth is now
  `max(eps_NN, eps_gen + idr_scale * eps_BT)`: a sequence-independent generic term
  `eps_gen_kj` (new optional `disordered:` key) added to the sequence-dependent BT
  sidechain–sidechain energy. `eps_gen` sets the overall level of collapse and is the
  knob to turn first when tuning how compact a disordered chain is; `eps_BT`
  modulates it by residue-pair type.

### Changed
- **IDR defaults are now SAXS-calibrated (behavior change).** `eps_gen_kj` defaults
  to **2.25 kJ/mol** and `idr_scale` to **1.0** (was `0.03`, i.e. the
  sequence-dependent BT term now enters at full strength). The pair was calibrated
  against SAXS radii of gyration for 24 fully-disordered proteins (90 ns each): RMS
  fractional Rg error
  **33%**, Pearson *r* **0.79**, OLS slope **0.71** — comparable to the HPS-Urry
  force field (29% / 0.69 / 0.54) on the same benchmark. See the new *Validation
  against SAXS* section of `usage/disordered_regions` for the plot and protocol.
  **A run with a `disordered:` section that omits these keys will produce
  substantially more compact ensembles than before**; pin `eps_gen_kj`/`idr_scale`
  explicitly to reproduce older results, and do not mix old and new IDR trajectories
  in one analysis.
- **`idr_scale: 0` alone is no longer a self-avoiding chain.** Because `eps_gen_kj`
  now defaults to 2.25, a pure self-avoiding chain requires **both** `idr_scale: 0`
  and `eps_gen_kj: 0`.
- **IDR well positions use a transferable per-AA van der Waals radius.** Disordered
  residues now carry `rvdw` (~11% tighter than the bare `Rmin/2` they used before),
  and pair wells are the plain sum of per-bead radii. A folded bead keeps its native
  Karanicolas–Brooks radius in folded–IDR cross pairs (previously it was also
  shrunk). Folded–folded pairs are unchanged. The same radius now reaches both the
  intra-chain and nascent↔ribosome excluded-volume channels consistently. See
  `usage/disordered_regions` for how `rvdw` is derived.

- **TOPO is now described as "a unified coarse-grained model for globular and
  disordered proteins."** The old tagline ("TOPOlogy-based coarse-grained model for
  folded prOteins") predated IDR support and undersold what the model covers; the
  acronym expansion is dropped and TOPO is simply the name. Applied across the
  README, docs, `pyproject.toml` description, module docstrings, and `CITATION.cff`
  (title + abstract, so the next Zenodo release deposits under the new name).

### Fixed
- Minimizer no longer spins forever when the max force cannot reach the 10 kJ/mol/nm
  target: the tolerance loop now exits at the floor and reports progress per
  iteration.
- `model_parameters` uses standard average residue masses.

### Docs
- **`usage/model_theory` per-residue table corrected.** The masses were stale and
  several were plain wrong (ARG 114 → 156.19, HIS 114 → 137.14, PRO 114 → 97.12,
  CYS 114 → 103.14, ASP, GLU, GLN); the `Radii (nm)` column was neither the right
  header nor the right numbers — it is `Rmin_2` and the values were ~1.8× too large
  (ALA 0.504 → 0.2862). All 20 rows now match `model_parameters` exactly.
- **`Rmin_2` is documented as a transferable parameter.** New section explaining
  that the table is per residue *type* and is used only where a bead has no native
  structure to measure (ribosome scenery, IDR beads, the nascent-chain fallback) —
  a folded domain never uses it, taking per-residue Karanicolas–Brooks radii from
  the input structure instead.
- **The contact-potential pair-class table covers disorder.** It previously listed
  only native and non-native pairs; IDR–ordered (excluded volume, folded bead keeps
  its K–B radius) and IDR–IDR (`rvdw` sum, `max(eps_NN, eps_gen + s_IDR * eps_BT)`)
  are now spelled out alongside them.
- **IDP/IDR support is stated up front** in the README, `index`, `overview`,
  `modules/introduction`, and `usage/model_theory`, which all described TOPO as
  folded-proteins-only. Removed the note in `modules/introduction` claiming COSMO
  was the more appropriate tool for disordered proteins — TOPO handles them.
- Removed the stale "Heads-up for readers of older docs" note about a
  Gaussian-quartic dihedral, and rendered the 12-10-6 equations with proper
  fractions.

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
