# File structure & procedure notes — `topo/csp/`

A map of what lives in this folder and the procedures we implement ourselves.
See [`DESIGN.md`](DESIGN.md) for the model. We **adopt the O'Brien-lab procedure**
(theory in `theory/`) but **implement everything ourselves** — we do not use their
code.


## This folder (`topo/csp/`)

| Path | What it is |
|------|------------|
| `DESIGN.md` | Design/spec for the protein-synthesis (co-translational folding) model on TOPO. |
| `FILES.md` | This file — folder map + procedure notes. |
| `theory/ja302305u_si_001.pdf` | SI of O'Brien et al., *JACS* 2012 — the RNC force-field theory we adapt. |
| `cg_ribosome.py` | Coarse-grains an all-atom protein+RNA structure to the TOPO convention (protein → Cα; RNA → P/R/BR beads). Writes a CG PDB. |
| `truncate_ribosome.py` | Truncates a CG ribosome around the exit tunnel (cylinder + exit half-space; per-residue keep). Writes a truncated CG PDB. |
| `core.py` | **Core nascent-chain MD engine (a library, not a runner).** Build-once-subset contacts, per-length build/seed/restrain/run (`run_length`, with the per-stage stability guard), new-residue placement at the A-anchor, harmonic C-terminus restraint, optional rigid ribosome (build step v2), `RunParams`. Reused by `protocol.py` (CSP) and `cylinder.py` (the analytic-tunnel model). See the build-step notes below. |
| `protocol.py` | **The CSP runner.** O'Brien's per-codon, 3-stage continuous-synthesis loop + INI: `RunParams` / `CSPConfig`, `read_csp_config`, `run_continuous_synthesis`, `csp()` CLI (`topo-csp` / `python -m topo.csp`). Calls `core.run_length` three times per residue. |
| `cylinder.py` | **The cylinder runner** (parallel to `protocol.py`). The **cylinder** ribosome model: an analytic exit tunnel (`add_tunnel_cylinder`, a hole in an infinite wall) instead of explicit beads. `CylinderParams` / `CylinderConfig`, `read_cylinder_config`, `run_cylinder_synthesis`, `cylinder()` CLI (`topo-cylinder`). Same O'Brien codon kinetics as CSP, but **one MD segment per residue** (no A→P sub-stages). Calls `core.run_length` once per residue. |
| `kinetics.py` | Pure timing core (no OpenMM): codon tables, FPT sampling, `scale_factor`→steps, the 3-stage split. |
| `ribosome.py` | **Rigid ribosome scenery (build step v2).** Loads the truncated CG ribosome as fixed (mass-0) beads and wires the ribosome↔nascent excluded-volume + electrostatics onto a built nascent model (`load_ribosome`, `append_ribosome`). |
| `movie.py` | **Movie stitcher** (CLI `topo-csp-movie`). Shared padding/stitching core (`stitch_segments` + the VMD `.tcl`) plus two layout adapters: `stitch_movie` for the CSP per-stage layout (`L_<L>/stage_<s>/`) and `stitch_length_movie` for the flat per-length layout (`L_<L>/`, used by the Tutorial-9 cylinder runner). The CLI auto-detects the layout. |
| *(tutorial)* | The ready-to-run CSP example (`csp.ini` + inputs) lives in `tutorials/12_auto/` and `tutorials/13_validate_claude_fix12/` (Sphinx pages: `docs/usage/continuous_synthesis.md`, `docs/topo.rst` autodoc). |
| `__init__.py` | Marks `topo.csp` as a package; re-exports the CSP public API (`run_continuous_synthesis` / `RunParams` / `CSPConfig` / `read_csp_config` / `csp` / `kinetics`). |
| `__main__.py` | `python -m topo.csp` → the `csp` (CSP protocol) CLI. |
| `data/` | Bundled package data. `ecoli_trans_times_310K.txt` — the default Fluitt *E. coli* (310 K) codon-time table (organism-universal; used when an INI gives no `trans_times`). |
| `structures/` | Input ribosome structures + their CG outputs (see below). |

### `structures/`

| File | What it is |
|------|------------|
| `4v9d_50S_tRNA.pdb` | *E. coli* ribosome **large subunit (50S)** with tRNA (from PDB 4V9D). Original orientation. **Not used** — kept for provenance. |
| `4v9d_50S_PtR_5jte_AtR_model.pdb` | **The working structure.** Same system with **P-site tRNA (PtR)** and **A-site tRNA (AtR)**, **re-oriented so the exit-tunnel central line lies along the X-axis**. We only work on this file. |
| `4v9d_50S_PtR_5jte_AtR_model_cg.pdb` | Coarse-grained output of the working structure, produced by `cg_ribosome.py` (generated; not raw input). |
| `4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb` | Truncated CG ribosome (exit-tunnel shell), produced by `truncate_ribosome.py` (generated). |

> **Biological-molecule labels (segID).** The input carries each molecule's
> identity in the PDB **segID** column (73–76): `AtR` (A-site tRNA), `PtR` (P-site
> tRNA), `23S` / `5S` (rRNAs), and `L2`…`L36` (large-subunit ribosomal proteins).
> `cg_ribosome.py` **preserves these** in the CG output, so you can select any
> molecule downstream (e.g. MDAnalysis `select_atoms("segid PtR")` for the P-site
> tRNA the nascent chain attaches to, or `segid L23`/`L24` for the tunnel
> proteins).


## Coarse-graining the ribosome (`cg_ribosome.py`)

Done **first**, on the full all-atom structure (before truncation). Mapping:

- **Protein** → one bead per residue at the **Cα** atom; **residue name kept**
  (ALA, GLY, …); bead atom name `CA`. Parameter lookup by residue name (as in
  TOPO's Cα model).
- **RNA** → **3 beads (pyrimidine C/U)** or **4 beads (purine A/G)** per
  nucleotide; each nucleotide stays **one residue** (original chain + residue
  number + A/U/G/C name preserved), with bead atoms:
  - `P`   — phosphate, placed at the `P` atom (absent on a 5′-terminal nucleotide),
  - `R`   — ribose-ring centroid (`C1' C2' C3' C4' O4'`),
  - `BR1` — base ring centroid (pyrimidine 6-ring; purine 6-ring `N1 C2 N3 C4 C5 C6`),
  - `BR2` — purine 5-ring centroid (`C4 C5 N7 C8 N9`), purines only.
  The **bead type** for parameter lookup (P/R/BR, defined in
  `topo.parameters.model_parameters`) is the atom name with trailing digits
  stripped (`BR1`/`BR2` → `BR`).

Non-protein/non-RNA residues (ions, water, ligands) are skipped. Run:

```bash
python cg_ribosome.py -i structures/4v9d_50S_PtR_5jte_AtR_model.pdb \
                      -o structures/4v9d_50S_PtR_5jte_AtR_model_cg.pdb
```

For `4v9d_50S_PtR_5jte_AtR_model.pdb` this yields **14,662 beads** (3,363 protein
Cα; 3,163 nucleotides → P 3,163, R 3,163, BR 4,973), across all 33 chains. Bead
coordinates are verified against the all-atom centroids.

> Naming note: the filename suggests P-site tRNA (`PtR`) from 4V9D and an A-site
> tRNA (`AtR`) modeled from PDB 5JTE. **(FILL IN / confirm)** the exact provenance
> if it matters downstream.

**Why the X-axis alignment matters:** with the tunnel central line on the X-axis,
the truncation and tunnel geometry can use the **X-axis directly** as the tunnel
center line — no separate tunnel-detection file (MOLEonline/CAVER) is needed. The
distance to the tunnel line is then just the radial distance in the Y–Z plane,
within an X-range.


## Why the ribosome can be truncated

The ribosome is **held rigid** (§2.2 of `DESIGN.md`) and the non-bonded
interactions are **cut off at 2 nm**. A ribosomal site farther than the cutoff
from every nascent-chain bead it could ever reach exerts **zero** force on the
chain. So removing ribosomal sites away from the exit-tunnel / emergence region is
**not an approximation for the nascent-chain dynamics** — it is exact — and it
cuts the particle count (and cost) dramatically.


## Pipeline order: coarse-grain first, then truncate

**CG the full ribosome, then truncate the CG beads** — not the other way around.
In the reference, `truncate_ribosome.py` takes the **CG** ribosome (`.psf`/`.cor`
from `create_cg_ribosome_model.py`) and keeps a subset of beads; it never touches
all-atom coordinates. Do the same:

```
all-atom ribosome (structures/*.pdb)  ->  CG the FULL ribosome  ->  truncate CG beads
```

Rationale: CG bead definitions (RNA P/R/BR centroids, protein Cα) are computed
from the *complete* structure, so every bead is well-defined; truncation is then
just bead selection. Truncating all-atom first would leave boundary residues
missing atoms needed for their centroids, giving ill-defined/shifted beads.


## Truncation procedure (we implement this ourselves)

Implemented in `truncate_ribosome.py`. The tunnel center line is the **X-axis**
(our structure is X-aligned), so radial distance `d = sqrt(y² + z²)`. Decide
**per residue** — keep the **whole residue/nucleotide if *any* of its beads
qualifies** (so each P/R/BR unit stays intact; matches the reference, and is safe
since the rigid ribosome is just excluded volume). A residue is kept if:

1. **Tunnel cylinder:** `d ≤ r_cyl` **AND** `x_lo ≤ x ≤ x_exit`, or
2. **Exit region:** `x ≥ x_exit` — kept **regardless of `d`**, or
3. its segID is in `--keep-segids` (optional; e.g. keep the P-site tRNA whole).

**Geometry of our structure (determined from the CG model):** the PTC (23S/A2602)
sits at **x ≈ 0** on the axis; the exit is toward **+x** (the tunnel/exit proteins
L4, L22, L23, L24, L29 are all at +x), while the tRNAs and 5S rRNA are on the −x
interface side. So PTC ≈ 0 is the low-x end and the exit is the high-x end —
directly analogous to the reference (PTC ≈ 0, `x_lo = −8`).

**Default parameters (follow O'Brien et al.):** `r_cyl = 30 Å`, `x_lo = −8 Å`,
`x_exit = 58 Å`. Our X-aligned structure shares the reference's coordinate frame
(PTC at x ≈ 0, exit plane at 58), so these reproduce the reference crop — see the
bead-count match below. (Larger `x_exit` removes more, since rule 2 keeps the
whole `x ≥ x_exit` cap.)

```bash
python truncate_ribosome.py \
    -i structures/4v9d_50S_PtR_5jte_AtR_model_cg.pdb \
    -o structures/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb \
    --r-cyl 30 --x-lo -8 --x-exit 58   # [--keep-segids PtR,AtR]
```

The PTC region needs no separate "20 Å around A2602" rule: with the PTC at x ≈ 0,
a 20 Å sphere around it lies inside the 30 Å cylinder, so rules 1–2 already cover
it.

### Truncation result

With the defaults above: **14,662 → 4,576 beads (68.8 % removed)**, 6,526 → 1,803
residues.

**Validation — this matches the reference.** O'Brien et al. report ≈ **4,577**
beads for the cropped *E. coli* ribosome; we get **4,576** — a near-exact match,
confirming both the parameters and that our structure is in the reference's
coordinate frame.

Sanity check — the retained ribosomal proteins are exactly the **tunnel/exit
machinery**: **L22 (110/110)**, **L23 (93/93)**, **L24 (102/102)**, **L29
(63/63)** kept whole, plus L4 (138/201), L17 (99/120), L32 (51/56), L34 (21/46);
the interface/far-side proteins (L5, L6, L9, L11, …) and **5S rRNA (0/418)** are
fully removed; 23S keeps 3,823/10,359 (tunnel wall + exit cap).

**tRNAs — acceptor ends retained (matches O'Brien).** The geometric crop keeps
only the PTC-proximal tRNA beads (PtR 17/267, AtR 14/255), but those are exactly
the **acceptor (CCA-3′) ends**: `AtR` residues 73–76 and `PtR` residues 72–76
survive — including **residue 76 of both**. This mirrors the reference, where the
ribosome that is truncated *includes* the tRNAs (`50S_tRNA_cg.psf` →
`..._truncated.psf`) and the elongation cycle uses the **A-site tRNA residue-76
`R` bead** to place each new amino acid (`new_AA = AtR76_R + 4.27 Å·[cos α, sin α,
0]`) while the P-site tRNA holds the nascent chain. So our truncated structure
already retains the functional attachment/placement points. If you want the
**whole tRNAs** (e.g. for visualization or a different scheme), add
`--keep-segids PtR,AtR`.


## What we implement (vs. reuse from TOPO)

- **Reuse from TOPO:** the whole force field (bonds, angles, torsions, Yukawa,
  native 12-10-6 contacts, non-native `(σ/r)¹²` excluded volume, BT potential,
  `nscale`) and the P/R/BR RNA bead parameters — already present.
- **Implement here:**
  1. **Coarse-grain the ribosome** from `structures/` (Cα protein beads + P/R/BR
     RNA beads) and hold it **rigid**.
  2. **Truncate** it around the exit tunnel (procedure above), using the X-axis as
     the tunnel center line.
  3. **Elongation driver** (the `topo.csp` runner) that grows the nascent
     chain N→C, reusing `topo.engine` (see `DESIGN.md` §4).
- **Do not use:** any external (O'Brien-lab) code — procedure only.


## Build step v1 — elongation loop (`elongate.py`)

Implements the mechanics-only loop of `PROMPT.md` (no ribosome forces yet): the
System is the **nascent chain only**; the truncated ribosome supplies just the two
fixed anchor points (`read_anchor` reads the `PtR`/`AtR` residue-76 `R` beads).
Reuses `topo.engine` (`setup_simulation` / `attach_reporters` /
`finalize_simulation`) for each length; the contact matrices are built **once** on
the full native PDB (`precompute_contacts`) and each length injects the top-left
`L×L` block through a **dedicated build path** (`build_length_model`) that mirrors
`topo.models.buildCoarseGrainModel` but supplies the precomputed matrices instead
of re-running STRIDE — so the stock builder and the rest of `topo.core` are left
untouched.

Two implementation decisions worth recording (both validated by the v1 run on
P0CX28, `L0=5`, which grows one residue per step, writes per-length outputs, and
lands each new C-terminus on the P-anchor to within ~0.005 nm):

- **Flexible (harmonic) bonds, not the package-default rigid `AllBonds`.** The new
  residue is seeded *at the A-anchor* while restrained to the P-anchor (~0.9–1.1 nm
  away), so the new bond starts far from equilibrium. A rigid distance constraint
  cannot be seeded that far off (the constraint solver / minimizer diverges); a
  harmonic bond absorbs the stretch, the per-step minimization relaxes it, and the
  restraint then translocates the residue **A→P** — exactly the `DESIGN.md` §2.5
  mechanism. (Override with `--constraints AllBonds` if seeding closer.)
- **A small zig-zag on the cold-start layout.** A perfectly collinear extended
  chain makes every bond angle 180°, where the angle-force gradient is singular
  (→ NaN on the first minimization). The cold-start beads are offset by ±0.03 nm
  perpendicular to the tunnel axis to break the exact collinearity; this is tiny
  relative to the 0.381 nm bond, so the chain is still "extended along the axis".


## Rigid ribosome forces (`ribosome.py` + `core.py`)

The supplied ribosome PDB is **always** loaded as rigid scenery (there is no
on/off flag). After the nascent model is built, the truncated ribosome is appended
to the System (`ribosome.append_ribosome`):

- **Rigid scenery.** ~4,576 ribosome beads added with **mass = 0** at indices
  `L..N-1`, coordinates as-is. OpenMM does not integrate mass-0 particles, so they
  stay fixed (verified: ~1e-15 nm drift). The nascent chain uses flexible bonds and
  the ribosome has none, so no constraint ever involves a mass-0 particle — the
  `rm_cons_0_mass` caution from `PROMPT.md` does not bite.
- **Contact force** restricted to interaction group `{nascent}×{nascent}`;
  ribosome beads get a dummy `id=0` so the `L×L` table stays nascent-only.
- **Ribosome–NC excluded volume:** a separate `CustomNonbondedForce`, pure
  `ε·(σ_ij/r)¹²` (`ε = 0.000132 kcal/mol`, σ from `model_parameters['radii']`),
  group `{nascent}×{ribosome}`.
- **Yukawa** extended with ribosome charges, groups `{nascent}×{nascent}` +
  `{nascent}×{ribosome}` (no intra-ribosome electrostatics).

Implementation notes (validated on P0CX28, `L0=5`, stable to full length; chain
stays around the tunnel and extrudes +x):

- **Identical exclusions.** OpenMM's CPU platform requires every
  `CustomNonbondedForce` to share one exclusion list, so the new ribosome–NC force
  copies the nascent bonded (1-2,1-3) exclusions from the contact force.
- **C-terminus tether (`trna_tether`, default on).** The C-terminus is tethered to
  the P-site tRNA `R` bead the O'Brien way: a harmonic bond `CA(L)–tRNA:R` plus an
  orienting angle `CA(L-1)–CA(L)–tRNA:R` (reusing the double-Gaussian backbone-angle
  force — its 91.7°/130° basins are exactly O'Brien's tether-angle parameters,
  `Continuous_synthesis_protocol/continuous_synthesis_v6.py`). The frozen tRNA bead
  + bond hold the C-terminus at the PTC; the angle is intended to **aim the chain
  down the tunnel** (constrains the terminal CA–CA vector's bend). (`trna_tether = no`
  falls back to the `ptc_offset` position restraint, whose auto-offset clears the
  P-anchor bead.)
  > **Not yet validated.** Quick single-trajectory diagnostics at the short test
  > dwell (1000 steps/residue) are too noisy and in the wrong regime (the PTC
  > collapse appears at long dwells) to show whether the tether improves extrusion.
  > Needs independent replicas at production-length dwell with a robust metric —
  > see `review/TODO.md` (§B).
- **Tunnel wall (`tunnel_wall`, default on; `add_tunnel_wall`).** O'Brien's one-sided
  planar restraint `U = k·min(x − x0, 0)²` (`k = 8368 kJ/mol/nm²` = 20 kcal/mol/Å²)
  on every nascent bead, keeping the chain at `x ≥ x0` so it can only **extrude
  forward** (+x, toward the exit) and cannot fold back past the synthesis point into
  the void where the 50S was truncated. **`x0` is auto-derived from the ribosome
  structure** (not a config knob) by the CSP runner: the lower / P-site C-terminus
  hold plane, `x0 = min(P-anchor.x, A-anchor.x) + ptc_offset` (≈ `1.05 nm` for the
  tutorial `ribosome_trunc.pdb`). Recomputed per structure, so it never goes stale.
  (O'Brien quote 58 Å, but that is their coordinate frame.) Applied throughout
  synthesis + post-elongation. **Units:** OpenMM defaults
  throughout — nm, ps, kJ/mol, K, and force constants in kJ/mol/nm².
- **PSF for the combined system** is written via parmed (`_dump_topology_psf`)
  rather than the model's nascent-only `dumpTopology`.
