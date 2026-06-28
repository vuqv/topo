# `topo.translation` — protein synthesis (elongation runner)

Simulate a protein being synthesized **vectorially** (N→C, one residue at a time)
on the ribosome, instead of refolding the full-length chain in bulk. This folder
holds the model spec, the ribosome coarse-graining/truncation tools, and the
**elongation runner** (`elongate.py`).

- **Model / design:** [`DESIGN.md`](DESIGN.md)
- **File map + procedures:** [`FILES.md`](FILES.md)
- **Build plan (v1 → v2):** [`PROMPT.md`](PROMPT.md)

> **Build steps.**
> - **v1 (default).** The simulated system is the **nascent chain only**. The
>   truncated ribosome supplies just two fixed points — the **P-anchor** (P-site
>   tRNA residue-76 `R` bead) and the **A-anchor** (A-site tRNA residue-76 `R`
>   bead) — used for placing the new residue and as the C-terminus restraint
>   target. No ribosome forces. Validates the elongation loop machinery.
> - **v2 (`rigid_ribosome = yes`).** The truncated ribosome is appended as fixed
>   (mass-0) scenery and the ribosome↔nascent **excluded volume + electrostatics**
>   are turned on (`ribosome.py`). The chain now feels the tunnel/ribosome surface.


## How it works

Synthesis is modeled as a **sequence of standalone simulations** (no live system
resizing). For each nascent-chain length `L = L0 .. N_full`:

1. **Build-once-subset contacts.** The native contact map is computed **once** on
   the full PDB (`R_full`, `eps_full`); length `L` uses the top-left `L×L` block.
   STRIDE / heavy-atom analysis are never re-run per length.
2. **Build** the length-`L` TOPO model on native residues `1..L` (bonds, angles,
   torsions, Yukawa, contacts), reusing `topo.engine`.
3. **Seed coordinates.** `L0`: an extended chain laid along the tunnel axis (+x)
   from the P-anchor. `L>L0`: residues `1..L-1` from the previous step's final
   structure, plus the new residue `L` placed at the A-anchor.
4. **Restrain only residue `L`** (the current C-terminus) to the P-anchor
   (harmonic, `k = 83680 kJ/mol/nm²` = 200 kcal/mol/Å²). The new residue therefore
   migrates **A→P**.
5. **Minimize → run `n_steps` → save.** Each length writes its own output folder
   and its final structure seeds the next length.


## Quick start

A ready-to-run example lives in `tutorials/07_protein_synthesis/`
(P0CX28 + the truncated *E. coli* ribosome). From that folder:

```bash
cd tutorials/07_protein_synthesis
topo-elongate -f elongate.ini
# equivalently:
python -m topo.translation -f elongate.ini
```

Two handy CLI overrides (everything else lives in the INI):

```bash
topo-elongate -f elongate.ini -o my_run --device GPU
```

Paths in the INI are resolved relative to the current working directory (as for
`md.ini`).


## Control file (`elongate.ini`)

A single `[OPTIONS]` section. **Required:** `pdb_file`, `ribosome`, `L0`. All
other keys are optional (defaults shown).

> **Units:** OpenMM defaults throughout — length **nm**, time **ps**, energy
> **kJ/mol**, temperature **K**, angle **rad**, force constants **kJ/mol/nm²**.

| key | default | meaning |
|-----|---------|---------|
| `pdb_file` | — | full native PDB of the target protein (the nascent chain) |
| `ribosome` | — | truncated CG ribosome PDB (source of the P-/A-anchors) |
| `L0` | — | starting nascent-chain length (cold-start layout) |
| `L_max` | full length | final length (blank → full residue count `N_full`) |
| `domain_def` | none | domain YAML for contact `n_scale` (one-time precompute) |
| `stride_output_file` | run STRIDE | precomputed STRIDE for the full PDB; if blank, STRIDE is run once (needs `stride` on PATH) and cached |
| `n_steps` | 1000 | integration steps per residue (constant schedule) |
| `dt` | 0.015 | time step (ps) |
| `ref_t` | 300 | temperature (K) |
| `tau_t` | 0.01 | Langevin friction / coupling constant (1/ps) |
| `nstout` | 50 | trajectory / log / checkpoint write frequency |
| `constraints` | `None` | `None` = flexible harmonic bonds (default); `AllBonds` = rigid |
| `restraint_k` | 83680 | C-terminus → P-anchor restraint force constant (kJ/mol/nm² = 200 kcal/mol/Å²) |
| `buffer` | 0.4 | new-residue placement offset beyond the A-anchor (nm) |
| `minimize` | yes | per-step energy minimization (relaxes the placement) |
| `rigid_ribosome` | no | **v2**: append the rigid ribosome + its forces |
| `trna_tether` | yes | v2: tether C-terminus to the P-site tRNA (bond + CA–CA–tRNA orienting angle) vs. a position restraint |
| `tunnel_wall` | yes | v2: one-sided planar wall keeping the nascent chain at x ≥ `tunnel_wall_x0` (forward-only extrusion) |
| `tunnel_wall_x0` | ~1.05 nm | wall plane = the C-terminal-AA addition plane (PTC) ≈ P-anchor x + tether bond length |
| `tunnel_wall_k` | 8368 | wall force constant (kJ/mol/nm² = 20 kcal/mol/Å²) |
| `ptc_offset` | auto | how far into the tunnel (+x) to hold the C-terminus from the P-anchor bead (nm); blank → 0 in v1, 0.4 in v2 |
| `nascent_only_output` | yes | v2: write only the nascent chain to the trajectory/PSF/final |
| `post_elongation` | stallation | post-elongation phase: `ejection` release the tether (→ `ejection/`); `stallation` keep it (→ `stallation/`) |
| `post_elongation_steps` | 0 | steps for the post-elongation phase (0 = skip) |
| `device` | CPU | `CPU` or `GPU` |
| `ppn` | 1 | CPU threads (device = CPU) |
| `outdir` | `synth_out` | per-length outputs go to `<outdir>/L_<L>/` |

Inline `#`/`;` comments are ignored; `n_steps = 1_000` (underscores) is allowed.


## Outputs

One self-contained folder per length, `<outdir>/L_<L>/`:

| file | what it is |
|------|------------|
| `traj.dcd` | trajectory for this length |
| `traj.log` | energy / temperature log |
| `traj.psf` | topology (for visualization / analysis) |
| `traj.chk` | OpenMM checkpoint |
| `traj_final.pdb` | final structure (seeds the next length) |
| `traj_runinfo.log` | run provenance (versions, hardware, timing) |
| `native_1_<L>.pdb` | length-`L` native CA structure used for the build |
| `seed.pdb` | seeded starting coordinates for this length |


## Build step v2 — the rigid ribosome

Set `rigid_ribosome = yes` to append the truncated ribosome as **fixed (mass-0)
scenery** and turn on the two ribosome↔nascent cross-interactions ([ribosome.py](ribosome.py)):

- **Excluded volume** — a dedicated `CustomNonbondedForce`, pure `ε·(σ_ij/r)¹²`
  (`ε = 0.000132 kcal/mol`, σ from `model_parameters['radii']`), restricted to
  `{nascent}×{ribosome}`.
- **Electrostatics** — the existing Yukawa force is extended with the ribosome
  charges (rRNA phosphates −1e, charged residues), restricted to
  `{nascent}×{nascent}` + `{nascent}×{ribosome}`.
- The nascent contact `L×L` table is restricted to `{nascent}×{nascent}` (ribosome
  beads get a dummy `id=0`), so it stays nascent-only.

**C-terminus tether (`trna_tether`, default on).** Instead of a generic position
restraint, the C-terminus is tethered to the P-site tRNA `R` bead the O'Brien way:
a harmonic **bond** `CA(L)–tRNA:R` plus an **orienting angle** `CA(L-1)–CA(L)–tRNA:R`
(added to the double-Gaussian backbone-angle force, whose 91.7°/130° basins are
exactly O'Brien's tether parameters). The frozen tRNA bead + bond hold the
C-terminus at the PTC; the angle is intended to **aim the chain down the tunnel**.
`trna_tether = no` falls back to the position restraint.

**Tunnel wall (`tunnel_wall`, default on).** O'Brien's one-sided planar restraint
`U = k·min(x − x0, 0)²` (`k = 8368 kJ/mol/nm²` = 20 kcal/mol/Å²) on every nascent
bead — keeps the chain at `x ≥ x0` so it can only **extrude forward** (+x, toward
the exit) and cannot
fold back past the synthesis point into the ribosome interior. **`x0` is the
C-terminal-AA addition plane** — the PTC / P-site where each new residue is placed and
tethered, ≈ `P-anchor x + tether bond length ≈ 1.05 nm` (`tunnel_wall_x0`). (O'Brien
quote 58 Å, but that is *their* coordinate frame.) Applied throughout synthesis and
the post-elongation phase.

> **Status: tether implemented, not yet validated.** Whether the tether actually
> improves extrusion is **not** established — quick single-run diagnostics at the
> short test dwell are too noisy and don't exercise the long-dwell PTC collapse.
> A proper test (independent replicas at production dwell + a robust metric) is
> pending (`TODO.md`). The tunnel wall is verified non-disruptive (in-tunnel chain
> undisturbed; nothing leaks below the cut).

The ribosome is held rigid — **no intra-ribosome forces are ever computed**, and
mass-0 freezes those beads (verified: zero drift). Because the P-anchor is itself a
ribosome bead, the C-terminus is held `ptc_offset` (default 0.4 nm) into the tunnel
from it, so the restraint doesn't pull the C-terminus onto that bead (which would
explode). Validated on P0CX28 (`L0=5`): stable runs, ribosome frozen, the nascent
chain stays around the tunnel and extrudes toward +x as it grows.

> **Output size:** by default v2 writes **nascent-only** trajectories / PSF /
> final structures (`nascent_only_output = yes`) — the rigid ribosome is static, so
> it is not written every frame (the checkpoint still holds the full system). Set
> `nascent_only_output = no` to write the whole system instead.


## Post-elongation: ejection vs. stallation

Once the chain reaches its final length, an optional phase continues the
final-length system (from its last structure) for `post_elongation_steps` steps:

- **`post_elongation = ejection`** — the C-terminus tether is **released** (no
  restraint), so the completed protein is free to move and leave the ribosome.
  Written to `<outdir>/ejection/`.
- **`post_elongation = stallation`** — the C-terminus restraint is **kept**, so the
  chain stays **stalled** on the ribosome (a control / ribosome-bound ensemble).
  Written to `<outdir>/stallation/`.

The phase runs only if `post_elongation_steps > 0`. It reuses the same machinery
(same system, including the rigid ribosome in v2), so its output folder has the
usual `traj.dcd` / `.log` / `.psf` / `_final.pdb`. Example:

```ini
post_elongation       = ejection
post_elongation_steps = 50_000
```

(Validated: with the restraint off the C-terminus drifts off the PTC; with it on
the C-terminus stays pinned there.)


## Visualize the elongation in VMD

Each length is a separate trajectory with a **different bead count** (length `L`
has `L` CA beads), and a single VMD molecule needs a constant atom count across
frames — so you can't just load all the per-length DCDs into one molecule. Stitch
them first into one fixed-width movie (padded to the longest length, with
not-yet-synthesized beads parked off-screen), which also writes a ready VMD script:

```bash
topo-elongate-movie -o synth_out
#   -> synth_out/movie.psf, synth_out/movie.dcd, synth_out/movie.tcl
vmd -e synth_out/movie.tcl
```

The movie plays as *F* frames of `L=L0`, then *F* frames of `L=L0+1`, … (e.g. 20
frames of L=5, then 20 of L=6, …). If a **post-elongation** phase was run, its
frames (`ejection/` or `stallation/`, at the final length) are appended **after**
the growth — so the movie continues into the protein's release/ejection (or stalled
wiggling). `movie.tcl` colors the chain by residue id and **hides the parked beads
each frame** (via a per-frame `selupdate` selection), so the chain visibly grows
N→C. Press Play (or run `animate forward`).

Options: `--park cterm` stacks future beads on the C-terminus instead of parking
them far away (no per-frame hiding needed, but leaves a small bead cluster at the
growing tip); `--prefix NAME` changes the output basename. Programmatically:
`from topo.translation.make_movie import stitch_movie`.

The stitcher works for **both v1 and v2** runs, because v2 writes nascent-only
per-length trajectories by default. For v2, pass `--ribosome <trunc.pdb>` to also
load the static ribosome as scenery the chain grows inside:

```bash
topo-elongate-movie -o synth_out \
    --ribosome topo/translation/structures/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
vmd -e synth_out/movie.tcl
```

> Running `vmd -dispdev text` may print `OptiXRenderer ... File not found`
> warnings — those are just the GPU ray-tracer being unavailable headless and are
> harmless.


## Notes / gotchas

- **Flexible bonds are the default** (`constraints = None`), *not* the package
  default rigid `AllBonds`. The new residue is seeded at the A-anchor while
  restrained to the P-anchor (~1 nm away), so the new bond starts far from
  equilibrium; a rigid distance constraint cannot be seeded that far off (the
  solver diverges). A harmonic bond absorbs the stretch and the restraint then
  translocates the residue A→P. (See the "Revision list" in `TODO.md` for the
  planned fix that would restore rigid bonds.)
- **No ribosome confinement in v1**, so the short nascent chain folds/collapses
  freely rather than extruding straight down the tunnel — expected until v2 adds
  the rigid ribosome forces.
- **Programmatic use:** `from topo.translation import run_elongation,
  ElongationParams` (or `read_elongate_config`) to drive the loop from Python.
