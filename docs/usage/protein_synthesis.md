# Protein synthesis (`topo.translation`)

The `topo.translation` module synthesizes a protein **vectorially on the ribosome**
(N→C, one residue at a time) and optionally **ejects** it, instead of refolding the
full-length chain in bulk. It is the third runner in the package, alongside
`topo.mdrun` and `topo.optimize`.

- **CLI:** `topo-elongate -f elongate.ini` (or `python -m topo.translation -f elongate.ini`)
- **Movie tool:** `topo-elongate-movie -o <out_root> [--ribosome ribo.pdb]`
- **Walkthrough:** Tutorial 7.
- **Design/spec & implementation notes:** `topo/translation/DESIGN.md`, `FILES.md`.

## Theory

### The biological picture

Proteins are synthesized **vectorially** — N-terminus first, one residue at a time
— by the ribosome. The nascent chain emerges through the ribosomal **exit tunnel**
and can begin to fold **co-translationally**, often before the full sequence is
available. The resulting folding pathways, intermediates, and yields can differ
from refolding of the full-length chain in bulk (what `topo.mdrun`'s equilibrium /
annealing runs simulate). The physical features the model aims to capture are:

- **Vectorial growth.** Chain length increases over time, N→C; at any moment only
  the already-synthesized residues `1..L` exist.
- **Tethering at the PTC.** The C-terminal residue is held at the peptidyl
  transferase center / tunnel exit; the chain extrudes from there.
- **Exit tunnel + ribosome surface.** Excluded volume and electrostatic
  confinement from the ribosome shape early folding.
- **Elongation kinetics.** Residues are added at a finite, codon-dependent rate;
  the competition between **elongation time and folding time** is the central
  observable (e.g. domain-folding order, misfolding, entanglement during
  synthesis).

### The force field (RNC Hamiltonian)

The ribosome–nascent-chain (RNC) force field is the **same model family as the rest
of TOPO** — the Best–Chen–Hummer variant of the Karanicolas–Brooks Cα Gō model,
following O'Brien et al. (2011/2012). Every term already exists in the engine; no
new physics is introduced for synthesis:

| Term | Functional form | TOPO force |
|------|-----------------|------------|
| Cα–Cα bonds | `k_b (r − r₀)²` | bonds |
| Bond angles | double-Gaussian `−(1/γ) ln[e^… + e^…]` | Gaussian angle |
| Dihedrals | `Σ_j k_φ[1 + cos(j φ − δ)]` | periodic torsion |
| Electrostatics | Debye–Hückel `q_i q_j /(4πε₀ε_r r)·e^(−r/l_D)` | Yukawa |
| Native contacts (intra-chain) | `ε_ij[13(σ/r)¹² − 18(σ/r)¹⁰ + 4(σ/r)⁶]`, `ε_ij ∝ BT` | structure-based contacts |
| Non-native / excluded volume | `ε_ij(σ_ij/r)¹²`, `ε = 0.000132 kcal/mol` | non-native excluded volume |

The key synthesis-specific consequence is that **native contacts to
not-yet-synthesized residues are simply absent**: the force field for a length-`L`
chain is the full native model *restricted to residues `1..L`*, so contacts
involving residues `> L` do not exist — no masking is required, and going `L→L+1`
only **reveals** residue `L+1`'s contacts. The contact map is therefore built
**once** on the full PDB and subset to the top-left `L×L` block each step (no
per-length STRIDE; this also keeps the `1..L` map identical at every length).

### Representing the ribosome (v2)

Following O'Brien et al. (2012), the ribosome is an **explicit rigid** coarse-grained
bead cloud — frozen (mass-0, not integrated), with coordinates used *as-is* so the
absolute tunnel/exit geometry is preserved. rRNA uses TOPO's existing P/R/BR beads
(phosphate `q = −1e`; ribose and base-ring centroids). Because it is rigid, **no
intra-ribosome forces are computed** — only ribosome↔nascent-chain interactions
matter, and those are **excluded volume + electrostatics only** (no attractive /
native contacts to the ribosome):

- a dedicated repulsive `ε·(σ_ij/r)¹²` force over `{nascent}×{ribosome}` (σ from
  the per-residue collision diameters in `model_parameters['radii']`, which match
  O'Brien's values), and
- the existing Yukawa electrostatics extended over the ribosome charges.

The source structure is the *E. coli* 50S subunit + P-/A-site tRNAs, re-oriented
with the exit-tunnel axis on **+x**, then coarse-grained (`cg_ribosome.py`) and
truncated to the residues near the tunnel (`truncate_ribosome.py`). The two tRNA
acceptor ends survive truncation and provide the **P-anchor** and **A-anchor**
used for residue placement and tethering (see Protocol). In **v1**
(`rigid_ribosome = no`) the bead cloud is omitted entirely and only those two
anchor points are read from the ribosome PDB — useful for validating the
elongation loop without the cost of the ~4,600-particle rigid scenery.

> **References.** The RNC force field and rigid-ribosome treatment follow O'Brien
> *et al.*, *J. Am. Chem. Soc.* **2012**, 134 (and the P/R/BR rRNA representation +
> truncation procedure from *JACS* **2011**, 133); the base Cα Gō model is
> Karanicolas & Brooks, *Protein Sci.* **2002**, the transferable variant is
> Best, Chen & Hummer, *Structure* **2005**, and native well depths use the
> Betancourt–Thirumalai potential (*Protein Sci.* **1999**). TOPO adopts the
> procedure/theory but implements everything itself.

## Protocol

Synthesis is a **sequence of standalone simulations**, one per nascent-chain
length `L` ("rebuild-and-continue"):

```
for L = L0 .. L_max:
    build the length-L CG system: nascent residues 1..L (+ rigid ribosome in v2)
    seed: residues 1..L-1 from step (L-1)'s final structure; new residue L at the A-site
    tether residue L to the P-site tRNA; run n_steps; save the final structure
after L_max (optional): ejection (release the tether) or stallation (keep it)
```

The native contact map is built **once** on the full PDB and each length uses the
top-left `L×L` block (no per-length STRIDE; contacts to not-yet-made residues are
absent). Coordinates carry across steps; velocities are re-drawn each step.

## Build steps

- **v1** (`rigid_ribosome = no`) — nascent chain only; the truncated ribosome
  supplies just the **P-anchor** (P-site tRNA res-76 `R` bead) and **A-anchor**
  (A-site tRNA res-76 `R` bead) coordinates for placement and the restraint target.
- **v2** (`rigid_ribosome = yes`) — append the rigid (mass-0) ribosome and wire:
  - **ribosome↔chain excluded volume** — a separate `ε·(σ_ij/r)¹²` force
    (`ε = 0.000132 kcal/mol`, σ from `model_parameters['radii']`), group
    `{nascent}×{ribosome}`;
  - **electrostatics** — the existing Yukawa force extended over the ribosome
    charges (rRNA phosphates −1e, charged residues), groups `{nascent}×{nascent}` +
    `{nascent}×{ribosome}`;
  - the nascent contact `L×L` table restricted to `{nascent}×{nascent}` (ribosome
    beads carry a dummy `id=0`, so the table stays nascent-only).

## C-terminus tether (v2)

With `trna_tether = yes` (default) the C-terminus is tethered to the P-site tRNA
`R` bead the O'Brien way: a harmonic **bond** `CA(L)–tRNA:R` plus an **orienting
angle** `CA(L-1)–CA(L)–tRNA:R` (added to the double-Gaussian backbone-angle force,
whose 91.7°/130° basins are O'Brien's tether parameters). The frozen tRNA bead +
bond hold the C-terminus at the PTC; the angle aims the chain down the tunnel.
`trna_tether = no` falls back to a harmonic position restraint of the C-terminus to
the P-anchor (`restraint_k`, with `ptc_offset` clearing the bead).

## Tunnel wall (v2)

With `tunnel_wall = yes` (default) a one-sided planar restraint
`U = k·min(x − x0, 0)²` is applied to every nascent bead, keeping the chain at
`x ≥ x0` so it can only **extrude forward** (+x, toward the exit) and cannot fold
back past the synthesis point into the ribosome interior. **`x0` is the
C-terminal-AA addition plane** — the PTC / P-site where each new residue is placed
and tethered, ≈ `P-anchor x + tether bond length ≈ 1.05 nm` (`tunnel_wall_x0`).
(O'Brien quote 58 Å, but that is their coordinate frame.)

## Post-elongation phase

After the final length, if `post_elongation_steps > 0`, the final-length system is
continued for that many steps:

- `post_elongation = ejection` — the C-terminus tether is **released**; the protein
  is free to move/leave the ribosome. Written to `<outdir>/ejection/`.
- `post_elongation = stallation` — the tether is **kept**; the chain stays stalled
  on the ribosome. Written to `<outdir>/stallation/`.

## Workflow

End-to-end, a synthesis study is: prepare inputs → configure → run → (post-elongate)
→ stitch a movie → analyze. The ready-to-run example lives in
`tutorials/07_protein_synthesis/` (P0CX28 + the truncated *E. coli* ribosome).

**1. Prepare inputs.**

- `pdb_file` — the **full** native structure of the target protein, cleaned (one
  model, Cα-resolvable). This defines the contact map and bonded terms for *all*
  lengths; the chain is its N-terminal residues `1..L`.
- `ribosome` — the truncated CG ribosome PDB. The bundled `ribosome_trunc.pdb`
  already carries the P-/A-site tRNA anchors; to build your own from a full
  structure use `cg_ribosome.py` then `truncate_ribosome.py` (see `FILES.md`).
- *(optional)* `domain_def` — a `domain.yaml` for per-domain contact `n_scale`
  (one-time precompute), and `stride_output_file` — a precomputed STRIDE file so
  STRIDE is not run at startup.

**2. Configure `elongate.ini`.** Copy the tutorial INI and edit the
inputs/schedule (full key reference below). The two decisions that matter most are
the **length window** (`L0`, `L_max`) and the **per-residue dwell** (`n_steps`).

**3. Choose the length window.** `L0` is the cold-start length — residues `1..L0`
are laid out extended along the tunnel axis from the P-anchor (pick a small value
with a few residues already in the tunnel, e.g. `L0 = 5`). `L_max` is the final
length; leave it blank to synthesize the whole protein, or cap it for a quick test.

**4. Run.**

```bash
cd tutorials/07_protein_synthesis
topo-elongate -f elongate.ini
# equivalently: python -m topo.translation -f elongate.ini
# handy overrides: topo-elongate -f elongate.ini -o my_run --device GPU
```

Each length `L` writes its own folder `<outdir>/L_<L>/`; the `traj_final.pdb` of
length `L` seeds length `L+1`. The run is restartable both *within* a length
(resume from its checkpoint) and *across* lengths (completed lengths are skipped).

**5. Post-elongation (optional).** With `post_elongation_steps > 0`, after `L_max`
the final system is continued either as **`ejection`** (tether released — the
protein can leave the ribosome, written to `<outdir>/ejection/`) or
**`stallation`** (tether kept — the chain stays stalled, written to
`<outdir>/stallation/`).

**6. Build the movie.**

```bash
topo-elongate-movie -o synth_out --ribosome ribosome_trunc.pdb
vmd -e synth_out/movie.tcl
```

This stitches the per-length trajectories (and any post-elongation phase) into one
growing-chain `movie.dcd`/`movie.psf` with a ready `movie.tcl`.

**7. Analyze.** Treat each length as a standalone trajectory: e.g. the
native-contact `Q` score per domain as residues appear (`topo.analysis`), the
folding order across `L`, and entanglement lifetimes. Synthesis is **stochastic**
— run many independent replicas and analyze the ensemble, not a single trajectory.

**Scaling to production.** The tutorial settings are deliberately tiny. For
production: set `n_steps` to a realistic dwell (see [Choosing
parameters](#choosing-parameters) and [Production notes](#production-notes)), run
on `device = GPU` (the v2 system is ~4,600 particles), and launch many replicas.

## `elongate.ini` options

A single `[OPTIONS]` section. **Required:** `pdb_file`, `ribosome`, `L0`.

```{note}
**Units:** OpenMM defaults throughout — length nm, time ps, energy kJ/mol,
temperature K, angle rad, force constants kJ/mol/nm².
```

| key | default | meaning |
|-----|---------|---------|
| `pdb_file` | — | full native PDB of the target protein (the nascent chain) |
| `ribosome` | — | truncated CG ribosome PDB (P-/A-anchors; rigid scenery in v2) |
| `L0` | — | starting nascent-chain length (cold-start layout) |
| `L_max` | full length | final length (blank → full residue count) |
| `domain_def` | none | domain YAML for contact `n_scale` (one-time precompute) |
| `stride_output_file` | run STRIDE | precomputed STRIDE for the full PDB |
| `n_steps` | 1000 | integration steps per residue (constant schedule) |
| `dt` | 0.015 | time step (ps) |
| `ref_t` | 300 | temperature (K) |
| `tau_t` | 0.01 | Langevin friction / coupling constant (1/ps) |
| `nstout` | 50 | trajectory / log / checkpoint write frequency |
| `constraints` | None | `None` (flexible bonds; required for A→P seeding) or `AllBonds` |
| `restraint_k` | 83680 | C-terminus position-restraint constant (kJ/mol/nm² = 200 kcal/mol/Å²; tether-off mode) |
| `buffer` | 0.4 | new-residue placement offset beyond the A-anchor (nm) |
| `minimize` | yes | per-step energy minimization |
| `rigid_ribosome` | no | **v2**: append the rigid ribosome + its forces |
| `trna_tether` | yes | v2: tether C-terminus to the P-site tRNA (bond + orienting angle) |
| `tunnel_wall` | yes | v2: keep the chain at `x ≥ tunnel_wall_x0` (forward-only extrusion) |
| `tunnel_wall_x0` | ~1.05 | wall plane (nm) = the C-terminal-AA addition plane (PTC) ≈ P-anchor x + tether bond length |
| `tunnel_wall_k` | 8368 | wall force constant (kJ/mol/nm² = 20 kcal/mol/Å²) |
| `ptc_offset` | auto | hold/seed the C-terminus this far (+x) from the P-anchor bead (nm) |
| `nascent_only_output` | yes | v2: write only the nascent chain to trajectory/PSF/final |
| `post_elongation` | stallation | `ejection` or `stallation` (if `post_elongation_steps > 0`) |
| `post_elongation_steps` | 0 | steps for the post-elongation phase (0 = skip) |
| `device` | CPU | `CPU` or `GPU` |
| `ppn` | 1 | CPU threads (device = CPU) |
| `outdir` | `synth_out` | per-length outputs go to `<outdir>/L_<L>/` |

(choosing-parameters)=
## Choosing parameters

The table above is the quick reference; this section explains the physics behind
the parameters you actually tune and gives test-vs-production guidance.

**Inputs and contacts.** `pdb_file` is the source of *all* force-field geometry
(native distances, bonded terms, contact map) — coordinates are decoupled from it
(seeding comes from the previous length). `domain_def` applies per-domain `n_scale`
to contact strengths exactly as in `topo.mdrun`; `stride_output_file` only saves
the one-time STRIDE call. None of these change between lengths.

**Elongation schedule and timescale.** `n_steps` is the per-residue **dwell** — how
long each length runs before the next residue is added — and is the single most
important physical knob, because the model's central observable is the competition
between dwell and folding time. Keep it tiny for traceable first runs
(`n_steps = 1000–2000`); for production, map it to real elongation kinetics. *E.
coli* runs ~20 aa/s (~0.05 s/residue); O'Brien map this to a mean
**12.6 ns = 840,000 steps × 0.015 ps**, so production uses e.g.
`n_steps = 840_000` (a constant schedule here; a per-codon variable schedule is a
planned extension). `dt = 0.015` ps is the standard CG step; `L0`/`L_max` set the
window (see Workflow).

**Integrator.** `ref_t`, `tau_t`, and `nstout` behave as in any TOPO run.
`constraints` should stay `None`: flexible harmonic bonds are **required** for the
A→P seeding (a new bead is placed at the A-anchor but restrained to the P-anchor,
so the connecting bond must stretch during relaxation). Velocities are re-drawn
from the Boltzmann distribution at `ref_t` each length (the particle count changes,
so they cannot carry over from a checkpoint).

**Elongation mechanics.** `restraint_k` (default 83680 kJ/mol/nm² = 200 kcal/mol/Å²)
is the stiffness holding the C-terminus at the PTC in the tether-off fallback;
`buffer` (0.4 nm) is how far past the A-anchor the new bead is placed. The buffer
**must clear excluded volume** — too small and the new bead overlaps its
neighbours, giving a huge `(σ/r)¹²` force that ejects it and crashes the step; the
per-step `minimize` relaxes residual clashes. If a step explodes, increase
`buffer` first.

**Ribosome (v2).** `rigid_ribosome = yes` turns on the full RNC environment
(rigid scenery + excluded volume + electrostatics); `no` keeps only the two
anchors (v1, for validating loop mechanics). With the ribosome on:

- `trna_tether = yes` (recommended) tethers the C-terminus to the P-site tRNA
  O'Brien-style — a harmonic `CA(L)–tRNA:R` bond plus an orienting
  `CA(L-1)–CA(L)–tRNA:R` angle (whose 91.7°/130° basins are O'Brien's tether
  parameters) — aiming the chain down the tunnel. `no` falls back to the plain
  position restraint (`restraint_k`, with `ptc_offset` clearing the bead).
- `tunnel_wall = yes` (recommended) applies a one-sided plane
  `U = k·min(x − x0, 0)²` so every nascent bead stays at `x ≥ x0` and can only
  **extrude forward** (+x) rather than folding back into the ribosome interior.
  `tunnel_wall_x0` (~1.05 nm) is the C-terminal-AA addition plane (PTC ≈ P-anchor
  x + tether bond length); `tunnel_wall_k` (8368 kJ/mol/nm² = 20 kcal/mol/Å²) is
  its stiffness.
- `nascent_only_output = yes` writes just the nascent chain to the
  trajectory/PSF/final (the ribosome is static); the checkpoint still holds the
  full system.

**Post-elongation.** `post_elongation_steps > 0` enables a final continuation;
`post_elongation` chooses `ejection` (release tether) vs `stallation` (keep it).

**Hardware.** `device = CPU` with `ppn` threads suits tutorials; the v2 system's
~4,600 particles make `device = GPU` the right choice for production.

## Outputs

One folder per length, `<outdir>/L_<L>/` (and `ejection/` or `stallation/`):
`traj.dcd`, `traj.log`, `traj.psf`, `traj.chk` (full system), `traj_final.pdb`
(nascent; seeds the next length), `traj_runinfo.log`, `native_1_<L>.pdb`. In v2 with
`nascent_only_output = yes` the trajectory/PSF/final hold only the nascent chain
(the checkpoint still holds the full system).

## Movie

`topo-elongate-movie -o <outdir>` stitches the per-length trajectories — and any
post-elongation phase — into one growing-chain `movie.dcd` + `movie.psf` + a ready
`movie.tcl` (colors by residue id, hides not-yet-synthesized beads each frame).
Pass `--ribosome <trunc.pdb>` to overlay the static ribosome. Then `vmd -e
<outdir>/movie.tcl`.

(production-notes)=
## Production notes

- **Dwell time.** *E. coli* ≈ 20 aa/s → ~0.05 s/residue; O'Brien map this to a mean
  of **12.6 ns = 840,000 steps × 0.015 ps** (a per-codon exponential; here a
  constant `n_steps`). So production uses, e.g., `n_steps = 840_000`.
- **Replicas.** Synthesis is stochastic — run many independent trajectories and
  analyze the ensemble.
- **GPU.** The v2 system has ~4,600 particles; use `device = GPU` for production.

## Status / caveats

The runner and the v2 force field are implemented and run stably. The C-terminus
**tether's quantitative effect is not yet validated** (it is physically motivated
and taken from O'Brien's protocol, but whether it measurably improves extrusion vs.
a plain restraint needs replicas at production dwell). A per-codon variable schedule
and the full 3-stage elongation cycle are planned extensions; see `TODO.md`.
