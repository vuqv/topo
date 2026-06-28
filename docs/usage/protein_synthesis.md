# Protein synthesis (`topo.translation`)

The `topo.translation` module synthesizes a protein **vectorially on the ribosome**
(N→C, one residue at a time) and optionally **ejects** it, instead of refolding the
full-length chain in bulk. It is the third runner in the package, alongside
`topo.mdrun` and `topo.optimize`.

- **CLI:** `topo-elongate -f elongate.ini` (or `python -m topo.translation -f elongate.ini`)
- **Movie tool:** `topo-elongate-movie -o <out_root> [--ribosome ribo.pdb]`
- **Walkthrough:** Tutorial 7.
- **Design/spec & implementation notes:** `topo/translation/DESIGN.md`, `FILES.md`.

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
