# Synthesis in cylinder ribosome model

The **cylinder** runner is the analytic-tunnel variant of protein synthesis: it
replaces the explicit coarse-grained ribosome of the
{doc}`coarse-grained ribosome model <continuous_synthesis>` with a **cylindrical
bore drilled through an infinite wall** (a "hole in an infinite wall"). There are
**no ribosome beads** — the simulated system is the nascent chain only — so it is
fast, never jams on ribosome excluded volume, and the analytic tunnel keeps the
in-tunnel segment extended so the chain threads out the exit and folds
**co-translationally** once it clears the bore.

- **CLI:** `topo-cylinder -f cylinder.ini` (or `python -m topo.csp.cylinder -f cylinder.ini`)
- **Worked example:** Tutorial 7 (`tutorials/07_translation_cylinder/`, the 106-residue P0CX28).
- **Module:** [`topo.csp.cylinder`](../topo) — a parallel module to the explicit-bead
  `topo.csp.protocol`. It **reuses** the shared low-level engine `topo.csp.core`
  (one-time contact precompute, build-once-subset length model, seed / restrain /
  output path) and the timing core `topo.csp.kinetics`, adding only the one analytic
  tunnel force (`add_tunnel_cylinder`) and a nascent-only synthesis loop.

```{tip}
Same **codon kinetics** as the coarse-grained-ribosome runner (`topo-csp`) — each
residue's MD length comes from its codon dwell time — but the cylinder runs **one MD
segment per residue** rather than three sub-stages, because the analytic tunnel has no
A→P translocation to model. For the full kinetics derivation (codon → seconds →
integration steps) see the {doc}`coarse-grained ribosome page <continuous_synthesis>`;
this page describes the tunnel model and the `cylinder.ini` options.
```

---

## Quick start

All paths in the INI are relative to the working directory; run from the tutorial
folder.

```bash
cd tutorials/07_translation_cylinder
topo-cylinder -f cylinder.ini          # -> synth_out/

# stitch the per-length trajectories into a movie AND draw the analytic tunnel
topo-csp-movie -o synth_out --tunnel cylinder.ini
cd synth_out && vmd -e movie.tcl        # movie.tcl loads its files by basename
```

`topo-cylinder` writes, per residue `L`, a standalone trajectory under
`<outdir>/L_<L>/`, optional post-synthesis phases (`stall/`, `ejection/`),
and a per-residue dwell-time log `<outdir>/dwell_times.dat`.

---

## The model: analytic exit tunnel

The tunnel is a cylindrical **bore** of radius `r` along the X-axis, drilled through
an **infinite wall** that spans from the closed PTC end (`x_lo`) to the exit face
(`x_exit = x_lo + tunnel_length`):

```text
   d                              (cytosol: free, any d)
   ^   |##### solid ribosome S #####|
 r |···|············ bore ··········|··············>  allowed past exit
   +---|----------------------------|----------------> x
     x_lo (PTC)                  x_exit
       |##### solid ribosome S #####|
                              ^ infinite exit-face wall (d > r)
```

A single `CustomExternalForce` over every nascent bead penalises its **penetration
depth into the solid region `S`** — everything outside the bore up to the exit face,
plus the closed PTC end:

```text
S   = { x < x_exit AND d > r } ∪ { x < x_lo },   d = |(y,z) − (y0,z0)|
U   = k·max(0, pen)²  +  k·min(0, x − x_lo)²
pen = (rounded) min( x_exit − x , d − r )        # 0 outside S; > 0 inside S
```

A bead escapes `S` via whichever face is nearer — the **bore wall** (`d − r`, a radial
inward push that keeps the in-tunnel chain extended) or the **exit face** (`x_exit − x`,
a +x push), so a bead in the cytosol can only re-enter the tunnel **through the bore**,
never off-axis. The 90° inner corner at the mouth is rounded by a fillet of radius
`rho = tunnel_mouth_round` so the potential stays continuous and the MD stays stable.

The C-terminus is **position-restrained on the tunnel axis at the PTC** `(x_lo, y0, z0)`
(stiffness `restraint_k`). Each new residue is *seeded* one equilibrium peptide bond
(0.381 nm) **deeper** than that rest point — toward the closed PTC end — so the new
`L-1↔L` bond starts at its equilibrium length instead of collapsed onto the previous
C-terminus; the restraint and the closed-end wall then lift the new residue back onto the
PTC while the older chain ratchets forward. (This is the analytic-tunnel analogue of the
explicit protocol's A/P-site offset, without the ribosome.) There is no A/P tRNA tether
and no translocation switch — the chain simply extrudes forward as it grows. Because the
chain carries TOPO's native Gō contacts, it folds co-translationally once residues clear
the bore.

```{note}
**How it differs from the {doc}`coarse-grained ribosome model <continuous_synthesis>`.**
(1) No `ribosome` PDB — the tunnel is analytic, its geometry set by the `tunnel_*` keys.
(2) **One MD segment per residue.** There are no peptidyl-transfer / translocation /
tRNA-binding sub-stages, so the whole codon dwell — the full per-codon time
`intrinsic[L]` from the codon-time table — is run as a **single** MD segment. This
gives the **same total simulated time per residue** as the explicit runner (whose
three sub-stages, by construction, sum in the mean to that same `intrinsic[L]`); the
cylinder just does not slice it into three. Consequently `time_stage_1` / `time_stage_2`
are **inert** here (they only *partition* the codon dwell in the explicit model, they
do not change its total), and the **`ribosome_traffic` / `initiation_rate`** keys are
**also not applied** — the cylinder runs the un-corrected `intrinsic[L]`, so no
inter-ribosome traffic correction is added. Set those keys in `cylinder.ini` and they
are silently ignored. (3) The ribosome-specific knobs (`nascent_ev_radii`,
`trna_tether`, `tunnel_wall`) do not apply.
The post-synthesis phases use the **same** `stall_steps` (hold at the PTC, restraint on)
then `ejection_steps` (free run) keys as the CSP runner.
```

### Physical scope — what the bore captures, and what it omits

The analytic tunnel reproduces **geometric confinement only**. Because there are no
ribosome beads, the cylinder model **omits** two interactions that the
{doc}`explicit coarse-grained ribosome model <continuous_synthesis>` includes:

- **Ribosome↔nascent-chain electrostatics.** There is no Debye–Hückel term between
  the ribosome surface (negatively charged rRNA phosphates, charged ribosomal-protein
  residues) and the charged nascent residues. The real exit tunnel is strongly
  electronegative; that field is entirely absent here.
- **Ribosome-surface excluded volume.** The wall is a smooth bore of one fixed radius,
  not the real ribosome's per-residue `(σ/r)¹²` surface. There is **no constriction
  site and no vestibule** — the *E. coli* tunnel is curved, ~8 nm long, and varies
  ~0.5–1 nm in radius, none of which a straight uniform cylinder reproduces.

It also has **no A→P translocation** — the C-terminus stays position-restrained at the
PTC, whereas the explicit model migrates it one register per residue.

```{warning}
**The cylinder and explicit-ribosome runners are comparable only in the *mean*
per-residue dwell time, not in confinement chemistry.** Use the cylinder model for
fast exploration of how *tunnel geometry + codon kinetics* shape co-translational
folding; use the {doc}`explicit ribosome <continuous_synthesis>` when tunnel-wall
charge, tunnel shape (constriction / vestibule), or translocation-coupled forces
matter. Do **not** compare folding observables (folding order, Q-vs-length, radius of
gyration) between the two models without accounting for these missing terms.
```

---

## Configuration

The cylinder reads a single INI control file with one `[OPTIONS]` section
(`topo.csp.cylinder.read_cylinder_config`). **Every control key is documented in one
place:** {doc}`synthesis_control` — the *shared* keys (inputs, kinetics, integrator,
post-synthesis, `resume`) plus the **Cylinder-runner-only** `tunnel_*` geometry keys
(`tunnel_radius`, `tunnel_length`, `tunnel_x_lo`, `tunnel_center`, `tunnel_k`,
`tunnel_mouth_round`) that define the analytic bore. A runnable `cylinder.ini` example is
on that page and in `tutorials/07_translation_cylinder/`.

Two cylinder-specific points: there is **no `ribosome` PDB** (the tunnel is analytic), and
`time_stage_1` / `time_stage_2` are accepted but have **no effect** — with a single MD
segment per residue the whole codon dwell `τ` is one segment, not a three-way split.

---

## Outputs

```text
<outdir>/
├── L_<L>/                  # one folder per residue L (single MD segment)
│   ├── traj.dcd            # (nascent-only) trajectory for that length
│   ├── traj_final.pdb      # last conformation (seeds the next residue)
│   ├── traj.log            # energies
│   └── ...
├── stall/                  # post-synthesis hold at PTC, restraint on (if stall_steps > 0)
├── ejection/               # post-synthesis free run (if ejection_steps > 0)
├── dwell_times.dat         # per-residue: codon, sampled dwell (s), ns, integration steps
└── progress.log            # append-only DONE/RUNNING resume status
```

**Resuming.** Like CSP, an interrupted cylinder run continues from the last completed
residue when re-invoked (`resume = auto`, on by default) — the schedule is re-read from
`dwell_times.dat` and the seed reloaded from the last `L_<L>/traj_final.pdb`, tracked by
`progress.log`. See {doc}`synthesis_resume`.

**Movie.** `topo-csp-movie -o synth_out --tunnel cylinder.ini` stitches the per-length
trajectories **and** draws the analytic tunnel — the bore tube, the closed PTC cap, and the
infinite exit-face wall as an annulus whose hole is the bore — reading the geometry from the
same `cylinder.ini` (the in-folder `make_movie_cylinder.py` is an equivalent wrapper). Plain
`topo-csp-movie -o synth_out` also works; it just omits the tunnel. For the interactive vs.
headless recipes and the `cd <outdir>` gotcha, see {doc}`synthesis_visualization`.

---

## See also

- {doc}`continuous_synthesis` — the explicit coarse-grained ribosome runner (`topo-csp`):
  the same codon kinetics with three MD sub-stages per residue and an explicit-bead exit
  tunnel.
- {doc}`synthesis_control` — the shared kinetics / MD control options in full.
- {doc}`model_theory` — the TOPO Gō-model force field the nascent chain carries.
