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
python make_movie_cylinder.py -o synth_out -f cylinder.ini
vmd -e synth_out/movie.tcl
```

`topo-cylinder` writes, per residue `L`, a standalone trajectory under
`<outdir>/L_<L>/`, optional post-synthesis free runs (`ejection/` and `dissociation/`),
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

The C-terminus is **seeded and position-restrained on the tunnel axis at the PTC**
`(x_lo, y0, z0)` (stiffness `restraint_k`); each new residue is seeded there. There is
no A/P tRNA tether and no translocation switch — the chain simply extrudes forward as
it grows. Because the chain carries TOPO's native Gō contacts, it folds
co-translationally once residues clear the bore.

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
The post-synthesis phases use the **same** `ejection_steps` / `dissociation_steps` keys
as the CSP runner.
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

## Configuration reference (`cylinder.ini`)

The cylinder reads a single INI control file with one `[OPTIONS]` section
(`topo.csp.cylinder.read_cylinder_config`). **Units are OpenMM defaults** — nm, ps,
kJ/mol, K, kJ/mol/nm² — and **dwell times are in seconds**. Integers may use `_` digit
separators. The kinetics and MD keys share the semantics of the
{doc}`synthesis control options <synthesis_control>` page; the analytic-tunnel geometry
and the post-elongation phase are specific to this runner.

Example `cylinder.ini` (Tutorial 7):

```ini
[OPTIONS]
; --- inputs (no `ribosome` PDB: the tunnel is analytic) ---
pdb_file           = P0CX28_clean.pdb           ; full native PDB of the nascent chain
domain_def         = domain.yaml                ; contact nscale (one-time precompute)
stride_output_file = P0CX28_clean_stride.dat    ; optional; else STRIDE is run for you

; --- length schedule ---
L0    = 1            ; starting nascent length (required)
L_max =              ; final length (blank -> full residue count)

; --- kinetics (same O'Brien codon timing as topo-csp) ---
mrna         = P0CX28_mrna.txt   ; one codon per residue (required for per-codon timing)
; codon_times omitted -> bundled E. coli 310 K table (or set a path, or a number of s = uniform)
scale_factor = 216564650         ; in-vivo s -> in-silico ns compression (larger = faster)
random_seed  = 20240629
max_steps_per_stage = 2000       ; TEST CLAMP (delete for production)
min_steps_per_stage = 50

; --- mechanics / integrator ---
constraints = None   ; flexible harmonic bonds (required for the seeding)
restraint_k = 83680  ; C-terminus -> PTC restraint (kJ/mol/nm^2)
minimize    = yes
dt     = 0.015
ref_t  = 300
tau_t  = 0.05
nstout = 100

; --- analytic exit tunnel ---
tunnel_radius      = 0.9        ; bore radius r (nm); ~3 CG beads wide
tunnel_length      = 10.0       ; bore length (nm); x_exit = x_lo + length
tunnel_x_lo        = 0.0        ; PTC / closed end (nm); C-terminus seeded on-axis here
tunnel_center      = 0.0, 0.0   ; tunnel axis (y0, z0) (nm); axis = X
tunnel_k           = 8368       ; wall stiffness (kJ/mol/nm^2 = 20 kcal/mol/A^2)
tunnel_mouth_round = 0.2        ; mouth-corner fillet radius rho (nm)

; --- post-synthesis free runs (after the chain reaches full length) ---
ejection_steps     = 300_000    ; release the C-terminus restraint; protein diffuses out (0 -> skip)
dissociation_steps = 0          ; continued free run; protein drifts off the ribosome (0 -> skip)

; --- hardware / output ---
device = GPU
ppn    = 4
outdir = synth_out
```

### Inputs & length schedule

| Key | Required | Default | Meaning |
|-----|----------|---------|---------|
| `pdb_file` | **yes** | — | Full native PDB of the target protein; the CG model is built from it. |
| `domain_def` | **yes** | — | `domain.yaml` — the protein's per-domain/per-interface contact-nscale definition. See {doc}`domain_definition`. |
| `stride_output_file` | no | — | Precomputed STRIDE file (skips re-running STRIDE). |
| `L0` | **yes** | — | Starting nascent-chain length (cold-start layout). |
| `L_max` | no | full length | Final nascent length (blank = whole chain). Must satisfy `1 ≤ L0 ≤ L_max ≤ N_full`. |
| `mrna` | for per-codon timing | — | mRNA file (one codon per residue); required unless `codon_times` is a number. |
| `codon_times` | no | bundled E. coli 310 K | Codon-timing key: a **table path** = per-codon; a **positive number of seconds** = uniform codon time (no `mrna` needed); omit = bundled Fluitt *E. coli* table. A table filename must **not** be a bare number. |
| `outdir` | no | `synth_out` | Output root; each residue writes `L_<L>/`. |

### Analytic tunnel geometry

| Key | Default | Meaning |
|-----|---------|---------|
| `tunnel_radius` | `0.9` | Bore radius `r` (nm); ~3 CG beads wide. |
| `tunnel_length` | `10.0` | Bore length (nm). The exit face sits at `x_exit = tunnel_x_lo + tunnel_length` (derived, not a key). |
| `tunnel_x_lo` | `0.0` | PTC / closed end of the bore (nm); the C-terminus is seeded on-axis here. |
| `tunnel_center` | `0.0, 0.0` | Tunnel axis `(y0, z0)` (nm). The axis runs along +x. |
| `tunnel_k` | `8368` | Wall stiffness (kJ/mol/nm² = 20 kcal/mol/Å²). |
| `tunnel_mouth_round` | `0.2` | Mouth-corner fillet radius `rho` (nm); rounds the 90° inner corner so the potential is continuous. |

### Post-synthesis free runs

Same keys as the CSP runner. Both phases run at full length and **release** the
C-terminus restraint so the finished protein diffuses out the exit (+x) and folds in the
cytosol; the analytic tunnel stays on throughout, so the only way out is the exit face.

| Key | Default | Meaning |
|-----|---------|---------|
| `ejection_steps` | `0` | Steps of the first free run (restraint OFF); `0` = skip. Use a **long** run so the protein can clear the tunnel. Writes `ejection/`. |
| `dissociation_steps` | `0` | Steps of a second, continued free run (restraint OFF); `0` = skip. Writes `dissociation/`. |

### Shared kinetics & MD keys

These behave exactly as documented on the {doc}`synthesis control options <synthesis_control>`
page (they are inherited from the shared `RunParams`). The clamps
`max_steps_per_stage` / `min_steps_per_stage` are **testing-only** and clamp the
per-residue step count.

| Key | Default | Meaning |
|-----|---------|---------|
| `scale_factor` | `4331293` | In-vivo-seconds → in-silico-ns compression (larger = fewer steps = faster). |
| `codon_times` | bundled E. coli 310 K | Table path = per-codon timing; a positive number of seconds = uniform codon time (no `mrna` needed). See Inputs above. |
| `random_seed` | — | Seed for the first-passage-time sampler (reproducible schedule). |
| `max_steps_per_stage` | — (uncapped) | **Testing only** — upper clamp on the per-residue MD step count. Leave unset in production. |
| `min_steps_per_stage` | `1` | **Testing only** — lower clamp on the per-residue MD step count. |
| `constraints` | `None` | Bond treatment; the cylinder needs flexible bonds — leave `None`. |
| `restraint_k` | `83680` | C-terminus → PTC harmonic restraint constant (kJ/mol/nm²). |
| `minimize` | `yes` | Energy-minimize the seeded structure before each residue's MD. |
| `dt` / `ref_t` / `tau_t` / `nstout` | `0.015` / `310` / `0.05` / `5000` | Timestep (ps), temperature (K), Langevin friction (1/ps), output interval (steps). |
| `device` / `ppn` | `CPU` / `1` | Compute platform and CPU thread count. |

```{warning}
`time_stage_1` / `time_stage_2` are accepted (inherited from the shared parameters) but
have **no effect** in the cylinder runner: with a single MD segment per residue, each
residue's step count comes from its **whole** codon dwell, not a three-way split. They
matter only for the {doc}`coarse-grained ribosome runner <continuous_synthesis>`.
```

```{note}
Boolean options accept `yes`/`no`, `true`/`false`, `1`/`0`.
```

---

## Outputs

```text
<outdir>/
├── L_<L>/                  # one folder per residue L (single MD segment)
│   ├── traj.dcd            # (nascent-only) trajectory for that length
│   ├── traj_final.pdb      # last conformation (seeds the next residue)
│   ├── traj.log            # energies
│   └── ...
├── ejection/               # post-synthesis free run (if ejection_steps > 0)
├── dissociation/           # continued free run (if dissociation_steps > 0)
└── dwell_times.dat         # per-residue: codon, sampled dwell (s), ns, integration steps
```

**Movie.** `make_movie_cylinder.py` stitches the per-length trajectories (reusing the
shared stitcher in `topo.csp.movie`) **and** draws the analytic tunnel — the bore tube,
the closed PTC cap, and the infinite exit-face wall as an annulus whose hole is the bore
— reading the geometry from the same `cylinder.ini`. Plain `topo-csp-movie -o synth_out`
also works; it just omits the tunnel.

---

## See also

- {doc}`continuous_synthesis` — the explicit coarse-grained ribosome runner (`topo-csp`):
  the same codon kinetics with three MD sub-stages per residue and an explicit-bead exit
  tunnel.
- {doc}`synthesis_control` — the shared kinetics / MD control options in full.
- {doc}`model_theory` — the TOPO Gō-model force field the nascent chain carries.
