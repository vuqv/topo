# Co-translational synthesis through an analytic exit tunnel (topo)

This is the **structure-based (topo) twin** of cosmo's tutorial 09. The full CSP
runner (Tutorials 12/13) builds the ribosome from **explicit beads** and threads the
nascent chain through the coarse-grained exit tunnel. Here the tunnel is modelled
**analytically** instead: a cylindrical **bore** of radius `r` along the X-axis
drilled through an **infinite wall** at `x_exit` (a "hole in an infinite wall").
There are **no ribosome beads** — the System is the nascent chain only, so it is
fast and never jams.

The chain is a **folded protein** built with topo's structure-based Gō contacts, so
it can **fold co-translationally** as it extrudes and once it clears the bore.

## Why a separate script (`cylinder.py`)

The analytic tunnel is a different *physics of confinement* than the explicit-bead
ribosome, so it lives here as a self-contained tutorial and does **not** modify the
shipped `topo.csp` package. It **reuses** that package's tested, unchanged low-level
machinery from `topo.csp.core` (the one-time contact precompute, the
build-once-subset length model, the seed / restrain / output path) and adds only the
one new force, `add_tunnel_cylinder`, plus its own minimal nascent-only elongation
loop. (Like the explicit-bead path it grows the chain at a fixed step count — fine
for this confinement demo; for codon-resolved kinetics use CSP, `topo-csp`.)

## The model (forbidden region `S`)

```
   d                              (cytosol: free, any d)
   ^   |##### solid ribosome S #####|
 r |···|············ bore ··········|··············>  allowed past exit
   +---|----------------------------|----------------> x
     x_lo (PTC)                  x_exit
       |##### solid ribosome S #####|
                              ^ infinite exit-face wall (d > r)
```

A bead is penalised by its **penetration depth into `S`** (everything outside the
bore up to the exit face, plus the closed PTC end), escaping via whichever face is
nearer — the bore wall (radial inward push → keeps the chain extended) or the exit
face (`+x` push → a cytosol bead can re-enter the tunnel **only through the bore**).
The 90° mouth corner is rounded by a fillet (radius `rho`) so the potential is
continuous and the MD is stable. The C-terminus is seeded and position-restrained
on the tunnel axis at the PTC `(x_lo, 0, 0)`; new residues are seeded there.

## Run

```bash
cd tutorials/09_translation_cylinder
python cylinder.py -f elongate.ini
```

All parameters live in [`elongate.ini`](elongate.ini) (`[OPTIONS]` section). The
nascent chain is the 106-residue P0CX28 (same as tutorial 07), with `domain.yaml` +
the precomputed STRIDE for the contact map. Tunnel defaults: bore radius **0.9 nm**,
length **10.0 nm** (`x_lo=0`, `x_exit=10`), axis on X, mouth fillet **0.2 nm**,
wall stiffness **8368 kJ/mol/nm²**.

## Post-elongation: ejection

Once the chain reaches its final length, an optional **post-elongation** phase
continues the same system from the finished structure (written to
`<outdir>/<phase>/`):

- **`ejection`** — releases the C-terminus restraint, so the completed protein is
  free to diffuse. The analytic tunnel stays on (bore + closed PTC end + exit
  wall), so the only way out is +x through the exit. This tests whether the
  nascent protein **diffuses out of the tunnel and folds in the cytosol**.
- **`stallation`** — keeps the restraint, so the chain stays threaded/stalled.

```ini
post_elongation       = ejection
post_elongation_steps = 300_000   # use a LONG run so the protein can clear the tunnel
```

## Visualize / validate

Use the tutorial's movie tool — it stitches the per-length trajectories (reusing the
shared stitcher in `topo.csp.movie`) **and** draws the analytic tunnel (bore
tube, closed PTC cap, and the infinite exit-face wall as an annulus whose hole is the
bore), reading the geometry from the same `elongate.ini`:

```bash
python make_movie_cylinder.py -o synth_out -f elongate.ini
vmd -e synth_out/movie.tcl
```

You then see the chain thread the (blue, transparent) bore, the red PTC end cap it
grows away from, and the grey exit wall it emerges through — then fold once it
clears the tunnel. (Plain `topo-csp-movie -o synth_out` also works
— it just omits the tunnel.)
