# Visualizing the synthesis process

A synthesis run — whether an interrupted-and-{doc}`resumed <synthesis_resume>` production
run or a short demo — writes the growing chain as **many short trajectories**, one per
residue (and, for {doc}`continuous_synthesis`, three MD sub-stages per residue). Different
lengths have different bead counts, so the segments cannot simply be concatenated. The
package stitches them into a single VMD-playable movie that *grows* the chain N→C, and
adds the scenery that makes the process legible: the **analytic exit tunnel** for the
cylinder model, or the **coarse-grained ribosome** for the explicit-ribosome model.

This page covers turning a finished `<outdir>` into an interactive VMD movie you can rotate
and replay, for **both** synthesis runners.

```{note}
Everything below works the same on a resumed run: `topo-csp-movie` reads whatever
`L_<L>/` segments are present under `<outdir>`, so stitch after the final residue completes.
```

## Stitch the movie (`topo-csp-movie`)

`topo-csp-movie` (`= python -m topo.csp.movie`) auto-detects the on-disk layout — the
3-stage CSP layout (`<outdir>/L_<L>/traj_s<1,2,3>.dcd`) or the flat per-length layout
(`<outdir>/L_<L>/traj.dcd`, used by the cylinder runner) — plus any `ejection/`
phase, and writes:

```text
<outdir>/movie.psf     # single topology (padded to the final length)
<outdir>/movie.dcd     # all segments back-to-back, in synthesis order
<outdir>/movie.tcl     # a VMD script that loads the movie and grows the chain
```

The not-yet-synthesized beads are parked at a far sentinel each frame; the generated
`movie.tcl` hides them per frame (the default `--park sentinel` scheme).

```{important}
`movie.tcl` loads `movie.psf` / `movie.dcd` (and, for the ribosome model,
`movie_ribosome.pdb`) **by basename**, so open it *from inside* `<outdir>`:

    cd <outdir> && vmd -e movie.tcl

Running `vmd -e <outdir>/movie.tcl` from the parent folder loads the scenery but fails to
find the trajectory.
```

Add the model-specific scenery with one flag, described next.

## Cylinder tunnel model (`topo-cylinder`)

The cylinder run has **no ribosome beads** — the exit tunnel is a pure force. Pass the same
`cylinder.ini` that drove the run so the drawn tunnel (bore tube, closed PTC end cap, and
the infinite exit-face wall as an annulus whose hole is the bore) matches the forces the
chain actually felt.

```bash
# Interactive: stitch + draw the tunnel into movie.tcl, then view from the run folder.
topo-csp-movie -o synth_out --tunnel cylinder.ini      # or: python -m topo.csp.movie -o synth_out --tunnel cylinder.ini
cd synth_out && vmd -e movie.tcl
```

You see the chain thread the (blue, transparent) bore, the red PTC cap it grows away from,
and the grey exit wall it emerges through — then fold once it clears the tunnel.
`--wall-outer NM` sets the drawn exit-wall radius (default: bore radius + 3 nm). The
in-folder `make_movie_cylinder.py` is an equivalent self-contained wrapper; both draw the
same geometry via `topo.csp.cylinder.tunnel_tcl`.

## Coarse-grained ribosome model (`topo-csp`)

The explicit-ribosome run synthesizes on TOPO's own CG ribosome. Overlay it as static
scenery by passing the truncated ribosome PDB:

```bash
# Interactive: stitch + overlay the ribosome, then view from the run folder.
topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
cd synth_out && vmd -e movie.tcl        # also uses movie_ribosome.pdb, written alongside
```

The ribosome loads as a separate, static molecule (a translucent silver bead cloud) that
the chain grows through; `topo-csp-movie` copies it to `<outdir>/movie_ribosome.pdb` so the
movie is self-contained.

## Notes & tips

- **Plain movie.** Omit `--tunnel` / `--ribosome` for a bare growing-chain movie (no
  scenery) — useful as a quick check.
- **`--prefix`** renames the output set (`<prefix>.{psf,dcd,tcl}`); **`--outname`** must
  match the per-segment basename the runner used (default `traj`).
- **Parking scheme.** `--park cterm` stacks future beads on the C-terminus instead of the
  far sentinel (no per-frame hiding needed, but leaves a small bead cluster); `sentinel`
  (default) is cleaner.
- **Post-synthesis phase.** The `ejection/` free run is appended automatically when
  present, so the movie continues past full length into release.

## See also

- {doc}`cylinder_synthesis` — the analytic-tunnel runner (`topo-cylinder`).
- {doc}`continuous_synthesis` — the explicit CG-ribosome runner (`topo-csp`).
- {doc}`synthesis_resume` — operating and resuming a long production run before you stitch.
