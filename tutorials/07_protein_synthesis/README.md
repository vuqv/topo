# Tutorial 7 — Protein synthesis (nascent-chain elongation)

**Goal:** simulate a protein being **synthesized vectorially on the ribosome**
(N-terminus first, one residue at a time) and then **ejected**, instead of
refolding the whole chain in bulk. You will run the `topo.translation` elongation
driver, understand its inputs and outputs, and turn the result into a movie.

**Time:** the demo settings finish in **~1 minute** on a CPU (16 chain lengths +
an ejection run). Production runs are far longer — see *Production settings* below.

> This tutorial is the hands-on companion to the design/spec in
> [`topo/translation/DESIGN.md`](https://github.com/vuqv/topo/blob/main/topo/translation/DESIGN.md), the file map in
> [`FILES.md`](https://github.com/vuqv/topo/blob/main/topo/translation/FILES.md), and the module reference in
> [`topo/translation/README.md`](https://github.com/vuqv/topo/blob/main/topo/translation/README.md).

---

## The idea in one minute

In a cell, the ribosome builds a protein **N→C, one residue at a time**, and the
chain emerges through the **exit tunnel** and starts folding **before it is even
finished**. That can give different folding pathways and intermediates than
refolding the full-length chain (what Tutorials 1–6 do).

TOPO models this as a **sequence of standalone simulations**, one per nascent-chain
length `L` (the "rebuild-and-continue" protocol):

```
for L = L0, L0+1, ..., L_max:
    1. BUILD the length-L CG system: nascent residues 1..L (+ the rigid ribosome)
    2. SEED: residues 1..L-1 from step (L-1)'s final structure; new residue L at the A-site
    3. TETHER residue L to the P-site tRNA; RUN n_steps; SAVE
    4. L -> L+1  (this step's final structure seeds the next)
after L_max: optional EJECTION (release the tether) or STALLATION (keep it)
```

The chain grows at its **C-terminus** (held at the peptidyl-transferase center,
PTC) while the **N-terminus extrudes** down the tunnel and folds outside.

---

## Files in this folder

| File | Role |
|------|------|
| `P0CX28_clean.pdb` | The target protein (106-residue single domain) = the nascent chain. |
| `domain.yaml` | Domain definition + calibrated contact strength (same as Tutorial 1). |
| `P0CX28_clean_stride.dat` | Precomputed STRIDE (backbone H-bonds); skips running STRIDE. |
| `ribosome_trunc.pdb` | **Truncated coarse-grained 50S ribosome** (4,576 beads) — the exit-tunnel shell + tRNAs, X-aligned. Source of the P-/A-anchors and (in v2) the rigid scenery. |
| `elongate.ini` | The control file — everything below is configured here. |

> The ribosome was coarse-grained (`cg_ribosome.py`) and truncated around the exit
> tunnel (`truncate_ribosome.py`); see [`FILES.md`](https://github.com/vuqv/topo/blob/main/topo/translation/FILES.md).
> Beyond the 2 nm force cutoff it exerts no force on the chain, so truncating is
> exact, not an approximation.

---

## Key concepts

- **Two anchors from the ribosome.** The **P-anchor** = P-site tRNA residue-76 `R`
  bead (the held C-terminus / PTC); the **A-anchor** = A-site tRNA residue-76 `R`
  bead (where each new residue is delivered).
- **Build-once-subset contacts.** The native contact map is computed **once** on
  the full PDB; length `L` uses the top-left `L×L` block. STRIDE/heavy-atom analysis
  never re-run per length, and contacts to not-yet-made residues are simply absent.
- **One step = one residue.** Each step places the new residue at the A-site,
  tethers it to the P-site, and runs `n_steps` of Langevin dynamics — collapsing
  O'Brien's 3-stage elongation cycle into a single step.
- **Build steps v1 vs v2.**
  - **v1** (`rigid_ribosome = no`): nascent chain only; the ribosome supplies just
    the anchor *coordinates*. Validates the loop machinery.
  - **v2** (`rigid_ribosome = yes`, this tutorial): the rigid (mass-0) ribosome is
    added with ribosome↔chain **excluded volume + electrostatics**, the **tRNA
    tether** (bond + `CA–CA–tRNA` orienting angle), and a **planar tunnel wall**.
- **Tunnel wall.** A one-sided restraint `k·min(x − x0, 0)²` keeping the chain at
  `x ≥ x0` so it can only **extrude forward** (+x, toward the exit). Here
  `x0 = 1.05 nm` = the **C-terminal-AA addition plane** (the PTC, where each new
  residue is placed and tethered ≈ P-anchor x + tether bond length). *(O'Brien quote
  58 Å, but that is their coordinate frame.)*

---

## Run it

From **this folder** (paths in `elongate.ini` are relative to the working dir):

```bash
cd tutorials/07_protein_synthesis
topo-elongate -f elongate.ini
# equivalently:  python -m topo.translation -f elongate.ini
```

You'll see the run echo its configuration, then grow the chain length by length:

```
ribosome forces (v2): on (rigid); tether: tRNA bond+angle; tunnel wall: x>=1.05 nm; output: nascent-only
Rigid ribosome: 4576 beads ... (appended as mass-0 scenery; ribosome<->nascent forces on).
Elongating P0CX28_clean.pdb: L = 5 .. 20 (N_full = 106), 2000 steps/residue.
# Nascent length L = 5  (+ rigid ribosome)
...
=== Post-elongation: ejection (L = 20, 2000 steps, C-terminus restraint OFF) -> synth_out/ejection/ ===
Done. Ejection written to synth_out/ejection/
```

> **Prerequisites:** `topo` + OpenMM importable, and (if you don't pass
> `stride_output_file`) `stride` on your `PATH`. See the top-level
> [tutorials README](https://github.com/vuqv/topo/blob/main/tutorials/README.md).

---

## What you get

One self-contained folder **per length**, `synth_out/L_<L>/`, plus a
post-elongation folder (`ejection/` or `stallation/`):

| file | what it is |
|------|------------|
| `traj.dcd` | trajectory for this length (**nascent chain only** by default) |
| `traj.log` | energy / temperature log |
| `traj.psf` | topology (nascent only) for visualization/analysis |
| `traj.chk` | OpenMM checkpoint (the **full** system, for restart) |
| `traj_final.pdb` | final structure — **seeds the next length** |
| `traj_runinfo.log` | provenance (versions, hardware, timing) |
| `native_1_<L>.pdb` | the length-`L` native CA structure used to build this step |
| `seed.pdb` | the seeded starting coordinates (v1 only) |

Quick check that the chain grew one residue per step:

```bash
for L in 05 10 20; do echo -n "L=$L atoms: "; grep -c " CA " synth_out/L_0$L/traj_final.pdb; done
```

---

## Make a movie

Stitch the per-length trajectories (and the ejection) into one VMD-playable movie
that grows the chain N→C, overlaying the static ribosome:

```bash
topo-elongate-movie -o synth_out --ribosome ribosome_trunc.pdb
vmd -e synth_out/movie.tcl
```

It plays as *F* frames of L=5, then L=6, …, then the ejection frames. The script
colors the chain by residue id and hides the not-yet-synthesized beads each frame,
so you watch the chain emerge and (after the tether is released) leave the tunnel.
Press **Play** in VMD (or run `animate forward`).

---

## The control file, explained

The important `elongate.ini` options (full reference: the
[module README](https://github.com/vuqv/topo/blob/main/topo/translation/README.md) and
[Protein synthesis options](https://vuqv.github.io/topo/usage/protein_synthesis.html)).

> **Units:** OpenMM defaults — length **nm**, time **ps**, energy **kJ/mol**,
> temperature **K**, force constants **kJ/mol/nm²** (e.g. `restraint_k`,
> `tunnel_wall_k`).

| key | here | meaning |
|-----|------|---------|
| `L0`, `L_max` | 5, 20 | grow from length 5 to 20 (use blank `L_max` for the full 106) |
| `n_steps` | 2000 | steps per residue — **tiny for the demo** (production ≈ 840,000) |
| `rigid_ribosome` | yes | build step v2: append the rigid ribosome + its forces |
| `trna_tether` | yes | tether the C-terminus to the P-site tRNA (bond + orienting angle) |
| `tunnel_wall` | yes | keep the chain at `x ≥ tunnel_wall_x0` (forward-only extrusion) |
| `tunnel_wall_x0` | 1.05 | wall plane (nm) = the C-terminal-AA addition plane (PTC) |
| `nascent_only_output` | yes | don't write the static ribosome every frame (save storage) |
| `post_elongation` | ejection | after full length: release the tether (or `stallation` to keep it) |
| `post_elongation_steps` | 2000 | length of the ejection/stallation run (0 = skip) |

---

## Production settings

The demo uses `n_steps = 2000` so it finishes fast. For realistic kinetics:

- **Dwell time.** *E. coli* adds ~20 aa/s → ~0.05 s per residue. O'Brien map this
  to a mean of **12.6 ns = 840,000 steps × 0.015 ps** (a per-codon exponential
  distribution; here we use a constant `n_steps`). So set, e.g.:
  ```ini
  n_steps = 840_000
  ```
- **Full length:** leave `L_max` blank to synthesize all 106 residues.
- **Replicas.** Folding is stochastic — run **many independent trajectories** (vary
  the run; O'Brien use 50 per protein) and analyze the ensemble, not one run.
- **GPU.** Set `device = GPU` — the v2 system has ~4,600 particles.

A constant `n_steps` is the simplest schedule; a per-codon (variable) schedule is a
planned extension.

---

## Caveats (read before drawing conclusions)

- **The tutorial L-range is short and the dwell is tiny**, purely to run fast. At
  `n_steps = 2000` the chain barely equilibrates per length — do **not** read
  folding pathways off the demo.
- **The C-terminus tether's quantitative effect is not yet validated.** It is
  implemented from O'Brien's protocol (and is physically motivated — it aims the
  chain down the tunnel), but whether it measurably improves extrusion vs. a plain
  restraint has not been established with proper replicas at production dwell. See
  the open item in [`TODO.md`](https://github.com/vuqv/topo/blob/main/topo/translation/TODO.md).
- **Single-domain caveat.** `P0CX28` is one intertwined domain; its native contacts
  are N-terminally enriched, so it is a fine machinery test but not a showcase of
  multidomain co-translational folding.

---

## Where to go next

- [`DESIGN.md`](https://github.com/vuqv/topo/blob/main/topo/translation/DESIGN.md) — the full model/spec (force field,
  protocol decisions, references to O'Brien et al.).
- [`topo/translation/README.md`](https://github.com/vuqv/topo/blob/main/topo/translation/README.md) — every option,
  the v2 forces, the tether, the wall, the movie tool.
- [`FILES.md`](https://github.com/vuqv/topo/blob/main/topo/translation/FILES.md) — how the ribosome was
  coarse-grained/truncated and the build-step implementation notes.
