# Tutorial 10 — Continuous Synthesis Protocol (O'Brien CSP)

**Goal:** synthesize a protein on the ribosome with **realistic, codon-resolved
translation kinetics** — each residue is timed from its mRNA codon and added through
the **three sub-stages of the elongation cycle** (peptidyl transfer → translocation →
tRNA binding), exactly as in O'Brien's *Continuous Synthesis Protocol*
(`continuous_synthesis_v6.py`). This is the kinetic upgrade of Tutorial 7.

**Time:** the demo settings finish in **~30 s** on a CPU (8 residues × 3 stages + an
ejection run). Production runs are far longer — see *Production settings*.

> This tutorial is the hands-on companion to [`PLAN.md`](PLAN.md) (the porting design
> + what is done / remaining) and [`TASK.md`](TASK.md) (the task checklist), and to
> the module reference in
> [`topo/csp/README.md`](https://github.com/vuqv/topo/blob/main/topo/csp/README.md).

---

## What's new vs. Tutorial 7

Tutorial 7 (`topo-elongate`) grows the chain one residue per step at a **fixed
`n_steps`**, collapsing the elongation cycle into a single MD segment. Its own README
flags the gap: *"a per-codon (variable) schedule is a planned extension."* **This
tutorial is that extension.** `topo-csp` adds the O'Brien kinetics:

1. **Per-codon translation times.** The mRNA is split into codons; a `trans_times`
   table gives each residue its own mean in-vivo translation time.
2. **`scale_factor`.** Maps in-vivo seconds → in-silico nanoseconds → integration
   steps, so simulated dwell times track real translation speed.
3. **Three sub-stages per residue**, each run for an **exponentially-sampled** dwell
   time (first-passage-time sampling):

   | stage | biology | dwell-time mean | restraint |
   |-------|---------|-----------------|-----------|
   | 1 | peptidyl transfer (A-site delivery) | `time_stage_1` | C-terminus → **A-site** |
   | 2 | translocation begins | `time_stage_2` (+ traffic) | C-terminus → **A-site** |
   | 3 | tRNA binding / waiting | codon total − stage 1 − stage 2 | C-terminus → **P-site** |

   Switching the restraint **A→P** between stages 2 and 3 reproduces translocation.
4. **Ribosome traffic** (optional): an upstream-queue correction to the stage-2 time
   (needs O'Brien's external `ribosome_traffic` binary; falls back to no traffic).
5. **Ejection / dissociation** end phases.

Everything else — the coarse-grained model, the rigid ribosome scenery, the tunnel
wall, build-once-subset contacts — is **reused verbatim** from `topo.translation`.

---

## Files in this folder

| File | Role |
|------|------|
| `P0CX28_clean.pdb` | The target protein (106-residue single domain) = the nascent chain. |
| `domain.yaml` | Domain definition + calibrated contact strength (same as Tutorial 1/7). |
| `P0CX28_clean_stride.dat` | Precomputed STRIDE (backbone H-bonds); skips running STRIDE. |
| `ribosome_trunc.pdb` | Truncated CG 50S ribosome (4,576 beads) — exit-tunnel shell + tRNAs, X-aligned. Source of the P-/A-anchors and the rigid scenery. |
| `trans_times.txt` | **Codon → mean in-vivo translation time (s)** — the Fluitt *E. coli* table (310 K), copied from the O'Brien example. Generic genetic-code table, reusable for any protein. |
| `P0CX28_mrna.txt` | **mRNA for P0CX28** — the PDB sequence back-translated to one representative synonymous codon per residue (+ 1 stop). Exercises the codon-resolved kinetics. |
| `csp.ini` | The control file — everything below is configured here. |

---

## Run it

From **this folder** (paths in `csp.ini` are relative to the working dir):

```bash
cd tutorials/10_csp_obrien
topo-csp -f csp.ini
# equivalently:  python -m topo.csp -f csp.ini
```

You'll see it echo the config, load the ribosome, precompute contacts **once**, then
print a per-residue **kinetic schedule** and grow the chain stage by stage:

```
[ O'Brien continuous synthesis -- kinetic schedule ]
  timing mode: per-codon (mRNA); scale_factor=4.33129e+06; dt=0.015 ps
  TEST CLAMP: <= 667 steps/stage (~2001 steps/residue). Remove for production.

# Residue L = 5  codon CCG  (total in-vivo dwell 0.03809 s)
#   stage 1 peptidyl transfer : 0.0002043 s -> 667 steps
#   stage 2 translocation     : 0.005106 s -> 667 steps
#   stage 3 tRNA binding/wait : 0.04863 s -> 667 steps
...
=== Ejection (L = 12, 667 steps, restraint OFF) -> synth_out/ejection/ ===
```

> **Prerequisites:** `topo` + OpenMM importable, and (if you don't pass
> `stride_output_file`) `stride` on your `PATH`.

---

## What you get

One self-contained folder **per residue**, with a sub-folder **per stage**, plus an
`ejection/` (and optional `dissociation/`) folder:

```
synth_out/
├── L_005/
│   ├── stage_1/   peptidyl transfer   (traj.dcd/.log/.psf/.chk, traj_final.pdb, ...)
│   ├── stage_2/   translocation begins
│   └── stage_3/   tRNA binding  ← traj_final.pdb seeds the NEXT residue
├── L_006/ ... L_012/
└── ejection/      restraint released; chain leaves the tunnel
```

Each stage folder is a standalone topo run (same files as Tutorial 7): `traj.dcd`
(nascent-only by default), `traj.log`, `traj.psf`, `traj.chk` (full system, for
restart), `traj_final.pdb` (seeds the next stage/residue), `traj_runinfo.log`,
`native_1_<L>.pdb`.

Quick check that the chain grew one residue per step:

```bash
for L in 005 008 012; do echo -n "L=$L: "; grep -c " CA " synth_out/L_$L/stage_3/traj_final.pdb; done
# -> 5, 8, 12
```

---

## Make a movie (VMD)

Stitch the per-residue/-stage trajectories (and the ejection) into one VMD-playable
movie that grows the chain **stage by stage**, overlaying the static ribosome:

```bash
topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
vmd -e synth_out/movie.tcl
```

It discovers every segment in synthesis order — `L=5` stage 1, 2, 3, then `L=6`
stage 1, 2, 3, … then `ejection` — pads each frame up to the final length (parking
the not-yet-made beads far away and hiding them), and writes `movie.psf` /
`movie.dcd` / `movie.tcl`. So you watch each new residue appear at the A-site (stage
1), settle (stage 2) and translocate to the P-site (stage 3) for every codon, then
the finished chain leave the tunnel. The chain is colored by residue id; press
**Play** in VMD (or run `animate forward`).

> Use `--park cterm` to stack future beads on the C-terminus instead of hiding them,
> or drop `--ribosome` to omit the static scenery. This is the CSP analog of
> Tutorial 7's `topo-elongate-movie` (it shares the same stitching core, but reads
> the nested `L_<L>/stage_<s>/` layout).

---

## The control file, explained

`csp.ini` is a **superset of `elongate.ini`** (Tutorial 7) — same structure / MD /
ribosome keys, plus the O'Brien kinetic keys.

> **Units:** OpenMM defaults — length **nm**, time **ps**, energy **kJ/mol**,
> temperature **K**, force constants **kJ/mol/nm²**. Dwell times are in **seconds**.

| key | here | meaning |
|-----|------|---------|
| `mrna` | `P0CX28_mrna.txt` | mRNA: one codon per residue (+ 1 stop) |
| `trans_times` | `trans_times.txt` | codon → mean in-vivo time (s) |
| `scale_factor` | 4331293 | in-vivo seconds → in-silico ns compressor |
| `time_stage_1` | 0.00034 | mean peptidyl-transfer dwell (s) |
| `time_stage_2` | 0.004201 | mean translocation dwell (s) |
| `ribosome_traffic` | no | apply the external traffic correction (off → no traffic) |
| `random_seed` | 20240628 | reproducible first-passage-time sampling |
| `max_steps_per_stage` | 667 | **TEST CLAMP** — cap each stage (~2000 steps/residue) |
| `min_steps_per_stage` | 50 | floor each stage's step count |
| `L0`, `L_max` | 5, 12 | grow from length 5 to 12 (blank `L_max` for the full 106) |
| `rigid_ribosome` | yes | build step v2: rigid ribosome + its forces (as Tutorial 7) |
| `tunnel_wall` | yes | keep the chain at `x ≥ 1.05 nm` (forward-only extrusion) |
| `ejection_steps` | 667 | post-synthesis free run (release the restraint) |
| `dissociation_steps` | 0 | optional second free run (drift off the ribosome) |

For the uniform-time mode (no mRNA), set `uniform_ta = yes` and `uniform_mfpt =
<seconds>` instead of `mrna` / `trans_times` (O'Brien's `uniform_ta = 1`).

---

## Production settings

The demo **caps** every stage at 667 steps (`max_steps_per_stage`) so it finishes in
seconds. For realistic kinetics, **delete the two `*_steps_per_stage` lines** so the
step counts follow the codon timing:

- A typical *E. coli* residue (~0.05 s) maps to ~12.6 ns ≈ **840,000 steps** at
  `dt = 0.015 ps` (split across the three stages by their dwell-time fractions).
- **Full length:** leave `L_max` blank to synthesize all 106 residues.
- **Replicas:** folding is stochastic — run **many independent trajectories** (vary
  `random_seed`; O'Brien use ~50 per protein) and analyze the ensemble.
- **GPU:** set `device = GPU` — the v2 system has ~4,600 particles.
- **Ribosome traffic:** set `ribosome_traffic = yes` if O'Brien's `ribosome_traffic`
  binary is on your `PATH`; otherwise it logs a warning and runs with no traffic.

---

## Caveats (read before drawing conclusions)

- **The demo is clamped and short**, purely to run fast. At ≤667 steps/stage the
  chain barely equilibrates — do **not** read folding pathways off the demo.
- **The synthetic mRNA is representative, not biological.** `P0CX28_mrna.txt` is the
  PDB sequence back-translated to one fixed codon per amino acid, so the *kinetics
  machinery* is exercised but the specific per-residue times are illustrative. Supply
  the real transcript for a real study.
- **3-stage mechanics are mapped, not literal.** The A→P translocation is reproduced
  by switching the C-terminus position-restraint target; the peptide bond is present
  in the bonded model from stage 1 (not toggled mid-residue), and explicit A/P tRNA
  bonded geometry is not modeled. See [`PLAN.md`](PLAN.md) *Remaining*.
- **Single-domain caveat.** `P0CX28` is one intertwined domain — a good machinery
  test, not a showcase of multidomain co-translational folding.

---

## Where to go next

- [`PLAN.md`](PLAN.md) — the porting design and the **done / remaining** matrix.
- [`TASK.md`](TASK.md) — the task checklist.
- [`topo/csp/README.md`](https://github.com/vuqv/topo/blob/main/topo/csp/README.md) —
  every option, the kinetics API, the 3-stage mapping.
- Tutorial 7 — the fixed-`n_steps` elongation driver this builds on.
