# Tutorial 13 — Full-length validation of the Tutorial-12 stability fix

This tutorial is an **independent, full-length stress test** of the per-stage
**stability guard** introduced for [Tutorial 12](../12_auto/) (`12_auto`). Tutorial 12
reproduced O'Brien's continuous-synthesis protocol on **4c5c** at the reference length
(L = 1 → 10) and fixed the energy blow-up that broke `tutorials/10_csp_obrien` (see
[`12_auto/WHY_10_FAILS.md`](../12_auto/WHY_10_FAILS.md)). The open question was whether
that fix holds **at scale** — for the **entire 306-residue chain**, which is exactly
the regime where the original Tutorial 10 blew up (≈5 of 306 stages reaching
PotE ≈ 10¹³ kJ/mol).

**It works as expected.** Synthesizing all 306 residues with `topo-csp`:

| check | result |
|-------|--------|
| residues synthesized | **1 → 306** (all, + ejection) |
| stage trajectories written | **919** (≈306 × 3 stages + ejection) |
| stages with PotE > 10⁹ kJ/mol (blow-ups) | **0** |
| worst max\|PotE\| over the whole run | **1755 kJ/mol** (the full-length ejection — finite and sane) |
| per-residue dwell-time log | `synth_out/dwell_times.dat`, **306 rows** |
| ejection | completed; chain leaves the tunnel |

So the dt-halving stability guard in `topo.translation.elongate.run_length`
(re-running any diverging stage at a halved timestep + proportionally more steps, which
preserves the physical dwell time) eliminates the Tutorial-10 blow-up across the
**full** chain, not just the short L=10 demo.

---

## How this differs from Tutorial 12

| | Tutorial 12 (`12_auto`) | Tutorial 13 (here) |
|--|--|--|
| length | L = 1 → **10** (the reference length) | L = 1 → **306** (full chain) |
| purpose | quantitative reproduction vs the reference run | **stability** of the fix at full scale |
| `scale_factor` | `4331293` (production) | `216564650` (= 50× — debug speed, so 306 residues finish in feasible time) |
| step caps | `max/min = 200000 / 100` | `max/min = 10000 / 400` (tighter, for tractability) |
| reference data | bundles `continuous_synthesis/` for the D6 match | none — this is the energy/robustness test, not the dwell match |

Because the kinetics here run at the 50× debug `scale_factor` with tight step caps,
the **physical dwell times** are unchanged (still written to `dwell_times.dat`) but the
**integrated MD per stage** is short. Use this run to confirm the chain grows the full
length with finite energies and a clean ejection; use **Tutorial 12** for the
quantitative dwell-time / R_g comparison against the reference.

---

## Files

| File | Role |
|------|------|
| `4c5c_model_clean.pdb`, `4c5c_model_clean_stride.dat` | 4c5c structure (306 aa) + precomputed STRIDE. |
| `domain.yaml` | 3-domain map + Go-scale strengths (as in Tutorial 12). |
| `4c5c_mrna.txt`, `trans_times.txt` | mRNA (one codon/residue) + Fluitt *E. coli* codon-time table. |
| `ribosome_trunc.pdb` | Truncated CG 50S + tRNAs (X-aligned); P-/A-anchors + rigid scenery. |
| `csp.ini` | **Debug** profile (L = 1 → 20) → `synth_out_debug/`. |
| `csp_val.ini` | **Full-length** profile (L = 1 → 306) → `synth_out/`. |
| `make_movie.py` | Stitch the per-stage trajectories into one VMD movie. |
| `analyze_validation.py`, `eject_demo.py` | Analysis helpers (copied from Tutorial 12). |
| `synth_out/`, `synth_out_debug/` | Outputs (per-residue `L_<L>/stage_<s>/`, `ejection/`, `dwell_times.dat`, `movie.*`). |

---

## How to run

All commands run **from this folder** (paths in the `.ini` files are relative). A GPU
is recommended (the v2 system has ~4,600 beads).

```bash
cd tutorials/13_validate_claude_fix12

# A. Quick debug run (L = 1 → 20, a few minutes) — smoke-test the pipeline.
topo-csp -f csp.ini            # -> synth_out_debug/

# B. Full-length validation (L = 1 → 306) — the actual stress test.
#    Watch for "[stability] ..." lines: those are stages auto-stabilised at a
#    halved timestep (the fix in action).
topo-csp -f csp_val.ini        # -> synth_out/
```

### Verify it worked (no blow-ups)

The headline check is that **no stage diverges**. Scan every stage's potential-energy
log (column 3) and flag anything above 10⁹ kJ/mol:

```bash
worst=0
for f in $(find synth_out -name traj.log); do
  m=$(awk 'NR>1 && NF>=3{pe=$3<0?-$3:$3; if(pe>m)m=pe} END{print m+0}' "$f")
  awk -v m="$m" -v f="$f" 'BEGIN{if(m>1e9) printf "BLOWUP %s  %.3g\n", f, m}'
  awk -v m="$m" -v w="$worst" 'BEGIN{exit !(m>w)}' && worst=$m
done
echo "worst max|PotE| = $worst kJ/mol"      # expect ~1.8e3, and zero BLOWUP lines
```

`analyze_validation.py` (copied from Tutorial 12) also prints the D5 energy scan and
the D5b ejection check; its **D6** section needs the reference run, so for the
quantitative dwell/R_g comparison run it in `../12_auto/` instead.

---

## How to generate the movie

The runner writes a standalone trajectory per sub-stage of every residue (each length
has a different bead count, so they can't be concatenated directly). `make_movie.py`
stitches them — in synthesis order, padding every frame up to the final length and
overlaying the static ribosome — into one VMD-playable movie under `synth_out/`:

```bash
python make_movie.py
# writes:  synth_out/movie.psf   synth_out/movie.dcd   synth_out/movie.tcl
#          synth_out/movie_ribosome.pdb   (static scenery)
```

`make_movie.py` is a path-pinned wrapper around the `topo-csp-movie` console tool;
equivalently:

```bash
topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
```

(point `-o` at `synth_out_debug` to make a movie of the quick debug run instead).

---

## How to visualize

Open the generated `.tcl` in VMD — it loads `movie.psf` + `movie.dcd`, adds the static
ribosome, and sets a sensible representation:

```bash
vmd -e synth_out/movie.tcl
```

Then press **Play** (or run `animate forward` in the VMD console). You will watch the
nascent chain **grow one residue at a time** out of the ribosome exit tunnel — each new
residue appearing at the A-site (stage 1), settling (stage 2) and translocating to the
P-site (stage 3) — for all 306 residues, then the finished chain leaving the tunnel in
the ejection phase. The chain is colored by residue id; the ribosome is shown as static
scenery. Beads not yet synthesized are parked off to the side and only appear as the
chain reaches their length.

> Headless / no display? Render frames from the command line with
> `vmd -dispdev text -e synth_out/movie.tcl` plus VMD's `render`/`movie` commands, or
> copy `synth_out/movie.{psf,dcd}` to a workstation and load them there.
