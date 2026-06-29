# The three elongation stages (CSP), explained

Each amino acid is added through **one ribosome elongation cycle**, which O'Brien
split into three kinetic sub-steps. `topo-csp` runs **one MD segment per sub-step**.
For residue `L`, the codon's total in-vivo translation time (from `trans_times`) is
partitioned across the three stages:

```
codon_total  =  stage 1  +  stage 2  +  stage 3
                (fixed)     (fixed+traffic)  (remainder, codon-dependent)
```

## Stage 1 — peptidyl transfer (A-site delivery)
- **Biology:** the aminoacyl-tRNA carrying the new amino acid sits in the ribosomal
  **A site**, and the peptide bond forms (the existing chain is transferred onto the
  new residue).
- **Simulation:** the new bead `L` is placed at the **A-anchor** (A-site tRNA
  position); its C-terminus is restrained there while a short MD runs.
- **Mean dwell:** `time_stage_1` (default 0.00034 s).

## Stage 2 — translocation
- **Biology:** the ribosome begins moving one codon along the mRNA; the tRNAs start
  shifting A→P.
- **Simulation:** continues from stage 1, **still held at the A-site**; MD runs.
- **Mean dwell:** `time_stage_2` (default 0.004201 s) **+ ribosome-traffic
  correction** if `ribosome_traffic = yes` (an upstream-queue delay).

## Stage 3 — tRNA binding / waiting
- **Biology:** translocation completes (the new residue's tRNA is now in the **P
  site**), the A site is empty, and the ribosome **waits for the next aminoacyl-tRNA**
  to arrive. Usually the longest part and **codon-dependent** — rare codons wait
  longer.
- **Simulation:** the C-terminus restraint **switches from the A-anchor to the
  P-anchor** — this A→P move is what reproduces translocation — and MD runs. Stage
  3's final structure **seeds the next residue's stage 1**.
- **Mean dwell:** the remainder = `codon_total − time_stage_1 − time_stage_2`.

## The loop

```
deliver at A (stage 1) → hold at A while translocating (stage 2)
   → move to P and wait for next tRNA (stage 3) → repeat for residue L+1
```

## Connected points

- **Sampling → steps.** The three dwell times are drawn from **exponential
  distributions** about their means (first-passage-time sampling), then mapped to MD
  steps via `scale_factor`: `steps = sampled_s · 1e9 / scale_factor / dt_ns`.
- **The test clamp.** `max_steps_per_stage` caps all three stages (demo: 667 ≈ ~2000
  steps/residue). Remove it for production so a rare codon's stage 3 genuinely runs
  longer than a fast codon's. Note: the clamp applies **even when `L_max` is blank**
  (full length) — delete the clamp lines for true codon-driven step counts.
- **Why the chain doesn't fall back into the ribosome.** The C-terminus is always
  tethered at the peptidyl-transferase center, and `tunnel_wall` keeps beads moving
  forward (+x) out the exit tunnel.
- **After the last residue.** A separate **ejection** phase releases the restraint so
  the finished chain leaves the tunnel (`ejection_steps`); an optional
  **dissociation** phase lets it drift off (`dissociation_steps`).

## Fidelity caveat

The *timing* (three codon-resolved dwell times) is faithful to O'Brien. The per-stage
*mechanics* are simplified: A→P translocation is mapped via the **restraint-target
switch**, and the peptide bond is present in the bonded model from stage 1 rather than
toggled on mid-residue. Explicit A/P tRNA bonded geometry is not modeled. See
[`PLAN.md`](PLAN.md) → *Remaining*.

## Code pointers

- Timing core: [`topo/csp/kinetics.py`](../../topo/csp/kinetics.py) —
  `stage_dwell_times`, `stage_steps`, `seconds_to_steps`.
- 3-stage loop: [`topo/csp/csp.py`](../../topo/csp/csp.py) →
  `run_continuous_synthesis` (three `run_length` calls per residue, A→A→P restraint).
