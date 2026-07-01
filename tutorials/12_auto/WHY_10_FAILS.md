# WHY_10_FAILS — post-mortem on `tutorials/10_csp_obrien` (D8)

> Written after D0–D7 passed (see [`NOTES.md`](NOTES.md)). This explains, with
> evidence, **why the original `tutorials/10_csp_obrien` does not reproduce O'Brien**,
> what `12_auto` did differently to succeed, and the minimal fix for Tutorial 10.

## TL;DR

Tutorial 10 fails on two independent levels:

1. **Numerical (the real bug).** Its CSP trajectory contains genuine **energy
   blow-ups** — PotE → ~10¹³ kJ/mol on a handful of stages — so the run is **not
   physically sane** (it would fail Goal D5). Root cause: an **integrator
   instability**, not the seed-placement clash its `OBSERVATIONS.md` #1 guessed.
2. **Scope.** Even setting the bug aside, Tutorial 10 is a *machinery demo*, never a
   *reproduction*: a different protein, a synthetic mRNA, a step clamp, and **no
   validation against the reference**.

`12_auto` fixes (1) with a per-stage stability guard and addresses (2) by running the
actual reference inputs and validating quantitatively.

---

## 1. The numerical failure, and its true root cause

`10_csp_obrien/OBSERVATIONS.md` #1 documents the symptom precisely: in a full-length
(306-residue) GPU run, **5 of 306 stages** explode (`L_014`, `L_094`, `L_095`,
`L_101`, almost all **stage 1/2**), e.g.

```
L_093/stage_3 last frame:  PotE = -6.7 kJ/mol     (fine)
L_094/stage_1 all frames:  PotE = 3.4e13 kJ/mol   (EXPLODED)
L_094/stage_2 last frame:  PotE = 617 kJ/mol      (recovered)
```

It self-recovers (the next stage's `minimizeEnergy` pulls the bad bead back), so the
final per-residue structures look fine — but the **trajectory frames of those stages
are garbage**, and that violates D5 ("per-stage PotE stays finite, no PotE ≳ 10¹²").

**Tutorial 10's hypothesis was wrong.** OBSERVATIONS #1 attributed the blow-up to the
new bead being *seeded on top of an existing bead* → a near-singular `(σ/r)¹²`
excluded-volume term the minimizer can't escape. We **measured** this on the 4c5c
system (which reproduces the identical bug — see below) and it is **not** the cause:

- Building the unstable stage in isolation and minimizing, the seed relaxes **cleanly**
  to PotE ≈ 243 kJ/mol; the new bond relaxes to 0.374 nm; the minimum nascent–ribosome
  distance is 3.07 Å — **no near-singular overlap**. (`analyze`/diagnostic in this
  folder; reproduced in the conversation log.)
- The divergence is therefore in the **dynamics**, and it is **deterministic in the
  timestep, not random in the velocity draw**:

  | config | dt = 0.015 ps (3 velocity seeds) | dt = 0.0075 ps (3 seeds) |
  |---|---|---|
  | L = 9 (4c5c) | 516 / 516 / 515 — stable | — |
  | **L = 10 (4c5c)** | **7.8e12 / 7.4e9 / 6.5e12 — diverges every seed** | **257 / 264 / 263 — stable every seed** |

**Mechanism.** topo's CSP build uses **flexible (harmonic) bonds** (it must — the new
residue is seeded ~1 nm from its bond partner at the A-site, which a rigid distance
constraint can't represent). When a newly added residue forms a **stiff native (Go)
contact**, that contact's vibrational period drops below what a **15 fs** Langevin step
can integrate, and the stage diverges — regardless of velocities or extra
minimization. It happens on only a few residues because only a few introduce such a
stiff well; that is exactly the "5/306, mostly stage 1/2" pattern.

**Why O'Brien's reference does not have this bug.** `continuous_synthesis_v6.py` builds
its system with `constraints = AllBonds` (rigid bonds via SHAKE) and runs three
minimization rounds (`rnc_l*_min_1/2/3.pdb`). Rigid bonds remove the fast bond mode,
so 15 fs is stable. Tutorial 10 inherited topo's flexible-bond elongation path without
that safeguard.

---

## 2. The scope failure

Independent of the bug, Tutorial 10 cannot be a *reproduction* of O'Brien:

| aspect | Tutorial 10 | what a reproduction needs (12_auto) |
|---|---|---|
| protein | **P0CX28** (106 aa, unrelated) | **4c5c** (306 aa) — the reference's protein |
| mRNA | synthetic back-translation of the PDB sequence | **byte-identical** to `setup/4c5c_mrna_sequence_fast.txt` |
| kinetics timescale | demo `scale_factor` + `max_steps_per_stage = 667` clamp | production `scale_factor = 4331293`, cap only as a tractability bound |
| validation | none | quantitative vs `continuous_synthesis/output/` (length, dwell, Rg) |
| CG model | demo domain/nscales | mapped from `protein_cg_model/domain_def.dat` |

Tutorial 10's own README is explicit that it is a clamped, illustrative demo ("do
**not** read folding pathways off the demo"; "the synthetic mRNA is representative, not
biological"). It exercises the *kinetics machinery* correctly; it does not *reproduce*
a specific O'Brien run, and its `csp.ini` was never compared to the reference outputs.

---

## 3. What `12_auto` did to succeed

1. **Fixed the blow-up at the source** (`topo/csp/elongate.py`, `run_length`):
   a per-stage **stability guard** that detects a diverging stage (max |PotE| > 10⁹)
   and re-runs it with a **halved timestep and double the steps** — keeping the
   physical dwell time `n_steps · dt` identical. In the validation run exactly the two
   unstable L=10 stages auto-stabilised; **no stage exceeds 524 kJ/mol** (D5 PASS).
   Crucially the guard judges divergence from the **maximum** energy over the stage,
   not the final frame (a stage can cool back under threshold yet still have ruined
   frames — which is precisely how Tutorial 10's blow-ups "self-recover").
2. **Used the reference inputs** (4c5c, identical mRNA, Fluitt table, production
   `scale_factor`) and **validated quantitatively**: length 10 = 10, total in-vivo
   dwell 1.01× the reference, final Rg 1.06× (D6 PASS).
3. **Added the missing artifacts**: `synth_out/dwell_times.dat` (the topo analog of
   the reference `output/1.out`) and an ejection analysis showing clean +x egress
   with no tunnel-wall penetration or ribosome collapse (D5b PASS).

---

## 4. The minimal fix for Tutorial 10

**Tutorial 10 shares the exact code path that was patched** — `topo.csp` calls
`topo.csp.core.run_length`, which now carries the stability guard. So:

- **Minimal fix (already in place):** simply **re-run Tutorial 10 against the patched
  `topo`**. The guard will catch its 5/306 diverging stages and re-integrate them at a
  finer timestep, eliminating the ~10¹³ kJ/mol frames with no change to the kinetics.
  (Evidence it works: the same bug, reproduced and root-caused on 4c5c at L=10 here,
  is removed by the identical guard.)
- **Cleaner long-term fix:** integrate with rigid `AllBonds` constraints (as the
  reference does), removing only the newly added bond's constraint during A-site
  delivery — this removes the stiff fast mode entirely rather than working around it.
- **To make Tutorial 10 a *reproduction* (not just stable):** point it at the 4c5c
  reference inputs and add the quantitative comparison — i.e. exactly what `12_auto`
  is. Tutorial 10 is best kept as the fast *machinery demo* it advertises, with
  `12_auto` as the validated reproduction.
