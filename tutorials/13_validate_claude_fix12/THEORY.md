# THEORY — how the continuous-synthesis simulation works, and what it models

> 🇻🇳 Bản tiếng Việt: [`THEORY.vi.md`](THEORY.vi.md).

This document explains, in detail, the physics and biology behind the run in this
tutorial: how `topo-csp` grows a nascent protein chain on a ribosome, and **how the
three simulation sub-stages map onto the real biological elongation cycle**. It is the
conceptual companion to the practical [`README.md`](README.md), and it **consolidates
the per-stage mechanics** that previously lived in a separate `STAGES.md` (now folded
into §3 below), so this document is self-contained.

Tutorial 13 runs the full **306-residue 4c5c** chain end to end and is the validation
that the protocol is numerically stable over an entire protein (see §7).

---

## 1. What is being modeled: co-translational synthesis

In a living cell a protein is not folded from a pre-existing full-length chain. It is
**synthesized vectorially, N-terminus first**, by the ribosome, one amino acid at a
time, while the growing ("nascent") chain threads out through the ribosomal **exit
tunnel** (~80 Å long, ~10–20 Å wide) and begins to fold **co-translationally** — i.e.
*as it emerges*. The kinetics of synthesis therefore matter: how long the ribosome
dwells on each codon sets how much time each segment of the chain has to sample
conformations before the next residue is added. Rare codons (decoded slowly) can act
as "translational pauses" that change folding outcomes.

The simulation reproduces this process: it grows a coarse-grained protein bead-by-bead
out of a coarse-grained ribosome, **timing each residue from its mRNA codon**, so the
nascent chain folds under realistic, codon-resolved kinetics.

### The real elongation cycle (one amino acid added)

Bacterial translation elongation repeats a three-step biochemical cycle per codon:

1. **Aminoacyl-tRNA selection / decoding.** A ternary complex (EF-Tu·GTP·aminoacyl-tRNA)
   delivers an aa-tRNA to the ribosomal **A site**. Correct codon–anticodon pairing
   triggers GTP hydrolysis and accommodation of the aa-tRNA. This is the **codon-dependent,
   highly variable, usually rate-limiting** step — its duration depends on how abundant
   the cognate tRNA is (codon usage bias).
2. **Peptidyl transfer (peptide-bond formation).** The peptidyl-transferase center (PTC)
   transfers the nascent peptide from the **P-site** tRNA onto the A-site aminoacyl-tRNA.
   The chain is now one residue longer and attached to the A-site tRNA. Fast (~0.3 ms).
3. **Translocation.** EF-G·GTP ratchets the ribosome forward by one codon: the tRNAs
   move **A→P** and **P→E**, the deacylated tRNA leaves via the E site, and the A site
   is freed for the next aa-tRNA. ~few ms.

O'Brien's *Continuous Synthesis Protocol* (CSP; `continuous_synthesis_v6.py`) partitions
the per-codon dwell time into exactly these three pieces, and `topo.csp` reproduces them
as **three MD sub-stages per residue**.

---

## 2. The simulation model (one elongation step)

The simulated System is built from reused, tested machinery
(`topo.csp.core` + `topo.csp.ribosome`):

- **Nascent protein — a structure-based (Gō-like) coarse-grained chain.** One bead per
  residue at the Cα position. Native contacts (residue pairs close in the folded 4c5c
  crystal structure) get attractive wells; everything else is repulsive. **The native
  structure is therefore the energy minimum**, so the growing chain folds *toward* 4c5c's
  native fold — co-translationally, exactly as the real chain folds as it emerges. Bonds
  are **flexible harmonic** (not rigid constraints — see §7).
- **Ribosome — rigid scenery.** The truncated CG 50S + tRNAs (`ribosome_trunc.pdb`,
  ~4,600 mass-0 beads) is fixed in space. It does not move, but its **excluded-volume and
  electrostatic interactions with the nascent chain are on**, so the chain feels the tunnel
  walls and the ribosome surface. The tunnel axis is aligned with **+x** (the chain exits
  toward +x).
- **The PTC anchors.** Two ribosome beads are singled out as fixed reference points: the
  **P-anchor** (P-site tRNA residue-76 "R" bead) and the **A-anchor** (A-site tRNA
  residue-76 "R" bead). These stand in for where the peptidyl-tRNA (P) and incoming
  aminoacyl-tRNA (A) hold the chain's C-terminus.
- **C-terminus tether — a harmonic restraint.** The current C-terminal bead is restrained
  to one of the anchors with `U = k·|r − r₀|²`, `k = restraint_k = 83680 kJ/mol/nm²`
  (= 200 kcal/mol/Å²). This reproduces the covalent attachment of the nascent chain's
  C-terminus to the tRNA sitting in the A or P site. **Switching the restraint target
  from the A-anchor to the P-anchor is how translocation is reproduced.**
- **Tunnel wall — a one-sided plane.** A half-harmonic wall keeps every nascent bead at
  `x ≥ x₀` (`k = 8368 kJ/mol/nm²`). Because the 50S is **truncated** for speed, there are
  no ribosome beads below the PTC — the lower tunnel is an empty void; the wall supplies
  the missing floor so the chain can only extrude **forward** (out the tunnel, +x) and
  cannot slip back through the synthesis point into the truncated region. The plane `x₀` is
  **auto-derived from the ribosome structure** (not a config knob): the lower / P-site
  C-terminus hold plane, `x₀ = min(P-anchor.x, A-anchor.x) + ptc_offset` — for this
  `ribosome_trunc.pdb`, `0.5705 + 0.476 = 1.0465 nm` (≈ the old hardcoded 1.05).
- **Thermostat.** Langevin dynamics at `ref_t = 310 K`, friction `tau_t = 0.01 /ps`,
  timestep `dt = 0.015 ps`.

### Build-once-subset contacts (why this is efficient)

The contact map is computed **once** on the full 4c5c structure (`R_full`, `eps_full`,
306×306). For nascent length `L`, the model uses the top-left `L×L` block — STRIDE and
the heavy-atom contact analysis are never re-run per length. So residues `1..L` carry
exactly the native contacts they will have in the final fold, and the chain can form
native structure as soon as the relevant residues exist.

---

## 3. The three stages: biology ↔ simulation

Each amino acid is added through **one ribosome elongation cycle**, which O'Brien split
into three kinetic sub-steps; `topo-csp` runs **one MD segment per sub-step**. For
nascent length `L`, each sub-stage is a standalone short simulation (its own
`synth_out/L_<L>/stage_<s>/` folder); stage 3's final structure seeds the next residue's
stage 1. The codon's total in-vivo translation time is **partitioned across the three
stages**, and the loop is one line:

```
codon_total  =  stage 1   +   stage 2            +   stage 3
                (fixed)       (fixed + traffic)      (remainder, codon-dependent)

deliver at A (stage 1) → hold at A while translocating (stage 2)
   → move to P and wait for the next tRNA (stage 3) → repeat for residue L+1
```

| stage | real biological process | what the simulation does | C-terminus restrained to | mean dwell time |
|-------|-------------------------|--------------------------|--------------------------|-----------------|
| **1** | **Peptidyl transfer** — peptide bond forms; chain now sits on the A-site tRNA | The new bead `L` is **placed at the A-anchor** (+`buffer = 0.4 nm` into the tunnel), bonded to residue `L−1`; minimize; run MD | **A-anchor** | `time_stage_1 = 0.34 ms` |
| **2** | **Translocation (onset)** — EF-G begins ratcheting the ribosome forward | Continue from stage 1, **still held at the A-anchor**; run MD | **A-anchor** | `time_stage_2 = 4.20 ms` |
| **3** | **Translocation completes + wait for next aa-tRNA** — peptidyl-tRNA now in the P site; ribosome waits to decode the next codon | **Switch the restraint A→P** (this geometric move *is* the translocation), then run MD while the chain relaxes/folds | **P-anchor** | remainder = (next codon's total time) − stage 1 − stage 2 |

### Stage 1 — peptidyl transfer / A-site delivery
- **Biology.** The incoming amino acid's tRNA is accommodated in the A site and the PTC
  forms the new peptide bond. The nascent chain is now length `L`, covalently held by the
  A-site tRNA.
- **Simulation.** Bead `L` is seeded at the A-anchor (offset `buffer` along +x so it does
  not overlap the anchor bead). A flexible bond to residue `L−1` is created — it starts
  stretched (the A and P anchors are ~0.92 nm apart), so the structure is **energy-minimized**
  first to relax the placement, then MD runs for the stage-1 step count. The C-terminus is
  restrained to the A-anchor.
- **Mean time.** `time_stage_1 = 0.000340 s` (the experimental peptide-bond-formation dwell).

### Stage 2 — translocation begins
- **Biology.** EF-G·GTP engages and the ribosome starts moving one codon along the mRNA;
  the tRNAs begin their A→P / P→E shift.
- **Simulation.** The system continues from stage 1's final structure, **still restrained
  at the A-anchor**, and runs MD for the stage-2 step count. (The actual A→P geometric move
  is applied at the start of stage 3.)
- **Mean time.** `time_stage_2 = 0.004201 s` (experimental translocation dwell), plus an
  optional **ribosome-traffic** correction (off in this tutorial; see §5).

### Stage 3 — translocation to the P site + tRNA-binding wait
- **Biology.** Translocation completes (the peptidyl-tRNA is now in the **P site**, the A
  site is empty) and the ribosome **waits for the next cognate aminoacyl-tRNA** to arrive
  and be decoded. This wait is **codon-dependent** and usually the longest part of the
  cycle — it is what makes rare codons slow.
- **Simulation.** The C-terminus restraint **switches from the A-anchor to the P-anchor**
  — reproducing the one-codon A→P translocation as a discrete geometric move — and MD runs
  for the stage-3 step count. With the C-terminus pulled back to the P site and the tunnel
  wall preventing back-folding, the chain extrudes one register further toward the exit.
  Stage 3's final structure becomes the **seed for residue `L+1`'s stage 1**.
- **Mean time.** the **remainder**: (mean total in-vivo time of the *next* codon) − `time_stage_1`
  − `time_stage_2`. If a fast codon makes this non-positive it is floored to a tiny positive
  value. This is where the per-codon variability enters the schedule.

> **Mechanics vs. timing — an honest note.** The restraint switch (an instantaneous,
> geometric A→P move) happens at the **start of stage 3**, while the **duration** charged to
> translocation is **stage 2** and the duration charged to the decoding wait is **stage 3**.
> So the physical translocation move and the time labelled "translocation" are slightly
> decoupled — a deliberate simplification. Likewise the peptide bond is present in the bonded
> model from stage 1 rather than being toggled on mid-stage, and explicit A/P tRNA bonded
> geometry is not modelled. The **timing** (three codon-resolved dwell times per residue) is
> faithful to O'Brien; the per-stage **mechanics** are a reduced model (explicit A/P tRNA
> bonded geometry is not modelled).

---

## 4. From codon to MD steps (the kinetics)

The timing core is `topo.csp.kinetics` (pure, no OpenMM). For every residue it answers:
*how many integration steps does each sub-stage run?*

**(a) Per-codon mean translation time.** The mRNA (`4c5c_mrna.txt`) is split into codons;
a lookup table (`trans_times.txt`, the Fluitt *E. coli* table at 310 K) maps each codon to
its **mean in-vivo translation time** in seconds. This is the codon's intrinsic **mean
first-passage time (mFPT)** — dominated by cognate-tRNA availability (codon usage). Call it
`τ(codon)`.

**(b) First-passage-time sampling.** Real elongation is stochastic: each step is governed by
a rate-limiting molecular event, so its waiting time is **exponentially distributed**. Each
sub-stage's actual dwell time is drawn as `t = −mean · ln(U)`, `U ∼ Uniform(0,1)`
(`random.expovariate`). A fixed `random_seed` makes the schedule reproducible.

**(c) The three-stage split** for nascent length `L` (1-indexed):
```
t1  ~ Exp(  time_stage_1 )                              # peptidyl transfer (fixed mean)
t2  ~ Exp(  time_stage_2 + max(0, traffic_correction) ) # translocation (+ optional traffic)
t3  ~ Exp(  τ(next codon) − time_stage_1 − time_stage_2 )# wait to decode the next codon
```
so the **per-cycle clock** = fixed peptide bond + fixed translocation + variable decoding
wait. (Indexing note: stage 3 uses the *next* codon's mean — the ribosome, having just
incorporated residue `L`, now waits for residue `L+1`'s tRNA.)

**(d) In-vivo seconds → in-silico steps.** The coarse-grained model evolves on a far shorter
timescale than real translation, so a **time-compression factor** maps seconds to MD steps:
```
t_sim (ns)  =  t_s · 1e9 / scale_factor
n_steps     =  t_sim (ns) / dt(ns) ,   dt(ns) = dt_ps · 1e-3
```
A **larger `scale_factor` ⇒ fewer steps per residue ⇒ a faster run**, while preserving the
*relative* timing of fast vs. slow codons. This tutorial uses
`scale_factor = 216564650` (= **50×** the production value 4331293) so the whole 306-residue
chain finishes quickly; the per-codon *ratios* — the physics that matters for
co-translational folding — are unchanged. Step counts are additionally clamped to
`[min_steps_per_stage, max_steps_per_stage] = [400, 10000]` for tractability (a clamp on MD
steps only; the sampled dwell **times in seconds** are recorded untouched in
`synth_out/dwell_times.dat`).

---

## 5. Ribosome traffic, and the `intrinsic` vs `real` per-codon times

The kinetics carry **two** per-codon time lists:

- **`intrinsic[i]`** — the codon's bare mean translation time, looked up straight from the
  `trans_times` table (no queueing). This is the intrinsic mFPT of §4.
- **`real[i]`** — the same time **with a ribosome-traffic correction**. In a real cell,
  ribosomes queue on one mRNA: a slow ribosome ahead delays those behind it, lengthening the
  effective dwell. O'Brien model this with an external `ribosome_traffic` binary that — given
  the mRNA, the intrinsic times and the `initiation_rate` — returns a traffic-corrected time
  per codon.

Only **stage 2** uses the correction:
```
correction      = real[L−1] − intrinsic[L−1]
mean(stage 2)   = time_stage_2 + correction   if correction > 0,   else time_stage_2
```

**When `ribosome_traffic = no`, `real == intrinsic`.** `build_mfpt_lists` initialises
`real = list(intrinsic)` and overwrites it **only** when traffic is both requested *and* the
external binary actually runs; otherwise `real` is left equal to `intrinsic`. So `real ==
intrinsic` whenever:

- `ribosome_traffic = no` (the `if ribosome_traffic:` branch is skipped entirely), **or**
- `ribosome_traffic = yes` but the `ribosome_traffic` binary is missing / fails (it returns
  `None` and the code falls back to `real = intrinsic`, with a printed warning).

In all of these `correction = 0`, so **stage 2's mean is simply `time_stage_2`** (no stretch).
This tutorial sets `ribosome_traffic = no`, so `real == intrinsic` everywhere and the schedule
is driven purely by the bare per-codon times; the traffic effect at these lengths is small.

Code: `topo.csp.kinetics.build_mfpt_lists` / `ribosome_traffic_times` / `stage_dwell_times`.

---

## 6. After the last residue: ejection (and dissociation)

Once residue 306 is added, the protein is complete. The simulation then runs a
**post-synthesis ejection phase** (`ejection_steps = 50000`): the C-terminus restraint is
**released** (the tether is cut), while the rigid ribosome and the one-sided tunnel wall
remain. Biologically this is **termination** — release factors hydrolyze the peptidyl-tRNA
bond and the finished protein is freed. With the tether gone, the chain **diffuses out of
the tunnel along +x** (the one-sided wall acts as a reflecting barrier that biases motion
forward) and clears the ribosome, free to complete folding. An optional **dissociation**
phase (`dissociation_steps = 0` here) would continue the free protein away from the ribosome.

For a longer, dedicated egress demonstration see `eject_demo.py` (releases the final chain
and runs it far enough to fully clear the tunnel).

---

## 7. Numerical integration and the stability guard (why Tutorial 13 exists)

The chain is integrated with **flexible harmonic bonds** at `dt = 0.015 ps`. Flexible (rather
than rigid `AllBonds`) bonds are required because stage 1 seeds the new bead ~1 nm from its
bond partner (A-site delivery), which a rigid distance constraint cannot represent — a
harmonic bond absorbs the stretch and the minimizer relaxes it.

The cost: at 15 fs, the integration is only **marginally stable** for some configurations.
When a newly added residue forms a **stiff native (Gō) contact**, that contact's vibrational
period drops below what a 15 fs step can integrate, and the dynamics **diverge** (potential
energy → ~10¹³ kJ/mol), corrupting that stage's frames. This was the latent bug in Tutorial
10 (`OBSERVATIONS.md` #1): ~5 of 306 stages blew up in a full run. It is **deterministic in
the timestep, not random** — measured on 4c5c, length 10 diverges at 15 fs on every velocity
seed but is stable at 7.5 fs on every seed. (O'Brien's reference avoids the problem entirely
by using rigid `AllBonds` constraints, which remove the fast bond mode.)

**The fix** (in `topo.csp.core.run_length`, the per-stage stability guard): each
stage is run in chunks while tracking the **maximum** |PotE|; if a stage diverges
(max |PotE| > 10⁹ kJ/mol) it is transparently **re-run with the timestep halved and the step
count doubled**. Because the physical dwell time is `n_steps · dt`, halving `dt` and doubling
`n_steps` **leaves the dwell time exactly unchanged** while making the integration stable
(up to 6 halvings). The common case runs once at 15 fs.

**Tutorial 13 is the validation of that fix at full scale.** Synthesizing all 306 residues
(919 stages), **zero** stages blow up and the worst max |PotE| over the entire run is
~1.8×10³ kJ/mol — confirming the guard holds across the whole chain, in the exact regime
where the earlier unguarded full-length run failed.

---

## 8. The 4c5c system and this tutorial's parameters

- **Protein:** 4c5c, **306 residues**, three domains (from `domain.yaml`, mapped from the
  O'Brien CG model): **A** = 1–84, **B** = 85–110 + 184–306 (discontiguous), **C** = 111–183,
  with per-domain and per-interface Gō-contact strength scaling (the structure-based analog of
  O'Brien's `nscal`).
- **mRNA:** `4c5c_mrna.txt` — one codon per residue (+ stop), **byte-identical** to the O'Brien
  reference fast mRNA, so the codon schedule is exactly the reference's.
- **Key numbers** (`csp_val.ini`): `dt = 0.015 ps`, `ref_t = 310 K`, `scale_factor = 216564650`
  (50× production), `time_stage_1 = 0.34 ms`, `time_stage_2 = 4.20 ms`,
  steps clamped to `[400, 10000]`, `restraint_k = 83680`, tunnel wall at `x ≥ 1.05 nm`
  (`k = 8368`), `ejection_steps = 50000`.

> **Config caveat.** The header comments at the top of `csp_val.ini` are stale (copied from
> Tutorial 12 — they mention "length 1 → 10" and `scale_factor = 4331293`). The *active*
> settings are the ones above (`L_max = 306`, `scale_factor = 216564650`); trust the keys, not
> the comment block.

---

## 9. Summary: the correspondence at a glance

```
REAL ELONGATION CYCLE (per codon)          SIMULATION (per residue L, 3 MD sub-stages)
─────────────────────────────────          ───────────────────────────────────────────
peptide bond forms (PTC)            →  stage 1: place bead L at A-anchor, hold at A, MD
                                              dwell mean = time_stage_1 (0.34 ms)
EF-G translocation begins           →  stage 2: continue, still held at A, MD
                                              dwell mean = time_stage_2 (4.20 ms)
translocation completes (→ P site)  →  stage 3: switch restraint A→P (the move), MD
+ wait to decode next aa-tRNA                 dwell mean = next-codon time − stage1 − stage2
(codon-dependent, rate-limiting)              (this is the variable, codon-resolved part)

termination / release               →  ejection: release tether, chain diffuses out +x
```

Codon-resolved dwell times (exponentially sampled, compressed by `scale_factor`) set how
long each segment relaxes before the next residue arrives — so the chain folds **under the
same relative kinetics as real translation**, while the ribosome scenery + tunnel wall keep
it threading forward out of the exit tunnel.

---

## Code & document pointers

- Kinetics: [`topo/csp/kinetics.py`](../../topo/csp/kinetics.py) — `codon_mfpt_list`,
  `sample_fpt`, `stage_dwell_times`, `seconds_to_steps`, `stage_steps`.
- 3-stage loop: [`topo/csp/protocol.py`](../../topo/csp/protocol.py) — `run_continuous_synthesis`
  (three `run_length` calls per residue; A→A→P restraint switch; ejection phase).
- Per-stage MD + the stability guard: [`topo/csp/elongate.py`](../../topo/csp/elongate.py)
  — `run_length`.
- The three-stage mechanics are in §3 of this document (consolidated from the former
  `STAGES.md`); the numerical fix and its rationale are in §7; validation numbers are in
  [`README.md`](README.md).
