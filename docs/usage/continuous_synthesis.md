# Synthesis in coarse-grained ribosome model

The **Continuous Synthesis Protocol (CSP)** is the **codon-resolved, kinetic** runner
for protein synthesis on an **explicit coarse-grained ribosome**. It times
**every residue from its mRNA codon** and splits each elongation cycle into **three
kinetic sub-stages** — the published continuous-synthesis protocol, in
topo style. (Codon-resolved kinetics are what make the model physically meaningful; a
fixed per-residue step count is not.) For the analytic-tunnel variant — the same codon
kinetics with the explicit ribosome replaced by a cylindrical bore — see
{doc}`cylinder_synthesis`.

- **CLI:** `topo-csp -f csp.ini` (or `python -m topo.csp -f csp.ini`)
- **Movie tool:** `topo-csp-movie -o <out_root> [--ribosome ribo.pdb]`
- **Worked example:** Tutorial 8 (`tutorials/08_ribosome_synthesis/`) — a smoke run
  (`csp_debug.ini`, L = 1 → 8) and a full-length validation (`csp_val.ini`, L = 1 → 306)
  on 4c5c, plus a second protein (P0CX28).
- **Architecture:** CSP is a thin outer loop. The per-length MD work — building the
  length-`L` model, seeding coordinates, restraints, running one stage under the
  stability guard, build-once-subset contacts — lives in the shared low-level engine
  `topo.csp.core` (`run_length`, `RunParams`); the rigid-ribosome scenery and
  tunnel wall live in `topo.csp.ribosome`; the timing lives in `topo.csp.kinetics`.
  CSP adds only the kinetics and the three-`run_length`-calls-per-residue loop. The
  force field and ribosome representation are described in the Theory section below.

---

## Quick start

All paths in the INI are relative to the working directory; run from the tutorial
folder. A GPU is recommended (the v2 system has ~4,600 ribosome beads).

```bash
# Run the CSP on 4c5c (Tutorial 8)
cd tutorials/08_ribosome_synthesis/4c5c
topo-csp -f csp_debug.ini       # smoke run -> synth_out_debug/  (or csp_val.ini -> synth_out/ for full length)

# Stitch the per-stage trajectories into one VMD movie
topo-csp-movie -o synth_out_debug --ribosome ribosome_trunc.pdb
```

`topo-csp` writes one folder per residue `L` under `<outdir>/L_<L>/` (a shared topology
plus a per-stage trajectory `traj_s<s>.dcd`), an optional `ejection/` (and
`dissociation/`) phase, and a per-residue dwell-time log `<outdir>/dwell_times.dat`.

---

## Theory

### 1. What is being modeled: protein synthesis

In a living cell a protein is **synthesized vectorially, N-terminus first**, by the
ribosome, one amino acid at a time, while the growing ("nascent") chain threads out
through the ribosomal **exit tunnel** (~80 Å long, ~10–20 Å wide) and begins to fold
**as it emerges**. The kinetics of synthesis matter: how long the
ribosome dwells on each codon sets how much time each segment of the chain has to sample
conformations before the next residue is added. Rare codons (decoded slowly) act as
"translational pauses" that can change folding outcomes.

CSP reproduces this: it grows a coarse-grained protein bead-by-bead out of a
coarse-grained ribosome, **timing each residue from its mRNA codon**, so the nascent
chain folds under realistic, codon-resolved kinetics.

#### The real elongation cycle (one amino acid added)

Bacterial translation elongation repeats a three-step biochemical cycle per codon:

1. **Aminoacyl-tRNA selection / decoding.** A ternary complex
   (EF-Tu·GTP·aminoacyl-tRNA) delivers an aa-tRNA to the ribosomal **A site**. Correct
   codon–anticodon pairing triggers GTP hydrolysis and accommodation. This is the
   **codon-dependent, highly variable, usually rate-limiting** step — its duration
   depends on cognate-tRNA abundance (codon-usage bias).
2. **Peptidyl transfer (peptide-bond formation).** The peptidyl-transferase center (PTC)
   transfers the nascent peptide from the **P-site** tRNA onto the A-site aminoacyl-tRNA.
   The chain is now one residue longer and attached to the A-site tRNA. Fast (~0.3 ms).
3. **Translocation.** EF-G·GTP ratchets the ribosome forward by one codon: the tRNAs move
   **A→P** and **P→E**, the deacylated tRNA leaves via the E site, and the A site is freed
   for the next aa-tRNA. ~few ms.

CSP partitions the per-codon dwell time into exactly these three pieces and reproduces
them as **three MD sub-stages per residue**.

### 2. The simulation model (one elongation step)

- **Nascent protein — a structure-based (Gō-like) coarse-grained chain.** One bead per
  residue at the Cα position. Native contacts (residue pairs close in the folded crystal
  structure) get attractive wells; everything else is repulsive — so **the native
  structure is the energy minimum** and the growing chain folds *toward* its native fold.
  Bonds are **flexible harmonic** (not rigid constraints — see §7).
- **Ribosome — rigid scenery.** The truncated CG 50S + tRNAs (~4,600 mass-0 beads) is
  fixed in space, but its **excluded-volume and electrostatic interactions with the
  nascent chain are on**, so the chain feels the tunnel walls and ribosome surface. The
  tunnel axis is aligned with **+x** (the chain exits toward +x).
- **The PTC anchors.** Two ribosome beads are singled out as fixed reference points: the
  **P-anchor** (P-site tRNA residue-76 "R" bead) and the **A-anchor** (A-site tRNA
  residue-76 "R" bead). They stand in for where the peptidyl-tRNA (P) and incoming
  aminoacyl-tRNA (A) hold the chain's C-terminus.
- **C-terminus tether — a harmonic restraint.** The current C-terminal bead is restrained
  to one of the anchors with `U = k·|r − r₀|²`, `k = restraint_k = 83680 kJ/mol/nm²`
  (= 200 kcal/mol/Å²), reproducing the covalent attachment of the C-terminus to the tRNA
  in the A or P site. **Switching the restraint target A→P is how translocation is
  reproduced.**
- **Tunnel wall — a one-sided plane on the nascent beads.** The rigid ribosome we
  supply is **truncated** to a shell around the exit tunnel (§8), so the model has no
  density on the **−x (synthesis-interface) side of the PTC**. In the real ribosome that
  side is *not* empty — it is occupied by the small subunit, the mRNA and the A-/P-site
  tRNAs, so the nascent chain can never move there — but truncation leaves an **unreal
  void** that the flexible nascent chain could drift *backward* (−x) into and tangle. A
  single **one-sided half-harmonic restraint on every nascent bead**,
  `U = k·min(x − x₀, 0)²` (`k = 8368 kJ/mol/nm²` = 20 kcal/mol/Å², a fixed model
  constant), supplies the missing barrier: it pushes any bead at `x < x₀` back to
  `x ≥ x₀`, so the chain can only extrude **forward** (+x, toward the cytosolic exit) and
  cannot slip below the synthesis point into the truncated region. The plane `x₀` is
  **auto-derived from the ribosome structure**, not a config knob: it sits at the lower
  (smaller-x) C-terminus hold plane, `x₀ = min(A-target.x, P-target.x)` —
  the P-site, where the C-terminus is tethered — so the held C-terminus sits on the plane
  and the chain grows away from it. It is recomputed for whatever ribosome PDB you supply,
  so it never goes stale; for Tutorial 8's `ribosome_trunc.pdb` it is **x₀ ≈ 1.05 nm**.
- **Thermostat.** Langevin dynamics at `ref_t = 310 K`, friction `tau_t = 0.05 /ps`,
  timestep `dt = 0.015 ps`.
- **Optional flexible exit-tunnel loop.** The ribosome is rigid by default, but
  `ribo_free_mask` can free a portion of the exit-tunnel protein (E. coli L24 /
  eukaryotic L26) as a mobile structure-based Go loop — the β-hairpin lining the tunnel
  can then breathe and contact the nascent chain. See
  {ref}`the flexible-loop notes <flexible-loop>`.

**Build-once-subset contacts.** The contact map is computed **once** on the full native
structure (`R_full`, `eps_full`, `N×N`). For nascent length `L`, the model uses the
top-left `L×L` block — STRIDE and the heavy-atom contact analysis are never re-run per
length. So residues `1..L` carry exactly the native contacts they will have in the final
fold, and the chain can form native structure as soon as the relevant residues exist.

### 3. The three stages: biology ↔ simulation

Each amino acid is added through **one ribosome elongation cycle**, split into three
kinetic sub-steps; CSP runs **one MD segment per sub-step**. For nascent length `L`, each
sub-stage is a short simulation writing its own `traj_s<s>.dcd` in the shared `L_<L>/`
folder; stage 3's final structure seeds the next residue's stage 1.

| stage | real biological process | what the simulation does | C-terminus restrained to | mean dwell time |
|-------|-------------------------|--------------------------|--------------------------|-----------------|
| **1** | **Peptidyl transfer** — peptide bond forms; chain now sits on the A-site tRNA | new bead `L` **placed at the optimal A-site target** (one equilibrium peptide bond from `L−1`), bonded to `L−1`; minimize; run MD | **A-target** | `time_stage_1 = 0.34 ms` |
| **2** | **Translocation (onset)** — EF-G begins ratcheting forward | continue from stage 1, **still held at the A-anchor**; run MD | **A-anchor** | `time_stage_2 = 4.20 ms` |
| **3** | **Translocation completes + wait for next aa-tRNA** | **switch the restraint A→P** (this geometric move *is* the translocation), then run MD while the chain relaxes/folds | **P-anchor** | remainder = (next codon's total) − stage 1 − stage 2 |

```{note}
**Mechanics vs. timing.** The restraint switch (an instantaneous A→P geometric move)
happens at the **start of stage 3**, while the **duration** charged to translocation is
**stage 2** and the duration charged to the decoding wait is **stage 3** — so the physical
move and the time labelled "translocation" are slightly decoupled. Likewise the peptide
bond is present in the bonded model from stage 1 rather than toggled on mid-stage, and
explicit A/P tRNA bonded geometry is not modelled. The **timing** (three codon-resolved
dwell times per residue) is faithful to the reference protocol; the per-stage
**mechanics** are a reduced model.
```

### 4. From codon to MD steps (the kinetics)

The timing core is `topo.csp.kinetics` (pure Python, no OpenMM). For every residue it
answers: *how many integration steps does each sub-stage run?*

**(a) Per-codon mean translation time.** The mRNA is split into codons; a lookup table
maps each codon to its **mean in-vivo translation time** in seconds — the codon's
intrinsic **mean first-passage time (mFPT)**. Call it `τ(codon)`.

```{note}
**What the codon-time table is.** `τ(codon)` is the **mean total per-codon translation
time** — the full elongation-cycle mean first-passage time — and is **dominated by, but
not equal to, the codon-dependent aminoacyl-tRNA decoding (selection) step** (cognate-tRNA
availability / codon-usage bias). It is *not* only the decoding time: in CSP, `τ` is the
whole codon dwell, and stage 3 (`τ − time_stage_1 − time_stage_2`) is the remainder left
after the codon-*independent* peptidyl transfer and translocation — effectively the
decoding/tRNA wait.

**The table is an organism + temperature property** (not of the protein). A per-codon run
therefore needs an explicit `codon_times` table path -- there is no bundled default. The
shared library under `assets/csp/codon_dwell_times/` holds tables per organism, including
the **Fluitt *E. coli* table at 310 K** (Fluitt, Pienaar & Viljoen, *Comput. Biol. Chem.*
2007; 61 sense + 3 stop codons, mean ≈ 0.068 s ≈ 15 aa/s) at
`assets/csp/codon_dwell_times/ecoli/ecoli_codon_dwell_times_310K.txt`. Both `mrna` (the
protein-specific sequence) and `codon_times` (the table) are mandatory for per-codon timing.
```

(fastest-slowest-synonymous-codon-mrna)=
#### Fastest / slowest / median (synonymous-codon) mRNA

Instead of supplying an mRNA file, you can have TOPO build one that keeps the **protein
sequence fixed** but reassigns each residue's **codon** — a controlled walk along the
synonymous-mutation (codon-optimization) axis. Because the protein, its native structure,
and its contacts are identical across these mRNAs and only the elongation *timing* changes,
any difference in the co-translational folding they produce is attributable to codon
kinetics alone. Set `mrna` to a keyword instead of a filename:

- **`fastest`** — the shortest-`τ` synonymous codon at every residue → the fastest
  translation the protein's codons allow (minimises the mean total dwell time).
- **`slowest`** — the longest-`τ` synonymous codon at every residue → the slowest
  translation (maximises the mean total dwell); the slow codons act as translational pauses
  that give each emerging segment more time to fold.
- **`median`** — the middle-`τ` synonymous codon at every residue → a neutral reference
  sitting between the two extremes. When an amino acid has an even number of synonymous
  codons there is no single middle, so TOPO takes the **faster (shorter-`τ`) of the two
  central** codons; the pick is deterministic.

Here `τ` is the codon's mean per-codon dwell time (introduced above), read from the
`codon_times` **table** — whose third column also supplies the codon→amino-acid map. The
runner reads the amino-acid sequence from `pdb_file`, picks one codon per residue as above,
appends the matching stop codon, and writes the raw-nucleotide mRNA next to the PDB as
`mrna_fastest.txt` / `mrna_slowest.txt` / `mrna_median.txt`; that file then feeds the normal
per-codon path. Because the pick is made per amino acid, the mode needs a codon-time
**table** — a single uniform `codon_times` number carries no synonymous choices and is
rejected. To pre-generate the file standalone (e.g. to inspect or version it):

```bash
topo-make-mrna --pdb protein.pdb \
    --codon-times assets/csp/codon_dwell_times/ecoli/ecoli_codon_dwell_times_310K.txt \
    --mode fastest
```

**(b) First-passage-time sampling.** Real elongation is stochastic: each stage is gated by
a single rate-limiting molecular event, so its waiting time is **exponentially
distributed**. Each sub-stage's actual dwell time is therefore **drawn at random** as

```text
t = −mean · ln(U) ,   U ∼ Uniform(0, 1)
```

(the inverse-transform draw of an exponential; the code uses the equivalent
`random.expovariate(1/mean)`). A fixed `random_seed` makes the whole schedule reproducible.

```{important}
**`time_stage_1` and `time_stage_2` (and the stage-3 remainder) are *means*, not the
per-residue dwell.** The values you set in `csp.ini` are the *averages* each stage
converges to over the chain; the **actual** time spent on each residue is a fresh random
draw from an exponential with that mean — so individual residues get shorter or longer
dwells (the exponential has a long tail: many short waits, occasional long ones), and only
their average equals the configured value. In the notation below, `Exp(m)` means "an
exponential whose **mean is `m`**" (the argument is the mean, *not* a rate). For example,
with `time_stage_1 = 0.00034 s`, draws of `t1` scatter around 0.00034 s (≈0.00006 s when
`U = 0.84`, ≈0.00046 s when `U = 0.26`) and average to 0.00034 s.
```

**(c) The three-stage split** for nascent length `L` (1-indexed). Draw three independent
`U₁, U₂, U₃ ∼ Uniform(0, 1)` and compute the **actual** per-residue dwells:

```text
t1(L) = −time_stage_1                                  · ln(U₁)   # peptidyl transfer
t2(L) = −time_stage_2                                  · ln(U₂)   # translocation
t3(L) = −(τ(next codon) − time_stage_1 − time_stage_2) · ln(U₃)   # wait to decode next codon
```

Each line is one **sample** from an exponential whose **mean** is the bracketed quantity —
i.e. `t1(L) = −time_stage_1·ln(U₁)` is the same statement as "`t1 ∼ Exp(mean =
time_stage_1)`". The bracketed quantities (`time_stage_1`, `time_stage_2`, and the stage-3
remainder) are the **means**; the `t(L)` values are the per-residue draws that scatter
around them. *Example:* residue L with `U₁ = 0.84` gets
`t1(L) = −0.000340·ln(0.84) ≈ 0.000057 s` (this residue), while another residue with a
different `U₁` gets a different `t1` — and they average to `time_stage_1 = 0.00034 s`.

So the **per-cycle clock** = fixed-mean peptide bond + fixed-mean translocation + a
variable-mean decoding wait. (Indexing: stage 3 uses the *next* codon's mean — having just
incorporated residue `L`, the ribosome now waits for residue `L+1`'s tRNA. This is why
`build_codon_time_lists` needs the codon list to extend to `L_max + 1`.)

**(d) In-vivo seconds → in-silico steps.** The coarse-grained model evolves far faster
than real translation, so a **time-compression factor** maps seconds to MD steps:

```text
t_sim (ns)  =  t_s · 1e9 / scale_factor
n_steps     =  t_sim (ns) / dt(ns) ,   dt(ns) = dt_ps · 1e-3
```

A **larger `scale_factor` ⇒ fewer steps per residue ⇒ a faster run**, while preserving the
*relative* timing of fast vs. slow codons (the physics that matters for folding during
synthesis). Step counts may additionally be clamped to
`[min_steps_per_stage, max_steps_per_stage]` for tractability — a clamp on **MD steps
only**; the sampled dwell **times in seconds** are recorded untouched in `dwell_times.dat`.

**(e) Worked example (real 4c5c mRNA, Tutorial 8).** The first residues of
`4c5c_mrna.txt` and their `τ` from the default E. coli 310 K table:

| residue | codon | τ (s) |
|--------:|:-----:|------:|
| 1 | AUG | 0.133321 |
| 2 | ACU | 0.027566 |
| 3 | GAU | 0.038593 |
| 5 | AUU | 0.048617 |
| 6 | GCU | 0.019547 |

Take **residue L = 1** (the ribosome has just made residue 1, codon AUG). Stage 3 looks
ahead to **residue 2's** codon, **ACU** (τ = 0.027566 s), with `time_stage_1 = 0.00034 s`
and `time_stage_2 = 0.004201 s`:

```text
mean(t3) = τ(ACU) − time_stage_1 − time_stage_2
         = 0.027566 − 0.00034 − 0.004201
         = 0.023025 s          # mean of the exponential; sampled t3 = −mean·ln(U)
```

Then → MD steps at Tutorial 8's `scale_factor = 216564650`, `dt = 0.015 ps`:

```text
t_sim = 0.023025 s × 1e9 / 216564650 = 0.10632 ns
steps = 0.10632 ns / (0.015e-3 ns)   ≈ 7088 steps   (for the mean; then clamped to [400, 10000])
```

Stage 3 is where **per-codon variability** enters: stages 1 and 2 are fixed means, but
stage 3's mean swings with the next codon's `τ` (e.g. L = 5 → next codon GCU →
`0.019547 − 0.00034 − 0.004201 = 0.015006 s`). Edge case: if a very fast next codon makes
`τ(next) − t1 − t2 ≤ 0`, the mean is floored to `1e-9 s`.

### 5. After the last residue: ejection (and dissociation)

Once the final residue is added, the protein is complete. The simulation then runs a
**post-synthesis ejection phase** (`ejection_steps`): the C-terminus restraint is
**released** (the tether is cut), while the rigid ribosome and one-sided tunnel wall
remain. Biologically this is **termination** — release factors free the finished protein.
With the tether gone, the chain **diffuses out of the tunnel along +x** (the one-sided
wall biases motion forward) and clears the ribosome. An optional **dissociation** phase
(`dissociation_steps`) continues the free protein away from the ribosome. For a longer,
dedicated egress demonstration, raise `ejection_steps` (and `dissociation_steps`) in the
Tutorial 8 `csp_val.ini`.

### 6. Numerical integration and the stability guard

The chain is integrated with **flexible harmonic bonds** at `dt = 0.015 ps`. Flexible
(rather than rigid `AllBonds`) bonds are required because stage 1 seeds the new bead
~1 nm from its bond partner (A-site delivery), which a rigid distance constraint cannot
represent — a harmonic bond absorbs the stretch and the minimizer relaxes it.

The cost: at 15 fs the integration is only **marginally stable** for some configurations.
When a newly added residue forms a **stiff native (Gō) contact**, that contact's
vibrational period drops below what a 15 fs step can integrate and the dynamics
**diverge** (potential energy → ~10¹³ kJ/mol), corrupting that stage's frames. This is
**deterministic in the timestep, not random** (the reference protocol avoids it entirely by
using rigid `AllBonds`, which remove the fast bond mode).

**The fix** (`topo.csp.core.run_length`, the per-stage *stability guard*): each stage
is run in chunks while tracking the **maximum** |PotE|; if a stage diverges
(max |PotE| > 10⁹ kJ/mol) it is transparently **re-run with the timestep halved and the
step count doubled**. Because the physical dwell time is `n_steps · dt`, halving `dt` and
doubling `n_steps` **leaves the dwell time exactly unchanged** while stabilising the
integration (up to 6 halvings). The common case runs once at 15 fs. Watch for
`[stability] ...` lines in the log: those are stages auto-stabilised at a halved timestep.
Tutorial 8's full-length run (`csp_val.ini`, L = 1 → 306) validates this across the whole
chain (919 stages, zero blow-ups).

---

## Configuration reference (`csp.ini`)

```{note}
The tables below cover the keys used on this page. For the **complete, canonical
per-key reference** (both `csp.ini` and `cylinder.ini`, grouped by concern) see
{doc}`synthesis_control`.
```

CSP reads a single INI control file with one `[OPTIONS]` section
(`topo.csp.protocol.read_csp_config`). **Units are OpenMM defaults** — nm, ps, kJ/mol, K,
kJ/mol/nm² — and **dwell times are in seconds**. Integers may use `_` digit separators
(e.g. `200_000`).

```{tip}
For a compact, single-page tabular reference of every `csp.ini` option (grouped
by role, with types and defaults) — the synthesis analogue of the single-protein
{doc}`simulation_control` page — see {doc}`synthesis_control`. The tables below
repeat the same options with more inline commentary on the physics.
```

### Inputs & schedule

| Key | Required | Default | Meaning |
|-----|----------|---------|---------|
| `pdb_file` | **yes** | — | All-atom native PDB; topo builds the CG model from it. |
| `ribosome` | **yes** | — | Truncated CG ribosome PDB (P-/A-anchors + rigid scenery). |
| `L0` | no | `1` | Start nascent-chain length (omit/blank = start from a single residue). |
| `L_max` | no | full length | Final nascent length (omit/blank = synthesize the whole chain). |
| `mrna` | cond. | — | mRNA file (one codon per residue), **or** `fastest`/`slowest`/`median` to auto-build a synonymous-codon mRNA (see [Fastest / slowest / median mRNA](#fastest-slowest-synonymous-codon-mrna)). Required for per-codon timing (unless `codon_times` is a number). A real filename must not be `fastest`/`slowest`/`median`. |
| `codon_times` | cond. | — | Codon-timing key. A **table path** = per-codon timing (required, no bundled default -- pick one under `assets/csp/codon_dwell_times/`); a **positive number of seconds** = uniform codon time (every codon, no `mrna` needed). A table filename must **not** be a bare number. |
| `domain_def` | **yes** | — | `domain.yaml` — the protein's **contact-nscale definition** (per-domain / per-interface Gō well-depth scaling, the structure-based analog of the reference model's `nscal`). |
| `stride_output_file` | no | — | Precomputed STRIDE file (skips re-running STRIDE). |
| `outdir` | no | `synth_out` | Output root. |

### Codon kinetics

| Key | Default | Meaning |
|-----|---------|---------|
| `scale_factor` | `4331293` | In-vivo-seconds → in-silico-ns compression (larger = fewer steps = faster). |
| `time_stage_1` | `0.00034` | Mean peptidyl-transfer (peptide-bond) dwell, **seconds**. |
| `time_stage_2` | `0.004201` | Mean translocation dwell, **seconds**. |
| `random_seed` | — | Seed for the FPT sampler (reproducible schedules). |
| `max_steps_per_stage` | — (uncapped) | **Testing only** — upper clamp on each stage's MD step count (tutorials use a small value for speed). See note below. |
| `min_steps_per_stage` | `1` | **Testing only** — lower clamp on each stage's MD step count. See note below. |
| `ejection_steps` | `0` | Post-synthesis ejection-phase length (steps); `0` = skip. |
| `dissociation_steps` | `0` | Post-synthesis dissociation-phase length (steps); `0` = skip. |

```{warning}
**`max_steps_per_stage` / `min_steps_per_stage` are testing-only knobs and will be
removed in production.** They clamp the MD step count per stage so the tutorials finish
quickly, which **breaks the physical timescale mapping** — the integrated MD per stage no
longer matches the sampled dwell time. In a production run, leave them **unset** so the
step count is driven entirely by the kinetics (`scale_factor`, the codon times, and `dt`).
The sampled dwell **times in seconds** are always written to `dwell_times.dat` regardless
of the clamp.
```

### MD / ribosome mechanics (`RunParams` fields)

| Key | Default | Meaning |
|-----|---------|---------|
| `dt` | `0.015` | Timestep, ps. |
| `ref_t` | `310` | Temperature, K. |
| `tau_t` | `0.05` | Langevin friction, 1/ps. |
| `nstout` | `5000` | Trajectory/log output interval (steps). |
| `device` | `CPU` | `GPU` / `CPU`. |
| `ppn` | `1` | CPU threads (CPU platform). |
| `constraints` | `AllBonds` | Bond constraints. Equilibrium PTC seeding (always on, below) keeps rigid bonds stable; set `None` for flexible harmonic bonds instead. |
| `restraint_k` | `83680` | C-terminus harmonic restraint constant, kJ/mol/nm². |
| `minimize` | `yes` | Energy-minimize the seeded structure before each stage's MD. |
| `tunnel_wall` | `yes` | Apply the one-sided tunnel wall (floor below the synthesis point); plane auto-placed. |

There is **no `rigid_ribosome` key**: the `ribosome` PDB is required, and supplying it *is*
the signal to load it as **rigid (mass-0) scenery** (excluded volume + electrostatics on) —
so it is always treated as rigid, and the tunnel wall defaults on with it. `trna_tether` is
**forced off** by the CSP runner — CSP needs the switchable A↔P position restraint.

**PTC geometry is always optimized.** Each new residue is seeded at the optimal A-site
target — one equilibrium peptide bond (0.381 nm) from the previous C-terminus, clear of
the ribosome excluded volume (`optimal_ptc_targets`) — and those optimized A-/P-site
points are the restraint targets and the tunnel-wall plane. That equilibrium seeding is
why `constraints = AllBonds` (rigid) is the default and stable.

```{note}
**tRNA presence / naming.** The P-/A-anchors and the `optimal_ptc_targets` solve read the
ribosome's tRNA beads under fixed names (segids `PtR`/`AtR`, resid 76, beads
`R`/`P`/`BR2`; the acceptor must be a purine A, which carries the `BR2` bead). A ribosome
PDB with **no tRNA**, or with **differently-named** tRNA segments, currently fails with a
generic "expected exactly one bead" error. Handling missing/renamed tRNA (a clear
up-front error, and/or configurable segid/resid/bead names) is a tracked TODO
(`review/TODO.md`).
```

---

## Outputs

```text
<outdir>/
├── L_<L>/                  # ONE folder per residue L (consolidated layout)
│   ├── traj.psf            # nascent topology (shared across the 3 stages; f(L) only)
│   ├── native_1_<L>.pdb    # length-L native structure (shared; f(L) only)
│   ├── traj_s1.dcd         # (nascent-only) trajectory, stage 1  (s2/s3 likewise)
│   ├── traj_s1.log         # energies, stage 1; column 3 = potential energy (kJ/mol)
│   ├── traj_runinfo.log    # folded run-info: one [run:...]/[result:...] per stage
│   └── traj_final.pdb      # stage-3 final — seeds L+1 and is the resume-reload target
├── ejection/               # post-synthesis ejection phase (if ejection_steps > 0)
├── dissociation/           # post-synthesis free run (if dissociation_steps > 0)
├── dwell_times.dat         # per-residue dwell-time log / schedule (see below)
└── progress.log            # append-only DONE/RUNNING resume status (see below)
```

Each residue's three kinetic sub-stages share **one** `L_<L>/` directory: `traj.psf` and
`native_1_<L>.pdb` depend only on `L` (the A/P differences are *forces*, not atoms) and are
written once; the trajectories stay split per stage (`traj_s{1,2,3}.dcd`) to preserve the
kinetic-dwell boundaries; and only stage 3 writes `traj_final.pdb`. There is no per-stage
`.chk` — per-residue resume reloads `traj_final.pdb`, never a checkpoint.

**`dwell_times.dat`** records, per residue, the codon, the three sampled dwell **times in
seconds** (`t1`/`t2`/`t3`), their nanosecond equivalents, and the integer MD step counts —
the physical schedule, independent of any step clamp. This is the file to compare against a
reference run for quantitative validation (Tutorial 8). It is drawn **once** before the run
and doubles as the **immutable plan** for resume (its `#PTC` header pins the restraint
geometry); see {doc}`synthesis_resume`.

**Resuming an interrupted run.** A production synthesis is hours to days of wall time.
Re-invoking `topo-csp` on an interrupted `outdir` **continues from the last completed
residue** rather than restarting — the schedule and PTC geometry are re-read from
`dwell_times.dat` and the seed is reloaded from the last residue's `traj_final.pdb`, tracked
by `progress.log`. Resume is on by default (`resume = auto`). See {doc}`synthesis_resume`
for the full mechanism, guarantees, and an HPC requeue pattern.

### Console progress log

`topo-csp` prints one compact, column-aligned line per residue followed by one line per
sub-stage, so a long synthesis stays readable. Each stage line reports the wall-clock time
and the **total system potential energy of the last integrated step** (nascent chain +
rigid ribosome + all cross-interactions), a quick per-stage health signal:

```text
L=  1    AUG  dwell   0.02757 s  steps  400/1136/2000
  L=  1  stage 1 peptidyl-transfer     400 steps    1.19 s  PE=  +4.9674e+01 kJ/mol
  L=  1  stage 2 translocation        1136 steps    0.42 s  PE=  +4.8717e+01 kJ/mol
  L=  1  stage 3 tRNA-binding         2000 steps    0.59 s  PE=  +4.8365e+01 kJ/mol
```

The residue line's `steps` field and the per-stage `steps` are the **configured** step
counts; a stage that trips the stability guard silently reruns at a halved timestep with
double the steps (the `[stability] ...` lines above), and always prints its concise summary
line afterwards. The post-synthesis `ejection` / `dissociation` phases print the same
summary line.

Set **`TOPO_CSP_VERBOSE=1`** to restore the full per-stage banners (build block, seeded-
structure minimization, run-metadata path, elapsed time) — useful when debugging a single
length:

```bash
TOPO_CSP_VERBOSE=1 topo-csp -f csp.ini
```

MDAnalysis emits cosmetic `UserWarning`s (missing `CRYST1` unit cell, absent
`formalcharges`) each time topo slices and writes a CA-only PDB — once per stage. `topo-csp`
(like `topo-mdrun` and `topo-optimize`) silences all MDAnalysis warnings for the run via a
process-local filter, so they never reach the console; other warnings are unaffected.

**Movie.** Each stage writes a standalone trajectory (different lengths have different bead
counts, so they cannot be concatenated directly). `topo-csp-movie` stitches them — in
synthesis order, padding every frame to the final length and overlaying the static
ribosome — into one VMD-playable movie:

```bash
topo-csp-movie -o <outdir> --ribosome ribosome_trunc.pdb
# writes <outdir>/movie.psf, movie.dcd, movie.tcl, movie_ribosome.pdb
vmd -e <outdir>/movie.tcl
```

---

## Python API

```python
from topo.csp.protocol import run_continuous_synthesis, read_csp_config

# (a) drive it from an INI, exactly like the CLI. read_csp_config returns a
#     CSPConfig dataclass; unpack its fields into the call:
cfg = read_csp_config("csp.ini")
run_continuous_synthesis(
    cfg.pdb_file, cfg.ribosome,
    L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
    mrna=cfg.mrna, codon_time_table_path=cfg.codon_time_table_path,
    domain_def=cfg.domain_def, stride_output_file=cfg.stride_output_file,
    params=cfg.params,
)

# (b) or construct parameters directly (the ribosome PDB is always rigid scenery;
#     the tunnel wall plane is auto-derived from it):
from topo.csp.core import RunParams
params = RunParams(tunnel_wall=True,
                   scale_factor=4331293.0, random_seed=20240629, ejection_steps=50000)
run_continuous_synthesis("4c5c_model_clean.pdb", "ribosome_trunc.pdb",
                         L0=1, L_max=10, mrna="4c5c_mrna.txt", params=params)
```

See the {doc}`API reference <../topo>` for the autodocumented `topo.csp.protocol` and
`topo.csp.kinetics` modules.

---

## See also

- {doc}`cylinder_synthesis` — the analytic-tunnel variant (`topo-cylinder`): the same
  codon kinetics with a nascent-only system and a single MD segment per residue, with the
  explicit-bead ribosome replaced by a cylindrical bore through an infinite wall.
- {doc}`synthesis_control` — the concise `csp.ini` control-options reference.
- {doc}`../usage/model_theory` — the TOPO Gō-model force field in full (the RNC
  Hamiltonian CSP uses, restricted to the synthesized residues).
- {doc}`API reference <../topo>` — the shared low-level engine `topo.csp.core`
  (`run_length`, `RunParams`), `topo.csp.ribosome`, and `topo.csp.kinetics`.
- Tutorial 8 (`tutorials/08_ribosome_synthesis/`) — runnable, validated CSP examples on
  4c5c (smoke `csp_debug.ini` + full-length `csp_val.ini`) and P0CX28.
