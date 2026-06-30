# Design: protein synthesis (co-translational folding) in the TOPO model

> **Status:** draft / idea capture. This document defines the *problem* and the
> *design space* before any code is written. Sections marked **(FILL IN)** are for
> you to draft. The goal is a shared, reviewable spec; implementation comes later.

Branch: `translation`. Future module: `topo.csp` (a sibling of
`topo.mdrun` / `topo.optimize`), reusing `topo.engine`'s build → setup → finalize
helpers and supplying its own driver loop (see the note in `topo/engine.py`).


## 1. The biological problem

Proteins are synthesized **vectorially**, N-terminus first, one residue at a time,
by the ribosome. The nascent chain emerges through the ribosomal **exit tunnel**
and begins to fold **co-translationally** — often before the full sequence is
available. This can give folding pathways, intermediates, and yields that differ
from refolding of the full-length chain in bulk (the process simulated by the
current equilibrium/annealing runners).

Key physical features to capture (or deliberately approximate):

- **Vectorial growth.** The chain length increases over time, N→C; only the
  already-synthesized residues exist at a given moment.
- **Tethering at the PTC.** The C-terminal residue is anchored at the peptidyl
  transferase center / tunnel exit; the chain dangles/extrudes from there.
- **The exit tunnel + ribosome surface.** Excluded volume and geometric
  confinement from the ribosome shape early folding.
- **Elongation kinetics.** Residues are added at a finite codon-dependent rate;
  the competition between elongation time and folding time is the central
  observable.

What specifically do you want to study? - All of these: co-translational
folding pathway of a specific multidomain protein; effect of elongation rate on
domain-folding order; misfolding/entanglement during synthesis.


## 2. The modeling problem (decisions)

How the biology maps onto a Cα structure-based model. Each item below is a design
decision with options to weigh.

### 2.1 Growing the chain — consecutive rebuild-and-continue (AGREED)

On-the-fly modification of a live OpenMM `System`/`Context` is hard and brittle,
so we **do not** grow the chain inside one running simulation. Instead, synthesis
is modeled as a **sequence of independent simulations**, each rebuilding a
one-residue-longer system and continuing from the previous step's final state:

```
for L = L0, L0+1, ..., N_full:                  # L = current nascent-chain length
    1. BUILD a CG system: nascent chain = native residues 1..L  (+ rigid ribosome)
    2. SEED coordinates: residues 1..L-1 from step (L-1)'s final structure;
       residue L (the new C-terminus) placed at the PTC
    3. RUN N_L steps (Langevin), WRITE this step's trajectory/log/final structure
    4. L -> L+1   (the final structure of step L seeds step L+1)
```

Rationale and consequences:

- **No live resizing.** Each length is a clean, standalone TOPO run; the chain
  "grows" only between runs. This sidesteps `Context` resizing entirely.
- **N→C growth.** Synthesis is vectorial N→C, so a length-`L` chain is the
  **N-terminal residues 1..L** of the target protein, and each step **appends a
  new C-terminal residue** `L+1` (the new C-terminus sits at the PTC). The chain
  is anchored at its C-terminus (§2.2) and extrudes N-terminus-first through the
  tunnel.
- **Native contacts to future residues are simply absent.** The force field for a
  length-`L` chain is built from the native structure **restricted to residues
  1..L**, so contacts involving not-yet-synthesized residues (`> L`) do not exist
  — no masking needed. (This resolves the earlier open question.)
- **Continuity between steps.** The L−1 already-synthesized residues carry their
  **coordinates** across the rebuild (via `init_position`); the new residue and
  the velocities are initialized fresh (see §2.4 / §5). The rigid ribosome is
  identical in every build.
- **Per-length outputs.** Each step writes its own trajectory/log/final structure
  (e.g. into a per-length output folder), giving one trajectory per chain length
  and a clear final structure to seed the next length.

**New-residue placement and the per-step mechanism are specified in §2.5**
(simplified elongation step): place the new residue at the A-site anchor, restrain
the current C-terminus to the P-site anchor, run. Both tRNA acceptor ends survive
truncation (`FILES.md`), providing the anchors.

The starting length `L0` (cold-start layout is agreed — §2.5).

### 2.2 Representing the ribosome — explicit rigid CG ribosome (AGREED)

**Decision (following O'Brien et al. 2012, see §3 + References):** an **explicit
rigid** CG ribosome — a fixed (frozen) bead cloud from a ribosome structure that
interacts with the nascent chain by excluded volume (+ electrostatics). The
ribosome is **held fixed**, so no intra-/inter-ribosome forces are computed; only
ribosome ↔ nascent-chain interactions matter, and those are **electrostatics +
excluded volume only** (no native/attractive contacts to the ribosome).
Concretely:

- **Rigid body.** Freeze ribosome particles (zero mass / not integrated, or an
  OpenMM `CMMotionRemover`-free fixed set). Their coordinates come from a
  ribosome structure and are used **as-is** — this is exactly why `init_position`
  was changed to stop shifting coordinates (absolute tunnel/exit geometry must be
  preserved).
- **RNA representation** (already parameterized in TOPO: `P`, `R`, `BR`):
  pyrimidines = 3 sites, purines = 4 sites — phosphate (**q = −1e**), ribose-ring
  centroid, and one centroid per conjugated base ring.
- **Ribosome ↔ NC interaction = a separate repulsive force**, pure
  `ε·(σ_ij/r)^12` (`ε = 0.000132 kcal/mol`, `σ_ij = (σ_i + σ_j)/2` from
  per-particle `radii`), restricted to {nascent}×{ribosome}. Electrostatics are
  folded into the existing Yukawa force over all charges. Full spec + the
  particle-ordering invariant: §3.2.

**Source structure.** `structures/4v9d_50S_PtR_5jte_AtR_model.pdb` — the *E. coli*
50S large subunit + tRNA, **re-oriented with the exit-tunnel central line on the
X-axis**. It is coarse-grained (`cg_ribosome.py`) and **truncated**
(`truncate_ribosome.py`, x_exit=58 → 14,662 → 4,576 beads, matching the O'Brien
crop). Procedure, parameters, and the structure inventory are in
[`FILES.md`](FILES.md); we implement both ourselves (no external code). The
C-terminal tether/PTC anchoring is specified in §2.5 (restrain the current
C-terminus to the P-anchor).

### 2.3 Elongation schedule
With the consecutive-build protocol (§2.1), the schedule is just the **per-step
run length `N_L`** — how many integration steps each length is simulated before
the next residue is added:

- **Constant rate:** the same `N_L = N` for every length (one dwell time per
  residue).
- **Variable rate:** a per-residue `N_L` (e.g. a codon-/residue-dependent
  dwell-time table), to model fast vs. slow codons.
- Map MD steps ↔ real elongation time (a few residues/s in vivo) — an
  effective-time mapping for a CG model; document the assumptions.

**Decision (AGREED, testing):** a single **constant `n_steps_per_residue`** — the
same `N_L = N` for every length. Keep it **small for traceable first runs**, e.g.
`n_steps_per_residue = 1000`. (A variable per-residue/per-codon schedule is a
later option, not needed now.)

### 2.4 Force field per step, and step-to-step continuity
*(Full term-by-term RNC force field — and what to keep vs. omit — is in §3.)*

Within each step the model is just the ordinary TOPO force field built for the
length-`L` chain plus the rigid ribosome — nothing special "during" a step,
because the chain does not change mid-run. The synthesis-specific handling is all
**at the rebuild boundary** between steps:

- **Each build is a normal TOPO build** of native residues 1..L (bonds, Gaussian
  angle, periodic torsion, Yukawa, structure-based contacts among 1..L) — the new
  C-terminal residue's bonded/non-bonded terms simply exist in the new system.
- **Coordinates carry over; velocities do not.** Residues 1..L−1 take their
  coordinates from step (L−1)'s final structure via `init_position` (which now
  uses coordinates as-is). Because the particle count changed, velocities cannot
  be loaded from a checkpoint — they are re-drawn from the Boltzmann distribution
  at `ref_t` each step (a short re-equilibration at the start of each step).
- **C-terminal tether.** The C-terminal bead is restrained near the PTC (it is
  covalently attached to the P-site tRNA); see §2.2.
- Domain scaling (`n_scale`) unchanged from `domain.yaml` (restricted to the
  residues present).

**Decisions (AGREED):**
- **No separate tunnel-confinement term** — the electrostatics + excluded volume
  (ribosome↔nascent) already confine the chain in the tunnel.
- **No per-step re-equilibration discard** — keep all frames of each step's run.
- **Test protein: P0CX28** (106-residue single domain,
  `tutorials/01_single_domain_quickstart/P0CX28_clean.pdb`) as the nascent chain
  for the first runs.

(C-terminal restraint is specified in §2.5.)


### 2.5 Simplified elongation step (AGREED)

A single-stage simplification of O'Brien's 3-stage cycle: their A-site binding →
peptidyl-transfer → translocation are **collapsed into one step**.

**Anchors** (fixed points taken from the rigid truncated ribosome):
- **P-anchor** = `PtR` resid-76 `R` bead + buffer — the held C-terminus site.
- **A-anchor** = `AtR` resid-76 `R` bead + buffer — where a new residue is delivered.

**Per step (L−1 → L):**
1. **Place** the new residue L at the **A-anchor**.
2. **Build** the length-L TOPO model on native residues **1..L**. This
   *automatically* includes the L−1↔L bond, all bonded terms, and the native
   contact map among 1..L — no manual bond needed.
3. **Restrain only the current C-terminus:** harmonic restraint **L → P-anchor**
   (k ≈ 200 kcal/mol/Å²). The hand-off (drop L−1's restraint, add L's) is
   automatic because each rebuilt step restrains only L.
4. **Run N steps** (small/hard-coded for tests); save final → seed L+1.

**Behavior.** L *starts* at A but is *restrained to* P, so the MD pulls it
**A→P**; occupying the fixed P-anchor pushes earlier residues up the tunnel (+x).
So one step folds **peptidyl-transfer + translocation** and the chain **advances
one register per cycle**; the C-terminus is always held at the *same fixed
P-anchor* while the chain extrudes past it.

**Placement buffer must clear excluded volume.** The buffer (the offset when
placing the new bead at the A-anchor, and the cold-start spacing) must be **large
enough that the new bead does not overlap neighboring beads** (the current
C-terminus, tRNA/ribosome beads). An overlap gives a large `(σ/r)¹²` repulsive
force that ejects the bead and crashes the run. So set the buffer ≥ the
excluded-volume contact distance of the surrounding beads; a short minimization at
the start of each step also relaxes residual clashes.

**Cold start (L0).** Lay residues 1..L0 **extended along the tunnel axis** from
the P-anchor — C-terminus (residue L0) at the P-anchor, N-terminus toward the exit
(+x) — at one CG bond length spacing, so no stretched bonds.

**v1 scope.** Loop mechanics only — nascent chain + the two anchor points (for
placement and the restraint target). **No ribosome excluded-volume forces yet**
(v2 adds the rigid ribosome, §2.2/§3.2).


## 3. Force field for the RNC system

The RNC force field is the **same model family as TOPO** (the Best–Chen–Hummer
variant of the Karanicolas–Brooks Cα Gō model; cf. O'Brien et al. 2011/2012). This
section records the exact energy terms and how each maps onto existing TOPO
machinery.

### 3.1 The Hamiltonian

| Term | Functional form | TOPO equivalent |
|------|-----------------|-----------------|
| Cα–Cα bonds | `k_b (r − r₀)²` | bonds |
| Dihedrals | `Σ_j k_φ[1 + cos(j φ − δ)]` | periodic torsion |
| Bond angles | double-Gaussian `−(1/γ) ln[e^… + e^…]` | Gaussian angle |
| Electrostatics | Debye–Hückel `q_i q_j /(4πε₀ ε_r r) · e^(−r/lᴅ)` | Yukawa |
| Native LJ, **intra-NC** | `ε_ij [13(σ/r)¹² − 18(σ/r)¹⁰ + 4(σ/r)⁶]`, `ε_ij ∝ BT` | structure-based contacts |
| **Non-native** (intra-, inter-protein, protein–RNA) | `ε_ij (σ_ij / r)¹²`, `ε_ij = 0.000132 kcal/mol` | non-native excluded volume |

So the RNC Hamiltonian = bonds + dihedrals + angles + Debye–Hückel + intra-NC
native 12-10-6 + non-native/excluded-volume `(σ/r)¹²` — every term already exists
in TOPO.

### 3.2 Ribosome ↔ nascent-chain interactions

No attractive (native) contacts to the ribosome. Ribosome↔nascent is a **separate
force** from the nascent-chain contacts (see the ordering/force invariant below):

- **Repulsive excluded volume only:** a dedicated `CustomNonbondedForce` with the
  pure form `ε·(σ_ij/r)¹²` (no attractive well), `ε = 0.000132 kcal/mol`,
  `σ_ij = (σ_i + σ_j)/2` where `σ_i` is the **per-particle `radii` from
  `topo.parameters.model_parameters`** (amino acids; P/R/BR = 0.710 nm). Cutoff
  2.0 nm, switch 1.8 nm. Restricted to **{nascent} × {ribosome}** via an
  interaction group.
- This is the paper's written form exactly. (The reference *code* used the
  12-10-6 with `ε = 0.000132` — a well ~5000× shallower than kT, i.e. numerically
  the same — but pure `(σ/r)¹²` is cleaner: strict excluded volume, no well.)
- **Electrostatics:** folded into the **existing Yukawa** `CustomNonbondedForce`,
  extended over all charged particles (nascent + ribosome rRNA phosphates −1e and
  charged residues). No separate electrostatic term is needed.

#### Particle ordering & force-separation invariant

- **Nascent chain occupies indices `0 … L-1`; the ribosome is appended at
  `L … N-1`** (as in O'Brien: `prot[1:L] + ribo`). This ordering is **required**,
  not cosmetic.
- TOPO's **native/non-native contact force is particle-indexed**: its tables are
  `L×L` and particle `i` carries `id = i` (`eps_table(id1,id2)` are
  structure-specific per-residue-pair values). The nascent block must be `0..L-1`
  for the table to align.
- **Table size vs. particle count are independent.** OpenMM requires one
  `addParticle` per *System* particle in every `CustomNonbondedForce` (so the
  contact force must list all `N` particles), **but the table stays `L×L`** — its
  size is set by the `id` values that are actually looked up, not by `N`. The
  contact force is restricted to **{nascent}×{nascent}** by an interaction group;
  ribosome beads get an `addParticle` entry with a **dummy in-range `id = 0`** that
  is **never read** (their pairs are never evaluated by this force), so the `L×L`
  table never needs ribosome rows/columns.

  Example (L=3 nascent, 2 ribosome, N=5):
  ```
  contact force:  Discrete2DFunction(3,3);  ids=[0,1,2,0,0];  group=([0,1,2],[0,1,2])
  (σ/r)^12 force: per-particle σ for all 5;                    group=([0,1,2],[3,4])
  ```
- **Ribosome↔nascent excluded volume is the separate `(σ/r)¹²` force above** — it
  does not touch the `L×L` table, so the table stays nascent-only and aligned.

### 3.3 Mapping to existing TOPO defaults

| Quantity | Reference paper | TOPO default | Status |
|----------|-----------------|--------------|--------|
| Debye length `lᴅ` | 10 Å | 1.0 nm | ✅ identical |
| Dielectric `ε_r` | 78.5 | 78.5 | ✅ identical |
| Charges | K/R +1, E/D −1 | ARG/LYS +1, ASP/GLU −1 | ✅ identical |
| Native well depths | ∝ Betancourt–Thirumalai | BT potential | ✅ same source |
| Non-native `ε` | 0.000132 kcal/mol | 0.000132 kcal/mol | ✅ identical |
| Native contact def. | side-chain heavy atoms < 4.5 Å | 4.5 Å (SS/BS) | ✅ same |
| LJ stability scaling | `λ_BT` (≈1.475–1.5) | `n_scale` / `strength` (+ optimizer) | ✅ same concept |
| RNA beads | P (q=−1), R, BR (3/4 sites) | `P`, `R`, `BR` already defined (q_P=−1, radii 0.71 nm) | ✅ already present |

Native well position `R_ij` = the native Cα–Cα distance (the 12-10-6 minimum sits
at the crystal distance), consistent with both TOPO and the reference model.

### 3.4 What actually needs building (vs. reused)

- **Reuse unchanged:** bonds, angles, torsions, Yukawa electrostatics, native
  12-10-6 contacts, non-native excluded volume, BT potential, `n_scale`, and the
  P/R/BR RNA parameters — all already in TOPO.
- **Add:** (1) the **build-once-subset** nascent-chain nonbonded (§3.5); (2) a
  **rigid ribosome** (load RNA P/R/BR + ribosomal-protein beads as fixed,
  non-integrated particles; ribosome↔NC via the existing Yukawa + the separate
  `(σ/r)¹²` force of §3.2; never compute intra-ribosome forces); (3) the
  elongation driver (§2.5).


### 3.5 Building the nascent-chain nonbonded (build once, subset to L×L) (AGREED)

The nascent-chain native/non-native contacts are built **once** from the full
native reference (`pdb_file`, all `N_full` residues) and **subset** to the current
length each step.

1. **Once at startup:** run the contact builder on the **full** native structure
   → `R_full`, `eps_full` (`N_full × N_full`), with STRIDE and `n_scale` applied;
   cache (including the STRIDE output).
2. **Per length L:** use the top-left block `R_full[:L, :L]`, `eps_full[:L, :L]`
   for the `L×L` tabulated `CustomNonbondedForce` (particle `i` → `id = i` =
   native residue `i+1`). Native contacts among 1..L are exactly the full-structure
   contacts restricted to 1..L; contacts to residues `> L` are absent; going
   L→L+1 only **adds residue L+1's row/column** (contacts "revealed").

**Why build-once-subset (not rebuild per L):** efficiency (STRIDE + heavy-atom
analysis run once); consistency (the contact map for 1..L is identical at every
length — a truncated-PDB rebuild would shift STRIDE/σ near the new C-terminus, so
existing contacts could silently change); and it matches O'Brien (`prot[1:L]`).

**Excluded-volume σ — two sources, by interaction class (matches O'Brien):**

| Interaction | σ source | distance `R_ij` |
|---|---|---|
| nascent native contact | full native structure | native Cα–Cα distance |
| nascent–nascent non-native | TOPO **structure-derived** σ on the **full** structure, **fixed across L** | `0.5(σ_i+σ_j)` |
| nascent–ribosome excluded vol. | **per-particle `radii` from `model_parameters`** | `0.5(σ_i+σ_j)` |
| ribosome–ribosome | — (rigid; not computed) | — |

So a nascent bead carries a structure-derived σ for intra-chain non-native and its
`model_parameters` radius for ribosome contacts — the same split O'Brien uses
(structure-based σ intra-protein; per-residue collision diameters
inter-protein / protein–RNA).

> **`model_parameters['radii']` = the O'Brien collision diameters, validated
> against both the code and the paper.** `radii × 10` (Å) equals the O'Brien
> **code's** `2 × rvdw` (per-residue sidechain vdw radius in
> `create_cg_protein_model.py`) to within **≤0.01 Å** (3-decimal-nm rounding); the
> SI **Table S2** is the same numbers rounded to 0.1 Å (RNA `P/R/BR = 7.1 Å`
> exact). So our `radii` come from the code at full precision — *not* a different
> source or an updated parameter set. `radii` holds the **collision diameter σ**
> (= 2×rvdw, not a radius); since the code's inter-molecular distance is
> `rvdw_i + rvdw_j`, our combining rule `σ_ij = (radii_i + radii_j)/2` is exactly
> O'Brien's — use `radii` directly as σ in the ribosome–NC `(σ/r)¹²` force.

**Coordinates are decoupled from the force field:** the native `pdb_file` sets only
the contact distances (well minima) and bonded terms; initial coordinates for
length L = previous (L−1) final structure + new residue placed at the A-anchor
(§2.5).


## 4. Reuse of existing infrastructure

What already exists and can be leveraged (avoid reinventing):

- **`topo.engine`** — `build_system`, `setup_simulation`, `finalize_simulation`.
  The translation runner is an **outer loop over chain length `L`** that, each
  iteration, calls these helpers to build + run a standalone length-`L` simulation
  (mirrors how `topo.mdrun` wraps a single run; here we wrap one run *per length*).
- **`topo.mdrun`** — each length is essentially one ordinary `mdrun` invocation;
  the translation driver mostly assembles each step's inputs (length-`L`
  structure, seeded coordinates) and calls the same build→run→finalize path.
- **`init_position` as-is coordinates** — used both to keep the rigid ribosome
  geometry and to **seed each step** from the previous length's final structure
  (no coordinate shifting).
- **`topo.reporter` / runinfo** — logging + provenance, reuse directly (one set
  per length).
- **Restart machinery** — restart applies *within* a length (resume an
  interrupted length's run from its checkpoint) and *across* lengths (skip lengths
  already completed, resume at the next length). No live resizing is involved, so
  the size-change concern is gone (each length is its own fixed-size system).

**(FILL IN)** anything else (analysis: `topo.analysis` Q per domain as residues
appear; `scripts/entangled_lifetime` for entanglements).


## 5. Open questions / risks

**Resolved:** live system resizing (each length is a standalone run, §2.1); native
contacts to future residues (absent — only 1..L built, §2.1); new-residue
placement (A-anchor, §2.5) and cold-start layout (extended along tunnel, §2.5);
velocity handling / re-equilibration (re-drawn Boltzmann each step, **no discard**,
§2.4); C-terminal tether (restrain current C-terminus to P-anchor, k≈200, §2.5);
performance (build-once-subset, no per-step STRIDE, §3.5). Remaining:

- **Starting length `L0`.** The numeric value to start synthesis from (layout is
  agreed). For the P0CX28 test, pick a small L0 (e.g. a few residues already in
  the tunnel).
- **Effective timescale mapping (deferred).** CG dynamics ↔ real elongation rate.
  Testing uses a constant `n_steps_per_residue = 1000` (§2.3); the in-vivo
  dwell-time mapping is a later concern, not needed for the first runs.


## 6. Proposed phases (rough)

The protocol is the consecutive rebuild-and-continue loop (§2.1); the target
environment is the rigid truncated ribosome (§2.2) with the §3 force field. Build
up to it:

1. **Single-step build.** Get one length-`L` run working: build native residues
   1..L as a TOPO model, run, write outputs. (Pure reuse of `topo.engine`.)
2. **The elongation loop, tether-only.** Wrap step 1 in the `L → L+1` loop:
   seed each step from the previous final structure, append residue `L+1` at the
   PTC, re-draw velocities, constant `N_L`. Use a simple PTC tether (no explicit
   ribosome yet) on one small protein. Goal: end-to-end N→C synthesis producing
   one trajectory per length.
3. **Rigid explicit ribosome.** Coarse-grain (`cg_ribosome.py`) + truncate the
   ribosome (`FILES.md`), add it as fixed particles; ribosome↔NC = Debye–Hückel +
   non-native `(σ/r)¹²` (§3.2). Optionally a variable elongation schedule.
4. **Analysis + a tutorial.** Once it runs, add a tutorial (Tutorial 7?) and
   `usage/` docs, with Q-vs-length and folding-order observables.

**(FILL IN)** adjust scope/order to your actual study.


## 7. References / prior work

- O'Brien, E. P.; Christodoulou, J.; Vendruscolo, M.; Dobson, C. M. **J. Am.
  Chem. Soc.** 2012, 134, ja302305u. — the ribosome–nascent-chain force field this
  design adapts (SI PDF kept under `topo/csp/theory/`).
- Karanicolas, J.; Brooks, C. L. **Protein Sci.** 2002, 11, 2351. — the base Cα
  Gō model.
- Best, R. B.; Chen, Y. G.; Hummer, G. **Structure** 2005, 13, 1755. — the
  transferable variant used by the reference RNC model.
- Betancourt, M. R.; Thirumalai, D. **Protein Sci.** 1999, 8, 361. — the BT
  statistical potential used for native well depths.
- O'Brien, E. P.; Christodoulou, J.; Vendruscolo, M.; Dobson, C. M. **J. Am.
  Chem. Soc.** 2011, 133, 513. — the CG ribosomal-RNA representation (P/R/BR) and
  the ribosome-truncation procedure we adopt.

> We adopt the **procedure/theory** from these references but **implement
> everything ourselves** — no external (O'Brien-lab) code is used.

**(FILL IN)** any additional group methods/papers to match.
