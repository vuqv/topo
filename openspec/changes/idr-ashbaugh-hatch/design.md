# Design

## The invariant everything else follows from

> **A force's interaction groups are set exactly once, by whoever creates it.
> No other code may add to them.**

OpenMM takes the **union** of a force's interaction groups. So any downstream code
that asserts a broader group over a force it did not create silently re-admits pairs
that force no longer owns, and those pairs are then evaluated twice — once by each
force — with no error raised.

Three sites in the codebase assert interaction groups. All three must respect the
rule:

| site | what it must do |
|---|---|
| the model builder | set each force's domain at construction |
| `csp.ribosome.append_ribosome` | create the nascent↔ribosome force with its own group; set a group on the contact force **only if it has none** |
| `utils.multichain.replicate_system_intra_only` | replicate the template's groups with a per-copy offset; fall back to a blanket `{copy}×{copy}` **only if the template had none** |

The cheapest check that the rule holds is energetic: N non-interacting copies must
have exactly N × the energy of one copy. Any double-count breaks that identity.

## The potential

With `L(r) = 4[(σ/r)¹² − (σ/r)⁶]` and `σ = 2^(−1/6)·R_ij`:

```
r ≤ R_ij :  U = eps_EV·L(r) + (eps_EV − eps_ij)
r >  R_ij :  U = eps_ij·L(r)
```

Properties, all verifiable by evaluating the shipped expression:

- minimum at `r = R_ij`, depth exactly `−eps_ij`, for any `eps_EV`;
- `dU/dr = 0` on both sides of `R_ij`, so the join is C¹;
- `eps_ij = 0` ⇒ `U(r > R) ≡ 0` — a clean WCA core of physical size;
- `eps_ij = eps_EV` ⇒ the plain 12-6.

The core position (`U = kT`) is set by `eps_EV` and moves only 0.846 R → 0.819 R as
`eps_ij` runs 0 → 2 kJ/mol — a 3.2 % residual coupling from the continuity shift,
against 58 % in the coupled 12-10-6.

### Why the functional form has to change: two separate defects

Moving IDR-involving pairs off the 12-10-6 bundles **two independent changes**. Both
are necessary; neither is sufficient alone.

#### Defect 1 — one `eps` scales both the wall and the well

In `U = eps·[13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶]`, the same `eps` multiplies the
repulsive `r^-12` term and the attractive well. There is no way to ask for
"more attraction at the same bead size". Measuring the core as the radius where
`U = +kT`:

| cell | eps_ij (kJ/mol) | core | well depth |
|---|---|---|---|
| `idr_scale = 0` (the `NON_NATIVE` floor) | 0.00055 | **0.583 R** | 0.000 kT |
| `eps_gen = 1.5` | 1.50 | 0.903 R | 0.60 kT |
| `idr_scale = 1`, mean BT | 2.42 | **0.912 R** | 0.97 kT |

Switching on *any* attraction moves the core from 0.58 R to ~0.91 R — a **56 %
increase in radius, ~3.7× in excluded volume** — to buy a well shallower than `kT`.
The excluded-volume gain dominates and the chain expands.

This also means the apparently attractive "zero attraction" limit is not a
self-avoiding chain of physical thickness: it is a chain whose beads are *thin
because they have no energy*. Those two things should be independent and are not.

The Ashbaugh–Hatch split separates them: `eps_EV` sets the core, `eps_ij` sets the
depth, and the core moves only 0.846 R → 0.819 R (3.2 %) as `eps_ij` runs 0 → 2
kJ/mol. The residual comes from the continuity shift `(eps_EV − eps_ij)`, and it runs
the *right* way — more attraction slightly softens the core rather than hardening it.

#### Defect 2 — the desolvation barrier

The Karanicolas–Brooks 12-10-6 is not a Lennard-Jones potential with a different
exponent. Its coefficients (13, −18, +4) place a **repulsive bump beyond the well**:

```
F(r) = 13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶
F crosses zero at r = 1.25 R,  peaks at r = 1.45 R  with F = +0.143
```

so the pair must climb `0.143·eps_ij` before it can reach the well. That barrier is
deliberate and physically motivated **for a native contact**: forming one requires
expelling the solvent between two residues, which is a real free-energy barrier, and
it is what gives Gō models their two-state contact kinetics.

Between two disordered residues there is no native contact and no desolvation event
to represent. Worse, the barrier actively fights the attraction, because it sits at
larger `r` than the well and therefore carries more `4πr²` weight in

```
B₂ = −½ ∫ (e^(−U/kT) − 1) 4πr² dr
```

`B₂ > 0` means good solvent (the chain swells), `B₂ = 0` is the θ point, `B₂ < 0` is
collapse. Scaling the tail amplifies barrier and well together, and over most of the
usable range the barrier wins:

| eps_att (kJ/mol) | B₂/hs — 12-10-6 split | B₂/hs — 12-6 AH |
|---|---|---|
| 0.0 | 0.832 | 0.744 |
| 0.5 | 0.911 | 0.340 |
| 1.0 | 0.967 | −0.109 |
| 2.0 | 0.993 | −1.169 |
| 3.0 | 0.869 | −2.516 |
| **θ point (B₂ = 0)** | **4.84 kJ/mol (1.94 kT)** | **0.88 kJ/mol (0.35 kT)** |

With the barrier, `B₂` *rises* up to `eps_att ≈ 1.5` — adding attraction makes the
chain **more** swollen — and θ is not reached until ~2 kT, an implausibly deep
contact for a generic IDP interaction. Without it, `B₂` falls monotonically from the
first increment and θ arrives at 0.35 kT, which the BT table reaches at
`idr_scale ≈ 0.32`.

So decoupling the core alone would **not** have fixed the model: it removes the
excluded-volume inflation but leaves a potential whose net effect is still repulsive
over the range anyone would use.

#### Why Ashbaugh–Hatch specifically

A plain 12-6 with a single `eps` would re-introduce Defect 1. The Ashbaugh–Hatch
construction splits the potential at its minimum, holding the repulsive branch fixed
and scaling only the attractive branch, which is exactly the separation required. It
also yields a clean WCA limit at `eps_ij = 0` — a purely repulsive bead of physical
size — so "no attraction" and "no excluded volume" finally become independent
statements.

This is the standard form in dedicated IDP force fields (HPS/Dignon, Ashbaugh–Hatch),
and it is *why* the hydropathy parameter behaves as a solvent-quality dial there and
did not here.

#### Empirical confirmation

The prediction is directional and was tested. Under the coupled 12-10-6, `ν` **rose**
with `idr_scale` (0.605 → 0.723 across 0 → 1). Under the Ashbaugh–Hatch 12-6, `ν`
**falls** (0.637 → 0.276 across 0 → 0.30), sweeping through the experimental 0.551
and out the other side into collapse. Same benchmark, same 18 proteins, opposite
sign — the functional form, not the parameter value, was the defect.

#### What is preserved

The folded model is untouched. Folded–folded pairs keep the 12-10-6 *with* its
desolvation barrier, because there the barrier is doing the job it was designed for.
The change is scoped precisely to pairs where a native contact does not exist.

## Parameters

| symbol | key | role | default |
|---|---|---|---|
| `R_ij` | (from the radius table) | length scale | `Rmin/2_i + Rmin/2_j` |
| `eps_ij` | `idr_scale` | well depth → ν | 0.10 |
| `eps_EV` | `eps_ev_kj` | repulsive core → bead size | 0.8368 |

`eps_ev_kj = 0.8368 kJ/mol` is the HPS/Dignon value (0.2 kcal/mol). It is
**compatible by construction**: `σ = (Rmin/2_i + Rmin/2_j)/2^(1/6)` reproduces the
published HPS per-residue σ to a mean −0.1 % (max 1.5 %, glycine exact).

`eps_EV` is a weak handle on size — the core goes as `eps_EV^(1/12)`, so a 100×
change moves the bead only ~46 %. Use the radius table if bead size must change.

## Decisions

**`eps_gen_kj` is removed, not deprecated.** A legacy `eps_gen_kj: 0.0` is accepted
(it was a no-op, and the calibration runs carry it); any non-zero value raises with
the `idr_scale` equivalent, rather than silently changing the physics of an old
config.

**Radii use `Rmin/2`, not a σ-radius.** The `R` slot of both 12-10-6 and the AH form
is an Rmin. The per-AA table value is already an `Rmin/2`, and
`calculate_rmin_2_values` produces the same quantity for folded beads, so both
populations carry one quantity and the sum rule is coherent. Passing
`R_ij` through `σ = 2^(−1/6)R_ij` recovers `rvdw_i + rvdw_j`, the van der Waals sum
— which is why the σ values match HPS.

**CSP: the nascent↔ribosome force stays 12-10-6 for every nascent bead**, disordered
or not. It is a separate force created by `append_ribosome`, so this holds
structurally with no `if disordered` branch anywhere. `RIBO_NC_EPS_KJ` is unchanged,
preserving the validated 4c5c reproduction. The consequence — a disordered bead is
on a different footing toward the ribosome than toward the chain — is deliberate and
should be recorded in a comment so it is not later "fixed".

**Declaring an IDR deletes cross-boundary native contacts.** This is intended (a
disordered residue's crystal contacts are artifacts) but is not free: for 4c5c with
residues 1–40 disordered, 88 of 819 contacts are removed and 34 folded residues lose
at least one contact. The folded remainder is a less stabilized fold than the full
Gō model.

## Risks

Each risk below is stated as: what would go wrong, how you would (or would not)
notice, and what prevents it. All three failure modes are **silent** — the
simulation runs to completion and produces plausible-looking numbers.

### 1. A pair evaluated by both contact forces

**What goes wrong.** OpenMM sums forces independently, and takes the *union* of a
force's interaction groups. If any code adds a broader group to a force after
construction, pairs that force no longer owns are re-admitted — and they are then
computed by both forces and added together.

Concretely, an IDR pair would receive the Ashbaugh–Hatch well **and** the coupled
12-10-6, which reinstates precisely the excluded-volume coupling this whole change
exists to remove. A folded native contact would receive its Gō energy **and** a
spurious AH well on top.

**How you would notice.** You would not. No exception, no warning; energies are
merely wrong, by an amount that depends on how many pairs overlap.

**What prevents it.** The ownership invariant at the top of this document, plus two
tests: explicit group-disjointness assertions, and the energy identity — N
non-interacting copies must have exactly N × the single-copy energy, which no
double-count can satisfy. Both tests must be shown to *fail* when the fix is
reverted; otherwise they are not testing what they claim.

### 2. Stale `idr_scale` values in existing configs

`idr_scale` scales the well depth `eps_ij = idr_scale · eps_BT`. Its calibrated value
is **0.10**; the θ point is at **~0.32**, above which the chain collapses.

The previous default was `1.0`, so existing `domain.yaml` files often carry
`idr_scale: 1.0` explicitly. That is legal YAML and a legal value, so it parses
without complaint — and under this force it is ~3× past θ, giving a collapsed
globule.

*The fix is trivial:* delete the key and take the new default, or set `0.10`. No
re-derivation is needed — the value is calibrated (see `proposal.md`).

*What cannot be carried over* is a **custom** `idr_scale` that someone fitted to
their own system against the previous functional form. That number has no meaning
here and would need re-fitting. This is a narrow case; the shipped default covers
everyone else.

*Prevented by:* nothing automatic — the number alone cannot say which model it was
tuned for. State the new value in the migration note, and note that anything above
~0.32 is past θ and should be treated as suspect. An implementer may choose to warn
above that threshold; left open.

### 3. CSP half-applying the disorder model

**What goes wrong.** `build_nonbonded_interaction` parses a `disordered:` section
whenever one is present in the `domain_def`, regardless of which caller asked. CSP
passes `domain_def` straight through. So if CSP does **not** also forward the IDR
handle to the force builder, `apply_disorder` still runs — the matrices get IDR radii
and IDR well depths — but the Ashbaugh–Hatch force is never created. Those pairs are
then evaluated by the **coupled 12-10-6** at the new calibrated `idr_scale = 0.10`,
a combination that has never been calibrated and puts the core near 0.7 R: a very
thin, nearly-phantom chain.

**How you would notice.** You would not. Synthesis completes with finite energies.

**What prevents it.** Task 5 wires CSP fully, so the case cannot arise. If that work
is deferred, CSP must instead **refuse** a `domain_def` carrying a `disordered:`
section, with a clear message. The one thing it must not do is half-apply.

### 4. The benchmark split is a judgement, not a datum

**What goes wrong.** The calibration excludes 6 of the 24 proteins (CspTm, R15, R17,
hCyp, Protein-L, sNase) as foldable globular proteins whose published `Rg` is a
folded-state value. That classification is made by protein identity and is not
recorded in the dataset. If it is wrong for any of them, the fitted `idr_scale`
shifts — including them moves the optimum measurably (the all-24 objective put
`eps_gen*` at 3.25 where the 18-IDP objective put it at 2.00).

**How you would notice.** Only by checking the provenance of the experimental values.

**What prevents it.** Nothing automatic. The split, its members and its rationale are
stated in the docs so the choice is visible and can be revisited, and the 6 are always
reported as a control rather than dropped silently.

### 5. Decoupling is very good, not perfect

**What goes wrong.** The continuity shift `(eps_EV − eps_ij)` means the repulsive
branch still moves slightly with the well depth: the core runs 0.846 R → 0.819 R as
`eps_ij` goes 0 → 2 kJ/mol. That is a 3.2 % residual coupling.

**How you would notice.** Only in a careful fit where bead size and attraction are
being separated to better than a few percent.

**What prevents it.** Nothing — it is accepted. Removing it would require either a
discontinuous potential at `r = R` or a separate WCA term with its own σ. The
residual runs in the benign direction (more attraction slightly *softens* the core)
and is 18× smaller than the 58 % coupling it replaces.
