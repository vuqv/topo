# Nscale Optimization Specification

## Purpose

Calibrate the non-bonded interaction strengths of a TOPO coarse-grained Gō-like
model so that every folded unit (each domain, and each domain–domain interface)
is stable at the simulation temperature, while perturbing the underlying
energetics as little as possible.

The optimizer (`topo.optimize`, console `topo-optimize`) searches for the
**smallest per-unit `nscale`** — a multiplier on the sidechain–sidechain (SS)
energy channel — on a discrete per-class ladder that keeps the native state
folded across independent trajectories. It owns the search logic and drives the
package tools (`topo.mdrun`, chain splitting, native-contact scoring) one round
at a time.

Scope notes and known limitations:

- `nscale` scales **only** the SS channel of native-contact well depths. Well
  positions (native Cα–Cα distances) and the hydrogen-bond (HB) and
  backbone–sidechain (BS) channels are never modified by this optimizer.
- The search is a greedy per-unit climb; it treats units as independent and
  assumes raising a unit's `nscale` only ever stabilizes it. Coupled
  destabilization is possible and is handled by reporting, not by re-search.
- The optimization is **not resumable**; each invocation starts fresh.
- Fully-disordered proteins have no SS native contacts to scale and are outside
  the optimizer's scope (see the short-circuit requirement); their residues are
  governed by `idr_scale`, not `nscale`.

## Requirements

### Requirement: nscale scales only the sidechain–sidechain channel

The per-pair `nscale` SHALL multiply only the sidechain–sidechain (SS) component
of a native contact's well depth. The total native-contact well depth is
`eps_ij = HB + BS + nscale * SS`, where the HB (hydrogen-bond, from STRIDE) and
BS (backbone–sidechain, distance-based) channels are independent of `nscale`.
The well **position** (`rmin`) SHALL remain the measured native Cα–Cα distance
regardless of `nscale`.

#### Scenario: Raising nscale deepens only the SS term

- **WHEN** a domain's `nscale` is increased from one ladder level to the next
- **THEN** the SS contribution to every intra-domain native contact scales by
  the same factor
- **AND** that contact's HB and BS contributions are unchanged
- **AND** the contact's well position (native Cα–Cα distance) is unchanged

#### Scenario: A contact with no SS energy is unaffected by nscale

- **WHEN** a native contact is formed solely through hydrogen bonds and/or
  backbone–sidechain contacts (no sidechain–sidechain contact)
- **THEN** its well depth is independent of `nscale`

### Requirement: Per-unit scope of nscale

Each domain SHALL carry its own intra-domain `nscale`, applied to all
native-contact pairs internal to that domain. Each domain–domain interface SHALL
carry its own `nscale`, applied to native-contact pairs that span the two
domains. Domains and interfaces are optimized as independent scoring units.

#### Scenario: Domain vs interface scaling

- **WHEN** the domain definition declares domains D1 and D2
- **THEN** intra-D1 and intra-D2 native contacts use their respective domain
  `nscale` values
- **AND** native contacts between a D1 residue and a D2 residue use the D1–D2
  interface `nscale`

### Requirement: Unassigned residues form an unoptimized domain X

Residues that are neither listed in any declared domain nor marked disordered
SHALL be collected, in the energy model, into a single synthetic domain named
`X` whose intra-domain `nscale` is fixed at 1.0, with every `X`–domain interface
also fixed at 1.0. Domain `X` retains its Gō native contacts (unlike disordered
residues, whose native contacts are masked out); it is merely left unscaled.

Domain `X` SHALL NOT be a scoring unit in the optimizer: the optimizer builds
scoring units only from declared domains and their pairwise interfaces, so `X`
is never scored for stability, never climbs the ladder, and never appears in the
per-round `domain.yaml` the optimizer writes. Unassigned residues SHALL be
surfaced as a warning, not silently absorbed.

This means the optimizer cannot rescue an unassigned region: because `nscale`
1.0 lies below every ladder's level 1, any genuinely structured domain
accidentally left undeclared is under-stabilized and the optimizer will not
detect or correct it — `X` is intended for flexible linkers and loops that
should not be rigidly folded.

#### Scenario: Leftover residues become X at nscale 1.0

- **WHEN** the domain definition declares domains for a strict subset of the
  protein's residues, leaving some residues unassigned and not disordered
- **THEN** the energy model places those residues in domain `X` with intra
  `nscale` 1.0 and all `X`–domain interfaces at `nscale` 1.0
- **AND** those residues' Gō native contacts are retained (not masked)

#### Scenario: X is invisible to the optimizer

- **WHEN** the optimizer runs on a domain definition with unassigned residues
- **THEN** it creates no scoring unit for domain `X`
- **AND** it never scores, climbs, or writes an `nscale` for `X`
- **AND** it emits a warning that those residues are unassigned

### Requirement: Discrete per-class nscale ladder

Candidate `nscale` values SHALL be drawn from a fixed discrete ladder keyed by
structural class. Each class provides five ordered levels plus a median
fallback. The classes and their levels (level 1 → level 5, then fallback) are:

- `alpha`:      1.1954, 1.4704, 1.7453, 2.0322, 2.5044  (fallback 1.7453)
- `beta`:       1.4732, 1.8120, 2.1508, 2.5044, 2.5044  (fallback 2.1508)
- `alpha_beta`: 1.1556, 1.4213, 1.6871, 1.9644, 2.5044  (fallback 1.6871)
- `interface`:  1.2747, 1.5679, 1.8611, 2.1670, 2.5044  (fallback 1.8611)

A domain's class SHALL be taken from its declaration in the domain definition
(`alpha`/`beta`/`alpha-beta`, with the aliases `a`/`b`/`c`). Every interface
SHALL use the `interface` class.

#### Scenario: Class determines the ladder

- **WHEN** a domain is declared with class `beta`
- **THEN** its `nscale` at level 1 is 1.4732 and climbs through the `beta` ladder

#### Scenario: Level beyond the ladder resolves to the median fallback

- **WHEN** a unit's ladder level is at or beyond the fallback index (past
  level 5)
- **THEN** its `nscale` is that class's median fallback value

### Requirement: Stability-gated greedy climb

The optimizer SHALL run in rounds, starting every non-frozen unit at ladder
level 1. Each round it writes the current per-unit `nscale` values into a
`domain.yaml`, runs one multi-copy MD producing `ntraj` independent
trajectories, scores each unit's stability, then:

- freezes every unit judged **stable** at its current level, and
- advances every **unstable** unit up one ladder level.

The loop SHALL stop when all units are stable (converged), when every remaining
unstable unit has reached the median fallback (accepted as-is), or when
`max_rounds` is reached. `max_rounds` defaults to 6 = 5 ladder levels plus the
fallback.

#### Scenario: Unstable unit climbs, stable unit freezes

- **WHEN** a round finds domain D1 stable and interface D1–D2 unstable
- **THEN** D1 is frozen at its current `nscale`
- **AND** D1–D2 advances to the next ladder level for the following round

#### Scenario: Convergence

- **WHEN** every non-frozen unit is stable in a round
- **THEN** the optimizer stops and reports convergence at that round

#### Scenario: Fallback exhaustion

- **WHEN** the only remaining unstable units are already at the median fallback
- **THEN** the optimizer stops and accepts those units at the fallback `nscale`

#### Scenario: max_rounds without convergence

- **WHEN** `max_rounds` is reached with one or more units still unstable
- **THEN** the optimizer emits a WARNING naming the unstable units
- **AND** still writes the calibrated model with those units left at their
  highest level / median fallback

### Requirement: Q-based per-unit stability criterion

A scoring unit SHALL be judged **stable** in a round only if its fraction of
formed native contacts `Q` exceeds `q_threshold` for at least `frame_fraction`
of frames in **every** one of the `ntraj` trajectories. Defaults: `ntraj` = 10,
`q_threshold` = 0.6688, `frame_fraction` = 0.98. A native contact counts as
formed when the instantaneous Cα–Cα distance is within a tolerance factor
(default 1.2) of its native distance. A unit with no native contacts (`Q` = NaN)
is treated as not-applicable and does not block stability.

#### Scenario: One trajectory unfolds

- **WHEN** a unit meets the folded-fraction threshold in 9 of 10 trajectories
  but not the tenth
- **THEN** the unit is judged unstable and climbs the ladder

### Requirement: Native-contact reference is built once and reused

The per-domain and per-interface native-contact lists (pair indices and
reference Cα–Cα distances) SHALL be built once from the all-atom reference
structure and reused to score every round, so the `Q` reference never drifts
between rounds. Native contacts are heavy-atom contacts within 4.5 Å with
sequence separation `|i − j|` greater than 3.

#### Scenario: Same reference across rounds

- **WHEN** the optimizer runs multiple rounds
- **THEN** each round's `Q` is measured against the identical native-contact
  reference derived from the input structure

### Requirement: STRIDE is computed once per optimization

The hydrogen-bond assignment (STRIDE) depends only on the fixed reference
structure. When the user does not supply a precomputed STRIDE file, the
optimizer SHALL run STRIDE exactly once into the optimization root and point
every round's `md.ini` at that single file. A user-supplied
`stride_output_file` SHALL be used as-is and SHALL be resolved to an absolute
path.

#### Scenario: STRIDE reused by every round

- **WHEN** no `stride_output_file` is given
- **THEN** STRIDE runs once at startup and its output is reused unchanged by all
  rounds

### Requirement: Weakly structured units are frozen, not optimized

A scoring unit whose native-contact count is below `min_contacts` SHALL be
pinned at ladder level 1 and excluded from the climb, because it is too weakly
structured to reach the folding threshold. `min_contacts` defaults to 0, which
disables the check.

#### Scenario: Below the contact floor

- **WHEN** `min_contacts` is 20 and a domain has 12 native contacts
- **THEN** that domain is pinned at level 1 and never advances the ladder
- **AND** it is reported as frozen in the log and final report

### Requirement: Fully-disordered input short-circuits the optimizer

When the total native-contact count across all scoring units is zero — the
protein is fully disordered, or every foldable contact is masked by the
`disordered:` section — there is no SS energy to scale. The optimizer SHALL exit
cleanly, leave the input domain definition unchanged, and report that there was
nothing to optimize, rather than reporting a vacuous convergence.

#### Scenario: No foldable contacts

- **WHEN** every scoring unit has zero native contacts
- **THEN** the optimizer returns the unchanged input domain definition and logs
  that nothing was optimized

### Requirement: Calibrated model output

On completion the optimizer SHALL write `domain_optimized.yaml` under the
optimization root, preserving each domain's residues and class from the input
and recording the final per-domain and per-interface `nscale` values. Per round
it SHALL also write, next to each trajectory's DCD, a `Q_<k>.csv` time series
(one row per frame, one column per scoring unit) for inspection of the values
behind each stability decision.

#### Scenario: Final artifact

- **WHEN** the optimizer finishes (converged or not)
- **THEN** `domain_optimized.yaml` is written with the final `nscale` values
- **AND** each round has one `Q_<k>.csv` per trajectory recording per-frame `Q`
