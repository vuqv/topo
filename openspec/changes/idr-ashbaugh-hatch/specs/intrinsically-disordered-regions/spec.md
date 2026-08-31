## MODIFIED Requirements

### Requirement: Three Pair Classes in the Energy Model

Residue pairs SHALL be partitioned across **two** `CustomNonbondedForce` objects
with **disjoint** interaction groups:

- **folded–folded** (neither residue disordered) — the existing Karanicolas–Brooks
  12-10-6, restricted to `{folded} × {folded}`. Native contacts, energies and radii
  SHALL be unchanged from a build with no `disordered:` section.
- **IDR–IDR** and **IDR–folded** — a second force using an Ashbaugh–Hatch split on a
  plain 12-6, on groups `{idr} × {idr}` and `{idr} × {folded}`.

The IDR force SHALL implement, with `L(r) = 4[(σ/r)¹² − (σ/r)⁶]` and
`σ = 2^(−1/6)·R_ij`:

```
r ≤ R_ij :  U = eps_EV·L(r) + (eps_EV − eps_ij)
r >  R_ij :  U = eps_ij·L(r)
```

The potential SHALL have its minimum at `r = R_ij` with depth exactly `−eps_ij` for
any `eps_EV`, and SHALL be C¹ at `r = R_ij`.

Well depths SHALL be `eps_ij = max(NON_NATIVE_KJ, idr_scale · eps_BT(i,j))` for
**every** IDR-involving pair — IDR–IDR and IDR–folded alike. The depth is a function
of the two residue types only; the folded/disordered status of the partner SHALL NOT
enter it. IDR–folded pairs SHALL NOT be restricted to the excluded-volume floor.

Note that `idr_scale` is calibrated on fully-IDP chains, in which no IDR–folded pair
occurs, so the cross-pair depth is an extrapolation of that calibration.

A model with **no** `disordered:` section SHALL build only the 12-10-6 force, with
no interaction group, taking the pre-existing code path unchanged.

#### Scenario: Folded–folded pairs are bit-identical
- **WHEN** a model is built with a `disordered:` section
- **THEN** every folded–folded well position and well depth SHALL equal those from an
  otherwise identical build without that section, to floating-point equality
- **AND** every folded residue's `Rmin/2` SHALL be unchanged

#### Scenario: The two forces partition the pairs
- **WHEN** a model with disordered residues is built
- **THEN** no residue pair SHALL appear in the interaction groups of both forces
- **AND** the 12-10-6 force's groups SHALL contain no disordered bead

#### Scenario: A cross pair is attracted like an IDR–IDR pair
- **WHEN** a pair has exactly one disordered residue and `idr_scale > 0`
- **THEN** its well depth SHALL equal `max(NON_NATIVE_KJ, idr_scale · eps_BT(i,j))`
- **AND** SHALL equal the depth an IDR–IDR pair of the same two residue types receives

#### Scenario: Zero attraction is a clean excluded-volume bead
- **WHEN** `idr_scale = 0`
- **THEN** the IDR force SHALL be exactly zero for `r > R_ij`
- **AND** SHALL remain repulsive for `r < R_ij`

#### Scenario: The core is set by eps_ev_kj alone
- **WHEN** `eps_ij` varies from 0 to 2 kJ/mol at fixed `eps_ev_kj`
- **THEN** the radius at which `U = kT` SHALL change by less than 5 %

### Requirement: Radii: Structural vs. Transferable

Disordered residues SHALL take the transferable per-amino-acid `Rmin/2` from the
parameter table directly. Folded residues SHALL keep their structure-derived
Karanicolas–Brooks `Rmin/2`. Both populations therefore carry the **same quantity**,
and every pair SHALL combine by the plain sum rule
`R_ij = Rmin/2_i + Rmin/2_j`, including IDR–folded cross pairs.

No σ-conversion (division by `2^(1/6)`) SHALL be applied when populating the
per-residue radius array. The conversion belongs inside the IDR force's own
expression, where `σ = 2^(−1/6)·R_ij` recovers the van der Waals sum.

#### Scenario: IDR radii are table Rmin/2 values
- **WHEN** a residue is marked disordered
- **THEN** its `Rmin/2` SHALL equal the parameter-table value for its residue type
- **AND** SHALL NOT be divided by `2^(1/6)`

#### Scenario: Sigma matches published HPS values
- **WHEN** `σ = (Rmin/2_i + Rmin/2_j)/2^(1/6)` is computed for `i = j`
- **THEN** it SHALL agree with the published HPS per-residue σ to within 2 % for all
  20 amino acids

## ADDED Requirements

### Requirement: Interaction-Group Ownership

A force's interaction groups SHALL be set exactly once, by the code that creates the
force. No other code SHALL add interaction groups to a force it did not create.

OpenMM takes the **union** of interaction groups, so a later broader assertion
silently re-admits pairs the force no longer owns and causes them to be evaluated by
two forces at once, with no error raised.

Code that must restrict a force it did not create SHALL do so **only when that force
has no interaction groups**, i.e. only when no other component has claimed ownership.

#### Scenario: Ribosome append does not widen the contact force
- **WHEN** a ribosome is appended to a nascent model that has disordered residues
- **THEN** the 12-10-6 force's interaction groups SHALL still contain no disordered bead
- **AND** no pair SHALL be present in both nascent-side forces

#### Scenario: Ribosome append still restricts an unclaimed force
- **WHEN** a ribosome is appended to a model with no disordered residues
- **THEN** the 12-10-6 force SHALL be restricted to `{nascent} × {nascent}`

#### Scenario: Replication preserves the partition
- **WHEN** a model with disordered residues is replicated into N non-interacting copies
- **THEN** each force's template groups SHALL be reproduced per copy with an index offset
- **AND** no pair SHALL be present in both forces
- **AND** no interaction group SHALL span two copies

#### Scenario: Replicated energy is exactly N times one copy
- **WHEN** a model is replicated into N non-interacting copies
- **THEN** the potential energy SHALL equal N × the single-copy energy to within
  floating-point tolerance

### Requirement: Excluded-Volume Strength Parameter

The `disordered:` section SHALL accept an optional `eps_ev_kj` (kJ/mol) setting the
repulsive-core strength of the IDR force, independent of the well depth. It SHALL
default to `0.8368` (0.2 kcal/mol, the HPS value).

#### Scenario: Default applied when omitted
- **WHEN** a `disordered:` section omits `eps_ev_kj`
- **THEN** the IDR force SHALL use `0.8368` kJ/mol

#### Scenario: Explicit value honoured
- **WHEN** `eps_ev_kj` is given
- **THEN** the IDR force's `epsEV` global parameter SHALL equal it

### Requirement: Continuous Synthesis with a Disordered Nascent Chain

The CSP path SHALL support a `disordered:` section in the nascent chain's domain
definition. The disorder mask SHALL be sliced to the emerged chain at each length,
using the existing correspondence that particle `i` (`0..L-1`) is native residue
`i+1` and the contact subset is the top-left `L × L` block.

The nascent↔ribosome interaction SHALL remain a 12-10-6 with `RIBO_NC_EPS_KJ` for
**every** nascent bead, disordered or not. It is a separate force and SHALL NOT be
switched by disorder.

#### Scenario: Fully disordered emerged chain
- **WHEN** the emerged chain length is within the disordered prefix
- **THEN** the 12-10-6 force's interaction group SHALL be empty and contribute zero
- **AND** the build SHALL succeed and produce finite energies

#### Scenario: Mixed chain
- **WHEN** the emerged chain contains both folded and disordered residues
- **THEN** the IDR force SHALL carry both `{idr}×{idr}` and `{idr}×{folded}` groups

#### Scenario: Every force sees every particle
- **WHEN** the ribosome is appended
- **THEN** every `CustomNonbondedForce` SHALL hold one parameter entry per system
  particle, and all SHALL have identical exclusion lists

## REMOVED Requirements

### Requirement: Two Attraction Parameters: idr_scale and eps_gen_kj

**Reason**: `eps_gen_kj` is degenerate with `idr_scale` — both feed one well depth
`eps_gen + idr_scale·eps_BT`, so any value of one trades for the other. Its only
non-degenerate effect is to raise the mean well depth without its spread, diluting
sequence contrast, which is already the model's weakest property.

**Migration**: use `idr_scale` alone. An equivalent mean depth is
`idr_scale + eps_gen_kj/⟨eps_BT⟩` with `⟨eps_BT⟩ ≈ 2.36 kJ/mol`. A config carrying
`eps_gen_kj: 0.0` SHALL be accepted unchanged (it was a no-op); a non-zero value
SHALL raise an error naming the `idr_scale` equivalent rather than being ignored.

The calibrated default for `idr_scale` becomes `0.10` (was `1.0`); the removed
`eps_gen_kj` default was `2.25`.
