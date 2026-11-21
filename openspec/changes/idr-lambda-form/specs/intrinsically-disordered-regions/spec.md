## ADDED Requirements

### Requirement: Model Selection
The `disordered:` section SHALL accept an optional `model:` key naming the parameter set
that supplies λ and ε. The supported models SHALL be `bt`, `hps_urry` and `calvados2`; `bt` SHALL be
the default. The selected model SHALL fill λ for every IDR-involving pair by the two
residue types; folded–folded entries SHALL NOT be touched.

#### Scenario: Default model is BT
- **WHEN** a `disordered:` section has no `model:` key
- **THEN** the build SHALL use the `bt` model

#### Scenario: Cross pairs use the same table
- **WHEN** a pair has exactly one disordered residue
- **THEN** its λ SHALL be the model's table value for the two residue types
- **AND** SHALL equal the λ an IDR–IDR pair of the same two residue types receives

#### Scenario: Unknown model
- **WHEN** `model:` is none of `bt`, `hps_urry` and `calvados2`
- **THEN** the build SHALL fail with an error listing the supported models

### Requirement: Shipped Pair-Table Parameter Sets
Every supported model SHALL be a shipped data file declaring a symmetric residue-pair
table `λ_ij` for the 20 standard amino acids, rounded to 4 decimals, together with its ε
(kJ/mol) and a citation. A model file SHALL NOT define bead sizes.

- `bt`: `λ_ij = 0.10·ε_BT,ij / ε` with `ε = 0.8368 kJ/mol`, where `ε_BT,ij` is the shifted
  Betancourt–Thirumalai energy topo already uses.
- `hps_urry`: `λ_ij = (λ_i + λ_j)/2` from the published per-residue λ (Regy et al.
  2021), with `ε = 0.8368 kJ/mol`. Published non-positive values SHALL be kept.
- `calvados2`: `λ_ij = (λ_i + λ_j)/2` from the published CALVADOS 2 per-residue λ
  (Tesei & Lindorff-Larsen 2022), with `ε = 0.8368 kJ/mol`. Only λ SHALL be taken from
  CALVADOS 2; its σ, charges and bond length SHALL NOT be used.

#### Scenario: Complete, symmetric coverage
- **WHEN** a model file is loaded
- **THEN** it SHALL define `λ_ij` for all 20×20 standard residue pairs with
  `λ_ij = λ_ji`, or the load SHALL fail naming what is missing

#### Scenario: BT table matches its source
- **WHEN** the shipped `bt` table is compared with `0.10·ε_BT,ij / 0.8368` recomputed
  from the shipped BT potential
- **THEN** every entry SHALL agree after rounding to 4 decimals

#### Scenario: HPS-Urry table matches the publication
- **WHEN** the shipped `hps_urry` table is loaded
- **THEN** selected entries SHALL equal `(λ_i + λ_j)/2` computed by hand from the cited
  per-residue values (pinned in a test)

#### Scenario: CALVADOS 2 table matches the publication
- **WHEN** the shipped `calvados2` table is loaded
- **THEN** selected entries SHALL equal `(λ_i + λ_j)/2` computed by hand from the cited
  per-residue values (pinned in a test)
- **AND** every entry SHALL be positive, the smallest being GLU–GLU

### Requirement: Model-Independent Bead Sizes
Bead sizes SHALL NOT depend on `model:`. An IDR residue's `Rmin/2` SHALL be its
per-residue value from the model parameters, a folded residue's SHALL be its
structure-derived Karanicolas–Brooks value, and every IDR-involving pair SHALL use
`R_ij = Rmin/2_i + Rmin/2_j`, `σ_ij = R_ij / 2^(1/6)`.

The per-residue table value SHALL be `rvdw·2^(1/6)`, where `rvdw` is the transferable
van der Waals radius of the residue type.

#### Scenario: Radii are the same under every model
- **WHEN** the same structure and `disordered:` residues are built with each supported
  model
- **THEN** the per-residue `Rmin/2` array and the well-position matrix SHALL be identical

#### Scenario: Self-pair sigma is twice the van der Waals radius
- **WHEN** `σ = (Rmin/2_i + Rmin/2_i)/2^(1/6)` is computed for a disordered residue
- **THEN** it SHALL equal `2·rvdw` of that residue type

### Requirement: Removed Tuning Keys
`idr_scale` and `eps_ev_kj` SHALL NOT be accepted in the `disordered:` section. Their
values are fixed by the selected model.

#### Scenario: idr_scale is rejected
- **WHEN** a `disordered:` section contains `idr_scale`
- **THEN** the build SHALL fail with an error naming `idr_scale` and pointing to `model:`

#### Scenario: eps_ev_kj is rejected
- **WHEN** a `disordered:` section contains `eps_ev_kj`
- **THEN** the build SHALL fail with an error naming `eps_ev_kj` and pointing to `model:`

### Requirement: Co-Translational Consistency
In continuous synthesis (ribosome present), the λ matrix SHALL be built once for the
full chain and sliced per nascent length exactly as the well-position matrix is. The
IDR force SHALL NOT evaluate any pair involving a ribosome bead, and NC↔ribosome
interactions SHALL remain the O'Brien excluded-volume and Yukawa terms for every
`model:`.

#### Scenario: Refactor leaves CSP energies unchanged
- **WHEN** the λ-form refactor is applied with λ computed at full precision
- **THEN** every force energy (including NC↔ribosome excluded volume and Yukawa) at each
  pinned nascent length SHALL match the pre-change snapshot

#### Scenario: Ribosome beads are outside the IDR force
- **WHEN** ribosome bead positions are changed while nascent positions are fixed
- **THEN** the IDR force energy SHALL be unchanged

#### Scenario: NC-ribosome radii do not depend on the model
- **WHEN** a co-translational system is built with each supported model
- **THEN** the radii passed to the NC↔ribosome force SHALL be identical

## MODIFIED Requirements

### Requirement: Three Pair Classes in the Energy Model

Residue pairs SHALL be partitioned across **two** `CustomNonbondedForce` objects
with **disjoint** interaction groups:

- **folded–folded** (neither residue disordered) — the existing Karanicolas–Brooks
  12-10-6, restricted to `{folded} × {folded}`. Native contacts, energies and radii
  SHALL be unchanged from a build with no `disordered:` section.
- **IDR–IDR** and **IDR–folded** — a second force using the Ashbaugh–Hatch potential
  on a plain 12-6, on groups `{idr} × {idr}` and `{idr} × {folded}`.

The IDR force SHALL implement the λ form, with `σ_ij = 2^(−1/6)·R_ij`:

```
r ≤ R_ij :  U = 4ε[(σ_ij/r)¹² − (σ_ij/r)⁶] + (1 − λ_ij)·ε
r >  R_ij :  U = λ_ij · 4ε[(σ_ij/r)¹² − (σ_ij/r)⁶]
```

The force SHALL take a per-pair λ table, the per-pair well-position matrix `R_ij`
(nm) and a scalar ε (kJ/mol), and SHALL contain no knowledge of how λ was derived.

`λ_ij` SHALL be the selected model's table value for the two residue types, for
**every** IDR-involving pair — IDR–IDR and IDR–folded alike. The folded/disordered
status of the partner SHALL NOT enter it. ε SHALL be the selected model's ε.

A model with **no** `disordered:` section SHALL build only the 12-10-6 force, with
no interaction group, taking the pre-existing code path unchanged.

#### Scenario: Folded–folded pairs unchanged
- **WHEN** a model is built with a `disordered:` section
- **THEN** every folded–folded well depth SHALL equal that from an otherwise identical
  build without that section, to floating-point equality
- **AND** every folded residue's `Rmin/2` SHALL be computed from the folded residues
  alone (disordered residues are not candidates for the nearest non-contact partner),
  so it is greater than or equal to its folded-only value
- **AND** every folded–folded well position not involving a residue whose `Rmin/2`
  changed SHALL equal its folded-only value

#### Scenario: The two forces partition the pairs
- **WHEN** a model with disordered residues is built
- **THEN** no residue pair SHALL appear in the interaction groups of both forces
- **AND** the 12-10-6 force's groups SHALL contain no disordered bead

#### Scenario: IDR–IDR pairs with attraction
- **WHEN** both residues in a pair are marked disordered
- **THEN** the pair's λ SHALL equal the model's table value for the two residue types
- **AND** the well SHALL be positioned at the sum of the two residues' per-AA `Rmin/2`

#### Scenario: IDR–folded pairs
- **WHEN** a pair has exactly one disordered residue
- **THEN** its λ SHALL be the model's table value for the two residue types
- **AND** SHALL equal the λ an IDR–IDR pair of the same two residue types receives

#### Scenario: Minimum and depth
- **WHEN** a single pair is evaluated with `λ > 0`
- **THEN** the minimum SHALL lie at `r = R_ij` with `U = −λε`, and U SHALL be C¹ there

#### Scenario: Zero attraction is a clean excluded-volume bead
- **WHEN** a pair has `λ_ij = 0`
- **THEN** the IDR force SHALL be exactly zero for `r > R_ij`
- **AND** SHALL remain repulsive for `r < R_ij`

#### Scenario: Non-positive lambda has no well
- **WHEN** a single pair is evaluated with `λ ≤ 0`
- **THEN** `U` SHALL be non-negative and non-increasing in `r` everywhere

#### Scenario: The core is set by eps_ev_kj alone
- **WHEN** `λ_ij·ε` varies from 0 to 2 kJ/mol at the model's fixed ε
- **THEN** the radius at which `U = kT` SHALL change by less than 5 %

#### Scenario: Independent oracle agreement
- **WHEN** the IDR force energy is computed by OpenMM for the A.7 test system
- **THEN** it SHALL equal a numpy evaluation of the formula above from the λ and R
  matrices, to relative 1e-9

### Requirement: YAML Configuration Parsing

The domain-definition YAML file SHALL accept an optional top-level `disordered:`
section containing residue IDs and an optional `model:` key. This section SHALL be
parsed by `read_yaml_config()` in `topo.utils.nonbonded`.

#### Scenario: disordered section with residue list only

- **WHEN** the YAML contains a `disordered:` section giving only residue IDs (1-indexed)
- **THEN** those residues are marked disordered
- **AND** the model defaults to `bt`

#### Scenario: disordered section with explicit parameters

- **WHEN** the YAML contains `disordered: { residues: [42, 43, ...], model: hps_urry }`
- **THEN** those residues are marked disordered and the named model is used

#### Scenario: Missing disordered section

- **WHEN** the YAML omits `disordered:` entirely
- **THEN** no residues are marked disordered
- **AND** apply_disorder() is not called (energy model is unaffected)
- **AND** build_native_contacts() is called without a `disorder` argument (analysis unaffected)

#### Scenario: Empty disordered section

- **WHEN** the YAML contains `disordered: []` or `disordered: { residues: [] }`
- **THEN** the disordered state is empty (no residues marked)
- **AND** apply_disorder() and analysis masking proceed as normal (no-op case)

#### Scenario: Invalid parameter types

- **WHEN** the YAML gives a `model:` that is not a supported model name, or any of the
  removed keys `idr_scale`, `eps_ev_kj`, or a non-zero `eps_gen_kj`
- **THEN** `read_yaml_config()` SHALL raise ValueError with a clear message

### Requirement: Native Contact Masking in the Energy Path

The energy path (`build_nonbonded_interaction()` → `apply_disorder()`) SHALL
drop all native-contact entries for pairs involving a disordered residue.
IDR–IDR and IDR–folded pairs are rebuilt from the selected model's λ table, and
their native contacts are permanently removed.

#### Scenario: Native contacts dropped for IDR-involving pairs

- **WHEN** `apply_disorder()` is called
- **THEN** every `eps_ij[i, j]` entry where `i` or `j` is in `dis_idx` is overwritten
  with `λ_ij · ε` from the selected model, by the two residue types
- **AND** the same rule applies whether one or both residues are disordered
- **AND** the folded–folded block is byte-identical to a folded-only run

#### Scenario: Apply disorder is called last

- **WHEN** the full non-bonded build is complete (HB, BS, SS all summed;
  domain scaling applied)
- **THEN** `apply_disorder()` is invoked, overwriting only IDR-involving pairs
- **AND** no subsequent build step modifies those entries

### Requirement: Relationship to nscale-Optimization

IDR pairs SHALL never be scaled by `nscale` and are **never part** of the
folding-stability ladder. The optimizer (`topo.optimize`) SHALL NOT tune IDR
energetics: they are fixed by the selected model.

**Design rationale:** `nscale` tunes folded-domain stability. IDR residues have
no folded reference; the ladder assumption (raising nscale → more stable) is
invalid. IDR energetics come from a calibrated parameter set instead.

#### Scenario: IDR pairs use raw (unscaled) SS energy

- **WHEN** `apply_disorder()` fills IDR-involving well depths
- **THEN** each is `λ_ij · ε` from the model's residue-pair table
- **AND** it is **not** multiplied by any domain `nscale` factor
- **AND** for the `bt` model the table derives from the raw BT energy, not from the
  domain-scaled sidechain–sidechain energy

#### Scenario: Optimizer must see disorder active

- **WHEN** the optimizer calibrates `nscale` values
- **THEN** it must run with the `disordered:` section active (energy via
  `apply_disorder()` applied, Q via masked `build_native_contacts()`)
- **AND** if the optimizer is run on a folded-only version (disorder removed),
  it will calibrate `nscale` values that are too low, leading to under-stability
  at production (when disorder is re-added)
- **AND** if Q masking is omitted (Q counts IDR contacts), the optimizer will
  calibrate nscale values that are too high, over-stabilizing the folded core

**Recommendation:** Always optimize with `disordered:` active. If domain
boundaries must be refined, do so in the folded model first, then re-optimize
with disorder.

#### Scenario: Domain X and disordered residues coexist

- **WHEN** the domain definition declares domains for a subset, leaving some
  unassigned AND some residues are marked disordered
- **THEN** unassigned residues form domain X (pinned at nscale 1.0)
- **AND** disordered residues are separate (not part of X, governed by the model)
- **AND** both paths (energy + analysis) respect both sets (see nscale-optimization
  spec for domain X details)

### Requirement: Behavioral Signatures in Simulation Trajectories

Simulations of mixed folded–IDR proteins SHALL exhibit characteristic RMSF
(root-mean-square fluctuation) and Rg (radius of gyration) signatures that
distinguish the folded and disordered regions.

#### Scenario: Mixed folded–IDR dynamics

- **WHEN** a protein is simulated with a folded core and disordered termini
- **THEN** the folded-core residues show RMSF ~0.7 Å (low fluctuation, stable)
- **AND** the IDR residues show RMSF ~17–26 Å (high fluctuation, diffusive)
- **AND** the overall Rg fluctuates around a mean that depends on the selected model

**Note:** This is an observation-based requirement, not a hard constraint. It is
validated by comparing to experimental SAXS or NMR data.

## REMOVED Requirements

### Requirement: Excluded-Volume Strength Parameter

**Reason**: ε is a property of the selected model, calibrated together with its λ
table, and is no longer a user input.

**Migration**: remove `eps_ev_kj` from the `disordered:` section. Every shipped model
uses `0.8368` kJ/mol, the previous default. A YAML that still
contains the key fails with an error pointing to `model:`.

