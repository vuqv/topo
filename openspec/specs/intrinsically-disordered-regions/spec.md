# Intrinsically Disordered Regions Specification

## Purpose

Enable structure-based Gō-like modeling of intrinsically disordered regions (IDRs)
alongside folded domains in the same protein. IDRs are sequence-specific but
structurally undefined; their energetics must decouple excluded-volume (repulsive)
from attraction (sequence-specific) and use transferable per-residue radii instead
of structure-derived ones.

The implementation lives in three files: config parsing and energy model
(`topo.utils.nonbonded`), native-contact analysis masking (`topo.analysis.native_contacts`),
and optimizer scoring (`topo.optimize.optimize`). The user-facing interface is
a `disordered:` YAML section in the domain definition file, specifying which
residues are disordered and tuning the two attraction parameters (`idr_scale`,
`eps_gen_kj`).

Scope notes and known limitations:

- IDR pairs always use `nscale = 1.0` and are **not part** of the folding-stability
  ladder; they are governed by `idr_scale` (and `eps_gen_kj` if specified).
  See Relationship to nscale-Optimization.
- Disordered residues marked in the YAML override any declared domain membership
  (if a residue appears in both, disorder wins; see the overlap scenario).
- Single-chain only: the `disordered:` section contains residue IDs only, with
  no chain qualifier. Multi-chain disorder requires separate sections per chain
  (not yet implemented).
- IDR↔folded cross-attraction (`cross_scale`/O'Brien `interact_scale`) is
  deliberately omitted; it is documented as an extension point only.
## Requirements
### Requirement: YAML Configuration Parsing

The domain-definition YAML file SHALL accept an optional top-level `disordered:`
section containing residue IDs and optional tuning keys. This section SHALL be
parsed by `read_yaml_config()` in `topo.utils.nonbonded`.

#### Scenario: disordered section with residue list only

- **WHEN** the YAML contains `disordered: [42, 43, 44, ...]` (residue IDs, 1-indexed)
- **THEN** those residues are marked disordered
- **AND** `idr_scale` defaults to `DEFAULT_IDR_SCALE = 1.0`
- **AND** `eps_gen_kj` defaults to `DEFAULT_EPS_GEN_KJ = 2.25`

#### Scenario: disordered section with explicit parameters

- **WHEN** the YAML contains `disordered: { residues: [42, 43, ...], idr_scale: 0.5, eps_gen_kj: 1.0 }`
- **THEN** those residues are marked disordered with the specified parameter values
- **AND** missing keys use their defaults

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

- **WHEN** the YAML specifies `idr_scale: "foo"` or `eps_gen_kj: -1.0`
- **THEN** `read_yaml_config()` SHALL raise ValueError with a clear message

---

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

Well depths SHALL be `eps_ij = idr_scale · eps_BT(i,j)` for
**every** IDR-involving pair — IDR–IDR and IDR–folded alike. The depth is a function
of the two residue types only; the folded/disordered status of the partner SHALL NOT
enter it. A negative `idr_scale` SHALL be rejected.

Note that `idr_scale` is calibrated on fully-IDP chains, in which no IDR–folded pair
occurs, so the cross-pair depth is an extrapolation of that calibration.

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
- **THEN** the pair's well depth SHALL equal `idr_scale · eps_BT(i,j)`, where `eps_BT`
  is the raw sidechain–sidechain BT energy (not domain-scaled)
- **AND** the well SHALL be positioned at the sum of the two residues' per-AA `Rmin/2`

#### Scenario: IDR–folded pairs
- **WHEN** a pair has exactly one disordered residue and `idr_scale > 0`
- **THEN** its well depth SHALL equal `idr_scale · eps_BT(i,j)`
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

#### Scenario: IDR radii override
- **WHEN** a residue is marked disordered
- **THEN** its `Rmin/2` SHALL equal the parameter-table value for its residue type
- **AND** SHALL NOT be divided by `2^(1/6)`

#### Scenario: CSP ribosome excluded-volume consistency
- **WHEN** CSP passes the overridden per-residue radius array to the ribosome–NC
  excluded-volume layer
- **THEN** that layer SHALL see the same IDR radii as the intra-chain model

#### Scenario: No disordered residues → radii unchanged
- **WHEN** the `disordered:` section is absent or empty
- **THEN** all residues SHALL keep their original radii

#### Scenario: Sigma matches published HPS values
- **WHEN** `σ = (Rmin/2_i + Rmin/2_j)/2^(1/6)` is computed for `i = j`
- **THEN** it SHALL agree with the published HPS per-residue σ to within 2 % for all
  20 amino acids

### Requirement: Native Contact Masking in the Energy Path

The energy path (`build_nonbonded_interaction()` → `apply_disorder()`) SHALL
drop all native-contact entries for pairs involving a disordered residue.
IDR–IDR and IDR–folded pairs are rebuilt as non-native (excluded-volume only or
with generic cohesion), and their native contacts are permanently removed.

#### Scenario: Native contacts dropped for IDR-involving pairs

- **WHEN** `apply_disorder()` is called
- **THEN** every `eps_ij[i, j]` entry where `i` or `j` is in `dis_idx` is overwritten
- **AND** if both are disordered: set to `max(NON_NATIVE_KJ, eps_gen_kj + idr_scale * SS)`
- **AND** if one is disordered: set to NON_NATIVE_KJ
- **AND** the folded–folded block is byte-identical to a folded-only run

#### Scenario: Apply disorder is called last

- **WHEN** the full non-bonded build is complete (HB, BS, SS all summed;
  domain scaling applied)
- **THEN** `apply_disorder()` is invoked, overwriting only IDR-involving pairs
- **AND** no subsequent build step modifies those entries

---

### Requirement: Native Contact Masking in the Analysis Path

The analysis path (`build_native_contacts()` in `topo.analysis.native_contacts`)
SHALL mask IDR-involving pairs from the native-contact count and Q metric,
separate from the energy masking. This ensures that Q values and stability
measures correctly exclude interactions that the simulation does not form
(disordered residues do not have native conformations).

**Design rationale:** The Q metric counts native contacts that remain in the
trajectory, compared to the reference structure. IDR residues have no reference
native state in the input structure; their contacts (if any) are coincidental
and do not reflect the disorder model. Including them inflates Q (denominator
is too small) and can make the folded core appear more unstable than it is.

#### Scenario: Disorder masking in Q calculation

- **WHEN** `build_native_contacts()` is called with `disorder=np.array([42, 43, ...])`
- **THEN** pairs involving any residue in the disorder array are excluded from
  the native-contact list
- **AND** the Q numerator (observed native contacts) and denominator (total native
  contacts) both exclude IDR-involving pairs
- **AND** effective domain membership becomes `declared_domain - disordered_residues`

#### Scenario: Optimizer threads disorder into Q calculation

- **WHEN** the optimizer runs with an active `disordered:` section
- **THEN** it passes the disorder array to `build_native_contacts(..., disorder=...)`
- **AND** Q values reflect folded-domain stability only

#### Scenario: No disorder → masking disabled

- **WHEN** no `disordered:` section is present in the domain definition
- **THEN** `build_native_contacts()` is called with `disorder=None` (default)
- **AND** all native contacts are included (no masking)

---

### Requirement: Relationship to nscale-Optimization

IDR pairs SHALL always use `nscale = 1.0` and are **never part** of the
folding-stability ladder. The optimizer (`topo.optimize`) SHALL handle IDR residues
separately from folded domains, using `idr_scale` and `eps_gen_kj` to tune
disorder energetics independently of the folding ladder.

**Design rationale:** `nscale` tunes folded-domain stability. IDR residues have
no folded reference; the ladder assumption (raising nscale → more stable) is
invalid. Instead, IDR residues are governed by `idr_scale` (and `eps_gen_kj`),
which tune the disorder model itself, not the climbing ladder.

#### Scenario: IDR pairs use raw (unscaled) SS energy

- **WHEN** `apply_disorder()` computes IDR–IDR well depths as
  `eps_gen_kj + idr_scale * SS`
- **THEN** `SS` is the raw sidechain–sidechain energy from the BT potential
- **AND** `SS` is **not** multiplied by any domain `nscale` factor
- **AND** the raw energy is accessed directly from `ss_interaction_energy`, not
  from `scaled_ss_interaction_energy`

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
- **AND** disordered residues are separate (not part of X, governed by idr_scale)
- **AND** both paths (energy + analysis) respect both sets (see nscale-optimization
  spec for domain X details)

---

### Requirement: Fully-Disordered Short-Circuit

When a protein is fully disordered (every residue is marked disordered, or zero
SS native contacts would result), the system SHALL recognize this and handle it
gracefully.

#### Scenario: Fully disordered input to optimizer

- **WHEN** all residues are marked disordered (or every SS contact is masked)
- **THEN** the optimizer detects zero native SS contacts to optimize
- **AND** it exits cleanly without searching (no ladder stages)
- **AND** it reports in the log that there was nothing to optimize
- **AND** the domain definition is left unchanged

#### Scenario: Fully disordered MD simulation

- **WHEN** the protein is fully disordered and MD is run (no optimization)
- **THEN** the simulation executes normally with IDR–IDR energies only
- **AND** the trajectory shows only disorder dynamics (no folded-domain motion)

---

### Requirement: Overlapping Domains and Disorder

If a residue is listed in both a declared domain and the `disordered:` section,
the disorder marking SHALL take precedence.

**Design rationale:** Disorder is a modeling choice that overrides structural
classification. A residue can be listed as part of a domain (e.g., from an
earlier fold prediction) but then marked disordered (e.g., from NMR or AlphaFold
confidence scores) in the same YAML.

#### Scenario: Residue in domain AND disordered

- **WHEN** residue 42 appears in `intra_domains: { D1: [40, 41, 42, ...] }`
  AND in `disordered: [42, 43, ...]`
- **THEN** residue 42 is treated as disordered
- **AND** in the energy path, residue 42 participates only in IDR–IDR and
  IDR–folded pairs
- **AND** in the analysis path, residue 42 is excluded from domain D1's native
  contacts (effective domain membership is D1 − {42})
- **AND** an info-level log line is emitted: "Residue 42 in both domain D1 and
  disordered; disorder takes precedence"

---

### Requirement: CSP Ribosome Integration

When the CSP (continuous synthesis pathway) is used with disordered residues,
the ribosome–nascent-chain excluded volume SHALL use the same per-residue radii
as the energy model. This ensures that a disordered nascent chain interacts
consistently with the ribosome.

**Design rationale:** The `rmin_2_nm` array is the single source of per-residue
radii. By overriding it in place, both energy and CSP consumers see the same
radii without duplicating the logic.

#### Scenario: Disordered nascent chain in CSP

- **WHEN** CSP calls `precompute_contacts()` with an active `disordered:` section
- **THEN** the `rmin_2_nm` array passed to `setParticlesRadii()` has already been
  overridden by `apply_disorder()`
- **AND** nascent-chain–ribosome excluded-volume interactions use the per-AA
  transferable radii for disordered residues
- **AND** the L24 free-loop (if used) is unaffected (it has no `disordered:` section)

---

### Requirement: Behavioral Signatures in Simulation Trajectories

Simulations of mixed folded–IDR proteins SHALL exhibit characteristic RMSF
(root-mean-square fluctuation) and Rg (radius of gyration) signatures that
distinguish the folded and disordered regions.

#### Scenario: Mixed folded–IDR dynamics

- **WHEN** a protein is simulated with a folded core and disordered termini
- **THEN** the folded-core residues show RMSF ~0.7 Å (low fluctuation, stable)
- **AND** the IDR residues show RMSF ~17–26 Å (high fluctuation, diffusive)
- **AND** the overall Rg fluctuates around an idr_scale-dependent mean
  (see Tune_idr_scale calibration)

**Note:** This is an observation-based requirement, not a hard constraint. It is
validated by comparing to experimental SAXS or NMR data.

---

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

## Scenarios for Testing

These scenarios correspond to test cases in `tests/test_idr.py`:

| Test ID | Scenario | Coverage |
|---------|----------|----------|
| 1 | Baseline: folded-only, no disordered section | No-op guarantee |
| 2 | Single IDR residue; IDR–IDR and IDR–folded pairs | Pair class masking |
| 3 | Multiple IDR regions; cross-domain interference | Overlaps + masking |
| 4 | Behavioral MD with mixed folded–IDR | Trajectory dynamics (manual) |
| 5 | Fully disordered input; short-circuit check | Edge case: zero SS contacts |
| 6 | idr_scale variations (0.0, 0.5, 1.0, 2.0) | Attraction tuning |
| 7 | eps_gen_kj variations | Generic cohesion tuning |
| 8 | Disorder + domain X coexistence | Multiple unassigned classes |
| 9 | Q masking: effective domain membership | Analysis path correctness |
| 10 | CSP integration: radii consistency | Ribosome excluded volume |
| 11 | Optimizer with disorder active vs. folded-only | Ladder calibration difference |

---

## Configuration Examples

### Minimal: disordered residue list only

```yaml
intra_domains:
  folded_core: [1, 2, 3, 50, 51, 52]

disordered: [10, 11, 12, 13, 14]
```

Uses defaults: `idr_scale = 1.0`, `eps_gen_kj = 2.25`.

### Full: explicit parameters

```yaml
intra_domains:
  folded_core: [1, 2, 3, 50, 51, 52]
  flexible_hinge: [25, 26, 27]

disordered:
  residues: [10, 11, 12, 13, 14]
  idr_scale: 0.7
  eps_gen_kj: 1.5
```

Tunes compaction (lower `idr_scale`) and reduces baseline cohesion.

### Edge case: fully disordered

```yaml
disordered: [1, 2, 3, ..., N]
```

No `intra_domains` declared; every residue is disordered. Optimizer will
short-circuit; MD will simulate disorder-only dynamics.

---

## References

- Implementation: `topo/utils/nonbonded.py` (config parsing, `apply_disorder()`)
- Analysis: `topo/analysis/native_contacts.py` (disorder masking in Q)
- Optimizer: `topo/optimize/optimize.py` (Scorer threads disorder)
- Tests: `tests/test_idr.py`, `tests/test_idr_baseline.py`
- User guide: `docs/usage/disordered_regions.rst`
- Design record: `review/DISORDER_IDR_SPEC.md` (closed design questions)
- CSP integration: `topo/csp/core.py` (ribosome excluded volume)

---

## Related Specifications

- **nscale-Optimization**: IDR pairs always use `nscale = 1.0` and do not climb
  the folding ladder. The optimizer must see disorder active to calibrate folded
  `nscale` correctly. See "Relationship to nscale-Optimization" above.
- **Isolated-Protein-Simulation**: The energy-model build (including disorder
  masking) is a dependency of the MD runner. See the runner's "Dependency:
  non-bonded force field" requirement.
