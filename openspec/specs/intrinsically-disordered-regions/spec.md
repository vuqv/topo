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

The non-bonded energy model SHALL enforce three distinct pair classes with
different interaction rules: folded–folded, IDR–IDR, and IDR–folded.
This is implemented by `apply_disorder()`, which runs **after** the folded model
is fully built.

**Design rationale:** Folded proteins have structure; the K-B excluded-volume
radius is structure-derived and physically meaningful. Disordered regions have
no fixed structure, so structure-derived radii are meaningless. Using transferable
per-residue radii (O'Brien's σ-radii per amino acid type) decouples the radii
from a particular conformation, matching the disorder premise.

#### Scenario: Folded–folded pairs unchanged

- **WHEN** neither residue is marked disordered
- **THEN** the pair's well depth (`HB + BS + nscale * SS`) and well position
  (native Cα–Cα distance) are unchanged from the folded build
- **AND** the folded residues keep their Karanicolas-Brooks Rmin/2 radii

#### Scenario: IDR–IDR pairs with attraction

- **WHEN** both residues in a pair are marked disordered
- **THEN** the pair's well depth is `max(NON_NATIVE_KJ, eps_gen_kj + idr_scale * SS)`
  where `SS` is the raw sidechain–sidechain energy from the BT potential
  (NOT domain-scaled; see Relationship to nscale-Optimization)
- **AND** the well is positioned at the bare sum of the two residues' per-AA
  transferable radii (O'Brien convention)
- **AND** the energy floor is NON_NATIVE_KJ, preventing the chain from passing
  through itself even at `idr_scale = 0`

**Design rationale:** The two-parameter model separates concerns. `idr_scale`
tunes the sequence-specific well depth (BT potential); `eps_gen_kj` adds a
sequence-independent attractive baseline that represents generic hydrophobic
cohesion. The floor prevents unphysical overlap.

#### Scenario: IDR–folded pairs (excluded-volume only)

- **WHEN** one residue is marked disordered and the other is folded
- **THEN** the pair has no native contact (no attraction)
- **AND** the well depth is NON_NATIVE_KJ (repulsive, excluded-volume only)
- **AND** the well is positioned at the sum of the IDR residue's per-AA radius
  and the folded residue's K-B Rmin/2
- **AND** the folded residue retains its native K-B radius (not shrunk)

**Design rationale:** IDR residues should not stabilize the folded core. Lacking
sequence information about cross-interactions, we use repulsive excluded-volume
only, keeping the folded bead on its native radius to avoid artificial
destabilization of core contacts.

---

### Requirement: Radii: Structural vs. Transferable

IDR residues SHALL use transferable per-amino-acid excluded-volume radii from
the model-parameters file, while folded residues keep their Karanicolas-Brooks
structure-derived radii. The per-residue `rmin_2_nm` array SHALL be overridden
in place by `apply_disorder()`, so both the energy model and any downstream
consumers (e.g., CSP ribosome excluded-volume) see consistent radii.

**Design rationale:** The K-B radius is derived from the native structure's
Cα–Cα contact distances and is meaningless for arbitrary coordinates. Transferable
σ-radii per amino acid type (O'Brien's values) capture the excluded-volume
stiffness without requiring a fixed structure.

#### Scenario: IDR radii override

- **WHEN** `apply_disorder()` is called with a `dis_idx` array of disordered residues
- **THEN** each residue `i` in `dis_idx` has its `rmin_2_nm[i]` set to
  `IDR_RMIN_2_NM[resname[i]] / 2^(1/6)`, the per-AA transferable σ-radius
- **AND** folded residues' `rmin_2_nm` remain exactly as before

#### Scenario: CSP ribosome excluded-volume consistency

- **WHEN** CSP calls `precompute_contacts()` → `setParticlesRadii()` with the
  overridden `rmin_2_nm` array
- **THEN** the ribosome–NC excluded-volume layer sees the same IDR σ-radii
- **AND** NC–NC and NC–ribosome excluded-volume models are internally consistent

#### Scenario: No disordered residues → radii unchanged

- **WHEN** the `disordered:` section is absent or empty
- **THEN** `apply_disorder()` is not called (or called with an empty `dis_idx`)
- **AND** all residues keep their original radii

---

### Requirement: Two Attraction Parameters: idr_scale and eps_gen_kj

The energy model SHALL support two independent attraction tuning parameters for
IDR–IDR pairs: `idr_scale` (sequence-specific) and `eps_gen_kj` (sequence-independent).

**Design rationale:** Physical decoupling. `idr_scale` scales the BT potential,
which encodes amino-acid pair specificity. `eps_gen_kj` is a flat hydrophobic
baseline, independent of sequence, representing non-specific attractive interactions
(e.g., hydrophobic effect in an aqueous environment). Together they allow
independent tuning of compaction (excluded-volume stiffness via `idr_scale`) and
cohesion (baseline attraction via `eps_gen_kj`).

#### Scenario: idr_scale tunes sequence-specific attraction

- **WHEN** `idr_scale = 0`
- **THEN** IDR–IDR pairs have only the NON_NATIVE_KJ floor (pure self-avoiding chain)
- **WHEN** `idr_scale = 1.0` (default)
- **THEN** IDR–IDR pairs use the full BT sequence-specific well depth, plus eps_gen_kj
- **WHEN** `idr_scale > 1.0`
- **THEN** IDR–IDR wells are deepened, favoring compaction

**Note:** The name "idr_scale" is historical. Its effect on compaction is
counterintuitive: raising `idr_scale` *stiffens* excluded-volume, not softens it,
because the BT well is shallow (~0.2 kJ/mol at default; raising α strengthens
excluded-volume repulsion relative to attraction).

#### Scenario: eps_gen_kj provides sequence-independent cohesion

- **WHEN** `eps_gen_kj = 0` and `idr_scale = 1.0`
- **THEN** IDR–IDR pairs use only the sequence-specific BT term
- **WHEN** `eps_gen_kj = 2.25` (default)
- **THEN** IDR–IDR pairs gain a flat 2.25 kJ/mol baseline attractive energy
- **AND** this baseline is added before the `idr_scale * SS` term, so both
  contribute to the total well depth

**Calibration note:** Default values (`idr_scale = 1.0`, `eps_gen_kj = 2.25`)
were determined by comparing simulated Rg to SAXS measurements on 24 fully-disordered
proteins (see `docs/usage/disordered_regions.rst`). At these defaults, the RMS
fractional Rg error is ~33% with a slope of 0.71 (Pearson r ≈ 0.79).

---

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
