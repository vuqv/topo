# Tasks

Branch from `dev`. Each numbered group is a self-contained commit with the suite
green at the end of it.

## 1. Radius convention — `Rmin/2`, not a σ-radius

- [ ] 1.1 In `apply_disorder`, stop dividing the per-AA table value by
      `RMIN_SCALE_FACTOR`; IDR beads take `IDR_RMIN_2_NM[aa]` directly.
      Keep `RMIN_SCALE_FACTOR` — `calculate_rmin_2_values` still needs it, where the
      σ→Rmin conversion is correct.
- [ ] 1.2 Update the assertion in the IDR radius test and drop the now-unused import.
- **Verify**: before editing the test, the suite shows **exactly one** failure, in
  the IDR radius assertion, and the folded-only baseline test stays green — proving
  the change touched only the IDR path. Then: full suite green; IDR radii lie inside
  the folded range rather than below every folded bead.

## 2. The IDR force

- [ ] 2.1 Add `DEFAULT_EPS_EV_KJ = 0.8368`; parse optional `eps_ev_kj` in the
      `disordered:` section.
- [ ] 2.2 `build_nonbonded_interaction` gains `return_idr=False`, returning
      `{'idx', 'eps_ev_kj'}` or `None`. **Keep it optional** so the three-value
      unpacks in `csp/core.py` and `csp/ribosome.py` continue to work.
- [ ] 2.3 `addCustomNonBondedForce(..., interaction_group=None)` — the 12-10-6 takes
      an opaque domain and knows nothing about disorder.
- [ ] 2.4 New sibling `addIDRNonBondedForce(..., idr_idx, folded_idx, eps_ev_kj)`
      building the Ashbaugh–Hatch 12-6, setting **its own** groups, and calling the
      same `createExclusionsFromBonds` so exclusion lists match.
- [ ] 2.5 Register the new force in `addSystemForces` under its own force group.
- [ ] 2.6 The caller computes the partition once and calls both constructors; with no
      disorder it calls only the first, with no interaction group.
- **Verify**: extract the energy expression **from the shipped source** and evaluate
  it in OpenMM — minimum at `R`, depth `−eps_ij`, `U ≡ 0` beyond `R` at
  `eps_ij = 0`, core moving < 5 % as `eps_ij` varies. Folded-only build bit-identical.

## 3. Remove `eps_gen_kj`, set the calibrated defaults

**Must ship together with task 2.** The force alone, with the previous
default of `idr_scale = 1.0`, puts a bare `disordered:` block ~3x past the
theta point -- a silently collapsed chain.

- [ ] 3.1 Delete the parameter, its default, and the additive term in `apply_disorder`.
- [ ] 3.2 Accept a legacy `eps_gen_kj: 0.0`; raise on any non-zero value with the
      `idr_scale` equivalent in the message.
- [ ] 3.3 `DEFAULT_IDR_SCALE = 0.10`; rewrite the calibration comment with the
      benchmark numbers from `proposal.md`.
- [ ] 3.4 Remove the `eps_gen` test and update the defaults test.
- **Verify**: bare `disordered:` block yields `idr_scale 0.10 / eps_ev_kj 0.8368`;
  `eps_gen_kj: 0.0` accepted; `eps_gen_kj: 2.25` refused with a usable message.

## 4. Honour group ownership downstream

- [ ] 4.1 `csp.ribosome.append_ribosome`: set `{nascent}×{nascent}` on the contact
      force **only if it has no groups**; pad the IDR force with the same `M` dummy
      particles the contact force gets.
- [ ] 4.2 `utils.multichain.replicate_system_intra_only`: replicate the template's
      interaction groups with a per-copy offset; fall back to a blanket
      `{copy}×{copy}` **only if the template had none**.
- **Verify**: N non-interacting copies have exactly N × the single-copy energy — the
  identity that any double-count breaks. Confirm each fix by reinstating the old line
  and watching the corresponding test fail.

## 5. CSP support

- [ ] 5.1 `precompute_contacts` returns the IDR handle; update the two unpack sites
      and the return annotation.
- [ ] 5.2 `build_length_model` accepts and forwards it.
- [ ] 5.3 `run_length` slices the mask to the emerged chain (`idx[idx < L]`); an empty
      slice yields `None` and the single-force path. Thread `idr_full` through the
      `run_length` call sites in `protocol.py` and `cylinder.py`.
- [ ] 5.4 Update the CSP test stubs that fake `precompute_contacts` to the new arity.
- **Verify**: a real CSP run with an **N-terminal** IDR — chosen so synthesis passes
  through fully-disordered, mixed, and ejection regimes — completes with finite
  energies and no NaN.

## 6. Tests

- [ ] 6.1 Potential shape, from the shipped expression (task 2 verify).
- [ ] 6.2 Force partition: no pair in both forces; the 12-10-6 sees no IDR bead;
      identical exclusion lists; the no-disorder path unchanged.
- [ ] 6.3 CSP: partition survives `append_ribosome`, in both the fully-disordered and
      mixed regimes.
- [ ] 6.4 Multichain: the N × energy identity, group disjointness, no cross-copy leak,
      and the no-IDR path unchanged.
- [ ] 6.5 A pinned aggregate regression for a fully-IDR build, in the style of the
      existing folded baseline test.
- **Verify**: every new test fails when its target fix is reverted. A test that
  passes both before and after is worthless.

## 7. Docs

- [ ] 7.1 `disordered_regions.rst`: rewrite the physics section for the two-force
      architecture, the AH form, why the desolvation barrier is dropped, and the
      `Rmin/2` convention. Replace the `eps_gen_kj` config rows with `eps_ev_kj`.
- [ ] 7.2 Rewrite the validation section: report the **18/6 split** and why the 6
      foldable globular proteins are excluded; give ν / R0 / RMS against experiment;
      and state the RMS against the right bar — the power law is fitted *to* the data
      and leaves a 9.5 % residual, so 12.0 % does not yet beat chain length alone.
- [ ] 7.3 `model_theory.rst` pair table and `optimization_control.rst` config example.
- [ ] 7.4 Record that a disordered region deletes cross-boundary native contacts, and
      that CSP's nascent↔ribosome force stays 12-10-6 for every bead by design.
- **Verify**: `cd docs && make html` succeeds; no stale `eps_gen` or `r_vdw`
  references outside the deliberate migration note.
