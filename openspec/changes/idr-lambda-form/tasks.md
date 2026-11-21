> All design questions are resolved. Start with the prerequisites in group 1.

## 1. Prerequisites

- [x] 1.1 Archive `idr-ashbaugh-hatch`; convert this change's deltas (REMOVED "Excluded-Volume Strength Parameter", MODIFIED "YAML Configuration Parsing" and "Three Pair Classes in the Energy Model")
- [x] 1.2 Confirm `tests/test_idr_ah_verification.py` is green on the starting commit (baseline)
- [x] 1.3 Extend the bundle with a co-translational case (toy ribosome, A.7 chain at two nascent lengths: before and after the first IDR emerges); pin every force energy incl. NC↔ribosome EV and Yukawa; commit before any refactor
- [x] 1.4 Sweep the remaining `idr_scale` / `eps_gen_kj` mentions in the main spec's requirements (native-contact masking, nscale relationship, behavioural signatures) into MODIFIED deltas
- [ ] 1.5 At archive time, update the main spec's non-requirement sections by hand (Purpose, the test table, Configuration Examples): deltas cannot reach them and they still show `idr_scale` / `eps_gen_kj`

## 2. λ-form refactor (D9 step 1 — no user-visible change, snapshot unchanged)

- [x] 2.1 Rewrite `addIDRNonBondedForce` to take `(lambda_matrix, rmin_matrix, epsilon)` with the `select()` expression (D1, D4); σ conversion stays inside the expression (D3)
- [x] 2.2 Compute `λ_ij = eps_ij/ε` in memory at full precision for IDR-involving pairs (temporary; replaced in group 3)
- [x] 2.3 Update callers and the per-length disorder handle (`models.py`, `csp/core.py`, `csp/cylinder.py`)
- [x] 2.4 Update `test_idr_force.py::test_idr_potential_*` to the new expression/globals
- [x] 2.5 Verify: full suite green AND `test_idr_ah_verification.py` (isolated and co-translational cases) green **without** regenerating the snapshot
- [x] 2.6 CSP: slice the λ matrix per nascent length exactly like `R_full`/`eps_full`; keep ribosome placeholder entries outside the IDR force's groups
- [x] 2.7 Resume: check `tests/csp/test_resume*.py` against the renamed global; keep the name stable or document the incompatibility for in-flight runs

## 3. Shipped `bt` table and YAML (D9 step 2 — BREAKING)

- [x] 3.1 Write the generator script that builds `topo/parameters/data/idr/bt.yaml` from `bt_potential.csv` (`λ_ij = 0.10·|BT_ij − 0.6|·4.184/0.8368`, 4 decimals)
- [x] 3.2 Add the loader: reads a model file, validates full 20×20 coverage and symmetry
- [x] 3.3 `apply_disorder` / `build_nonbonded_interaction` fill IDR-involving pairs from the model's λ table and return the model's ε; remove `DEFAULT_IDR_SCALE` / `DEFAULT_EPS_EV_KJ` from the run-time path
- [x] 3.4 YAML reader: accept `model:` (default `bt`); raise on `idr_scale`, on `eps_ev_kj`, and on an unknown model
- [x] 3.5 Rewrite the tests that set the removed keys (`test_idr.py`, `test_idr_force.py`, `test_idr_ah_verification.py`); drop the `idr_scale = 0` and custom-`eps_ev_kj` cases
- [x] 3.6 Regenerate the verification snapshot in a commit that changes nothing else; record the largest energy and force difference against the old snapshot

## 4. `hps_urry`

- [x] 4.1 Copy the 20 per-residue λ (`hps` field of `hps_urry` in cosmo's `cosmo/parameters/model_parameters.py`, Regy et al. 2021, shift already applied) into the generator script with the citation; keep `GLU = −0.08`
- [x] 4.2 Generate `topo/parameters/data/idr/hps_urry.yaml` (`λ_ij = (λ_i + λ_j)/2`, 4 decimals, `ε = 0.8368`)
- [x] 4.3 Pin a few hand-computed `λ_ij` in a test, including GLU–GLU = −0.08 and its purely repulsive pair energy

## 5. Tests

- [x] 5.1 Shipped tables equal their recomputation from source (both models)
- [x] 5.2 Extend the numpy oracle to take the λ matrix directly; run it for both models on A.7
- [x] 5.3 Removed-key and unknown-model error scenarios
- [x] 5.4 Partition tests (multichain, CSP append) under `hps_urry`
- [x] 5.5 Radii are model-independent: `rmin_2` and the NC↔ribosome radii are identical for `bt` and `hps_urry`; IDR force energy is independent of ribosome bead positions
- [x] 5.6 Short CSP smoke run (A.7 or B.2 tutorial settings) under both models: stages complete, no blow-up through the tunnel

## 6. Scientific validation

- [ ] 6.1 Under `hps_urry`, simulate a few IDPs from the 18-IDP benchmark and compare Rg with the published HPS-Urry values and with `bt`
- [ ] 6.2 Record the gap and attribute it (bonded terms, exclusions, electrostatics, cutoff) in the docs caveat for `hps_urry`
- [x] 6.3 Benchmark `bt`, `hps_urry` and `calvados2` against SAXS Rg for 22 IDPs at the experimental temperatures (`sandbox/2026_11_19_IDR_model_benchmark/REPORT.md`, 5 copies × 500 ns each). RMS relative error: `bt` 9.9%, `hps_urry` 17.5%, `calvados2` 10.2%. `bt` and `hps_urry` used the old `Rmin/2` table, `calvados2` the new one
- [ ] 6.4 Re-run `bt` and `hps_urry` of that benchmark with the current `Rmin/2` table, so all three models share one set of bead sizes
- [x] 6.5 Docs validation: `disordered_regions.rst` reports the 22-IDP, three-model benchmark of 6.3 (figures, summary and per-protein tables) in place of the 18-protein `bt` table; the "not benchmarked" statements are removed from the docs, Tutorial A.7 and the CHANGELOG. By decision the `Rmin/2` difference between the runs is not mentioned in the docs

## 7. Documentation

- [x] 7.1 `disordered_regions.rst`: λ form, the two models, the reduced YAML, removal of `idr_scale` / `eps_ev_kj`, the "HPS-Urry λ within topo's chain model" caveat
- [x] 7.2 `model_theory.rst`: equation in λ form; cross-pair conventions
- [x] 7.3 `optimization_control.rst`: remove references to the removed keys
- [x] 7.4 Tutorial A.7: revise `domain.yaml` and `README.md` to the new YAML
- [x] 7.5 CHANGELOG entry marked BREAKING; rebuild docs (`cd docs && make html`)

## 8. Per-residue `Rmin/2` table (D11; commit 70426f3)

- [x] 8.1 Replace the 20 per-residue `Rmin_2` values in `topo/parameters/model_parameters.py` with `rvdw·2^(1/6)` (nm) from the transferable residue vdW radii
- [x] 8.2 Update the comments that described the values as O'Brien's rounded radii (`model_parameters.py` docstring, `IDR_RMIN_2_NM` in `nonbonded.py`)
- [x] 8.3 Regenerate the verification snapshot and the two pinned radius sums in `test_idr_force.py::test_fully_idr_build_pinned`; well depths unchanged
- [x] 8.4 CHANGELOG entry marked BREAKING (IDR and ribosomal-protein radii shift; RNA radii unchanged)
- [ ] 8.5 Re-run the O'Brien 4c5c reproduction (Tutorial 12) with the new ribosomal-protein radii, or record that it was validated with the previous table

## 9. `calvados2` (D10; commit cacd9c6)

- [x] 9.1 Copy the 20 per-residue λ from `examples/single_RNA/residues_C2RNA.csv` (KULL-Centre/CALVADOS) into the generator script with the citation
- [x] 9.2 Generate `topo/parameters/data/idr/calvados2.yaml` (`λ_ij = (λ_i + λ_j)/2`, 4 decimals, `ε = 0.8368`)
- [x] 9.3 Add `calvados2` to `IDR_MODELS`
- [x] 9.4 Tests: shipped table equals its recomputation; hand-computed `λ_ij` pinned, all positive with GLU–GLU the smallest; the force's ε and the CSP-append partition under `calvados2`; unknown-model tests use a name that is still unsupported
- [x] 9.5 Docs: `calvados2` section and reference in `disordered_regions.rst`, with the "CALVADOS 2 λ within topo's chain model" caveat; `model_theory.rst`, `optimization_control.rst` and Tutorial A.7 list the model
- [x] 9.6 CHANGELOG entry; rebuild docs (`cd docs && make html`)
- [ ] 9.7 Confirm the citation against the paper: the values were recognised as the CALVADOS 2 set but not checked against the publication
