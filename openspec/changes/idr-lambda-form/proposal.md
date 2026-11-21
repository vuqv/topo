> **Status: revised 2026-11-17 with the author's decisions on the YAML surface, the
> parameter files and the supported models.** No design question is open (see
> design.md). Prerequisite before implementing: archive `idr-ashbaugh-hatch`.
>
> **Revised 2026-11-21, after implementation:** a third model, `calvados2`, was added
> (commit cacd9c6), and the per-residue `Rmin/2` table that sets IDR bead sizes was
> replaced by `rvdw·2^(1/6)` from the transferable residue vdW radii (commit 70426f3).

## Why

The IDR contact force is already an Ashbaugh–Hatch potential, but it is written in terms
of a pair well depth `eps_ij` and a core strength `epsEV`, and its parameters are built
at run time from two user knobs (`idr_scale`, `eps_ev_kj`). Published IDP force fields
(HPS-Urry and others) ship a hydrophobicity **λ** for the same functional form.
Rewriting the force in the canonical λ form, and shipping each parameter set as a
precomputed table, lets a set be selected by name without touching the force.

The rewrite is energy-identical for the current model: with `ε = epsEV` and
`λ_ij = eps_ij / epsEV`,

    U_ij(r) = 4ε[(σ/r)^12 − (σ/r)^6] + (1 − λ_ij) ε        r ≤ 2^(1/6) σ_ij
    U_ij(r) = λ_ij · 4ε[(σ/r)^12 − (σ/r)^6]                r > 2^(1/6) σ_ij

and the current `epsEV = 0.8368 kJ/mol` is already the HPS ε (0.2 kcal/mol).

`idr_scale` and ε are properties of a model, not free inputs: for the BT-derived set,
`idr_scale = 0.10` with `ε = 0.8368 kJ/mol` is the calibrated pair, and a λ set from
another force field comes with its own ε. They therefore move out of the user YAML and
into the shipped parameter set.

## What Changes

- The IDR force (`system.addIDRNonBondedForce`) is expressed in the λ form above, taking
  a per-pair **λ table**, the existing well-position matrix and a scalar **ε** — the
  functional form only, no knowledge of where the parameters come from.
- Each model is a shipped data file holding a precomputed **20×20 residue-pair table
  `λ_ij`** (4 decimals), its ε and a citation:
  - `bt` (default): `λ_ij = 0.10·ε_BT,ij / 0.8368`, i.e. today's `eps_ij/epsEV`.
  - `hps_urry`: `λ_ij = (λ_i + λ_j)/2` from the published per-residue λ.
  - `calvados2`: `λ_ij = (λ_i + λ_j)/2` from the published CALVADOS 2 per-residue λ.
- Only `bt`, `hps_urry` and `calvados2` are supported. There is no user-supplied table.
- Bead sizes do **not** come from the model file. As today, an IDR residue's size is its
  per-residue `Rmin/2` from `model_parameters`, a folded residue's is structure-derived
  (Karanicolas–Brooks), and `σ_ij = (Rmin/2_i + Rmin/2_j)/2^(1/6)`.
- **Numerical change (small, deliberate)**: the per-residue `Rmin/2` table in
  `model_parameters` is now `rvdw·2^(1/6)` from the full-precision transferable residue
  vdW radii, so `σ_ii = 2·rvdw`. The previous values were O'Brien's rounded radii; the
  new ones differ by at most 1.4%. The table also sets the ribosomal-protein bead
  radii, so those shift too. The verification snapshot is regenerated a second time.
- The `disordered:` YAML section accepts `residues` and an optional `model:` only.
- **BREAKING**: `idr_scale` and `eps_ev_kj` are no longer accepted in the `disordered:`
  section; a YAML containing either fails with an error pointing to `model:`.
- **Numerical change (small, deliberate)**: the shipped `bt` table is rounded to 4
  decimals, so `bt` energies differ slightly from the current full-precision values and
  the verification snapshot is regenerated once (see design.md, D9).

## Capabilities

### New Capabilities
<!-- none: the work extends the existing IDR capability -->

### Modified Capabilities
- `intrinsically-disordered-regions`: the IDR force is specified in λ form; parameters
  come from a named, shipped set selected by `model:`; `idr_scale` and `eps_ev_kj` are
  removed from the YAML. Requirements are written as ADDED here because the predecessor
  change `idr-ashbaugh-hatch` (which introduced the IDR force and `eps_ev_kj`) is not
  yet archived; once it is, its "Excluded-Volume Strength Parameter" requirement is
  REMOVED and "YAML Configuration Parsing" becomes MODIFIED.

## Impact

- `topo/core/system.py` — `addIDRNonBondedForce` signature and expression.
- `topo/utils/nonbonded.py` — the YAML reader (`model:`; rejects `idr_scale`,
  `eps_ev_kj`); `apply_disorder` / `build_nonbonded_interaction` fill IDR-involving
  pairs from the selected λ table; `DEFAULT_IDR_SCALE` / `DEFAULT_EPS_EV_KJ` leave the
  run-time path.
- `topo/core/models.py`, `topo/csp/core.py`, `topo/csp/cylinder.py` — callers of the IDR
  force and builders of the per-length disorder handle (`{'idx', 'eps_ev_kj'}` today).
  Co-translational runs (design D8): the λ table is sliced per nascent length like
  today's matrices; ribosome beads stay outside the IDR force; NC↔ribosome remains
  O'Brien excluded volume + Yukawa for every model; radii are unchanged. CSP checkpoints
  written before the change may not resume after it if global parameter names change.
- New data: `topo/parameters/data/idr/bt.yaml`, `topo/parameters/data/idr/hps_urry.yaml`,
  `topo/parameters/data/idr/calvados2.yaml`, and a generator script that rebuilds them
  from their sources.
- `topo/parameters/model_parameters.py` — the 20 per-residue `Rmin_2` values. This
  table is read by `IDR_RMIN_2_NM` (IDR beads) and by `OBRIEN_SC_RMIN_2_NM`
  (ribosomal-protein beads in `topo/csp/ribosome.py`).
- Tests: `test_idr.py`, `test_idr_force.py` and `test_idr_ah_verification.py` set
  `idr_scale` / `eps_ev_kj` in YAML and must be rewritten; the verification snapshot is
  regenerated once, in its own step (D9).
- Tutorial A.7 (`domain.yaml`, `README.md`) is revised to the new YAML.
- Docs: `docs/usage/disordered_regions.rst`, `model_theory.rst`,
  `optimization_control.rst`.
- Depends on: `idr-ashbaugh-hatch` (archive first).
