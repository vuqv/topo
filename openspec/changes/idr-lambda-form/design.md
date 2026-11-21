> **Status: revised 2026-11-17 with the author's decisions** (YAML surface, precomputed
> pair tables, no σ in the model file, `bt` and `hps_urry` only, topo's bonded terms,
> electrostatics and cutoff for every model). No design question is open; the
> remaining prerequisite is archiving `idr-ashbaugh-hatch`.
>
> **Revised 2026-11-21, after implementation:** `calvados2` added as a third model
> (D5, D7, D10); the per-residue `Rmin/2` table replaced by `rvdw·2^(1/6)` (D3, D11).

## Context

The IDR force (introduced by `idr-ashbaugh-hatch`) evaluates every IDR-involving pair with
an Ashbaugh–Hatch split 12-6 written as

    r ≤ R:  epsEV·LJ + (epsEV − eps_ij)      r > R:  eps_ij·LJ,   σ = 2^(−1/6) R

with `eps_ij = idr_scale·ε_BT,ij` and `R_ij = Rmin/2_i + Rmin/2_j`.
Substituting `ε = epsEV`, `λ_ij = eps_ij/ε` gives the canonical HPS λ form exactly. The
current `epsEV = 0.8368 kJ/mol` equals the HPS ε (0.2 kcal/mol).

Folded–folded pairs stay on the Karanicolas–Brooks 12-10-6; the two forces partition the
pairs via interaction groups. Folded Rmin/2 is computed from the folded residues only
(commit 35ea4b3). An IDR residue's Rmin/2 is the per-residue value in
`model_parameters` (`IDR_RMIN_2_NM`).

A verification bundle (`tests/test_idr_ah_verification.py`, Tutorial A.7 system P28566,
folded / mixed / fully-disordered builds at three configurations) pins every force's
energy and the contact forces' per-bead forces, plus a numpy λ-form oracle.

## Goals / Non-Goals

**Goals:**
- Express the IDR force in the λ form.
- Ship each supported model as a precomputed residue-pair λ table selected by name.
- Support three models: `bt` (default, the current model), `hps_urry` and `calvados2`.

**Non-Goals:**
- User-tunable `idr_scale` or ε, scale sweeps, or a pure-WCA (`λ = 0`) reference.
- User-supplied λ tables, and any model other than `bt`, `hps_urry` and `calvados2`.
- Force fields whose short-range term is not Ashbaugh–Hatch (Mpipi / Wang–Frenkel, etc.).
- Changing bead sizes: σ never comes from a model file (D3).
- Changing the folded-domain (12-10-6) model or the force partition.
- Importing anything from a force field other than its λ and ε. For every model the
  bonded terms and exclusions are topo's, electrostatics are the Yukawa force, and the
  cutoff is topo's 1.8–2.0 nm switch.

## Decisions

### D1. Force = functional form only; the model fills the table
The force takes `(lambda_matrix, rmin_matrix, epsilon)`. It contains no knowledge of how
λ was derived.
*Alternative:* one force per force field — rejected, duplicates the partition logic that
has already produced silent double-counting bugs.

### D2. Per-pair tables, in the force and in the file
BT is genuinely pairwise (`ε_BT,ij` is not a function of `ε_i`, `ε_j`), so a per-residue
λ with a combining rule cannot represent it. Both the force (n×n `Discrete2DFunction`)
and the shipped file (20×20 residue-pair table) are therefore pairwise. A per-residue
set such as HPS-Urry is expanded to a pair table once, when its file is generated.

### D3. Bead size is not part of a model
No model file carries σ. For every model:

- an IDR residue's size is its per-residue `Rmin/2` from `model_parameters`;
- a folded residue's size is its structure-derived Karanicolas–Brooks `Rmin/2`;
- `R_ij = Rmin/2_i + Rmin/2_j` and `σ_ij = R_ij / 2^(1/6)`.

This is what the code does today. The force keeps receiving the well-position matrix
`R` and converts inside its expression, so `apply_disorder`'s radius handling, the
per-residue `rmin_2` array and the NC↔ribosome radii are all untouched by this change.
For `hps_urry` and `calvados2` this means topo's radii are used, not the published
per-residue σ. The two sets share one σ table, and with the `Rmin/2` values of D11
topo's `σ_ii` agrees with it to a mean +0.04% and at most 0.20% (GLY).

### D4. Write the expression with `select()`
`eps*select(step(r - R), lambda*LJ, LJ + 1 - lambda)` mirrors the piecewise definition
directly, replacing the current branch-free identity. Must be checked against the
bundle (C¹ join, no double count at `r = R`).

### D5. Model files: precomputed pair λ tables
`topo/parameters/data/idr/<model>.yaml`, next to `bt_potential.csv`:

```yaml
name: bt
citation: "Betancourt & Thirumalai, Protein Sci 1999; scaled for topo IDRs"
epsilon_kj: 0.8368
source: "lambda_ij = 0.10 * |BT_ij - 0.6| * 4.184 / 0.8368"
lambda:            # 20 x 20, symmetric, 4 decimals
  ALA: {ALA: ..., ARG: ..., ...}
  # ... all 20
```

- `bt`: `λ_ij = 0.10·ε_BT,ij / 0.8368` — the calibrated `idr_scale` and ε are baked in.
- `hps_urry`: `λ_ij = (λ_i + λ_j)/2` from the per-residue λ of Regy et al. (Protein Sci
  2021), taken from the `hps` field of `hps_urry` in the cosmo package's
  `cosmo/parameters/model_parameters.py` (the 0.08 shift is already applied there);
  `ε = 0.8368`.
- `calvados2`: `λ_ij = (λ_i + λ_j)/2` from the per-residue λ of Tesei & Lindorff-Larsen
  (Open Res. Europe 2, 94, 2022), copied from the `lambdas` column of
  `examples/single_RNA/residues_C2RNA.csv` in KULL-Centre/CALVADOS (amino acids only;
  the RNA beads of that file are not used); `ε = 0.8368`, which is also the CALVADOS ε.
  Every value is positive; the smallest pair is GLU–GLU at 0.0007.
- Values are rounded to **4 decimals**.
- The files are **derived**: a generator script rebuilds each from its source
  (`bt_potential.csv`; the 20 HPS-Urry and the 20 CALVADOS 2 per-residue λ, copied
  into the script), and a
  test checks the shipped file against that recomputation. They are not hand-edited.

### D6. Cross pairs (IDR–folded)
For every model, λ is looked up by the two residue types in the model's pair table, for
IDR–IDR and IDR–folded pairs alike: a folded residue facing an IDR residue takes the λ
of its residue type. Sizes follow D3. Folded–folded
pairs are never touched.

### D7. YAML surface

```yaml
disordered:
  residues: [1-18, 218-281]
  model: bt            # bt (default) | hps_urry | calvados2
```

`idr_scale` and `eps_ev_kj` are **removed**. A YAML that contains either fails with an
error naming the key and pointing to `model:`; silently ignoring them would change the
physics without telling the user. An unknown `model:` fails listing the two supported
names.

### D8. Co-translational simulation (ribosome present)
How the CSP path (`topo/csp/core.py`, `topo/csp/cylinder.py`) interacts with the change.
Current behaviour, which the design must preserve:

- The full-chain matrices are built **once** (`build_nonbonded_interaction` on the full
  PDB) and each nascent length L takes the top-left `L×L` block, plus the slice of the
  IDR index among the emerged beads. λ depends only on the two residue types, so the
  λ matrix is sliced exactly like `R_full` / `eps_full` today.
- Ribosome beads get placeholder `(0,)` entries in the IDR force (every
  `CustomNonbondedForce` needs one entry per particle) but are in **none** of its
  interaction groups, so the IDR force never sees a ribosome pair.
- NC↔ribosome contact is a **separate** force: O'Brien's 12-10-6 pure excluded volume
  (`RIBO_NC_EPS_KJ`), sum rule over per-bead Rmin/2. It does **not** use λ: no published
  λ set defines protein–rRNA/r-protein λ, so NC↔ribosome stays pure excluded volume +
  Yukawa for every `model:`.
- **Radii are model-independent (D3)**, so the `rmin_2` array fed to
  `append_ribosome(nascent_rmin_2=…)` is the same for every model. An IDR bead
  keeps one size against the chain and against the tunnel wall with no extra work.
- The per-length disorder handle (`{'idx', 'eps_ev_kj'}` today, built in
  `csp/core.py` and `csp/cylinder.py`) carries the model's ε instead of a user value.
- Interaction-group ownership: `append_ribosome` only restricts forces that carry no
  groups; the IDR force's groups stay set in exactly one place
  (`addIDRNonBondedForce`). The refactor must not move that.

### D9. Two-step landing, one deliberate snapshot regeneration
Rounding λ to 4 decimals moves `bt` energies by more than the 1e-9 the verification
bundle is pinned to. To keep the refactor provable:

1. **Refactor step.** Rewrite the force in λ form with λ computed in memory at full
   precision (`eps_ij/ε`). The bundle must pass **without** regenerating its snapshot.
2. **Table step.** Switch to the shipped 4-decimal tables and remove the YAML keys.
   Regenerate the snapshot once, in a commit that changes nothing else, and record the
   largest energy and force difference against the old snapshot.

Any other snapshot regeneration in this change is a review red flag.

### D10. `calvados2` is a third table, not a new mechanism
`calvados2` reuses everything `hps_urry` uses: the same combination rule
`λ_ij = (λ_i + λ_j)/2`, the same ε, the same generator, loader and file format. Adding it
is one entry in the generator, one data file and one name in `IDR_MODELS`. Nothing from
CALVADOS other than the 20 λ enters topo: not its σ (D3), charges (HIS is neutral in
both), bond length, exclusions, electrostatics or cutoff.

Alternative considered: also adopt the CALVADOS σ for this model. Rejected, because D3
keeps bead sizes out of model files and the two σ sets already agree to 0.2% (D3, D11).

### D11. Per-residue `Rmin/2` is `rvdw·2^(1/6)`
The 20 per-residue `Rmin_2` values in `model_parameters` are computed from the
transferable residue vdW radii `rvdw` (full precision) as `Rmin/2 = rvdw·2^(1/6)`, i.e.
`σ = 2·rvdw` and `Rmin = 2^(1/6)·σ`. The previous values were O'Brien's `S<aa>` radii,
which are the same radii rounded; the new values differ from them by at most 1.4%
(ILE, LEU, MET, CYS +1.42%; GLN −1.38%).

This does not change D3: bead sizes still come from `model_parameters` and not from a
model file, and the conversion to σ still happens inside the force expression. It
changes the numbers only, so the verification snapshot and the pinned radius sums are
regenerated in the same commit, and it is recorded as BREAKING in the changelog because
the table is shared with the ribosomal-protein beads.

## Risks / Trade-offs

- **[`hps_urry` λ alone is not the published force field]** → the λ set was calibrated
  with HPS's bonded model: harmonic bonds only, no angles or torsions, and only
  directly bonded beads excluded (`bonded_exclusions_index: 1` in cosmo's `hps_urry`).
  topo gives IDR beads the same bonds, angles and torsions as folded beads and excludes
  1-2 **and** 1-3 pairs from every nonbonded force (`bonded_exclusions_index: 2`).
  Electrostatics (topo's Yukawa) and the cutoff also differ. All of this is by
  decision: `hps_urry` means "HPS-Urry λ within topo's chain model". Mitigation:
  document the model as "HPS-Urry λ within topo's chain model"; benchmark Rg against
  published values (task group 6) to measure the gap.
- **[`calvados2` λ alone is not the published force field]** → same situation as
  `hps_urry`. CALVADOS 2 has bonds only (0.38 nm), excludes only bonded neighbours, and
  has its own σ, charges, salt- and temperature-dependent electrostatics and a 2.0 nm
  Ashbaugh–Hatch cutoff; its λ were optimised together with those. By decision,
  `calvados2` means "CALVADOS 2 λ within topo's chain model": charges, bond length,
  bonded terms, exclusions, electrostatics and cutoff stay topo's. Mitigation: the same
  docs warning as `hps_urry`; benchmark Rg against experiment (task group 6).
- **[Benchmark mixes two radius tables]** → the 22-IDP benchmark (task 6.3) ran `bt` and
  `hps_urry` before D11 and `calvados2` after it. The radii differ by at most 1.4%,
  which is comparable to the `bt`–`calvados2` gap in RMS error (9.9% against 10.2%), so
  that benchmark does not rank those two. Mitigation: task 6.4 re-runs `bt` and
  `hps_urry` with the current table.
- **[Ribosomal-protein radii moved with D11]** → `OBRIEN_SC_RMIN_2_NM` reads the same
  table, so the ribosome scenery is no longer bit-identical to O'Brien's `.prm`. The
  O'Brien 4c5c reproduction (Tutorial 12) was validated with the old values and has not
  been re-run. RNA bead radii are unchanged.
- **[Negative λ]** → HPS-Urry has `λ_GLU = −0.08`, so exactly one pair type, GLU–GLU,
  has `λ_ij < 0` (−0.08) and is purely repulsive (no well); the next smallest is
  GLU–ASP at 0.067. The published value is kept and the repulsion is intended. Tests
  and docs must not assume a well at `−λε` for every pair.
- **[Cross-pair size mixes conventions]** → folded K-B Rmin/2 (up to ~5.4 Å) is larger
  than the per-residue IDR Rmin/2. Same situation as today; state it in docs.
- **[4-decimal rounding]** → small λ lose relative precision (a λ of 0.05 is stored to
  0.2%). Accepted; quantified once in the D9 table step.
- **[Derived files drift from their source]** → mitigated by the generator script and
  the recomputation test (D5).
- **[Existing inputs stop building]** → BREAKING by design (D7). Mitigation: a clear
  error message, a CHANGELOG entry, the revised tutorial.
- **[Expression rewrite breaks equivalence subtly]** → D9 step 1 must pass without
  regenerating the snapshot.
- **[CSP / multichain]** → the ribosome append and replicator must preserve the
  partition; existing `test_idr_force.py` CSP tests cover this and must stay green.
- **[CSP not covered by the verification bundle]** → the bundle pins isolated-chain
  builds only. Mitigation: extend it with a co-translational case (toy ribosome as in
  `test_idr_force.py`, two nascent lengths — one before and one after an IDR emerges)
  **before** the refactor, pinning every force energy including NC↔ribosome EV and
  Yukawa.
- **[Resume across the change]** → CSP stages resume from OpenMM checkpoints/state.
  Renaming the `epsEV` global can make a checkpoint written before the change fail to
  load, or load with stale parameters. Mitigation: check `tests/csp/test_resume*.py`;
  either keep the global name stable or document that in-flight CSP runs must finish on
  the old version.

## Migration Plan

1. Archive `idr-ashbaugh-hatch` first; convert this change's deltas (REMOVED
   "Excluded-Volume Strength Parameter", MODIFIED "YAML Configuration Parsing").
2. D9 step 1: λ-form refactor, full-precision λ, snapshot unchanged.
3. D9 step 2: shipped tables for `bt`; remove the YAML keys; regenerate the snapshot;
   revise tests, Tutorial A.7 and docs.
4. Add `hps_urry`.
5. Replace the per-residue `Rmin/2` table (D11); regenerate the snapshot and the pinned
   radius sums in the same commit (70426f3).
6. Add `calvados2` (D10; commit cacd9c6).

Rollback: revert the commits. Inputs that used `idr_scale` / `eps_ev_kj` need those keys
removed to build after step 3.

## Open Questions

None open. Resolved 2026-11-17:

- Bonded terms and exclusions are topo's for every model, `hps_urry` included.
- Electrostatics are the Yukawa force for every model.
- The cutoff is topo's for every model.
- A folded residue in an IDR–folded pair takes λ by residue type under every model.
- The HPS-Urry per-residue λ come from the cosmo package; GLU–GLU stays at −0.08
  (purely repulsive).

Resolved 2026-11-20:

- `calvados2` is supported as a third model, built like `hps_urry`: only the
  per-residue λ are taken; charges and bond length stay topo's (D10).
- The per-residue `Rmin/2` is `rvdw·2^(1/6)` from the transferable residue vdW radii
  (D11).
