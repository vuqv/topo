# PLAN — Additive generic-cohesion term `eps_gen` for the IDR potential (2-parameter model)

**Status:** CODE IMPLEMENTED + SCAN RUN 2026-07-19. **RESULT: INCONCLUSIVE — blocked by
a SAMPLING/EQUILIBRATION failure, NOT a model failure.** (An earlier note here claimed
the additive approach "fundamentally can't compact" — that was WRONG and is retracted;
see below.)

**What the scan showed (scan2d, idr_scale=1, eps_gen_kj ∈ {0,2.51,3.77,5.02} kJ, 24 IDPs):**
RMS frac dev vs experiment stayed 43–48% at every eps_gen; the Rg cloud stayed
over-expanded — no compaction. BUT the runs are **not equilibrated**, so this says
nothing about the model:
- **Implementation verified correct:** eps_gen deepens EVERY IDR-IDR well by exactly
  +eps_gen_kj (hCyp mean 2.37→7.39 kJ = ~1→~3 kBT at eps_gen=5.02, 100% of pairs).
- **Deep wells SHOULD compact at equilibrium** (deeper attraction = poorer solvent) —
  the physics is sound. The sims just don't reach equilibrium.
- **Kinetic trapping:** the deep 12-10-6 has a repulsive shoulder at ~1.4·R that scales
  with eps (~0.4 kBT/pair at eps_gen=5); collapse needs many contacts to form
  cooperatively across those shoulders → a real barrier. From an EXTENDED start in 40 ns
  the chain can't cross it: hCyp minimizes to 2.73 nm then MD DRIFTS UP to 4.05 nm
  (escapes into the entropic expanded basin). Start-dependent hysteresis (folded α-scan
  runs collapsed; extended eps_gen runs don't).

**Next step:** this is an EQUILIBRATION problem. Enhanced sampling (replica exchange /
simulated tempering) would be the textbook fix, but it is **OUT OF SCOPE (user decision
2026-07-19)**. So the equilibration limit is accepted: large-chain (N≳100) Rg from
canonical MD is qualitative only, and conclusions lean on the small chains + the
fast folded→coil regime. See `review/IDR_TUNING_STATUS.md` for the accepted-constraint
framing and the prioritized combining-rule (`(Rmin_2_i+Rmin_2_j)/2^(1/6)`) experiment.
(The minimize infinite-loop bug in `system.py checkLargeForces` was fixed en route;
folded starts work now.)

---
_Original plan (kept for context):_

**Owner decision captured (2026-07-19):**
- Keep the BT shift at **0.6** (unchanged O'Brien/topo convention — do *not* adopt
  SOP-IDP's 0.7).
- Add **one** new parameter: an additive, sequence-independent generic-cohesion term
  `eps_gen_kj`.
- `domain.yaml` accepts `eps_gen_kj` **already in kJ/mol** (no kcal→kJ conversion in
  code — it is added directly to the other kJ quantities).
- Starting value for the first 2-D scan: **`eps_gen_kj = 2.5104`** ( = 0.6 kcal/mol
  × `KCAL_TO_KJ` = 0.6 × 4.184).

---

## 1. Motivation

Two independent lines of evidence say the single-knob `idr_scale` model cannot
reproduce IDP SAXS Rg, and that the missing ingredient is a *generic* (sequence-
independent) cohesive channel:

1. **Expanded-α scan (empirical).** Rg(α) rises to a peak near α ≈ 0.3–0.5 then
   collapses (turnover at eps ~ 1 kBT, as predicted), so the model *can* compact —
   but the per-protein turnover points **cross**: at α = 2.0, hCyp is still +32 %
   (needs more attraction) while IN is already −22 % (over-collapsed). No single
   global α fits compact and expanded chains at once. (See
   `sandbox/Tune_idr_scale/scan/REPORT.md` once the scan completes.)

2. **SOP-IDP comparison (model architecture).** Baul, Chakraborty, Mugnai, Straub &
   Thirumalai, *JPCB* 2019, 123, 3462 (DOI 10.1021/acs.jpcb.9b02575) reproduce SAXS
   Rg for 24 IDPs with a **two-bead, three-channel** energy:
   - ε_BB = 0.12 kcal/mol = 0.20 kBT (backbone–backbone, generic)
   - ε_BS = 0.24 kcal/mol = 0.41 kBT (backbone–sidechain, generic)
   - ε_SS = 0.18 kcal/mol = 0.30 kBT (sidechain–sidechain, × |ε_ij − **0.7**|)

   topo (Cα-only, fully-IDP) currently has **only the SS-analog** channel, and at the
   O'Brien default α = 0.03 it is ~7× weaker than SOP-IDP's SS *and* entirely lacks
   the ~0.6 kBT generic BB+BS cohesion. Because topo places every interaction on the
   single Cα bead sharing one 12-10-6 radial form, SOP-IDP's BB + BS channels can be
   **lumped into a single additive constant** on that bead — that constant is
   `eps_gen`.

See memory notes `idr-scale-expansion-mechanism`, `sop-idp-vs-topo-energy-scale`.

---

## 2. Model change

Current IDR-IDR well depth (`topo/utils/nonbonded.py:998`):

```python
eps_ij[dd] = np.maximum(NON_NATIVE_KJ, idr_scale * ss_interaction_energy[dd])
```

New form (shift stays 0.6; `ss_interaction_energy` is already `KCAL_TO_KJ·|w−0.6|`):

```python
eps_ij[dd] = np.maximum(NON_NATIVE_KJ,
                        eps_gen_kj + idr_scale * ss_interaction_energy[dd])
```

- `eps_gen_kj` — additive generic cohesion in **kJ/mol**, applied equally to every
  IDR–IDR pair (sequence-independent). Default **0.0** → byte-identical to today.
- `idr_scale` — unchanged role: scales the sequence-specific BT term `|w−0.6|`.
- Well **position** unchanged (sum-rule collision distance); only the depth gains the
  additive floor. The `max(NON_NATIVE_KJ, …)` guard is retained (moot once
  `eps_gen_kj > NON_NATIVE_KJ`, but keeps `eps_gen_kj = 0, idr_scale = 0` at the EV
  floor).

**Physical reading:** `eps_gen_kj` is the solvent-quality dial (monotone-decreasing
Rg, well-behaved), and `idr_scale·|w−0.6|` adds sequence specificity on top — instead
of forcing one BT-weighted channel to set both baseline compaction and specificity
(the original design flaw that produced the backwards Rg(α)).

---

## 3. Exact code changes (surgical, backward-compatible)

All in `topo/utils/nonbonded.py`:

1. **Config parse (`:806-808`)** — read optional `eps_gen_kj`, default 0.0:
   ```python
   idr_scale  = float(dis_cfg.get('idr_scale', DEFAULT_IDR_SCALE))
   eps_gen_kj = float(dis_cfg.get('eps_gen_kj', 0.0))     # kJ/mol, additive, generic
   disorder = {'residues': dis_residues, 'idr_scale': idr_scale,
               'eps_gen_kj': eps_gen_kj}
   ```
   (No unit conversion: the value in `domain.yaml` is already kJ/mol.)

2. **`apply_disorder` signature (`:936-939`)** — add `eps_gen_kj: float = 0.0`
   (defaulted, so any other caller is unaffected).

3. **Energy line (`:998`)** — as in §2.

4. **Callsite (`:1162-1164`)** — pass `disorder['eps_gen_kj']`.

5. **Docstrings** — update `apply_disorder` (§2.3 pair-class description) and the
   `domain.yaml` example block (`:684, :708, :729`) to document `eps_gen_kj`.

**New symbol name:** `eps_gen_kj` (mirrors `NON_NATIVE_KJ` — the `_kj` suffix signals
the unit and that it is added directly to kJ quantities).

### Backward compatibility / tests
- `eps_gen_kj` absent ⇒ 0.0 ⇒ current behavior exactly. All existing
  `tests/test_idr.py` cases (esp. `test_eps_construction`, `test_default_idr_scale`,
  `test_contact_removal_self_avoiding`) remain valid unchanged.
- **New test** (add to `tests/test_idr.py`): with `eps_gen_kj: X` and `idr_scale: a`,
  every IDR–IDR pair equals `max(NON_NATIVE_KJ, X + a·ss)`; and `eps_gen_kj: 0` matches
  the old expectation.

### NOT in scope for this change (raise separately if needed)
- The optimizer (`topo/optimize`) is a folded-Go strength tuner, unrelated to IDR
  fitting — do not touch.
- Well-position multi-scale geometry (SOP-IDP's BB is short-range, SS long-range):
  lumping onto one Cα well loses this. Accepted approximation; documented, not fixed.
- 12-10-6 (topo) vs 12-6 (SOP-IDP) well width: minor; don't change the functional form.

---

## 4. Second simulation set — **1-D `eps_gen_kj` scan first** (idr_scale fixed)

**Decision (2026-07-19):** the "2-D scan" starts as a **1-D scan** — fix
**`idr_scale = 1.0`** and vary **`eps_gen_kj` only**, starting at **2.51 kJ/mol**
(0.6 kcal/mol). Rationale: `idr_scale = 1` uses the full-strength sequence-specific BT
term (mean SS depth `⟨|w−0.6|⟩·KCAL_TO_KJ` ≈ 2.79 kJ/mol ≈ 1.1 kBT), and `eps_gen_kj`
then adds the generic cohesion that the α-only scan showed is missing. Once this 1-D
curve is understood, promote to the full 2-D `(eps_gen_kj, idr_scale)` scan if needed.

A **new sibling** directory `sandbox/Tune_idr_scale/scan2d/` (leave the α-only `scan/`
tree intact and reproducible):

- Reuse `gen_scan.py` logic with `idr_scale = 1.0` fixed and an `EPS_GEN_KJ` list.
  Each run dir: `<idx>_<name>/g<eps_gen_kj>/`.
- `domain.yaml` per run: `idr_scale: 1.0` + `eps_gen_kj: <value>`.
- Analysis: 1-D Rg(eps_gen_kj) per protein → objective `J(eps_gen_kj)` → optimum,
  per-protein residuals, ν_sim vs ν_exp. (Generalize `analyze_rg.py`'s α axis to an
  `eps_gen_kj` axis; same PCHIP-interpolation machinery.)

**Reuse the α-scan baseline:** the point `(idr_scale = 1.0, eps_gen_kj = 0)` is
*exactly* the existing `a1.00` runs in `scan/` (24 runs, `idr_scale = 1`, `eps_gen`
absent = 0). Symlink/point the analysis at them rather than recomputing — the 1-D scan
then only needs the `eps_gen_kj > 0` points.

**Starting grid (APPROVED 2026-07-19):** `eps_gen_kj ∈ {0.0 (reused), 1.25, 2.51,
3.77, 5.02}` kJ/mol ≈ {0, 0.5, 1.0, 1.5, 2.0} kBT at 300 K, centered on 2.51. From the
α-scan, `(idr_scale=1, eps_gen=0)` is still too expanded for most chains (e.g. hCyp
3.80, sNase 4.18 nm), so Rg should fall monotonically as `eps_gen_kj` rises — the
expected well-behaved solvent-quality dial. Watch the polyampholytes for over-collapse
(see §5.2).

---

> **Note (idr_scale = 1 is ~5× SOP-IDP's SS prefactor).** SOP-IDP uses ε_SS = 0.18
> kcal/mol; topo at `idr_scale = 1` (treat-BT-as-kcal convention) corresponds to a
> prefactor of 1 kcal/mol on `|w−0.6|`, ~5.5× stronger sequence-specific weighting.
> This is **fine for the first 1-D scan** — we are isolating the `eps_gen_kj` axis with
> the SS term at full strength — but flag it as the first knob to **reduce later**
> (toward ~0.18, i.e. `idr_scale` ≈ 0.18–0.30) once `eps_gen_kj` is carrying baseline
> compaction, if per-protein residuals show the sequence-specific term is over-weighted.

> **⚠️ Deep-well minimization hang (learned from the α-scan, 2026-07-19).** The
> α = 3.0 runs (and most α = 2.0) **hung indefinitely in build/minimize** — a folded
> starting PDB has many beads closer than the collision distance, so at deep wells
> (idr_scale ≥ 2, ~2–5 kBT) the initial `(R/r)^12` forces are enormous and OpenMM's
> minimize/context step chokes before dynamics start (0 output for 100+ min; all
> α=3.0 produced nothing). **The `eps_gen` scan at `idr_scale=1` + `eps_gen_kj` up to
> ~5 kJ/mol gives ~3 kBT total wells and WILL hit the same hang.** Mitigation before
> submitting `scan2d/`: start these runs from an **extended chain** (not the folded
> PDB), and/or add a brief soft-core pre-minimization; also throttle concurrency
> (`%8`–`%12`) so each task gets a dedicated GPU. Verify one deep-well run completes
> before fanning out.

## 5. Open questions for the post-scan discussion
1. Final `eps_gen_kj` and `idr_scale` grid points (coarse bracket sizes).
2. Confirm **electrostatics are active** in the fully-IDP runs (topo Yukawa/DH force,
   `system.py:649`, ε_r = 78.5) — a ~1 kBT generic attraction without the
   compensating charge repulsion will over-collapse polyampholytes (ProTα, osteopontin).
3. Whether to also let `eps_gen_kj` vary per-pair by bead identity later (currently a
   single global constant — the minimal, intended first step).
4. Objective weighting (equal per protein, as in the α-scan) — keep or revisit.
