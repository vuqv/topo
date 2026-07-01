# Nascent-chain Rmin/2 in the NC↔ribosome excluded volume — Option A vs Option B

The NC↔ribosome force (O'Brien 12-10-6, sum rule `R_ij = Rmin/2_i + Rmin/2_j`, ε matched) needs a
per-bead Rmin/2 for the **nascent** side. Two choices were considered:

- **Option B (per-AA `SA..SY`)** — O'Brien's per-amino-acid *ribosomal-protein* sidechain radii
  (`OBRIEN_SC_RMIN2_NM`, ~0.25–0.38 nm). Simple 20-value table. **What was implemented first.**
- **Option A (per-residue Karanicolas–Brooks)** — the structure-derived collision diameter
  `Rmin/2_i = min_{j: non-native, |i−j|>2} d_native(Cα_i,Cα_j) · 2^(1/6)/2` (varies per residue,
  2.24–5.44 Å for 4c5c). **This is what O'Brien's nascent chain actually uses** (verified from
  `rnc.xml`/`rnc.prm`: nascent = per-position `A1..An`; NBFIX only A–A + B–B; nascent↔ribosome =
  LB of self-Rmin ⇒ nascent side = per-position A_i). See `TOPO_OBrien_NCribosome_nonbonded_compare.md`.

**Decision (user):** switch to **Option A** (faithful to O'Brien). This file records Option B's
results as the baseline so Option A can be compared after implementation.

---

## Option B baseline (committed state `0c8f7d0`; code: `OBRIEN_SC_RMIN2_NM[resname]` nascent side)

Runs: `synth_out/` (4c5c, `synth_out_softEV/` = the earlier pure-r¹² baseline), `P0CX28/synth_out/`.
Config: equil-PTC + AllBonds, O'Brien ribosome (`ribosome_obrien.cor`), NC-EV 12-10-6 + sum rule,
scale_factor 216564650, seed 20240629.

### 4c5c (L=1→306)
| metric | Option B value |
|--------|----------------|
| hard clashes (# nascent center <3 Å) | **0** |
| min nascent–ribosome center-dist (synth) | **3.32 Å** |
| min O'Brien-gap (center − ΣRmin/2) | −4.76 Å |
| # residues gap<0 | 63/306 |
| extrusion x-range | 10..**179 Å** |
| threading corr(res,x) | −0.863 |
| D5b (wall / clash / egress) | PASS / PASS / PASS |
| ejection min NC dist | 3.08 Å |
| dwell ratio vs ref | 1.01× |
| final R_g (L=10) vs ref | 0.96× (0.722 vs 0.750 nm) |
| worst max\|PotE\| | 1.9e6 kJ/mol (≪1e9) |
| dt-halving recoveries (blow-ups / retries) | **96 / 201** |

### P0CX28 (L=1→106)
| metric | Option B value |
|--------|----------------|
| hard clashes | **0** |
| min nascent–ribosome center-dist (synth) | **3.91 Å** |
| threading corr(res,x) | −0.969 |
| extrusion x-range | 10..**98 Å** |
| final R_g | 2.306 nm (1.09× native 2.111) |
| D5b (wall / clash / egress) | PASS / PASS / PASS |
| ejection: CoM-x net / min NC dist | +15.7 Å / 3.18 Å |
| worst max\|PotE\| | 2.75e3 kJ/mol |
| dt-halving recoveries | **22 / 47** |

**Key Option B observations to compare against:** clash eliminated on both; strong extrusion;
near-native P0CX28 R_g; BUT heavy dt-halving (96 / 22 recoveries) and high transient energy
(4c5c 1.9e6) because the seed sits inside the wall (`optimal_ptc_targets` uses the old soft-EV
metric). Option A uses *larger* per-residue radii on average (~0.385 vs ~0.32 nm), so without a
seeding fix the dt-halving count may rise — watch that metric.

---

## Seeding-consistency fix (Option A) — big dt-halving win

After adopting Option A, `optimal_ptc_targets` was changed to minimize the **same 12-10-6 sum-rule
NC↔ribosome EV** the simulation uses (was the old soft pure-`(σ/r)¹²` average-rule metric), so the
A/P seed is placed at the least-buried point of the *real* wall (a_target clearance 2.37 → 4.27 Å).
4c5c L=1→306:

| metric | Option A, no seeding fix (`synth_out_optA_noseed`) | Option A + seeding fix (`synth_out`) |
|--------|-----------------------------------------------------|--------------------------------------|
| **dt-halving (blow-ups / retries)** | 263 / 534 | **0 / 0** |
| worst max\|PotE\| | 1.94e4 kJ/mol | **1.02e3 kJ/mol** |
| hard clashes (<3 Å) | 0 | 0 |
| min nascent–ribosome center-dist | 4.12 Å | **4.79 Å** (better) |
| threading corr / x<0 (leak) | −0.852 / 0 | −0.773 / 0 |

**The seeding fix eliminates all blow-ups** (seed no longer born inside the wall) → no dt-halving
retries → faster + cleaner energies, while keeping 0 clashes and no leak (egress re-checked). This
also removes the main downside of Option A vs B (its extra dt-halving from larger radii).

## Option A results (fill after implementation)

### 4c5c (L=1→306) — DONE
| metric | Option A value | vs B |
|--------|----------------|------|
| hard clashes (<3 Å) | **0** | tie (both 0) |
| min center-dist (synth) | **4.12 Å** | better (B 3.32) — held further off ribosome |
| extrusion x-range | 10..173 Å | ~tie (B 10..179) |
| threading corr | −0.852 | ~tie (B −0.863) |
| D5b (wall/clash/egress) | **PASS / PASS / PASS** | tie (both pass) |
| ejection egress: CoM-x net / end dist-from-ribo | **+160.8 Å / 140.8 Å (fully clears)** | **A ≫ B** (B +22.5 / stayed near) |
| ejection min NC dist | 3.88 Å | ~tie (B 3.08) |
| dwell ratio | 1.01× | tie |
| final R_g (L=10) vs ref | 0.849 nm = 1.13× | B closer (0.96×); both within ~15% |
| worst max\|PotE\| | **1.94e4 kJ/mol** | **A cleaner** (B 1.9e6) |
| dt-halving recoveries (blow-ups / retries) | **263 / 534** | **A ~2.7× worse** (B 96 / 201) |
| K–B Rmin/2 reproduces O'Brien `.prm` A-values? | **yes** (76/306 exact, mean 0.06 nm; residual = topo native-map defn) | — |

**4c5c verdict:** Option A is **more faithful** (per-position K–B = O'Brien's actual nascent radii),
**excludes better** (4.12 vs 3.32 Å), has **cleaner per-stage energies** (1.94e4 vs 1.9e6), and shows
a **much more complete egress** (chain fully leaves, CoM +161 Å, ends 140 Å from the ribosome, vs
Option B's +22 Å). Both pass all D-checks. **Cost:** ~2.7× more dt-halving (larger radii ⇒ seeds sit
deeper in the wall — the seeding-consistency fix would reduce this).

### P0CX28 (L=1→106) — DONE
| metric | Option A value | vs B |
|--------|----------------|------|
| hard clashes (<3 Å) | 0 | tie (both 0) |
| min center-dist (synth) | 3.86 Å | ~tie (B 3.91) |
| threading corr | −0.969 | tie (B −0.969, identical) |
| extrusion x-range | 10..91 Å | ~tie (B 10..98) |
| final R_g | 2.219 nm | ~tie (B 2.306; both ≈ native 2.111) |
| D5b (wall/clash/egress) | PASS / PASS / PASS | tie (both pass) |
| ejection: CoM-x net / min NC dist | +21.1 Å / 3.01 Å | ~tie (B +15.7 / 3.18) |
| worst max\|PotE\| | 7.96e3 kJ/mol | B lower (2.75e3) |
| dt-halving recoveries | 61 / 127 | **A ~2.8× worse** (B 22 / 47) |

**P0CX28 verdict:** Option A and B are **essentially equivalent** in outcome (clash, threading,
extrusion, near-native R_g, D5b all pass); Option A just costs ~2.8× more dt-halving.

---

## Overall verdict & recommendation

| | Option B (per-AA `SA..SY`) | Option A (per-residue K–B) |
|--|--|--|
| Faithful to O'Brien's nascent radii | ✗ (those are ribosomal-protein radii) | ✅ (O'Brien's actual `A_i`) |
| 4c5c: clash / energy / egress | 0 / 1.9e6 / +22 Å | **0 / 1.94e4 / +161 Å (fully clears)** |
| P0CX28 | ✅ all pass | ✅ all pass (≈ B) |
| D-checks (both proteins) | all PASS | all PASS |
| dt-halving cost | 96 (4c5c) / 22 (P0CX28) | **263 / 61 (~2.7–2.8× more)** |

**Recommendation: adopt Option A as the default.** It is the O'Brien-faithful choice (per-position
Karanicolas–Brooks = his `A_i` nascent radii), is *at least as good* on every D-check for both
proteins, and is clearly better for 4c5c (better exclusion, cleaner energies, a far more complete
egress). Its only downside is ~2.7–2.8× more dt-halving, caused by the *larger* per-residue radii
seeding the new residue deeper into the (now correct) wall. **Next improvement:** make
`optimal_ptc_targets` place the A/P seed against the *same* 12-10-6 + sum-rule EV so the new residue
is born clear of the wall — that should cut the dt-halving overhead back down while keeping Option A's
fidelity. (Option B remains available as the `append_ribosome(..., nascent_rmin2=None)` fallback.)
