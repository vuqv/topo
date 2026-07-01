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

### P0CX28 (L=1→106)
| metric | Option A value | vs B |
|--------|----------------|------|
| hard clashes | | |
| min center-dist / corr / x-range / R_g | | |
| D5b | | |
| dt-halving recoveries | | |
