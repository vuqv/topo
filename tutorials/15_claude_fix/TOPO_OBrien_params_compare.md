# TOPO vs O'Brien CG model — parameter & functional-form comparison

Side-by-side comparison of the coarse-grained force field as implemented in **topo**
(`topo.csp` / `topo.core`) vs the **O'Brien lab** reference
(`cg_simtk_protein_folding`), written in the **same functional form** so the energy terms
can be compared directly. All values are in **OpenMM units** (length nm, energy kJ/mol,
angle rad) unless noted.

## Sources

| Side | Files |
|------|-------|
| **O'Brien** | protein `.prm` `4c5c_model_clean_nscal1_fnn1_go_bt.prm`; ribosome `.prm` `combine_ribo_L24_Yang.prm`; consolidated `rnc.xml`; generator `CG_protein_parameterization/parse_cg_prm.py` (defines the OpenMM functional forms + unit conversions). |
| **topo** | `topo/core/system.py` (force builders), `topo/parameters/model_parameters.py` (per-residue mass/radii/charge, bond constants), `topo/parameters/dihedral.py` (`dihedral_params.csv`), `topo/utils/nonbonded.py` (native-contact energies), `topo/csp/ribosome.py`. |

**Unit conversions** used by O'Brien's `parse_cg_prm.py` (for reading their CHARMM `.prm`):
energy ×4.184 (kcal→kJ), length ÷10 (Å→nm), angle ×π/180 (deg→rad); bond `k` additionally
×100 (Å⁻²→nm⁻²) **×2** (CHARMM `Kb(r-r₀)²` → OpenMM `½k(r-r₀)²`).

## Summary

| Term | Functional form | Constants/parameters | Verdict |
|------|-----------------|----------------------|---------|
| Bond | identical (`½k(r−r₀)²`, both via `HarmonicBondForce`) | **FIXED**: topo k 20920→**41840** (added the missing CHARMM→OpenMM ×2); now matches O'Brien | ✅ (2026-06-30) |
| Angle | identical (dual-basin log-sum-exp) | identical (106.4/91.7/26.3/130/0.1/4.3) | ✅ |
| Dihedral | identical (Σ₁⁴ periodic) | identical (topo's ×0.756 reconstructs the prm — verified ratio 1.000) | ✅ |
| Improper | none (Cα model) in both | — | ✅ |
| Electrostatics | identical (Debye–Hückel) | identical (f=138.935, εᵣ=78.5, l_D=1 nm) | ✅ |
| Native contacts | identical (12-10-6 Gō) | same scheme (HB/BS/BT); source differs (NBFIX vs computed) | ✅ form / ⚠️ source |
| Non-native EV | identical (12-10-6, tiny ε) | **radii source + combination rule differ** | ⚠️ |
| Charges | per-residue ±1 | identical | ✅ |
| Masses | per-residue | identical (standard residue masses) | ✅ |
| Ribosome beads | excluded volume + electrostatics | topo previously mis-mapped R (C5′) + fat radii — **now fixed** to O'Brien's | 🔧 fixed |

---

## 1. Bonds

Both use OpenMM `HarmonicBondForce`:

$$U_\text{bond}(r) = \tfrac{1}{2}\,k\,(r-r_0)^2$$

| | O'Brien | topo |
|--|---------|------|
| force object | `mm.HarmonicBondForce` (`E=½k(r−r₀)²`) | `mm.HarmonicBondForce` (`E=½k(r−r₀)²`) |
| r₀ | 3.81 Å = **0.381 nm** (prm `BOND`) | **0.381 nm** (`bond_length_protein`) |
| k (in the ½-form) | `50·4.184·100·`**`2`** = **41840 kJ/mol/nm²** | **20920 kJ/mol/nm²** (`bond_force_constant`) |
| effective Kb | **50 kcal/mol/Å²** (= CHARMM) | **25 kcal/mol/Å²** |
| E at +1 Å stretch | **209.2 kJ/mol** | **104.6 kJ/mol** (empirically evaluated) |

✅ **FIXED (2026-06-30):** `model_parameters['topo']['bond_force_constant']` changed **20920 →
41840**. Both codes use OpenMM `HarmonicBondForce` (the ½ is in both); the gap was in `k` —
`parse_cg_prm.py` multiplies CHARMM's `Kb` by **2** so `½k = Kb·conv` reproduces `Kb(r−r₀)²`,
whereas topo's old `20920 = 50·4.184·100` omitted that ×2 (the value for a *no-½* form dropped
into a *½*-form → 2× too soft). Now `41840` gives **E = 209.20 kJ/mol at a 1 Å stretch =
O'Brien/CHARMM** (was 104.6). Moot on tut-15's `AllBonds` path (bonds are rigid constraints, so
`HarmonicBondForce` isn't even created), but it corrects any flexible-bond run (legacy 12/13,
package default `buildCoarseGrainModel`).

## 2. Backbone angles

Both use a `CustomAngleForce` with the **dual-basin (double-Gaussian, log-sum-exp)** form:

$$U_\text{angle}(\theta) = -\frac{1}{\gamma}\ln\!\Big[
 e^{-\gamma\,[\,k_\alpha(\theta-\theta_\alpha)^2 + \epsilon_\alpha\,]}
+e^{-\gamma\,k_\beta(\theta-\theta_\beta)^2}\Big]$$

| param | O'Brien (`.prm`, converted) | topo (`addGaussianAngleForces`) |
|-------|-----------------------------|----------------------------------|
| k_α | 106.4 kcal → **445.18** kJ/mol/rad² | **445.1776** |
| θ_α | 91.7° → **1.60047** rad | **1.6004669** |
| k_β | 26.3 kcal → **110.04** | **110.0392** |
| θ_β | 130.0° → **2.26893** rad | **2.2689280** |
| γ | 0.1 /4.184 → **0.023901** mol/kJ | **0.0239006** |
| ε_α | 4.3 kcal → **17.991** kJ/mol | **17.9912** |

✅ **Identical.** O'Brien's protein `.prm` uses these **uniform** values for every angle
triplet (verified: a single unique tuple across all `ANGLE` lines), and topo hardcodes the
same as global parameters.

## 3. Backbone dihedrals

Both use OpenMM `PeriodicTorsionForce`, sequence-dependent, four periodicities:

$$U_\text{tors}(\varphi) = \sum_{n=1}^{4} k_{D,n}\,[\,1+\cos(n\varphi-\delta_n)\,]$$

| | O'Brien | topo |
|--|---------|------|
| source | per-quartet from `.prm` `DIHEDRAL` (303 unique quartets; `k`×4.184, phase×π/180) | `dihedral_params.csv` keyed by the **two central residues** |
| scaling | prm values used as-is | `dihedral_params.csv` `k_D` **×0.756** (`dihedral.py:23`) |

✅ **Identical in effect.** topo's CSV stores the force constants pre-divided, so the 0.756
"calibration factor" exactly **reconstructs** O'Brien's `.prm` values. Verified numerically on
the first 4c5c quartet (`A1-A2-A3-A4`, central THR-ASP): topo(×0.756) = O'Brien for all four
periodicities (1.773 / 3.860 / 0.450 / 0.595 kJ/mol, **ratio 1.000**). topo keys on the central
residue *pair*; O'Brien stores per atom quartet — equivalent for a Cα chain.

## 4. Impropers

✅ **None in either** for the Cα protein model. `parse_cg_prm.py` explicitly skips the
`IMPROPER` section ("No improper torsion energy term for Ca model"); topo adds no improper.
(The improper *form* `k·min(Δθ,2π−Δθ)²` is reused by the tutorial-15 **tRNA tether**, §7.)

## 5. Electrostatics (screened Coulomb / Debye–Hückel)

$$U^\text{el}_{ij}(r) = f\,\frac{q_i q_j}{\epsilon_r\, r}\, e^{-r/l_D}$$

| param | O'Brien | topo (`addYukawaForces`) |
|-------|---------|--------------------------|
| f (Coulomb) | √ke·√ke = **138.935485** | **138.935458** kJ·nm·mol⁻¹·e⁻² |
| ε_r | EPS = **78.5** (prm `CUTNB`) | **78.5** |
| l_D | √ld·√ld = **1.0 nm** (10 Å) | **1.0 nm** |
| charges | per-residue ±1 (psf) | per-residue ±1 (`model_parameters`) |

✅ **Identical** (the 7-digit f difference is round-off). O'Brien folds this term **into** the
single combined `CustomNonbondedForce`; topo keeps it as a **separate** `yukawaForce` — same
energy, same 2.0 nm cutoff / 1.8 nm switch, bonded (1-2,1-3) exclusions.

## 6. Native + non-native contacts (12-10-6 Gō)

Both use the Karanicolas–Brooks **12-10-6** Gō potential (well minimum −ε at r=R):

$$U^\text{nb}_{ij}(r) = \varepsilon_{ij}\Big[\,13\big(\tfrac{R_{ij}}{r}\big)^{12}
 -18\big(\tfrac{R_{ij}}{r}\big)^{10} + 4\big(\tfrac{R_{ij}}{r}\big)^{6}\Big]$$

O'Brien writes it as `a/r¹²+b/r¹⁰+c/r⁶` with `a=13εR¹², b=−18εR¹⁰, c=4εR⁶` (tabulated per
pair); topo evaluates the factored form directly with per-pair `eps_table`/`R_table`.
✅ **Form identical.**

**Native contacts** (ε_ij, R_ij):

| | O'Brien | topo (`build_nonbonded_interaction`) |
|--|---------|--------------------------------------|
| R_ij | native pair distance (prm `NBFIX`, 846 pairs) | native Cα–Cα distance from the structure |
| ε_ij | from `NBFIX` (the `go_bt` scheme baked into the prm) | **HB + BS + SS**, domain-scaled |
| HB | 0.75 kcal/H-bond (cap 1.5) | **0.75** kcal/mol (cap 1.5) ✅ |
| BS | 0.37 kcal/mol | **0.37** kcal/mol ✅ |
| SS | Betancourt–Thirumalai (`go_bt`), shift 0.6 | **BT**, `bt_shift` 0.6 ✅ |
| domain scale | `nscal` (=1 here; baked into NBFIX) | per-domain `nscale` (`domain.yaml`) |

✅ Same energy scheme (HB/BS/BT constants match). ⚠️ **Different provenance**: O'Brien ships
the final per-pair ε/R in the `.prm` `NBFIX`; topo **recomputes** them at run time from the
structure + STRIDE + `domain.yaml`. (The tut-15 `nscale` in `domain.yaml` is topo's analog
of O'Brien's `nscal`.)

**Non-native pairs** (soft excluded volume — tiny ε, repulsion-dominated):

| | O'Brien | topo |
|--|---------|------|
| ε | √(ε_i·ε_j), ε = 0.000132 kcal/mol → **5.5×10⁻⁴ kJ/mol** | negligible well depth (same magnitude) |
| R_ij | **R_min/2(i) + R_min/2(j)** (sum of per-**position** Rmin/2 from prm) | **½(σ_i + σ_j)** (per-**residue-type** σ from `model_parameters`) |

⚠️ Two differences: (a) **radii source** — O'Brien assigns a per-residue-*position* Rmin/2
in the prm (`A1`=4.566, `A2`=4.205, `A3`=3.711 Å, …); topo uses a per-amino-acid-*type* table
(ALA 0.504, GLU 0.592, LYS 0.636, GLY 0.45, TRP 0.678 nm, …). (b) **combination rule** —
O'Brien sums the two half-radii; topo averages the two `σ`. The soft ε (5.5×10⁻⁴ kJ/mol)
makes this a weak repulsion in both, but the excluded-volume *range* differs.

## 7. Per-residue properties (mass / radii / charge)

topo `model_parameters['topo']` (selected): mass / radii(nm) / charge —
ALA 71/0.504/0, ARG 114/0.656/+1, ASP 114/0.558/−1, GLU 128/0.592/−1, GLY 57/0.45/0,
LYS 128/0.636/+1, TRP 186/0.678/0, VAL 99/0.586/0 (all 20 defined). **Masses** = standard
residue masses (match O'Brien's psf). **Charges** = ±1 on D/E/K/R, 0 otherwise (match).
**Radii** = the per-type excluded-volume σ used in §6's non-native rule (see ⚠️ above).

## 8. Ribosome beads (tutorial-15 fix)

O'Brien's CG ribosome (`combine_ribo_L24_Yang.prm` NONBONDED + `.psf`/`.cor`):

| bead | O'Brien Rmin/2 (nm) | O'Brien charge | topo (old, `ribosome_trunc.pdb`) |
|------|---------------------|----------------|----------------------------------|
| P (phosphate) | **0.645** | −1 | 0.71, −1 |
| R (ribose, **C5′**) | **0.523** | 0 | 0.71, 0 — **position also off by up to 3.6 Å** |
| BR (base) | **0.534** | 0 | 0.71, 0 |
| protein (per type) | 0.25–0.52 | ±1/0 | per-AA-type table |

> **The NC↔ribosome excluded-volume *interaction* (form + combination rule), which is what
> shapes the tunnel and gates egress, is analyzed in depth in
> [`TOPO_OBrien_NCribosome_nonbonded_compare.md`](TOPO_OBrien_NCribosome_nonbonded_compare.md).
> Key result: topo's pure `ε(σ/r)¹²` + average-σ rule is ~1000× softer than O'Brien's 12-10-6 +
> sum-of-Rmin/2 rule — the mechanistic root of the residual clash. (Confirms "Rmin/2 = σ/2".)**

ε (ribosome–NC excluded volume) = **0.000132 kcal/mol = 5.5×10⁻⁴ kJ/mol** in both. 🔧 **Fixed
in tut 15**: topo had built its *own* CG ribosome whose **R (ribose) bead used a different atom
than O'Brien's C5′** (P beads matched exactly; R beads differed mean 0.49 Å, max 3.6 Å,
shifting the PtR:76 R **P-anchor by 3.6 Å**) and a **uniform 0.71 nm** radius for all RNA beads.
`topo.csp.ribosome.load_obrien_ribosome` now reads O'Brien's authentic
`50S_tRNA_cg_truncated.{cor,psf,prm}` verbatim — his C5′ positions, per-type Rmin/2 radii, and
psf charges — so the ribosome geometry and excluded volume match the reference exactly.

## 9. tRNA tether (tutorial-15 §1b ✅2, O'Brien `continuous_synthesis_v6.py`)

Not part of the base force field; reproduced for the tether feature: bond `½k(r−r₀)²`
(k=200 kcal/mol/Å² = 83680, r₀ 4.27 Å A-site / 4.76 Å P-site); two harmonic angles
`½k(θ−θ₀)²` (k=25 kcal/mol/rad² = 104.6; θ₀ A: 106°/127°, P: 117°/130°); improper
`k·min(Δθ,2π−Δθ)²` (k=25 kcal/mol/rad²; θ₀ A 128°, P −161°). ✅ Matches O'Brien's
`A_site_tRNA_binding` / `translocation_AtR` constants.

---

### Bottom line

The topo CG force field is a **faithful port** of the O'Brien model: the **angle**,
**dihedral**, **electrostatic**, and **contact (12-10-6 Gō)** functional forms **and** their
constants match exactly (dihedrals verified numerically, ratio 1.000), as do charges and masses.
The substantive differences are: (1) the harmonic-**bond k** was 2× soft in topo (`20920`,
missing the CHARMM→OpenMM ×2) — **now fixed to `41840`** to match O'Brien (moot under `AllBonds`,
corrects flexible-bond runs); (2) **non-native excluded-volume radii** use a
different source (per-AA-type table vs per-position prm) and combination rule (½(σ_i+σ_j) vs
Rmin/2_i+Rmin/2_j) — a weak repulsion either way given ε≈5.5×10⁻⁴ kJ/mol; and (3)
native-contact ε/R are **recomputed** by topo from structure+STRIDE+`domain.yaml` rather than
read from the prm `NBFIX` (same HB/BS/BT scheme + constants, so they agree by construction).
The **ribosome** bead positions/radii inconsistency (C5′ R bead + uniform fat radii) has been
**fixed** in tutorial 15 by loading O'Brien's authentic truncated ribosome
(`load_obrien_ribosome`).
