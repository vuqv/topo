# Step 2 — Optimal PTC geometry: A-site / P-site restraint points 3.81 Å apart

**Goal.** Place the new amino acid (A-site, residue **L**) and hold the previous one
(P-site, residue **L−1**) at restraint points that are **exactly one peptide-bond length
(3.81 Å) apart** and **clear of the ribosome excluded volume**, so the `AllBonds` peptide-bond
constraint is never stretched and the **dt-halving guard can be deleted** (step 1 goal).

This is the geometric prerequisite for step 1. It uses topo's CG ribosome **as-is** (O4′ ribose
convention — *not* switched; see [`AGENTS.md`](AGENTS.md) §CG-mapping) and reproduces O'Brien's
PTC resting geometry on it.

Source: [`step1_allbonds_no_dt_halving.md`](step1_allbonds_no_dt_halving.md),
[`../../topo/csp/core.py`](../../topo/csp/core.py) (`seed_positions`, `run_length`),
[`../../topo/csp/ribosome.py`](../../topo/csp/ribosome.py) (radii, exclusions),
O'Brien `continuous_synthesis_v6.py` (`A_site_tRNA_binding`, `create_elongation_system`).

---

## The problem being solved

Current topo seeds residue L at `AtR:76@R + 0.4 nm·x̂` and holds L−1 near
`PtR:76@R + offset`, which puts them **≈ 9.2 Å ≈ 1 nm apart** while the peptide-bond
equilibrium is **3.81 Å**. With the bond always present, a rigid `AllBonds` constraint can't be
seeded that far off → flexible bonds + the dt-halving stability guard are forced.

**Fix:** seed/hold the two residues at points that are 3.81 Å apart from the start.

---

## Result — optimal restraint points (topo CG ribosome)

> **These are fixed restraint *target points*, not bonds.** topo's CSP path holds the
> C-terminus with a `CustomExternalForce` position restraint to a **fixed point in space**
> ([`core.py:add_cterm_restraint`](../../topo/csp/core.py)) — it does **not** bond the residue
> to `AtR:76@R`/`PtR:76@R`. The "derived from" column below is only **how the point's location
> was computed** (a geometric construction using O'Brien's tRNA bond lengths so the target lands
> at the real PTC); once computed, each is just an absolute coordinate the residue is pinned to.

| | coordinate (nm) | coordinate (Å) | derived from (construction only) |
|--|--|--|--|
| **A-site** (new AA **L**, seed + stage-1/2 restraint) | `(0.9514, −0.1990, 0.0000)` | (9.514, −1.990, 0.000) | point 4.27 Å from `AtR:76@R` |
| **P-site** (prev AA **L−1**, stage-3 / carried restraint) | `(1.0280, 0.1671, −0.0725)` | (10.280, 1.671, −0.725) | point 4.76 Å from `PtR:76@R` |
| PTC midpoint | `(0.9897, −0.0159, −0.0362)` | (9.897, −0.159, −0.362) | — |

- **`|A − P| = 3.810 Å`** — exactly the peptide-bond equilibrium. ✅
- `|A − AtR:76@R| = 4.27 Å`, `|P − PtR:76@R| = 4.76 Å` — O'Brien aminoacyl / peptidyl-tRNA
  bond lengths, used **only to locate** the target points at the true PTC (no bonded term).
- **No ribosome overlap.** Nearest *non-excluded* beads: `23S:2506@R` at 5.73 Å (A),
  `PtR:76@P` at 5.49 Å (P). Excluded-volume energy ≈ 0.005–0.007 kcal/mol ≪ kT. ✅

Direction vectors (if re-deriving the *point* as `R_bead + d·û`; still a fixed point, not a bond):
- `û_A = (a_target − AtR:76@R)/4.27 ≈ (0.288, 0.902, 0.320)`
- `û_P = (p_target − PtR:76@R)/4.76 ≈ (0.961, −0.274, −0.029)`

---

## Method (why these are "optimal")

Minimize **O'Brien's own restraint energy** restricted to the new + previous AA:

```
E(A,P) = 200·(|A−P|−3.81)²            # peptide bond           (stiff, k=200 kcal/mol/Å²)
       + 200·(|A−AtR:R|−4.27)²        # A-site tRNA bond       (stiff)
       + 200·(|P−PtR:R|−4.76)²        # P-site tRNA bond       (stiff)
       +  25·(angle(A,AtR:R,AtR:P)−106°)² + 25·(angle(A,AtR:R,AtR:PU2)−127°)²   # soft orient
       +  25·(angle(P,PtR:R,PtR:P)−117°)² + 25·(angle(P,PtR:R,PtR:PU2)−130°)²
       +  25·(dih(A,AtR:R,AtR:P,AtR:PU2)−128°)² + 25·(dih(P,PtR:R,PtR:P,PtR:PU2)+161°)²
       + Σ_beads ε·(σ_ij/r)¹²         # ribosome excluded volume, ε=0.000132 kcal/mol,
                                       # σ_ij=½(σ_AA+σ_bead), σ_bead=0.71 nm, σ_AA=0.678 nm (TRP)
```

with **O'Brien's exclusions** (new AA↔`AtR:76@R`, prev AA↔`PtR:76@R`;
`continuous_synthesis_v6.py:759,762`) and a 60-point multistart over the residual rotational
freedom. Equivalent framing: solve the 3 hard distance constraints (4.27 / 4.76 / 3.81) and use
the remaining 3 DOF to minimize ribosome clash.

**Caveats on "do not overlap":**
- Nominal pair contact σ = ½(0.678+0.71) nm = **6.94 Å**, so 5.5–5.7 Å is *inside* σ — but with
  ε = 0.000132 kcal/mol the repulsion only reaches kT at **r ≈ 3.4 Å**. The soft excluded volume
  is genuinely clear; reaching 6.94 Å from *every* bead is impossible at the crowded PTC and is
  not required.
- The **only hard requirement** for `AllBonds` stability is `|a_target − p_target| = 3.81 Å`.
  The 4.27 / 4.76 tRNA distances merely anchor the pair at the PTC — relax them for more
  clearance freedom if ever needed.

Reproduce: `scipy.optimize.minimize` (SLSQP, 3 equality constraints) over
`ribosome_trunc.pdb` beads + `model_parameters['topo']` radii.

---

## How it removes dt-halving (link to step 1)

With A→P migration over the 3 stages, **every residue is born exactly 3.81 Å from the current
C-terminus** (the just-migrated previous residue sits at `p_target`; the new one is seeded at
`a_target`, 3.81 Å away). So the last peptide bond is never stretched → a rigid `AllBonds`
constraint seeds and minimizes cleanly at 15 fs → **the dt-halving guard is unnecessary.**

---

## Implementation plan (step 3 — next)

Behind a flag (default off; Tutorials 12 & 13 unchanged):

1. In `protocol.py` / `core.py`, derive `a_target` and `p_target` from the ribosome (O'Brien
   restraint-energy minimization, or store the precomputed offsets `û_A·4.27`, `û_P·4.76`).
2. **Seed** residue L at `a_target` (replace `seed_positions`' `a_anchor + buffer_nm·x̂`);
   restrain L to `a_target` in stages 1–2 and to `p_target` in stage 3 (migration).
3. Build with `constraints="AllBonds"` (the last peptide bond is now unstretched, so no
   exemption is even required — unlike step 1's fallback).
4. **Delete / bypass** the dt-halving stability guard in `run_length`.
5. Validate the Tutorial-13 way: short small-`L_max`, large `scale_factor` debug run — green on
   energies + ejection, dwell/geometry vs `../12_auto/.../output/`.

---

## Step 3 — IMPLEMENTED + TESTED ✅

**Code.** New flag `ElongationParams.optimize_ptc_geometry` (default **False**, so Tutorials
12 & 13 are byte-for-byte unchanged). When on, `protocol.py` calls
`core.optimal_ptc_targets(ribo)` for `(a_target, p_target)`, seeds residue L at `a_target`
(new `seed_point` arg threaded `protocol → run_length → seed_positions`), restrains L→`a_target`
(stages 1–2) / →`p_target` (stage 3), and rederives the tunnel-wall plane to
`min(a_target.x, p_target.x)`. INI key: `optimize_ptc_geometry = yes`. Files:
[`../../topo/csp/core.py`](../../topo/csp/core.py) (`optimal_ptc_targets`, `seed_positions`,
`run_length`, `ElongationParams`), [`../../topo/csp/protocol.py`](../../topo/csp/protocol.py).

**Implemented target points** (the in-package solver picks this *deterministic* member of the
3.81-Å solution family; it differs slightly from the §Result table above — both are valid,
3.81 Å apart, clash-clear):
- `a_target = (0.8696, −0.1713, −0.2373)` nm, `p_target = (1.0052, 0.1083, −0.0169)` nm,
  `|A−P| = 3.810 Å`. Clearance: nearest non-excluded beads 4.60 Å (A) / 5.58 Å (P),
  EV ≤ 0.024 kcal/mol ≪ kT.

**Test** (`csp_step2_allbonds.ini`: L=1→5, `constraints=AllBonds`, dt=15 fs, CPU →
`synth_out_step2/`):

| run | geometry | bonds | seed bond L−1↔L | max\|PotE\| | dt-halving |
|--|--|--|--|--|--|
| **step 2 (fix)** | equilibrium | **AllBonds** | **3.79–3.81 Å** | **30.5 kJ/mol** | **never fired** ✅ |
| control | old far-seed | AllBonds | 7.30–7.41 Å | 500.8 kJ/mol | (n/a) |
| baseline | old far-seed | flexible | 7.30–7.41 Å | 501.6 kJ/mol | not in 5-res run¹ |

- All 15 stages of the fix run completed with sane energies (global max \|PotE\| = 30.5 kJ/mol
  vs the 1e9 divergence limit); the `[stability]` dt-halving guard **never triggered**.
- Mechanism confirmed: new geometry seeds every peptide bond at **3.81 Å = the AllBonds
  constraint length**; old geometry seeds it at **~7.4 Å (1.9× stretched)** → the constraint
  solver must yank the new bead ~0.36 nm → ~16× higher peak strain (500 vs 30 kJ/mol).
- Residue L−1 rests **0.04–0.11 Å from `p_target`** after its stage-3 migration — the A→P
  hand-off lands it exactly where the next bond expects it.

¹ The short L=1→5 debug run is too small to hit the rare stiff-Go-well configs that trigger
dt-halving in production (longer chains); the seed-bond-length + strain-energy gap is the direct
evidence. **TODO:** a longer `L_max` baseline to capture an actual dt-halving event for a full
before/after, then delete the guard (step 1 item 4) once the equilibrium geometry is the default
for AllBonds runs.

### Step-3 follow-ups (units, robustness, exit-side constraint)

After the first implementation, `optimal_ptc_targets` was hardened:

1. **Consistent OpenMM units.** The solver now works entirely in **nm / kJ/mol / rad** (it used
   to compute in Å / kcal internally and convert at the end). Bond lengths
   `0.427 / 0.476 / 0.381 nm`; angle/improper `k = 104.6 kJ/mol/rad²` (= 25 kcal/mol/rad²);
   excluded-volume `eps = RIBO_NC_EPS_KJ = 5.52e-4 kJ/mol` (imported from `ribosome.py`).
2. **Exit-side inequality constraint.** Each target must satisfy `x ≥ x` of its own tRNA R bead
   (the +x exit convention), so the residue lies **between the tRNA and the exit port**, never
   buried back in the ribosome. Encoded as two SLSQP `ineq` constraints
   (`A.x ≥ AtR:76@R.x`, `P.x ≥ PtR:76@R.x`).
3. **Full-sphere multistart.** Starts now cover the sphere of new-residue offset directions
   (deterministic Fibonacci sphere) instead of a planar xy-ring. The planar version missed
   out-of-plane minima and — because SLSQP is scale-sensitive — returned a **different, worse,
   partially-buried** solution under nm vs Å (A.x = 0.517 nm < tRNA, EV = 5.7 kJ/mol). The
   sphere version is unit-invariant and converges to the true optimum at both n=60 and n=100.

4. **Hard peptide bond, soft tRNA bonds.** Only `|A−P| = 0.381 nm` is a **hard** equality
   constraint — it is topo's AllBonds constraint length, so it must hold exactly. The two
   tRNA distances (0.427 / 0.476 nm) are now **soft harmonic penalties** in the objective at
   O'Brien's own `k = RESTRAINT_K_KJ = 83680 kJ/mol/nm²` (his bonds are finite-k harmonic, not
   rigid). Benefits: (a) faithful to O'Brien's model; (b) the solve stays **feasible** even on a
   ribosome where 0.427/0.476 can't be met exactly alongside the 0.381 peptide bond (the old
   all-hard version would `RuntimeError`); (c) strain redistributes instead of dumping entirely
   on the angles. `|A−RA|` / `|P−RP|` now come out *near* 0.427/0.476 rather than exactly.

Final implemented targets (essentially unchanged; now robust + feasibility-safe):
`a_target = (0.8711, −0.1652, −0.2359) nm`, `p_target = (1.0069, 0.1149, −0.0162) nm`,
`|A−P| = 0.381 nm` (hard), `|A−AtR:76@R| = 0.433`, `|P−PtR:76@R| = 0.475 nm` (soft); exit-side
satisfied (A.x 0.871 ≥ 0.828; P.x 1.007 ≥ 0.571); EV ≤ 0.09 kJ/mol ≪ kT. Re-validated: L=1→5
AllBonds run completes, dt-halving never fires, max\|PotE\| = 85 kJ/mol. Cost ≈ 17 s (once per run).

### Step-4 — now the CSP DEFAULT

The equilibrium-bond PTC geometry + `AllBonds` is now the **default** for continuous
synthesis (`topo-csp`). Implemented via `protocol.csp_default_elong()`
(`optimize_ptc_geometry=True`, `constraints="AllBonds"`), used by `CSPParams.elong` and
`read_csp_config`. Key points:

- The **shared** `ElongationParams` dataclass defaults are left as the old far-seed +
  flexible-bond values, so the elongation / Tutorial-9 cylinder paths (which subclass it as
  `CylinderParams`) are **unaffected**. Only the CSP layer overrides them.
- `read_csp_config` now keeps the CSP default when `constraints` is **absent** (it used to
  force `None`); an explicit `constraints = None` still selects flexible bonds.
- **Legacy path stays reachable**: pin `optimize_ptc_geometry = no` + `constraints = None`.
  **Tutorials 12 & 13** (the validated O'Brien reproduction) are pinned this way and run
  exactly as before — verified: no new-geometry banner, `flexible (harmonic)` bonds.
- `scipy` added to `pyproject.toml` dependencies (the default path now needs it for the
  one-time target solve; imported lazily).
- Verified: keyless INI / bare `CSPParams()` → `equil=True, constraints=AllBonds`, runs clean
  (0 dt-halving, max\|PotE\| 24 kJ/mol); pinned INIs → `equil=False, constraints=None`, old path.
