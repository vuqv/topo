# Step 1 — Use AllBonds constraints and remove the dt-halving guard

**Goal.** Run the `topo.csp` continuous-synthesis at 15 fs with **rigid `AllBonds`
constraints** (O'Brien's setup) and **delete the dt-halving stability guard**, without
the seeded-far-residue instability that originally forced flexible bonds + dt-halving.

This note records the geometry analysis behind the fix. Source paths:

- topo placement/restraint: [`../../topo/csp/core.py`](../../topo/csp/core.py)
  (`seed_positions`, `build_length_model`, `run_length`, the dt-halving guard).
- O'Brien reference: `…/Continuous_synthesis_protocol/continuous_synthesis_v6.py`
  (`A_site_tRNA_binding`, `create_elongation_system`).
- Structures: topo `ribosome_trunc.pdb`; O'Brien
  `../12_auto/continuous_synthesis/input/setup/50S_tRNA_cg_truncated.cor`.

---

## 1. Anchor coordinates (Å) — the two structures are NOT fully identical

| Bead | O'Brien `.cor` | topo `.pdb` | \|Δ\| |
|--|--|--|--|
| AtR:76@P | (4.431, −5.260, 0.736) | (4.431, −5.260, 0.736) | 0.000 |
| AtR:76@PU2 | (6.414, −5.252, −4.221) | (6.414, −5.252, −4.221) | 0.001 |
| **AtR:76@R** | (8.242, −5.594, −0.957) | (8.284, −5.843, −1.368) | **0.48** |
| PtR:76@P | (7.516, 4.736, −4.340) | (7.516, 4.736, −4.340) | 0.000 |
| PtR:76@PU2 | (6.367, 4.466, 2.308) | (6.367, 4.466, 2.308) | 0.000 |
| **PtR:76@R** | (9.277, 2.419, −0.959) | (5.705, 2.977, −0.588) | **3.64** |

⚠️ **Finding 1 — the P-site ribose anchor is wrong by 3.6 Å in topo.**
Every P / PU1 / PU2 bead matches O'Brien to 3 decimals, but the **R (ribose) bead of
`PtR:76` is displaced 3.64 Å** in `ribosome_trunc.pdb`. `read_anchor(..., bead="R")`
uses exactly that bead for the P-anchor (restraint target *and* the tunnel-wall plane),
so this single bad bead inflates the A→P geometry:

- A→P anchor separation: **O'Brien 8.08 Å vs topo 9.22 Å** (+1.14 Å).

---

## 2. L and L−1 placement, computed from `AtR:76@R`

**New residue L (added to the A site):**

| | formula | position (Å) |
|--|--|--|
| O'Brien | `AtR_R + 4.27·[cos10°, sin10°, 0]` | (12.448, −4.852, −0.957) |
| topo | `AtR_R + 0.4 nm·[1,0,0]` (`buffer_nm`) | (12.284, −5.843, −1.368) |

The two L placements differ by only **1.09 Å** (4.0 Å vs 4.27 Å magnitude; pure +x vs a
10° xy-tilt). **The A-site placement is essentially fine and is NOT what forces dt-halving.**

**Previous residue L−1 (held at the P side):**

- O'Brien: harmonic **bond 4.76 Å to `PtR:76@R`** + 2 angles (106°→P, 130°→PU2) +
  improper → orients L−1 down the tunnel.
- topo: in the current build L−1 is *not* re-restrained; it keeps its previous-cycle
  final coords (held at `PtR_R + ptc_offset ≈ PtR_R + 4.0 Å·x̂`). With the misplaced
  P-anchor, that target is itself 3.6 Å off.

---

## 3. ⚠️ Finding 2 — the real cause of dt-halving

It is **not** the A-site placement distance. It is that **topo keeps the peptide bond
L−1↔L present (and `AllBonds`-constrainable) while L−1 and L sit at opposite tRNA sites
~8–9 Å apart**, whereas its equilibrium length is 3.81 Å:

```
initial |L − L−1| ≈ A→P anchor separation − (a-arm + p-arm geometry)
topo seeds this bond at ~1.6× its equilibrium length
```

A **rigid `AllBonds` constraint cannot be initialized that far from its length** — the
constraint solver diverges — which is why the code falls back to flexible harmonic bonds,
and then the stiff Go-well that forms on relaxation blows up at 15 fs → dt-halving.

**O'Brien avoids this without ever stretching a constraint:**

- **Stage 1:** `system.removeConstraint(L−2)` — **deletes the L−1↔L peptide-bond
  constraint**. New AA held only by the soft 4.27 Å bond to AtR; L−1 by the soft 4.76 Å
  bond to PtR.
- **Stage 2:** removes it again and re-adds it as a **harmonic** bond at 3.81 Å
  (never a constraint).
- Every *other* backbone bond stays a rigid **`AllBonds`** constraint — and none is ever
  stretched, because residues 1..L−1 keep their relaxed coordinates.

---

## 4. The fix — AllBonds everywhere except the last peptide bond

topo already seeds residues 1..L−1 from the previous relaxed structure, so **the only
over-stretched bond in the whole system is the last peptide bond L−1↔L.** Mirror O'Brien:

1. **Build with `constraints="AllBonds"`** but **exempt the single C-terminal peptide
   bond** (L−1↔L): in `build_length_model`, skip `addConstraint` for the last bond and add
   one `HarmonicBondForce` term for it instead. No constraint is then ever seeded far from
   its length → solver is happy → 15 fs is stable → **dt-halving guard can be deleted.**
2. **Fix the `PtR:76@R` bead** in `ribosome_trunc.pdb` → O'Brien's (9.277, 2.419, −0.959).
   Restores the 8.08 Å A→P separation and the correct restraint / wall target.
3. *Optional fidelity (already on the DIFFERENCES work list):* L placement
   `4.27 Å·[cos10°, sin10°, 0]` in `seed_positions`; freeze all but the C-terminal 15
   residues in `run_length`.

**Net:** `AllBonds` on the whole nascent chain **except the one freshly-formed peptide
bond, which stays harmonic** — that single exemption is what lets a 15 fs constrained
integrator stay stable and removes the dt-halving guard.

---

## 5. Implementation plan (proposed, not yet applied)

- Add a flag (e.g. `exempt_last_bond` / `peptide_bond_flexible`) to `build_length_model`
  and thread it through `run_length` / `ElongationParams`. **Default off** so Tutorials
  12 & 13 stay byte-for-byte reproducible.
- When on + `constraints="AllBonds"`: constrain bonds 0..L−3, leave bond L−2 (the
  L−1↔L peptide bond) harmonic.
- Patch `ribosome_trunc.pdb` `PtR:76@R` coordinate (keep a copy of the original).
- Validate the Tutorial-13 way: short small-`L_max`, large `scale_factor` debug run —
  green on energies + ejection, dwell/geometry vs `../12_auto/.../output/`.
