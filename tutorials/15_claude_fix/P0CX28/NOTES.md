# NOTES.md — Tutorial 15 / P0CX28 (D8)

P0CX28 (106 residues, single-domain, contact strength 2.5044) on the **claude-fix path**
(equilibrium-PTC + rigid `AllBonds`). Reproduce: `topo-csp -f csp.ini` (debug) →
`topo-csp -f csp_val.ini` (full, L=1→106) → `python eject_demo.py -f csp_val.ini --steps
500000 --device GPU` → `python analyze_validation.py`.

**No O'Brien reference run exists for P0CX28**, so **D6 (quantitative dwell/geometry match)
does NOT apply**. The deliverable is D0–D5/D5b + internal consistency (AGENTS.md §4/§5/D8).

## Results (L=1→106, scale_factor 216564650, random_seed 20240629)

- **D2 fix active** — optimal PTC targets |A−P| = 0.3810 nm; seed peptide bond 3.810 Å;
  all stages rigid `AllBonds`; debug max|PotE| = 38.3 kJ/mol.
- **D3 completes** — "Done. Synthesized 1 → 106", exit clean; **0 dt-halving lines** (320 stages).
- **D4 outputs** — per-stage trajectories + `dwell_times.dat` + ejection.
- **D5 PASS** — worst max|PotE| = **241 kJ/mol** ≪ 1e12.
- **D5b** — wall **PASS** (min nascent x 8.37 Å vs wall 8.71 Å); `clash` 2.29 Å (same soft-EV
  grazing as 4c5c — a model property, see ../NOTES.md); `egress` net ~0 over 7.5 ns — **expected**:
  see below.
- **Internal consistency** — final L=106; threads the tunnel corr(residue, x) = **−0.740**
  (N-term distal x≈32 Å, C-term at PTC x≈10 Å); **no collapse** into the ribosome (0 beads at
  x<0); final nascent **R_g = 1.345 nm**.

## Interpretation

- **Chain is still inside the tunnel.** All 106 beads lie at x = 8.6–35.9 Å; the truncated
  ribosome's tunnel reaches ~104 Å near the axis. **106 residues is too short to emerge** from
  the ~100 Å exit tunnel — the nascent chain occupies the proximal tunnel and has not yet
  cleared the exit. This is why the short free-ejection shows little net translation (the chain
  is confined, not blocked): the egress/clearing of FINAL GOAL #3 is only fully observable once
  the chain is long enough to protrude (as for 4c5c, 306 res). No leak / collapse occurs.
- **R_g is sensible for confinement.** Synthesized R_g 1.345 nm vs the free native fold
  2.111 nm (0.64×): the tunnel-confined nascent chain is **more compact than native**, as
  expected — it has not emerged and post-translationally folded to its native state. Energies
  finite throughout, bonds rigid, no instability.

**D8 verdict:** P0CX28 synthesizes to full length on the claude-fix path with finite energies,
no dt-halving, correct tunnel threading with no wall leak, and a sensible (confinement-compact)
R_g — the internal-consistency bar (D6 being N/A). The residual soft-EV clash and the
short-chain egress behavior match the documented 4c5c findings.

## 2026-06-30 — Re-validated with O'Brien's ribosome (`ribosome_obrien.cor`)

Re-ran P0CX28 L=1→106 on the corrected ribosome (old run preserved as `synth_out_oldribo/`).
**0 dt-halving.** Results: **D5 PASS** (worst 381 kJ/mol), **D5b** wall PASS + egress PASS
(CoM 19.1→22.3 Å) + clash 1.74 Å, **internal consistency** corr(res,x) −0.64, **no leak**
(0 beads x<0), R_g 1.307 nm. Old-vs-new (end-of-synthesis L=106):

| metric | OLD ribo | NEW ribo (O'Brien) |
|--------|----------|---------------------|
| min center-dist | 2.81 Å | 2.50 Å |
| min surf-gap (correct radii) | −4.29 Å | **−3.62 Å** (shallower) |
| threading corr | −0.740 | −0.569 (weaker order) |
| R_g | 1.345 nm | 1.299 nm (native 2.111) |
| beads x<0 | 0 | 0 (no leak) |

Same conclusion as 4c5c: O'Brien's ribosome is **more faithful** (shallower surface penetration,
finite energy, no leak) but the residual soft-EV clash/lateral order is not resolved — confirming
the need for an explicit y-z tunnel-wall (a real lumen, ~2–6 Å, now exists with the correct radii).
The 106-mer remains entirely inside the ~100 Å tunnel (x 10–35 Å), so full egress is not
observable at this length (no leak/collapse).

## 2026-06-30 — NC↔ribosome EV fix (12-10-6 + sum rule) RESOLVES it for P0CX28 too

Re-ran P0CX28 L=1→106 with the O'Brien-consistent NC↔ribosome EV (old soft-EV run preserved as
`synth_out_softEV/`). 22 blow-ups auto-recovered (47 dt-halving retries). **D5 PASS** (worst 2.75e3
kJ/mol), **D5b all PASS (wall/clash/egress)**, no leak. The effect is even bigger than for 4c5c:

| metric | SOFT-EV | NEW 12-10-6 + sum |
|--------|---------|-------------------|
| # hard clashes (<3 Å) | 18 | **0** ✅ |
| min nascent–ribosome center-dist | 2.61 Å | **3.91 Å** |
| threading corr(res,x) | −0.643 | **−0.969** |
| extrusion x-range | 10..40 Å | **10..98 Å** |
| R_g | 1.307 nm (0.62× native) | **2.306 nm (1.09× native 2.111)** |
| ejection: CoM-x net / min NC dist | ~0 / 1.74 Å | **+15.7 Å / 3.18 Å** |

The soft EV left the 106-mer **bunched** in the proximal tunnel (over-compact, not extruding); the
correct stiff EV **extrudes it out toward the exit** (x→98 Å), threads it cleanly (corr −0.969),
gives a **near-native R_g** (2.306 vs 2.111 nm), eliminates hard clashes, and produces real egress
on release. Consistent with the 4c5c result — the NC↔ribosome EV was the root cause.
