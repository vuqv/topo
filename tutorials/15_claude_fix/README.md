# Tutorial 15 — `claude_fix`: equilibrium-PTC + rigid `AllBonds` continuous synthesis

Validates **"the claude fix"** — the **equilibrium-bond PTC geometry + rigid `AllBonds`**
elongation path (Tutorial 14, steps 2–4) — by reproducing O'Brien's continuous-synthesis
results on **4c5c** (then extending to **P0CX28**), on the *unpinned* fix path (no
`equil_peptide_geometry = no`, no `constraints = None`), and closing the O'Brien-consistency
gaps selected in [`AGENTS.md`](AGENTS.md) §1b.

See [`AGENTS.md`](AGENTS.md) for the full directive, [`NOTES.md`](NOTES.md) for the decision
log + numbers, and [`TASK.md`](TASK.md) for the ticked checklist.

## The fix in one paragraph

Each new residue is seeded **exactly one peptide bond (3.81 Å) from the current C-terminus**
at the optimal PTC target points (`core.optimal_ptc_targets`), so the always-present peptide
bond is born at its equilibrium length and a **rigid `AllBonds`** build seeds/minimizes
cleanly at 15 fs — the **dt-halving stability guard never fires**. This is more physically
faithful to O'Brien (whose bonds are rigid constraints) than the legacy far-seed + flexible
+ dt-halving path of Tutorials 12/13.

## Reproduce

From this folder (env: `bioenv` — OpenMM 8.2, scipy, MDAnalysis; a CUDA GPU for the full run):

```bash
# 1. Debug smoke test (L=1->8, CPU, ~1 min): confirms the fix path is ACTIVE + stable.
topo-csp -f csp.ini                 # -> synth_out_debug/
#    expect: "[equil_peptide_geometry] optimal PTC restraint targets (|A-P| = 0.3810 nm)",
#            all stages "rigid (AllBonds)", seed peptide bond 3.810 A, ZERO "[stability]" lines.

# 2. Full-length validation (L=1->306, GPU, ~30 min) + post-synthesis ejection.
topo-csp -f csp_val.ini             # -> synth_out/

# 3. Extended free-MD ejection demo (egress along +x; ~4 min GPU).
python eject_demo.py -f csp_val.ini --steps 500000 --device GPU   # -> synth_out/ejection_long/

# 4. Validation report (D5 energies, D5b egress/clash, D6 vs the O'Brien reference).
python analyze_validation.py        # reference run read from ../12_auto/ (read-only, L=1->10)

# (optional) §1b feature A/B: the O'Brien tRNA tether + orientation (debug).
topo-csp -f csp_tether.ini          # -> synth_out_tether/   (trna_tether = yes)
```

Outputs (`synth_out*/`) are git-ignored and safe to regenerate. Raw inputs and the
`../12_auto/` reference are never written.

## Results (4c5c)

Full run: equilibrium-PTC + `AllBonds`, L=1→306, `scale_factor` 216564650 (50× production;
dwell times logged in seconds are scale-independent), random_seed 20240629.

| Check | Result |
|-------|--------|
| **D2** fix active | optimal PTC targets \|A−P\| = **0.3810 nm**; seed peptide bond **3.810 Å**; all stages rigid `AllBonds` |
| **D3** completes | "Done. Synthesized 1 → 306", exit clean; **0 dt-halving lines** over 919 stages |
| **D5** finite energy | worst max\|PotE\| = **1.48e3 kJ/mol** ≪ 1e12 |
| chain topology | corr(residue index, x) = **−0.926** — threads the tunnel, N-term extruded (x≈109 Å), C-term at PTC (x≈10 Å); **no leak** through the truncation (no bead x<0) |
| **D5b** wall | **PASS** — min nascent x 8.37 Å vs wall 8.71 Å |
| **D5b** egress | **PASS** — on release C-term +12.0 Å, CoM +22.8 Å along +x (7.5 ns); no collapse back |
| **D5b** clash | residual 2.2–2.9 Å nascent–23S interpenetration (soft-EV model property; finite energy) |
| **D6** dwell | per-codon in-vivo dwell matches ref; **total ratio 1.01×** (≪ 2× tol) |
| **D6** geometry | final R_g (L=10) **0.799 vs 0.750 nm = 1.06×** |

**Fix vs legacy (12/13).** The fix reproduces the reference at least as well (1.01× dwell,
1.06× R_g) while being more physically faithful — rigid bonds, **no stability guard ever
firing** (the legacy far-seed path peaked ~500 kJ/mol from a 1.9×-stretched seed bond and
relied on dt-halving).

## §1b O'Brien-consistency features (validate-first)

Of the five selected features, only the tRNA tether is applicable to topo's deliberate
simplifications (see [`NOTES.md`](NOTES.md) for the full assessment + numbers):

| # | Feature | Outcome |
|---|---------|---------|
| ✅2 | tRNA tether + orientation (bond + 2 angles + improper) | **Implemented** behind `trna_tether = yes` (A-site st1–2, P-site st3). Stable, kinetics identical; O'Brien-faithful but does not reduce the clash (holds the C-term at the tRNA, closer to the PTC). Position restraint stays the validated default. |
| ✅3 | restrain previous AA (L−1) | **Implemented** (stage-1 P-site tether; pairs with ✅2). |
| ✅1 | C-terminal mobility window | **Deferred** — mass-0 freeze is forbidden with `AllBonds`, and freezing the bulk breaks topo's *diffusion-driven* extrusion (no explicit register translocation). |
| ✅4 | ribosome L24 free loop | **Deferred** — topo's ribosome is scenery-only (no intra-ribosome FF to keep a freed loop physical); the clashes are 23S, not L24. |
| ✅5 | 10° off-axis tilt | **N/A** — superseded by equil-PTC seeding (the seed already sits at the optimal off-axis target). |

## The residual clash (model finding)

The lone unmet strict criterion is the **2.2–2.9 Å nascent–23S interpenetration** as the
chain threads the narrow truncated tunnel. It is a property of O'Brien's deliberately-soft
excluded volume (ε = 0.000132 kcal/mol, bead radii ~7 Å) at *full-length* tunnel packing —
energy stays finite and the chain otherwise respects the ribosome (threads the lumen, no leak
through the truncation faces, clean egress). At the only length with a reference (L=10), ours
(3.36 Å) and O'Brien (4.57 Å) are both clash-free; the sub-3 Å contacts appear only as the
chain fully packs the tunnel. The §1b feature that could target it (mobility window) is
incompatible with topo's extrusion route. → reported as a model property, not a bug.

## Files

- `csp.ini` — debug profile (L=1→8) → `synth_out_debug/`.
- `csp_val.ini` — full-length profile (L=1→306) → `synth_out/`.
- `csp_tether.ini` — §1b ✅2/✅3 A/B (debug, `trna_tether = yes`) → `synth_out_tether/`.
- `eject_demo.py` — extended free-MD ejection (egress demo) → `synth_out/ejection_long/`.
- `analyze_validation.py` — D5 / D5b / D6 report (reference read from `../12_auto/`).
- `P0CX28/` — the second target (L=1→106); D6 N/A (no O'Brien reference).
- Raw inputs (never overwritten): `4c5c_model_clean.pdb`, `*_stride.dat`, `4c5c_mrna.txt`,
  `trans_times.txt`, `ribosome_trunc.pdb`, `domain.yaml`.
