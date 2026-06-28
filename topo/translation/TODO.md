# topo.translation — TODO

**Status (implemented & working end-to-end):** v1 + v2 elongation runner
(`elongate.py`); build-once-subset contacts; rigid mass-0 ribosome with
ribosome↔chain excluded volume + Yukawa (`ribosome.py`); tRNA tether (bond +
`CA(L-1)-CA(L)-tRNA` orienting angle); planar tunnel wall; nascent-only output;
ejection/stallation post-phase; movie tool (`make_movie.py`); INI control file;
Tutorial 7 + Sphinx docs. Implementation details are in `DESIGN.md` / `FILES.md` /
`README.md` — this file tracks only what is **still open**.

---

## Open — validation & production

### 1. Validate the C-terminus tether  *(highest priority)*
Whether the tether (bond + `CA-CA-tRNA` angle) actually improves extrusion vs. a
plain position restraint is **not established**. (An earlier "N-term 2.7→4.1 nm"
claim was a single noisy snapshot and was retracted.)
- Run **independent replicas** (≥5–10, different seeds) for `trna_tether = no` vs
  `yes`.
- Use a **production-relevant dwell** — the PTC collapse only shows at long dwell, so
  the 1000-step demo can't test it. Use `n_steps` ≳ 50k (ideally 840k); truncate
  `L_max` (e.g. 40) to keep cost down.
- **Robust metric**, averaged over frames *and* replicas (not one final snapshot):
  per-residue mean radial distance from the tunnel (x-)axis and mean x; classify
  "in-tunnel-extended" vs "emerged-folded". Decide if the tether keeps the in-tunnel
  segment extended/aimed +x.
- Outcome: keep the tether on by default only if it measurably helps; else flip the
  default and document.

### 2. Realistic production parameters (incl. dwell time)
The example INIs use tiny demo values. Define/validate a **production** `elongate.ini`:
- **Dwell:** `n_steps = 840_000` (mean 12.6 ns at dt = 0.015 ps; from ~20 aa/s +
  O'Brien scaling). *(This is the "set realistic dwell time" item.)*
- `L_max` blank (full length), `device = GPU`.
- **Many independent trajectories** (O'Brien run 50/protein) — add a small replica
  driver/loop or document how to launch them.
- Benchmark wall-clock per length on GPU for the ~4,600-particle v2 system and
  document expected runtime.

### 3. Analysis layer  *(DESIGN §6 phase 4 — none of this exists yet)*
- **Q vs. length:** fraction of native contacts (reuse
  `topo.analysis.native_contacts`) per domain on each `L_<L>/traj.dcd`, plotted vs.
  L → co-translational folding curve / folding order.
- **Ejection:** does the released chain clear the tunnel (C-term x ≥ threshold) and
  fold? Define an ejection-time observable (steps from tether release to exit).
- Wire a worked example into Tutorial 7.

### 4. Full-length / threading validation run
One full `L0 → N_full` run (P0CX28) at a realistic dwell, checking the chain
**threads the tunnel** (radial distance vs. length stays small in-tunnel; the
N-terminal domain folds *outside*). End-to-end confidence check for v2 geometry.
*(Can be combined with #1/#2.)*

### 5. tunnel_wall move to 5.8nm in ejection phase
---

## Open — model extensions

### 5. Variable (per-codon) elongation schedule + timescale mapping
Replace the constant `n_steps` with a per-residue dwell drawn from an exponential
distribution whose mean is the codon's decoding time (Fluitt–Viljoen), scaled to a
12.6 ns overall mean. Needs a codon→dwell table input and per-length `n_steps`.

### 6. Restart / resume across lengths  *(DESIGN §4)*
Skip lengths already completed (their `L_<L>/traj_final.pdb` exists) and resume an
interrupted length from its checkpoint. Each length is a fixed-size system (no live
resizing) — mostly bookkeeping in `run_elongation`.



---

## Revision list

### Restore rigid `AllBonds` for the elongation runner
v1/v2 currently use flexible harmonic bonds (`constraints = None`) instead of the
package-default rigid `AllBonds`, because the new residue is seeded *at the A-anchor*
while restrained to the P-anchor (~0.9–1.1 nm away): a rigid distance constraint on
the new C-terminal bond can't be seeded that far off (constraint solver / minimizer
diverges, observed E→2.3e44). Flexible bonds absorb the stretch.
- **Proposed fix:** make the *placement* point and the *restraint/tether* target sit
  **one CG bond length (0.381 nm) apart** — derive placement from
  `P-anchor + 0.381 nm toward the A-anchor direction`, not the raw A-anchor
  coordinate. Then the new bond starts at equilibrium and rigid `AllBonds` can be
  seeded directly (matching the rest of TOPO, allowing the larger 15 fs step safely).
- Revisit *with* the v2 geometry (real ribosome beads + tether), since the tether
  bond length (0.476 nm) and excluded-volume clearance also constrain the placement.
