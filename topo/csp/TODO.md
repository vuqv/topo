# topo.csp — TODO

**Status (implemented & working end-to-end):** v1 + v2 elongation runner
(`elongate.py`); build-once-subset contacts; rigid mass-0 ribosome with
ribosome↔chain excluded volume + Yukawa (`ribosome.py`); tRNA tether (bond +
`CA(L-1)-CA(L)-tRNA` orienting angle); planar tunnel wall; nascent-only output;
ejection/stallation post-phase; movie tool (`movie.py`); INI control file;
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

### 7. `ejection_steps = auto` — stop ejection once the chain clears the ribosome
Let `ejection_steps` accept an **int** (current: fixed step count) **or the string
`'auto'`**. With `auto`, run the ejection phase and **periodically check** whether the
nascent chain has fully cleared the ribosome, stopping as soon as it has (instead of a
guessed fixed length):
- Every `check_interval` steps (e.g. a new `ejection_check_every`, default ~10k–50k),
  read the C-terminus x and `max(x)` over the ribosome beads.
- **Stop** when `x(C-terminal AA) − max(x_ribosome) > cutoff` (default **2.0 nm**) — the
  C-terminus has emerged past the far (+x) edge of the ribosome, so the chain is out of
  the tunnel.
- Implementation: chunked stepping loop in the ejection branch of
  `run_continuous_synthesis` (mirror the chunked loop already used by the stability
  guard in `core.run_length`); the ribosome beads are fixed (mass-0), so `max(x_ribosome)`
  is computed once. Add a sane safety cap (max steps) so `auto` can't run forever.
- Config: parse `ejection_steps` as int-or-`auto`; add `ejection_cutoff_nm` (2.0) and
  `ejection_check_every`. (Consider the same for `dissociation_steps`.)



---

## Revision list

### Remove the per-stage step clamps for production
`max_steps_per_stage` / `min_steps_per_stage` (`CSPParams`, parsed by `read_csp_config`)
are **testing-only** — they cap/floor each stage's MD step count so the tutorials run fast,
but they **break the physical timescale mapping** (the integrated MD per stage no longer
matches the sampled dwell time). For production the step count must come entirely from the
kinetics (`scale_factor` × codon time ÷ `dt`).
- **Decide:** drop both knobs entirely (remove the `CSPParams` fields + `read_csp_config`
  parsing + the `stage_steps` clamp args), or keep them but guard/loudly warn when set so
  a production run can't silently use them. Currently documented as testing-only in
  `docs/usage/continuous_synthesis.md` and the example INIs set them for speed.

### Ribosome-traffic correction — hidden, revisit later
The per-codon **ribosome-traffic** correction (`ribosome_traffic` / `initiation_rate`
in `CSPParams`; `topo.csp.kinetics.ribosome_traffic_times` + the intrinsic-vs-real
split in `build_mfpt_lists`; stage-2 stretch in `stage_dwell_times`) is **off by
default and now HIDDEN** from the user-facing docs and the example `csp.ini` files
(2026-06-30, user request) — the code remains and is still parsed if a key is present,
but it is not advertised because it is not validated and depends on an external
`ribosome_traffic` binary that is not bundled.
- **Decide:** either (a) finish it properly — vendor/locate the `ribosome_traffic`
  binary (or reimplement the upstream-queue correction natively), validate the
  intrinsic→real stretch against a reference run, then re-expose it in the docs +
  config; or (b) drop the feature entirely (remove the `CSPParams` fields, the
  `read_csp_config` parsing, and `kinetics.ribosome_traffic_times`).
- While hidden, `real == intrinsic` everywhere and stage-2's mean is exactly
  `time_stage_2` (no traffic stretch), so results are unaffected.
- Also unhide the runtime banner line in `protocol.py` (it still prints
  `ribosome_traffic=off`) when the feature is re-exposed or removed.

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
