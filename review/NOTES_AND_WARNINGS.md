# Review notes & warnings — tick when addressed

Every observation, warning, and recommendation I raised while reviewing `topo/csp/` against
the O'Brien reference and sweeping the tree for TODOs (2026-07-01, branch `tut15-claude-fix`,
HEAD `e9667c1`). Tick `[x]` as each is handled. Sources: this review + companions
[`DIFFERENCES.md`](DIFFERENCES.md), [`REVIEW_NOTES.md`](REVIEW_NOTES.md).

Severity: 📌 = decision needed · 📝 = doc/accuracy · 🔧 = suggestion.

---

## Decisions needed (📌 — "finish or delete" / pick one)

- [x] ✅ **DECIDED (2026-07-09) — flexible exit-tunnel loop is topo-only.** The L24/L26
  free-loop capability (`ribo_free_mask` / `ribo_free_pdb`) is implemented in `topo.csp`
  (`append_flexible_l24_loop`, commit `30d6162`), reproducing O'Brien's mobile-L24 setup
  with topo's own structure-based native-contact model (`build_nonbonded_interaction`).
  Deliberately **not** ported to **cosmo**: cosmo uses the `hps_kr` hydropathy-scale model
  with a different excluded-volume / non-bonded basis and has no `build_nonbonded_interaction`,
  so topo's structure-based Go loop would be physically inconsistent there. **topo-only by
  design — do not re-attempt the cosmo port.**

- [ ] 📌 **Ribosome-traffic correction: finish or delete.** Off by default and hidden since
  2026-06-30, but the code is still carried in `RunParams`, `read_csp_config`,
  `kinetics.ribosome_traffic_times`, and a `ribosome_traffic=off` banner. It depends on an
  external binary that is not bundled (exits 127). Decide: (a) vendor/reimplement +
  validate + re-expose, or (b) remove the fields/parsing/function/banner.

---

## Doc / accuracy notes (📝)

- [ ] 📝 **Tether item is closed only in the tether path.** "Restrain previous AA +
  orientation control" is done when `trna_tether = yes`, still absent in the default
  position-restraint path — the split verdict must be kept when quoting status.

---

## Suggestions for improvement (🔧 — beyond the existing TODOs)

Forward-looking ideas I noticed that aren't tracked anywhere yet. Grouped; tick if adopted
or consciously declined.

### Science
- [ ] 🔧 **Report a confidence interval on dwell/Rg from replicas.** Current D6 ratios
  (1.01×, 1.06×) are single-run point estimates; O'Brien runs ~50 trajectories/protein.
  Once the replica driver exists, quote mean ± CI so "matches reference" is statistically
  grounded.
