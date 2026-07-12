# Review notes — DIFFERENCES.md revision (2026-07-01)

Companion to [`DIFFERENCES.md`](DIFFERENCES.md). Records how the revision was produced,
the code evidence behind each verdict, and the delta from the original.

## Scope / method

- **Original**: `tutorials/14_obrien_topo_consistency/DIFFERENCES.md`, dated 2026-06-30 09:30.
  Copied verbatim to the repo root as `DIFFERENCES.md` (baseline for diffing).
- **Revised**: `review/DIFFERENCES.md` — re-checked against current source on branch
  `tut15-claude-fix`, HEAD `e9667c1`.
- **Evidence gathered** by reading both codebases in full:
  - O'Brien: `.../Continuous_synthesis_protocol/continuous_synthesis_v6.py` (1697 lines).
  - topo: `topo/csp/{core.py, protocol.py, ribosome.py, kinetics.py, cg_ribosome.py}` +
    `CHANGELOG.md` (standalone-ribosome migration log).

## Why the original is out of date

The original was written before the following commits landed (all `topo/csp/`):

| Commit | What it aligned |
|--|--|
| `040a064`, `a526da1` | NC↔ribosome EV → O'Brien 12-10-6 + sum rule (root cause of clashes) |
| `83dfd6a`, `40f54e3` | nascent Rmin/2 → per-residue Karanicolas–Brooks (Option A default) |
| `e69c229` + CHANGELOG | standalone O'Brien-faithful ribosome (C5′ ribose bead) |
| `3d7a7c6` | full O'Brien tRNA tether (bond + 2 angles + improper) behind `trna_tether` |
| `2c01dec`, `7bda1c6` | seeding consistency (`optimal_ptc_targets` minimizes the sim EV) |
| `75d4c88` | flexible-bond force constant 20920→41840 kJ/mol/nm² |
| `a2b8519` | stage-2 minimize skip (topo speedup, not an O'Brien alignment) |

## Verdict evidence (current source)



## Nuance worth flagging

- `RunParams.trna_tether` **defaults to `False`** (`core.py`), unified with the CSP
  **runner** (`read_csp_config`, `protocol.py`). Position restraint is therefore the
  default on every path — `topo-csp` and direct `run_length`/`RunParams` callers agree.
  (Was a footgun: the dataclass used to default `True` while the runner defaulted `False`;
  fixed 2026-07-09.)
- The original's "placement geometry" 🟡 item is now better described as *superseded*: topo's
  `optimize_ptc_geometry` solves for a least-buried seed rather than matching O'Brien's fixed
  4.27 Å + 10° tilt. Same intent, different mechanism — so "still different" but not a
  regression to close.

## Files produced

- `DIFFERENCES.md` (repo root) — verbatim copy of the 2026-06-30 original (baseline).
- `review/DIFFERENCES.md` — revised, current-state version.
- `review/REVIEW_NOTES.md` — this file.
