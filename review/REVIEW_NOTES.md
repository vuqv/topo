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

**Still different**

- **Mobility window** — O'Brien freezes all but C-terminal 15 (`id > L-15`, ref lines
  254-258). topo: grep for `setParticleMass`/freeze/window in `core.py` → none on the
  nascent chain; whole chain 1..L integrates. **Open.**
- **L24 free loop** — O'Brien `ribo_free_mask = L24 : 42-56` (ref line 37; original doc
  mis-stated 42–59). topo: ribosome appended mass-0 with no free mask
  (`core.py:1012`). **Open.**
- **Tunnel wall lifetime** — O'Brien adds the x-wall only at `stage == 5` (ref lines
  819-828). topo `add_tunnel_wall(..., range(L))` called every `run_length`
  (`core.py:1046-1053`), default `tunnel_wall = True` (`core.py:859`). **Open (toggleable).**
- **Chain chemistry** — O'Brien `removeConstraint` then `addBond(..., 3.81 Å, k=200)`
  (ref lines 770-772), `constraints=AllBonds` (ref line 644). topo builds full 1..L each
  stage; `constraints` default `None` = flexible (`core.py:810`). **KEEP by design.**
- **Prev-AA / orientation (default path)** — `add_cterm_restraint` restrains only `L-1`
  (`core.py:578-618`); no angles/improper. **Open in default path only.**

**Closed / newly aligned**

- **NC↔ribosome EV** — `ribosome.py:252-260` 12-10-6 + sum rule, matched ε; nascent K-B
  radii threaded `precompute_contacts → run_length → append_ribosome`. **Done.**
- **tRNA tether** — `ribosome.py:373-466` bond + 2 angles + improper; `protocol.py:198,
  354, 388` wire A-site st1-2 / P-site st3, prev-AA st1. Off by runner default
  (`protocol.py:654-655`). **Done (opt-in).**
- **Seeding** — `optimal_ptc_targets` (`core.py:111-241`) minimizes the same
  `ε[13x¹²−18x¹⁰+4x⁶]` with `R_pair = aa_rmin_2 + Rmin/2_ribo`. **Done.**
- **Ribosome geometry** — see `CHANGELOG.md §1` (C5′ set, median 0.0003 Å vs `.cor`). **Done.**

## Nuance worth flagging

- `ElongationParams.trna_tether` **dataclass default is `True`** (`core.py:848`), but the
  CSP **runner** (`read_csp_config`) defaults it to **False** (`protocol.py:654-655`). What
  `topo-csp` actually runs is therefore the position-restraint path. Direct callers of
  `run_length`/`ElongationParams` who don't go through `read_csp_config` would get the tether
  path unless they set it — a minor footgun worth a comment or a unified default.
- The original's "placement geometry" 🟡 item is now better described as *superseded*: topo's
  `equil_peptide_geometry` solves for a least-buried seed rather than matching O'Brien's fixed
  4.27 Å + 10° tilt. Same intent, different mechanism — so "still different" but not a
  regression to close.

## Files produced

- `DIFFERENCES.md` (repo root) — verbatim copy of the 2026-06-30 original (baseline).
- `review/DIFFERENCES.md` — revised, current-state version.
- `review/REVIEW_NOTES.md` — this file.
