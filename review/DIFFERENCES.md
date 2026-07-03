# O'Brien reference vs `topo.csp` — differences & consistency plan (REVISED)

> **Revision.** This is a revised copy of
> `tutorials/14_obrien_topo_consistency/DIFFERENCES.md` (original dated 2026-06-30),
> re-checked line-by-line against the **current** `topo/csp/` source (branch
> `tut15-claude-fix`, HEAD `e9667c1`) and the O'Brien reference
> `continuous_synthesis_v6.py` on **2026-07-01**. Several items the original listed as
> open have since been closed or newly aligned; the tables below are re-tagged and a
> new **"What is still different"** summary heads the document. Companion:
> [`REVIEW_NOTES.md`](REVIEW_NOTES.md) (evidence + what changed since the original).

This document lists the **remaining differences** between the O'Brien continuous-synthesis
reference and the `topo.csp` port (what `topo-csp` runs). **Full parity is not the goal** —
some topo simplifications are deliberate. Each difference is tagged:

- 🟢 **CLOSE** — worth aligning; meaningful science impact, tractable effort.
- 🟡 **CONSIDER** — align if cheap / if it matters for a given study.
- ⚪ **KEEP** — intentional simplification; document the deviation, don't "fix" it.
- ✅ **DONE** — aligned since the original doc (new this revision).

Reference (read-only):
`/storage/home/qzv5006/programs/cg_simtk_protein_folding/Continuous_synthesis_protocol/continuous_synthesis_v6.py`
Port: `topo/csp/` (`core.py`, `protocol.py`, `ribosome.py`, `kinetics.py`, `cg_ribosome.py`).

> **Which code path is compared.** The CSP runner (`read_csp_config`) defaults
> `trna_tether = False` ([`protocol.py:654`](../topo/csp/protocol.py)), so by default `topo-csp`
> uses the **fixed-point position-restraint** path with A→P migration. The full O'Brien
> **tRNA tether** (bond + 2 angles + improper) now *exists and is selectable* via
> `trna_tether = yes` (added since the original doc — see §"What changed"). Note the raw
> `RunParams` dataclass default is `True` ([`core.py:848`](../topo/csp/core.py)); the
> runner's INI default is what governs `topo-csp`, and that is **off**.

---

## What is STILL different (answer, at a glance)

Ordered by science impact. "Path" = which restraint path it applies to.

| # | Difference | O'Brien | topo (current) | Tag | Path |
|--|--|--|--|--|--|
| 1 | **Mobility window** | only C-terminal **15** residues integrate (`id > L−15`; rest mass-0) | **entire** nascent chain 1..L mobile — *no freeze mask exists* | 🟢 CLOSE | both |
| 2 | **Ribosome L24 free loop** | L24 loop residues **42–56** mobile (`ribo_free_mask`) | ribosome **entirely** frozen (mass-0); L24 radii also per-AA-approximated | 🟡 CONSIDER | both |
| 3 | **Tunnel wall lifetime** | one-sided x-wall only in **stage 5** (dissociation) | wall on **throughout** synthesis + post (now behind `tunnel_wall` toggle, default on) | 🟡 CONSIDER | both |
| 4 | **Chain chemistry** | peptide bond explicitly **deleted → formed at 3.81 Å**; **rigid `AllBonds`** | peptide bond **always present** (full 1..L each stage); **flexible** harmonic bonds by default | ⚪ KEEP | default |
| 5 | **Restrain previous AA / orientation** | new AA **and** prev AA restrained, with orienting angles + improper (R/P/PU2 beads) | only residue **L** restrained, **no orientation** | 🟢/⚪ | default path only |
| 6 | **PtR:76 A76 P-anchor** | O'Brien source coord | topo CG coord **3.45 Å** off (source-structure diff) | ⚪ KEEP | both |
| 7 | **Placement geometry (fixed heuristic)** | 4.27 Å offset + **10° xy tilt** | *superseded*: optimized least-buried seed (equil mode), or 0.4 nm +x (default, no tilt) | 🟡 CONSIDER | both |

Item **5** is **closed in the `trna_tether` path** (which now restrains prev-AA in stage 1
and adds the 2 angles + improper) but still open in the default position-restraint path.

Everything else the original doc tracked is now **✅ DONE** (see next section).

---

## What changed since the original (2026-06-30) — now ✅ DONE

These were the largest sources of divergence; all have been aligned and are **no longer
differences** (some were not even in the original tables):

- ✅ **NC↔ribosome excluded volume** — the biggest fix. topo previously used a pure
  `(σ/r)^12` term ~10³–10³× softer than O'Brien, the root cause of nascent–ribosome
  clashes. Now O'Brien-**consistent**: the **12-10-6** form
  `U = ε[13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶]` with the **sum rule**
  `R_ij = Rmin/2_i + Rmin/2_j` and matched ε ([`ribosome.py:252`](../topo/csp/ribosome.py)).
  Nascent Rmin/2 = **per-residue Karanicolas–Brooks** (Option A default, reproduces
  O'Brien's `.prm` A-values; per-AA `SA..SY` is the Option B fallback). Hard clashes
  36→0 at full length. *(commits `040a064`, `a526da1`, `83dfd6a`, `40f54e3`)*
- ✅ **Standalone O'Brien-faithful ribosome** — CG ribose `R` bead fixed to the 5-carbon
  centroid `C1'–C5'` (was `C1'–C4'+O4'`), reproducing O'Brien's `.cor` beads to **median
  0.0003 Å**; RNA per-type radii/charges **exact**. No `.cor/.psf/.prm` needed at runtime.
  *(commits `e69c229`, plus the standalone-migration `CHANGELOG.md`)*
- ✅ **Full O'Brien tRNA tether + orientation** — bond + **2 angles + improper**, A-site
  in stages 1–2, P-site in stage 3, previous-AA restrained in stage 1 — implemented behind
  the `trna_tether` flag ([`ribosome.py:373`](../topo/csp/ribosome.py),
  [`protocol.py:198`](../topo/csp/protocol.py)). This is the original worklist item #2.
  Position restraint stays the CSP default. *(commit `3d7a7c6`)*
- ✅ **Seeding consistency** — `optimal_ptc_targets` now minimizes the **same** 12-10-6
  sum-rule EV as the simulation, placing the A/P seed at the least-buried point of the
  real wall. A-target clearance 2.37→4.27 Å; 4c5c dt-halving 263/534→0/0.
  ([`core.py:111`](../topo/csp/core.py)) *(commits `2c01dec`, `7bda1c6`)*
- ✅ **Flexible-bond force constant** — `bond_force_constant` 20920→**41840** kJ/mol/nm²
  (= O'Brien 50 kcal/mol/Å²). Moot under `AllBonds`, corrects flexible-bond runs.
  *(commit `75d4c88`)*

Unrelated but noted: stage-2 minimize is now skipped (`minimize_override=False`,
[`protocol.py:379`](../topo/csp/protocol.py)) — a topo-side speedup, not an O'Brien
alignment; and CSP now runs on GPU where available (commit `18af082`).

---

## Already consistent (unchanged from original)

- **3 stages per residue**, **dwell-time kinetics** (`steps = t·1e9/scale_factor/dt`),
  **new-residue anchor** `AtR:76@R`, **restraint k = 200 kcal/mol/Å²** (`restraint_k = 83680`),
  **A→P migration intent**, **timestep 15 fs**, **post-synthesis ejection + dissociation**.

---

## 1. Remaining differences (detail)

### Mobility / what moves each stage — 🟢 CLOSE (STILL OPEN)

| | O'Brien | topo | Tag |
|--|--|--|--|
| Mobile nascent atoms | only C-terminal **15** residues (`id > L−15`, ref lines 254-258); rest mass-0 | **entire** chain 1..L mobile — no freeze mask ([`core.py`](../topo/csp/core.py), searched: none) | 🟢 CLOSE |
| Ribosome | frozen except L24 loop **42–56** (`ribo_free_mask`, ref line 37) | entirely frozen (mass-0, [`core.py:1012`](../topo/csp/core.py)) | 🟡 CONSIDER |

Highest value/effort ratio and **still not done**: O'Brien integrates only the C-terminal
15 residues per step (co-translational relaxation is localized + cheap); topo moves the
whole chain. Add an optional "freeze all but C-terminal *N* residues" mass-0 mask in
`run_length` (default *N*=15 to match; "all" keeps current behaviour). *(Original doc said
the L24 loop was 42–**59**; the reference is 42–**56** — line 37.)*

### C-terminus restraint — split verdict

| | O'Brien | topo default (position restraint) | topo `trna_tether=yes` | Tag |
|--|--|--|--|--|
| Mechanism | bond + 2 angles + improper to **tRNA beads** | isotropic spring `k·Σ(x−x0)²` to a **fixed point** | bond + 2 angles + improper (matches) | ✅ in tether path |
| Orientation control | yes | **none** | yes | 🟢 CLOSE (default only) |
| Residues restrained | new AA **and** prev AA | only residue **L** | new AA + prev AA (st1) | ✅ in tether path |
| Anchor sub-beads | `R`, `P` (R−1), `PU2` (R+2) | `R` only | `R` (+ orienting angle) | 🟡 CONSIDER |

The most physically meaningful divergence (orientation control) is **closed when the flag
is on**; the validated **default** path still omits it. Decision to keep position restraint
as default is intentional (A↔P migration is simpler and validated).

### Chain chemistry & integration — ⚪ KEEP

| | O'Brien | topo | Tag |
|--|--|--|--|
| Peptide bond L−1↔L | explicit: deleted (st1) → formed at **3.81 Å** (st2), k=200 | **always present** (full 1..L each stage) | ⚪ KEEP |
| Bond treatment | rigid `AllBonds` constraints | flexible harmonic bonds (default) | ⚪ KEEP |
| Stability at 15 fs | constraints absorb stiff Go wells | dt-halving guard (`run_length`, up to 6 halvings) | ⚪ KEEP |

Deliberate (flexible bonds + guard needed to seed the far-placed new residue) and validated
full-length. **New option** since the original: `optimize_ptc_geometry = yes` +
`constraints = AllBonds` in the INI gets much closer to O'Brien — rigid constraints with the
new equilibrium-bond seeding (0.381 nm targets) so `AllBonds` seeds cleanly with 0 dt-halving.
Still KEEP as **default**, but the rigid path is now available.

### Tunnel wall — 🟡 CONSIDER (STILL OPEN)

| | O'Brien | topo | Tag |
|--|--|--|--|
| One-sided x-wall | only in stage 5 (dissociation), k=20 kcal/mol/Å², x0=`x_eject−2` (ref lines 819-828) | on **throughout** synthesis + post; k=20 kcal/mol/Å²; plane auto-derived from structure ([`ribosome.py:468`](../topo/csp/ribosome.py), [`core.py:1046`](../topo/csp/core.py)) | 🟡 CONSIDER |

topo's persistent wall is a pragmatic anti-leak guard for the truncated ribosome (harmless
during synthesis — the chain only pushes +x). Now gated by the `tunnel_wall` INI toggle
(default on) and its plane is auto-derived (no stale knob). Set `tunnel_wall = no` if strict
stage-5-only fidelity is wanted.

### Placement geometry (stage 1) — 🟡 CONSIDER (superseded, not matched)

| | O'Brien | topo | Tag |
|--|--|--|--|
| Offset distance | 4.27 Å (ref lines 238-240) | 0.4 nm (4.0 Å) `buffer` +x (default), **or** optimized seed (equil mode) | 🟡 CONSIDER |
| Off-axis tilt | **10°** in xy (`[cos10°, sin10°, 0]`) | none (pure +x default); equil mode solves 3-D least-buried point | 🟡 CONSIDER |

topo's `optimize_ptc_geometry` mode replaces O'Brien's fixed 4.27 Å + 10° tilt **heuristic**
with a solved least-buried seed point on the real wall (same intent: clear excluded volume).
The plain default still uses 0.4 nm +x with no tilt.

### Accepted source deviations — ⚪ KEEP

- **PtR:76 A76 P-anchor** 3.45 Å from O'Brien (genuine source-structure difference; PTC
  seeding absorbs it). Splice O'Brien's single A76 R coord in if bit-exact P-site ever needed.
- **L24** radii approximated per-AA rather than O'Brien's per-residue B-types (mean protein
  Rmin/2 error 0.013 nm).

---

## Prioritized work list (best-effort, updated)

1. 🟢 **C-terminal mobility window** — the last high-value 🟢 not yet done. Optional mass-0
   freeze of all but the last *N*=15 nascent residues in `run_length`.
2. ✅ ~~Per-stage tRNA tether with orientation~~ — **done** (`trna_tether` flag).
3. 🟡 **Placement geometry** — largely superseded by `optimize_ptc_geometry`; only add the
   fixed 4.27 Å + 10° tilt if reproducing O'Brien's *exact* seed is required.
4. 🟡 **Restrain previous AA in the default path** (already done in the tether path).
5. 🟡 **Ribosome L24 free loop (42–56)** and traffic correction — only if a study needs them.
6. ⚪ **Keep + document**: flexible bonds + stability guard, always-bonded chain, persistent
   tunnel wall (now toggleable), A76 / L24 accepted deviations.

Validate each the Tutorial-13 way: short debug run (small `L_max`) green on energies +
ejection, then compare dwell/geometry to the reference bundle.
