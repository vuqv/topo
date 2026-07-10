# O'Brien reference vs `topo.csp` — differences & consistency plan (REVISED)

> **Revision.** This is a revised copy of
> `tutorials/14_obrien_topo_consistency/DIFFERENCES.md` (original dated 2026-06-30),
> re-checked line-by-line against the **current** `topo/csp/` source (branch
> `tut15-claude-fix`, HEAD `e9667c1`) and the O'Brien reference
> `continuous_synthesis_v6.py` on **2026-07-01**. Several items the original listed as
> open have since been closed or newly aligned; the tables below are re-tagged and a
> new **"What is still different"** summary heads the document. Companion:
> [`REVIEW_NOTES.md`](REVIEW_NOTES.md) (evidence + what changed since the original).
>
> **Update 2026-07-09.** The **ribosome L24 free loop is now implemented**
> (`ribo_free_mask` / `append_flexible_l24_loop`); the **C-terminal mobility window** is a
> **deliberate keep** (topo moves the whole chain — the window will not be adopted); and the
> `trna_tether` dataclass default is unified to `False` (position restraint on every path).

This document lists the **remaining differences** between the O'Brien continuous-synthesis
reference and the `topo.csp` port (what `topo-csp` runs). **Full parity is not the goal** —
some topo simplifications are deliberate. Each difference is tagged:

- 🟢 **CLOSE** — worth aligning; meaningful science impact, tractable effort.
- 🟡 **CONSIDER** — align if cheap / if it matters for a given study.
- ⚪ **KEEP** — intentional simplification; document the deviation, don't "fix" it.
- ✅ **DONE** — aligned since the original doc (new this revision).

Reference (read-only):
`/storage/home/qzv5006/group/programs/cg_simtk_protein_folding/Continuous_synthesis_protocol/continuous_synthesis_v6.py`
Port: `topo/csp/` (`core.py`, `protocol.py`, `ribosome.py`, `kinetics.py`, `cg_ribosome.py`).

> **Which code path is compared.** The CSP runner (`read_csp_config`) defaults
> `trna_tether = False` ([`protocol.py:654`](../topo/csp/protocol.py)), so by default `topo-csp`
> uses the **fixed-point position-restraint** path with A→P migration. The full O'Brien
> **tRNA tether** (bond + 2 angles + improper) now *exists and is selectable* via
> `trna_tether = yes`. The `RunParams` dataclass default is `False` (unified with the
> runner), so position restraint is the default on every path.

---

## What is STILL different (answer, at a glance)

Ordered by science impact. "Path" = which restraint path it applies to.

| # | Difference | O'Brien | topo (current) | Tag | Path |
|--|--|--|--|--|--|
| 1 | **Mobility window** | only C-terminal **15** residues integrate (`id > L−15`; rest mass-0) | **entire** nascent chain 1..L mobile — no freeze mask (deliberate) | ⚪ KEEP | both |
| 2 | **Tunnel wall lifetime** | one-sided x-wall only in **stage 5** (dissociation) | wall on **throughout** synthesis + post (now behind `tunnel_wall` toggle, default on) | 🟡 CONSIDER | both |
| 3 | **Chain chemistry** | peptide bond explicitly **deleted → formed at 3.81 Å**; **rigid `AllBonds`** | peptide bond **always present** (full 1..L each stage); **flexible** harmonic bonds by default | ⚪ KEEP | default |
| 4 | **Restrain previous AA / orientation** | new AA **and** prev AA restrained, with orienting angles + improper (R/P/PU2 beads) | only residue **L** restrained, **no orientation** | 🟢/⚪ | default path only |
| 5 | **PtR:76 A76 P-anchor** | O'Brien source coord | topo CG coord **3.45 Å** off (source-structure diff) | ⚪ KEEP | both |
| 6 | **Placement geometry (fixed heuristic)** | 4.27 Å offset + **10° xy tilt** | *superseded*: optimized least-buried seed (equil mode), or 0.4 nm +x (default, no tilt) | 🟡 CONSIDER | both |

Item **4** is **closed in the `trna_tether` path** (which now restrains prev-AA in stage 1
and adds the 2 angles + improper) but still open in the default position-restraint path.

---

## 1. Remaining differences (detail)

### Mobility / what moves each stage — ⚪ KEEP

| | O'Brien | topo | Tag |
|--|--|--|--|
| Mobile nascent atoms | only C-terminal **15** residues (`id > L−15`, ref lines 254-258); rest mass-0 | **entire** chain 1..L mobile — no freeze mask (deliberate) | ⚪ KEEP |

topo moves the whole nascent chain rather than only O'Brien's C-terminal 15 residues — a
**deliberate keep** (whole-chain mobility is topo's design; the C-terminal window will not be
adopted). *(The ribosome L24 free loop O'Brien also mobilizes is matched — see `ribo_free_mask`
/ `append_flexible_l24_loop`.)*

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
- **L24 (frozen scenery)** collision radii (Rmin/2) approximated per-AA rather than O'Brien's
  per-residue B-types (mean protein Rmin/2 error 0.013 nm). A loop freed with `ribo_free_mask`
  instead uses **structure-derived per-residue** Rmin/2 (same convention as the mobile nascent
  chain) for both its native contacts and its excluded volume — no per-AA approximation there.

---

## Prioritized work list (best-effort, updated)

1. 🟡 **Placement geometry** — superseded by `optimize_ptc_geometry` (settled); only add the
   fixed 4.27 Å + 10° tilt if reproducing O'Brien's *exact* seed is ever required.
2. 🟡 **Restrain previous AA in the default path** (already done in the tether path).
3. ⚪ **Deliberate keeps** (do not "fix"): whole-chain mobility (no C-terminal window),
   flexible bonds + stability guard, always-bonded chain, persistent tunnel wall (toggleable),
   A76 / L24-radii accepted deviations.

Validate each the Tutorial-13 way: short debug run (small `L_max`) green on energies +
ejection, then compare dwell/geometry to the reference bundle.
