# O'Brien reference vs `topo.csp` — differences & consistency plan

This document lists the **remaining differences** between the O'Brien continuous-synthesis
reference and the `topo.csp` port (what `topo-csp` runs, as in Tutorials 12 & 13), and a
**best-effort** plan to close them. **Full parity is not the goal** — some topo simplifications
are deliberate. Each difference is tagged:

- 🟢 **CLOSE** — worth aligning; meaningful science impact, tractable effort.
- 🟡 **CONSIDER** — align if cheap / if it matters for a given study.
- ⚪ **KEEP** — intentional simplification; document the deviation, don't "fix" it.

> **Which code path is compared.** `topo-csp` forces `trna_tether = False`
> ([`topo/csp/protocol.py:191`](../../topo/csp/protocol.py)), so it uses the **fixed-point
> position-restraint** path with A→P migration — *not* the `add_trna_tether` bond+angle code
> (which exists but is unused by CSP).

Reference (read-only):
`/storage/home/qzv5006/programs/cg_simtk_protein_folding/Continuous_synthesis_protocol/continuous_synthesis_v6.py`
Port: [`topo/csp/`](../../topo/csp/) (`core.py`, `protocol.py`, `ribosome.py`, `kinetics.py`).

---

## Already consistent (not repeated below)

These already match and are **dropped** from the difference tables:

- **3 stages per residue** (A-site/peptidyl-transfer → translocation → tRNA-binding).
- **Dwell-time kinetics** — FPT sampling, `steps = t·1e9/scale_factor/dt` (`kinetics.stage_steps`).
- **New-residue anchor bead** — `AtR:76@R` in both.
- **Restraint force constant** — k = 200 kcal/mol/Å² (`restraint_k = 83680`).
- **A→P migration intent** — equivalent (O'Brien switches the tRNA the bond attaches to;
  topo switches the restraint target point a→a→p).
- **Timestep** — 15 fs.
- **Post-synthesis phases** — ejection + dissociation (free runs) exist in both.

---

## The 3-stage cycle, side by side (adding residue L)

Both models split each elongation into 3 timed sub-stages with the **same dwell-time
partition**, but implement the per-stage geometry very differently. Indices are 0-based
(residue *n* = index *n*−1).

| | Stage 1 | Stage 2 | Stage 3 |
|--|--|--|--|
| **Biology** | A-site tRNA binding / peptidyl transfer | peptide-bond formation / translocation begins | translocation A→P / wait for next tRNA |
| **MD dwell (both)** | sampled `time_stage_1` (peptidyl-transfer) | sampled `time_stage_2` (translocation) | sampled remainder (next-tRNA-binding) |

### O'Brien `continuous_synthesis_v6.py`

| Aspect | Stage 1 `A_site_tRNA_binding` | Stage 2 `peptide_bond_formation` | Stage 3 `translocation_AtR` |
|--|--|--|--|
| New residue L | placed at `AtR:76@R + 4.27 Å·[cos10°,sin10°,0]`, inserted as C-terminus | (kept) | (kept) |
| L held to | **A-site tRNA**: bond 4.27 Å (k=200) + 2 angles + improper | **A-site tRNA**: bond 4.27 Å + angles + improper | **P-site tRNA**: bond 4.76 Å + angles + improper |
| L−1 held to | **P-site tRNA**: bond 4.76 Å + angles + improper | (released; via peptide bond) | (interior) |
| Peptide bond L−1↔L | **removed** (not yet bonded) | **formed at 3.81 Å** (k=200) | native chain bond |
| Mobile atoms | C-terminal 15 residues only | C-terminal 15 residues only | C-terminal 15 residues only |
| Procedure | minimize → MD `step_stage_1` | minimize → MD `step_stage_2` | minimize → MD `step_stage_3` → save `_stage_3_final.cor` (seeds next L) |

### topo `topo.csp` (position-restraint path, `trna_tether=False`)

| Aspect | Stage 1 (peptidyl transfer) | Stage 2 (translocation) | Stage 3 (tRNA binding) |
|--|--|--|--|
| New residue L | seeded at `a_anchor + buffer` (0.4 nm +x); residues 1..L−1 from prev stage's final | continue (`seed_override` = stage-1 final) | continue (`seed_override` = stage-2 final) |
| L held to | **position restraint** to `a_target` (= A-anchor + `ptc_offset`)¹ | **position restraint** to `a_target` | **position restraint** to `p_target` (= P-anchor + offset) → A→P migration |
| L−1 held to | not restrained | not restrained | not restrained |
| Peptide bond L−1↔L | **always present** (full 1..L model built each stage; no formation event) | always present | always present |
| Mobile atoms | **entire** nascent chain 1..L | entire chain | entire chain |
| Procedure | build 1..L + minimize → MD `s1` (dt-halving stability guard) | minimize → MD `s2` | minimize → MD `s3` → final coords seed next L |

¹ For the cold-start length (L0) stage 1 restrains to `p_target` instead. Each stage is a
separate `run_length` call; only the current C-terminus (residue L) is ever restrained.

**Net effect — same outcome, different route.** O'Brien walks the anchor
**P → (new at A) → peptide bond → A → P** by re-attaching tRNA bonds and explicitly forming
the peptide bond, moving only the last 15 residues. topo seeds the new residue at A and lets a
**moving position restraint** pull it A→A→P over the three stages, with the whole chain mobile
and the peptide bond always present. After one cycle, residue L occupies the P-site slot that
L−1 held, and the chain has extruded one register down the tunnel.

---

## 1. Remaining differences

### Placement geometry (stage 1)

| | O'Brien | topo | Tag |
|--|--|--|--|
| Offset distance | 4.27 Å | `buffer` = 0.4 nm (4.0 Å) | 🟡 CONSIDER |
| Off-axis tilt | 10° in xy (`[cos10°, sin10°, 0]`) | none (pure +x) | 🟡 CONSIDER |

Trivial change to `seed_positions` ([`core.py:322`](../../topo/csp/core.py)): use 0.427 nm + the 10° tilt.

### C-terminus restraint

| | O'Brien | topo | Tag |
|--|--|--|--|
| Mechanism | bond + 2 angles + improper to **tRNA beads** | isotropic position spring `k·Σ(x−x0)²` to a **fixed point** | 🟢 CLOSE (partial) |
| Orientation control | yes (angles + improper aim the chain down the tunnel) | none in the used path | 🟢 CLOSE |
| Anchor beads | `R`, `P` (R−1), `PU2` (R+2) | only the `R` coordinate | 🟡 CONSIDER |
| Restraint target distance | bond eq. A 4.27 Å, P 4.76 Å | target = anchor + `ptc_offset` (~4 Å +x) | 🟡 CONSIDER |
| Residues restrained | new AA **and** previous AA | only residue L | 🟡 CONSIDER |

🟢 `add_trna_tether` ([`ribosome.py:283`](../../topo/csp/ribosome.py)) already does a bond + one
orienting angle to `PtR:76@R` — it just isn't used by CSP. Closing this = make CSP use a
per-stage tRNA tether (A-site beads in stages 1–2, P-site in stage 3) + add the 2nd angle +
improper. Restores **orientation control**, the most physically meaningful divergence. Do it
behind a flag so the position-restraint path stays available.

### Chain chemistry & integration

| | O'Brien | topo | Tag |
|--|--|--|--|
| Peptide bond L−1↔L | explicit: deleted (st1) → formed at 3.81 Å (st2) | always present (full 1..L model each stage) | ⚪ KEEP |
| Bond treatment | rigid `AllBonds` constraints | flexible harmonic bonds | ⚪ KEEP |
| Stability at 15 fs | constraints absorb stiff Go wells | dt-halving stability guard (`run_length`) | ⚪ KEEP |

Flexible bonds + the guard are **deliberate** (needed to seed the far-placed new residue) and
validated full-length in Tutorial 13. Do **not** convert to rigid constraints — keep + document.

### Mobility / what moves each stage

| | O'Brien | topo | Tag |
|--|--|--|--|
| Mobile nascent atoms | only C-terminal **15** residues (`id > L−15`); rest mass-0 | **entire** chain 1..L mobile | 🟢 CLOSE |
| Ribosome | frozen except L24 loop 42–59 (`ribo_free_mask`) | entirely frozen | 🟡 CONSIDER |

🟢 Big, often-overlooked difference: O'Brien samples only the C-terminal 15 residues per step;
topo moves the whole chain. Changes both the physics (co-translational relaxation during
synthesis) and the cost. Add an optional "freeze all but C-terminal N residues" mask in
`run_length` (default N=15 to match; "all" keeps current behavior).

### Tunnel wall

| | O'Brien | topo | Tag |
|--|--|--|--|
| One-sided x-wall | only in stage 5 (dissociation), k=20 kcal/mol/Å², x0=`x_eject−2` | on throughout synthesis + post, k=20 kcal/mol/Å² | 🟡 CONSIDER |

topo's persistent wall is a pragmatic anti-leak guard for the truncated ribosome (harmless
during synthesis — the chain only pushes +x). Gate it to the post phase only if strict fidelity
is wanted.

---

## Prioritized work list (best-effort)

In order — stop whenever the rest are "intentional / not worth it":

1. 🟢 **C-terminal mobility window.** Optional mass-0 freeze of all but the last `N` nascent
   residues (default 15). Highest value/effort ratio.
2. 🟢 **Per-stage tRNA tether with orientation.** Switch CSP to the bond+angle(+improper)
   linkage (A-site beads st1–2, P-site st3); behind a flag.
3. 🟡 **Placement geometry.** 0.427 nm offset + 10° tilt in `seed_positions`. Trivial.
4. 🟡 **Restrain previous AA too** and use the `P`/`PU2` sub-beads for the angle/improper.
5. 🟡 **Ribosome L24 free loop** and **traffic correction** — only if a study needs them.
6. ⚪ **Keep + document**: flexible bonds + stability guard, always-bonded chain, persistent
   tunnel wall.

Validate each the Tutorial-13 way: short debug run (small `L_max`) green on energies +
ejection, then compare dwell/geometry to the reference bundle in
[`../12_auto/continuous_synthesis/output/`](../12_auto/continuous_synthesis/output/).
