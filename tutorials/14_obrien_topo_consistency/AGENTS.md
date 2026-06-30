# AGENTS.md — Tutorial 14 orientation for AI agents

> **Read me first.** This file tells an agent **what this tutorial is for**, **where the
> reference code and the topo port live**, and **how its files are organized**. It is the
> Tutorial-14 analogue of [`../13_validate_claude_fix12/AGENTS.md`](../13_validate_claude_fix12/AGENTS.md).

---

## 1. What this tutorial is

Tutorial 14 is a **consistency-improvement** task: bring the `topo.csp` continuous-synthesis
port **closer to** the O'Brien reference (`continuous_synthesis_v6.py`), **best-effort**.

**The goal is NOT full parity.** Several topo simplifications are deliberate and validated
(flexible bonds + the Tutorial-13 stability guard, the always-bonded chain, A→P migration via
a moving restraint target). The aim is to close the *meaningful* gaps and explicitly **keep +
document** the rest. The gap analysis and the prioritized plan are in
[`DIFFERENCES.md`](DIFFERENCES.md) — **start there.**

## 2. Where the code is

**O'Brien reference** (read-only — semantics only, never run as the deliverable, never modify):

| Path | What it is |
|------|-----------|
| `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` | O'Brien lab CG model + protocol suite (Jiang, Nissley, O'Brien). |
| `…/Continuous_synthesis_protocol/continuous_synthesis_v6.py` | **The protocol being matched.** Key fns: `run_elongation`, `elongation`, `A_site_tRNA_binding`, `peptide_bond_formation`, `translocation_AtR`, `create_elongation_system` (per-stage bonds/angles/impropers to tRNA beads). |
| `…/Continuous_synthesis_protocol/ribosome_traffic` | Upstream-queue delay binary. |
| `…/CG_protein_parameterization/` | How the reference CG protein `.psf/.top/.prm` were built. |
| `…/CG_ribosome_parameterization/` | **How the reference CG ribosome/tRNA was built** — the analogue of topo's [`../../topo/csp/cg_ribosome.py`](../../topo/csp/cg_ribosome.py). Key file: `create_cg_ribosome_model.py` (all-atom → CG bead mapping: protein→CA, RNA→P/R/BR; `gen_50S_pdb.py`, `truncate_ribosome.py`, `gen_ribosome_FF.py`). Compare against `cg_ribosome.py` — they should map an identical all-atom input to identical beads. **They don't (see below).** |

### CG ribosome mapping: O'Brien `create_cg_ribosome_model.py` vs topo `cg_ribosome.py`

Both map the **same all-atom ribosome** (protein→1 CA bead/residue; RNA→`P`/`R`/`BR`
beads). Bead-by-bead they agree **except the ribose `R` bead atom set**, which is the
tRNA anchor (`AtR:76@R`, `PtR:76@R`) the nascent-chain bond attaches to:

| Ribose `R` centroid | atoms used |
|--|--|
| O'Brien (`create_cg_ribosome_model.py:326`) | `C1' C2' C3' C4' C5'` (incl. exocyclic **C5'**, **no** ring O) |
| topo (`cg_ribosome.py:64`, `RIBOSE_RING`) | `C1' C2' C3' C4' O4'` (incl. ring **O4'**, **no** C5') |

(Protein CA, `P`, and all base-ring `BR` beads — purine 6-ring `N1 C2 N3 C4 C5 C6`,
5-ring `C4 C5 N7 C8 N9`, pyrimidine `N1 C2 N3 C4 C5 C6` — are **identical** in both.)

This single difference is the root cause of the anchor mismatch in
[`step1_allbonds_no_dt_halving.md`](step1_allbonds_no_dt_halving.md) §1:
- `AtR:76@R` differs by **0.48 Å** — exactly the O4'↔C5' swap (one of 5 centroid atoms
  moved ~2 Å → ~0.4–0.5 Å centroid shift). **Consistent with the atom-set difference alone.**
- `PtR:76@R` differs by **3.64 Å** — ~7× too large for the swap; **P-site residue 76**
  (3'-CCA terminal A, the peptidyl carrier) has an **additional** problem (missing/extra
  ribose atoms in the all-atom input, or a different residue selected). Verify residue 76's
  ribose atom completeness in the all-atom `PtR` input before trusting that bead.

**Fix:** align `RIBOSE_RING` in `cg_ribosome.py` to O'Brien's `C1' C2' C3' C4' C5'` (drop
`O4'`, add `C5'`) and regenerate `ribosome_trunc.pdb`; then re-check `PtR:76@R` specifically.

**topo port** (this is what you edit to improve consistency):

| Path | Role |
|------|------|
| [`../../topo/csp/core.py`](../../topo/csp/core.py) | Per-length/per-stage MD: `run_length`, `seed_positions` (placement), `add_cterm_restraint`, the dt-halving stability guard, mass/freeze logic. |
| [`../../topo/csp/protocol.py`](../../topo/csp/protocol.py) | The 3-stage outer loop; **forces `trna_tether=False`** (uses the position-restraint path). `.ini` → params. |
| [`../../topo/csp/ribosome.py`](../../topo/csp/ribosome.py) | Rigid ribosome scenery; `add_trna_tether` (bond+angle to `PtR:76@R` — exists but unused by CSP), `add_tunnel_wall`. |
| [`../../topo/csp/kinetics.py`](../../topo/csp/kinetics.py) | FPT/dwell-time sampling, ribosome-traffic correction. |
| [`../../topo/csp/DESIGN.md`](../../topo/csp/DESIGN.md) | Design rationale + invariants. |

**Reference *run* to validate against** (inputs + outputs) — in Tutorial 12:
[`../12_auto/continuous_synthesis/`](../12_auto/continuous_synthesis/)
(`input/cont_synth_ecoli.cntrl`, `input/setup/`, `output/info.log`, `output/output/1.out`).

## 3. This folder's file structure

| File | Role |
|------|------|
| `DIFFERENCES.md` | **The deliverable so far:** full O'Brien-vs-topo comparison + a tagged (🟢 CLOSE / 🟡 CONSIDER / ⚪ KEEP), prioritized consistency plan. |
| `AGENTS.md` | This orientation file. |

As the work proceeds, add (suggested): a `PLAN.md`/`NOTES.md` log of what was changed and
why, a `csp.ini` + inputs if you run validation here, and a `synth_out*/` for outputs (write
**only** under `synth_out*/`).

## 4. How to work this tutorial

1. Read [`DIFFERENCES.md`](DIFFERENCES.md); pick the **next unticked item** from its
   "Prioritized work list".
2. Make the smallest change to `topo/csp/*` that closes that gap, **behind a flag/default**
   so the existing (validated) behavior stays reachable.
3. Validate the Tutorial-13 way: a short debug run (small `L_max`, large `scale_factor`) green
   on energies + ejection, then compare dwell/geometry to the reference bundle in
   `../12_auto/continuous_synthesis/output/`.
4. Record the change + result in your notes; re-read the work list; continue or stop (the goal
   is best-effort, not exhaustive).

## 5. Safety rules

- The reference source in `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` is
  **read-only** — consult it; do not modify or run it as the deliverable.
- **Never overwrite** the reference data under `../12_auto/continuous_synthesis/`. Write run
  outputs only under a local `synth_out*/`.
- Changes to `topo/csp/*` are real package edits — keep them flagged/back-compatible and
  re-validate (Tutorials 12 & 13 must still pass).
