# CSP + ribosome-prep review — 2026-07-04

Scope: the continuous-synthesis subsystem `topo/csp/` (core / ribosome / kinetics /
protocol / cylinder / movie) + the new ribosome-prep asset bundle
`assets/csp/prepare_ribosome/` + the synthesis docs. Focus on **biological /
physical correctness** and flow.

Method: three focused deep passes (force-field physics; kinetics + runners; docs
+ asset scripts) plus independent spot-verification of the highest-stakes claims
against the actual source.

**How to use this doc.** Each finding has a **Decision:** line. Mark one of
`FIX` / `SKIP` / `DISCUSS` (and add notes). Hand it back and I'll implement the
`FIX`es. Nothing here is changed yet.

Confidence key: **✓ confirmed** = re-verified in source this session · **~ plausible**
= specific and credible from a deep pass, not independently re-verified.
"(my code)" = added this session in the asset bundle; everything else is
pre-existing CSP engine code.

---

## Suggested execution order (once you've marked decisions)

1. **#3** (my asset code; small, re-verify the 4 organisms) and the remaining
   **9f–9j** cleanups.
2. **Docs reorg** — as a separate pass; new pages + consolidation + cross-links.

Nothing is committed. I'll implement exactly the `FIX`-marked items and re-verify.

---

## Summary (open items, ranked)

Resolved & removed 2026-07-04: **#1** (½-convention — **Option A** applied: the two
`Harmonic{Bond,Angle}Force` tether terms now pass `2·k` so they realize the full
nominal O'Brien stiffness, matching the improper / wall / position-restraint; verified
the effective tether-bond energy doubled to nominal), **#2** (`resrange = 1-76` on
native A-site tRNAs — user decision: **no fix / won't-fix**; harmless for the current
4 organisms, revisit only if a future organism has a 77-nt A-site tRNA), **#4** (cylinder "ignores time_stage_1/2" — the full
codon dwell is run as one segment = same total time as the explicit model, so this is
correct; documented, incl. that traffic keys are inert), **#5** (cylinder omissions —
documented), **#6** (0.5 nm seed radius — not a bug; docstring + docs fixed), **9a**
(model_theory stale/incorrect claim), **9b** (stale deleted-doc code comments),
**9c** (RUNBOOK rewritten to a current maintainer guide), **9d** (dropped the duplicate
`out/` PROVENANCE copies), **9e** (fixed `view_trunc.tcl` default path). **#7** (cold-start singularity + unused `ptc_offset_nm` — resolved by the always-on PTC optimization: `ptc_offset` removed entirely, and the runner now always seeds/holds at the EV-clear `optimal_ptc_targets` point, never a raw tRNA bead), **9h** (AllBonds + far-seed foot-gun — moot: the far-seed path is removed, equilibrium seeding is always on).
Numbers below are kept stable for continuity, so there are gaps.

| # | Sev | Conf | One-liner | Where |
|---|-----|------|-----------|-------|
| 3 | Med | ✓ | `via_trna` graft matches acceptor arm by absolute resid on un-renumbered tRNAs (my code) | `gen_subunit.py:184-201` |
| 9 | Low | ~ | Grab-bag (remaining 9f,9g,9i,9j): malformed-tRNA handling, tether-angle borrowing, log/API polish, magic numbers | various |

**No high-severity crashing bugs.** Units (all kcal→kJ conversions correct), the
12-10-6 force forms, the mean-conserving 3-stage kinetics, the cylinder tunnel
potential (continuous, correct sign everywhere), the rigid-scenery interaction
groups, and the L→L+1 rebuild are all correct as written.

---

## Part 1 — Code findings

### 3. `via_trna` graft matches acceptor arm by absolute residue number (Med, ✓, my code)

**What.** `gen_subunit.superpose_via_trna` aligns the donor P-site tRNA onto the
target P-site tRNA using acceptor-arm atoms `{1..7} ∪ {66..76}`, matched by
**identical residue number**, on the **raw (un-renumbered) cif models**:

```python
accept = set(range(1, 8)) | set(range(66, 77))     # gen_subunit.py:~178
# ... bb() keeps atoms whose seqid.num is in `accept`, matched by number
```

**Why it matters.** For the yeast graft the donor is human `Pt`, whose acceptor is
at **77** in the deposition — so donor residues 66–76 are offset one nucleotide
from the true 3′ end and pair with **non-homologous** target residues across the
amino-acyl acceptor stem (the rigid region the graft must reproduce). This is why
the yeast graft fit is ~3.2 Å (markedly worse than the same-species LSU fits). It
still produces "a sane A-site anchor" (already caveated in `MANIFEST.md`), so the
shipped yeast structure is usable, but the method is fragile.

**How to fix — match by offset-from-the-3′-terminus.** In `bb()`, instead of keeping
atoms whose absolute `seqid.num` is in `{1..7} ∪ {66..76}`, select the acceptor region
**counted back from each tRNA's own 3′-terminal nucleotide** — e.g. the last ~11
residues (the CCA acceptor + the 3′ acceptor-stem strand, the part nearest the PTC
whose pose the graft must reproduce), matched by *position-from-3′-end* (0, 1, 2, …)
rather than residue number. This registers the acceptor ends of two tRNAs of **any
length**, without the artifacts of the numbering-based approaches:

- Do **not** use "renumber `max→76` then match by absolute number." A 77-nt tRNA's
  extra nucleotide sits in the body (variable/D arm), not at a terminus, so shifting
  the whole numbering to align the 3′ end pushes the 5′ acceptor strand (residues 1–7)
  out of register by one — it fixes one strand of the acceptor helix and breaks the
  other — and it introduces a residue-0 artifact. Offset-from-3′ avoids both.
- (Optional) also match the 5′ strand 1–7 counted from the 5′ start; but the 3′
  region is what actually places the A-site acceptor, so matching it alone is enough.

Then re-verify the yeast graft RMSD improves (target: well under the current ~3.2 Å).
Note: this is **independent of** the write-time `max→76` renumbering of the *output*
model (which stays — the CSP engine looks up the acceptor at `resid 76`).

**Decision:** ______

---

### 9. Grab-bag (Low)

Batchable cleanups; mark the ones you want. (9a, 9b resolved & removed 2026-07-04.)

- **9f** Inconsistent malformed-tRNA handling: `add_trna_tether` soft-skips missing
  P/BR2, `optimal_ptc_targets` hard-fails. Pick one. **Decision:** ______
- **9g** KB backbone double-Gaussian angle (basins 91.7°/130°) is reused for the
  CA–CA–ribose tether angle (`ribosome.py:482-484`) — a physically arbitrary
  borrowing of a protein-backbone potential onto a protein–RNA triple. Consider a
  dedicated angle term. **Decision:** ______
- **9i** Log/API polish: codon-label vs dwell off-by-one in the dwell logs
  (kinetics correct, O'Brien convention; the *logged* codon isn't the one governing
  the logged time); `t_invivo_total` logs the mean not the sampled sum; `L0`
  required in cylinder but optional in protocol; `p_anchor`/`a_anchor` param names in
  `run_length` are generic restraint-target/seed-anchor; euk `*_model.pdb` ATOM
  serials wrap >100k (cosmetic). **Decision:** ______
- **9j** Undocumented magic numbers: `scale_factor=4331293.0`, `time_stage_1=0.00034`,
  `time_stage_2=0.004201` (`core.py:907-909`) — the entire real↔sim ~4.3×10⁶
  timescale rests on `scale_factor` with no cited provenance. Add a source comment.
  **Decision:** ______

---

## Part 2 — Docs review + proposed layout

### What's wrong now
- **Overlap:** `usage/continuous_synthesis.md` §"Configuration reference (csp.ini)"
  and `usage/synthesis_control.rst` both document the same `csp.ini` keys → they
  drift. (Also: cylinder-model diagram duplicated in `cylinder_synthesis.md` and
  tutorial 07; the 3-stage explanation in 3 places; the "tunnel frame + renumber-to-76"
  in every PROVENANCE.)
- **Gap:** the Sphinx site never tells a reader **where the `ribosome_trunc.pdb`
  comes from** — the prepared structures + prep pipeline (`assets/csp/prepare_ribosome/`)
  are invisible to doc readers. Real break in the flow.
- **No synthesis landing page** routing cylinder-vs-explicit-vs-prep.
- **Split physics:** base force field in `model_theory.rst`, RNC additions in
  `continuous_synthesis.md` — acceptable, but the seam needs an explicit cross-link
  (the factual error there is already fixed).

### Proposed layout (organized bio → physics → prep → run → outputs)

```
The TOPO model
  usage/model_theory              base Cα Gō force field (add "→ extended
                                  for synthesis" pointer)

Co-translational synthesis
  1. usage/synthesis_overview     NEW, short — the biology (elongation cycle, exit
                                  tunnel, codon kinetics) + the fork: explicit-bead
                                  vs analytic-cylinder, WHEN to use which (incl.
                                  the cylinder-vs-explicit "not directly
                                  comparable" caveat — now in cylinder_synthesis.md)
  2. usage/synthesis_model        the PHYSICS: RNC additions (ribosome↔NC excluded
                                  volume + Debye–Hückel, tRNA tether/anchors, PTC
                                  geometry, 3-stage mechanics, tunnel wall, stability
                                  guard). Absorbs continuous_synthesis §Theory.
  3. usage/synthesis_kinetics     codon→MD-step mapping, Fluitt table, mFPT, 3-stage
                                  split, effective-timescale caveats (document 9j)
  4. usage/ribosome_preparation   NEW (the missing page) — the shipped structures
                                  (MANIFEST table + caveats) AND how to build/
                                  re-truncate one; surfaces the asset bundle
  5. usage/synthesis_control      ONE canonical csp.ini + cylinder.ini reference
                                  (merge continuous_synthesis §Config + synthesis_control)
  6. usage/synthesis_outputs      per-length/per-stage layout + movies

Tutorials → Part B
  07_translation_cylinder (analytic, simplest first)
  08_ribosome_synthesis   (explicit beads; uses a bundled structure)
```

Linear path: overview → model → kinetics → preparation → control → outputs → tutorials.

**Decision (docs):** ✅ DONE (2026-07-04) — flow-critical reorg implemented; Sphinx
builds clean (new pages render, no new warnings):
- **NEW `usage/synthesis_overview.md`** — landing page: the biology + the
  explicit-vs-cylinder fork ("which to use") with the not-directly-comparable warning.
- **NEW `usage/ribosome_preparation.md`** — surfaces the asset bundle (the missing
  "where does `ribosome_trunc.pdb` come from" page): shipped structures table +
  caveats, how to use one by path, how to re-truncate / build a new organism.
- **`index.rst`** — caption → "Co-translational synthesis"; toctree reordered:
  overview → ribosome_preparation → continuous_synthesis → cylinder_synthesis →
  synthesis_control.
- **`model_theory.rst`** — added the "→ extended for synthesis" seam pointer.
- **`tutorials/README.md`** — Part B intro links the overview + prep page.
- **Config overlap** — light-touch: `continuous_synthesis` §Config now points to
  `synthesis_control` as the canonical per-key reference.

**Deferred** (larger/riskier, not done): the full *split* of the 503-line
`continuous_synthesis.md` into separate model/kinetics pages, and the full *merge* of
its config section into `synthesis_control` (would risk destabilizing good working
docs — do as a focused follow-up if wanted).
