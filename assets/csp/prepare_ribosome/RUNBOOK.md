# Runbook — prepare standard oriented, truncated ribosomes for TOPO

**You are an autonomous coding agent. Read this runbook top to bottom, then work
end-to-end from it until every acceptance test in §7 passes — do not stop partway.**
It is self-contained: it gives you the goal, the exact pipeline, what already exists
in `topo`, the two steps you must build, the per-organism landmark definitions, and
the acceptance tests that decide when you are done.

Read `topo/csp/FILES.md` and `topo/csp/DESIGN.md` **in full first** — they define
the CG bead convention, the segID convention, and how the truncated ribosome is
consumed downstream. Do not re-litigate that design.

---

## 1. Goal

TOPO's co-translational synthesis engine (`topo.csp`) needs a **rigid,
truncated, coarse-grained ribosome** whose exit tunnel is placed on a known axis.
Today `topo` only ships one such structure — *E. coli* from PDB **4V9D** — and it
was oriented/truncated by hand outside the repo. We want a **reproducible,
scripted pipeline** that turns a raw ribosome deposition into a TOPO-ready
truncated CG ribosome, and we want **standard reference ribosomes for the most
studied organisms**: *E. coli*, *S. cerevisiae* (yeast), *N. crassa*, and
*H. sapiens* (human).

**Deliverable:** for each organism, a validated
`<tag>_cg_trunc.pdb` under `sandbox/prep_rib/out/<organism>/`, produced by a
scripted, re-runnable pipeline, plus the two new pipeline stages that `topo`
is currently missing (**gen** and **orient**), plus a short provenance note per
organism recording every choice you made (PDB IDs, chain map, landmark residues).

---

## 2. The full pipeline (5 stages)

For every organism the pipeline is the same five stages. **Stages 4 and 5 already
exist in `topo` — reuse them, do not reimplement.** You must build stages 2 and 3.

| # | Stage | Status | Module |
|---|-------|--------|--------|
| 1 | **Acquire** raw structure(s) (mmCIF) from the PDB | you script the download | — |
| 2 | **Gen** large subunit + A-/P-site tRNAs → one all-atom PDB with segIDs | **BUILD THIS** | `sandbox/prep_rib/gen_subunit.py` |
| 3 | **Orient** — PTC→origin, tunnel→+x, tRNA tails→y | **BUILD THIS** | `sandbox/prep_rib/fix_orient.py` |
| 4 | **Coarse-grain** to TOPO beads (protein→Cα; RNA→P/R/BR) | **DONE** | `python cg_ribosome.py` |
| 5 | **Truncate** around the exit tunnel (cylinder + exit half-space) | **DONE** | `python truncate_ribosome.py` |

The two existing stages assume the input is **already oriented** the way stage 3
produces it: **exit tunnel on the X-axis, PTC near x = 0, exit toward +x, radial
distance `d = sqrt(y² + z²)`**. Stage 3 is what makes stages 4–5 valid — get it
right.

> **Provenance of the current E. coli file.** The working structure
> `topo/csp/structures/4v9d_50S_PtR_5jte_AtR_model.pdb` = the **50S from 4V9D** +
> the **P-site tRNA from 4V9D (segID `PtR`)** + an **A-site tRNA grafted from PDB
> 5JTE (segID `AtR`)**, then oriented as in stage 3. Your *E. coli* run must
> **reproduce this file** (that is the regression test in §7). The A-site tRNA is
> grafted from 5JTE because 4V9D lacks a well-resolved A-site tRNA; keep that graft.

---

## 3. Reference procedure (O'Brien lab)

We **adopt** the O'Brien-lab procedure and **implement it ourselves** (same policy
as the rest of `topo.csp`). The two originals are, for reference only:

- `gen_50S_pdb.py` — control-file driven (`parmed`). Inputs: `cif_file`,
  `sub_unit` (space-separated segment names, e.g. `PtR 23S 5S L2 … L36`),
  `chain_id` (the source chain letters in the deposition), `new_chain_id` (the
  reassigned output chains). It fetches the components of the **large subunit**
  (E. coli 50S: 23S + 5S rRNA + proteins L2…L36; yeast 60S: 25S + 5.8S + 5S +
  L2…L43) plus the P-site tRNA, and writes one PDB.
- `fix_orein_50S_pdb.py` — aligns **PTC `23S:A2602@N6` to the origin**, then the
  **exit-tunnel vector `23S:A2602@N6 → uL24:51@N` along the x-axis**, then the
  **tRNA-tail vector `centroid(AtR:76 ribose ring) → centroid(PtR:76 ribose ring)`
  parallel to the y-axis**, in that order (`parmed` + `pdbfixer`).

You do **not** need their code and should not depend on `parmed`/`pdbfixer` unless
you want to; `MDAnalysis` or `Biopython` + `numpy` is enough and matches the rest
of the repo. Match their **geometry and segID convention**, not their code.

---

## 4. Stage 2 — `gen_subunit.py` (BUILD)

Write `sandbox/prep_rib/gen_subunit.py`: a small, config-driven script that
assembles the large subunit + tRNAs into **one all-atom PDB** with the TOPO
**segID convention** (cols 73–76): `23S`/`25S`/`5S`/`5.8S` for rRNAs, `L2`…`L##`
for large-subunit proteins, `PtR` for the P-site tRNA, `AtR` for the A-site tRNA.

Requirements:

1. **Input:** a per-organism config (INI/JSON/YAML — your choice, keep it human
   editable) giving the source `cif` file(s), and for each biological molecule to
   keep: its **source chain ID** in the deposition, the **output segID**, and
   (optionally) the **output chain letter**. See the table in §6.
2. **A-site tRNA graft:** the config may point the `AtR` molecule at a **different
   source structure** than the main assembly (e.g. E. coli `AtR` comes from 5JTE).
   Support grafting a chain from a second cif and tagging it `AtR`. Do **not**
   re-superpose it onto the tunnel here — orientation is stage 3's job, and stage
   3 uses the A-/P-tRNA tails as a landmark, so the graft only needs to be in a
   biologically sane A-site pose (as deposited in 5JTE relative to the shared
   subunit). If the deposition frames differ, superpose the donor's large subunit
   onto the target's before lifting `AtR` (document if you do).
3. **Keep** only standard protein + RNA residues; **drop** ions/water/ligands
   (stage 4 drops them anyway, but keep the intermediate clean).
4. **Keep ONLY the ribosome scenery** — the large-subunit rRNAs + large-subunit
   proteins + the A-/P-site tRNAs. Several depositions are captured mid-translation
   or as complexes and carry extra molecules that are **not** part of the rigid
   ribosome and must be **dropped**: a **bound nascent-chain / nascent-peptide
   segment** in the exit tunnel (e.g. 7R81 is a translating ribosome — remove the
   nascent chain), mRNA, elongation/release factors, antibiotics, and the small
   subunit if the entry is a full 70S/80S (extract only the 50S/60S). Whitelist by
   the segIDs you assign in the config; anything not whitelisted is dropped. A
   leftover nascent-chain bead inside the tunnel would corrupt every downstream
   synthesis run.
5. **Write** `out/<organism>/<tag>_model.pdb`. The segID column must be populated
   for every atom — stage 4 (`cg_ribosome.py`) and everything downstream key off
   it.

Keep it dependency-light and re-runnable. Print a per-segID atom/residue count so
a human can eyeball that nothing important was dropped.

---

## 5. Stage 3 — `fix_orient.py` (BUILD)

Write `sandbox/prep_rib/fix_orient.py`: read the all-atom model from stage 2 plus
the organism's **four landmarks**, apply a single rigid-body transform (translation
+ rotation, **no scaling, no reflection**), and write `out/<organism>/<tag>_oriented.pdb`.

Landmarks (all are single atoms or centroids; §6 gives the organism-specific IDs):

- `P_PTC`  = PTC atom (E. coli `23S:A2602@N6`)
- `P_exit` = tunnel-exit atom (E. coli `uL24:51@N`)
- `P_A`    = centroid of the A-site tRNA residue-76 ribose ring (`AtR:76`)
- `P_P`    = centroid of the P-site tRNA residue-76 ribose ring (`PtR:76`)

Transform, in order:

1. **Translate** so `P_PTC → (0,0,0)`.
2. **Rotate** so the tunnel vector `v_x = P_exit − P_PTC` lies along **+x**.
3. **Rotate about the now-fixed x-axis** so the tRNA-tail vector
   `v_y = P_P − P_A` lies in the **xy-plane pointing along +y** (project `v_y`
   onto the yz-plane and rotate that projection onto +y; the x-component of `v_y`
   is untouched so step 2 is preserved). z is then fixed by right-handedness.

Implementation notes:

- Build the rotation explicitly (e.g. `numpy`); verify `det(R) = +1` (a proper
  rotation, no reflection) and `RᵀR = I` before applying.
- The **ribose ring** for a nucleotide is the sugar carbons; match stage 4's
  convention — TOPO's `R` bead centroid is `C1' C2' C3' C4' O4'` (see
  `topo/csp/FILES.md`). Use the same atom set for the `:76` centroid so the
  landmark and the CG bead agree.
- After transforming, **assert**: `|P_PTC| ≈ 0`, `v_x` is `(+,0,0)` within
  tolerance, and `v_y` is `(≈0,+,≈0)` within tolerance. Fail loudly otherwise.

---

## 6. Per-organism landmark table

The *E. coli* and *S. cerevisiae* landmark atoms are **given and authoritative**
(E. coli also defines the regression test). The **N. crassa** and **H. sapiens**
rows carry **suggested** landmark atoms you **must verify** — ribosomal-protein
numbering and rRNA numbering differ between organisms, so the uL24-homolog residue
and the PTC adenine must be confirmed by structural homology, not assumed. Even for
the given rows, sanity-check that the named atoms exist in your chosen deposition
(chain/segID and numbering can vary by entry) and that the `PTC → exit` vector
threads the tunnel.

| Organism | Large subunit | Main cif | A-site tRNA donor | PTC atom | Tunnel-exit atom | Notes |
|----------|---------------|----------|-------------------|----------|------------------|-------|
| *E. coli* | 50S | **4V9D** | **5JTE** (`AtR`) | `23S:A2602@N6` | `uL24:51@N` | **Authoritative — reproduce the shipped file.** |
| *S. cerevisiae* | 60S | **6Q8Y** (O'Brien used this) | confirm/graft | `25S:A2971@N6` | `L26:91@N` | **Landmarks given.** Yeast uL24 homolog = ribosomal protein **L26** (segID `L26`); large subunit also has 5.8S rRNA. |
| *N. crassa* | 60S | **7R81** (2.7 Å 80S translating, cycloheximide-arrested) | **graft** — 7R81 has P-site tRNA but the A-site is likely empty (cycloheximide arrest) | homolog of A2971 on **25S** — verify by superposition (fungal; expect close to the yeast site) | **likely `L26:91@N`** (same as yeast; uL24 homolog = **L26**) — verify | Full 80S — extract the **60S** + P-site tRNA; graft A-site tRNA. Fungal 25S + 5.8S like yeast, so map landmarks from the *S. cerevisiae* row. **7R81 contains a bound nascent-chain segment in the tunnel — DROP it** (keep only 60S rRNAs/proteins + tRNAs); it is not part of the rigid ribosome. |
| *H. sapiens* | 60S | **8G61** (preferred: 2.94 Å 80S decoding complex, has A- **and** P-site tRNA + mRNA bound) or 4UG0 / 6OLE / 6EK0 — confirm | **likely none** — 8G61 already has A-/P-site tRNAs; graft only if a chosen deposition lacks one | **likely `28S:A4548@N6`** (homolog of A2602) — verify | uL24 homolog residue — verify | 8G61 is a full 80S; extract the **60S** large subunit + its bound tRNAs. Pick a complete, high-resolution deposition; record why. |

**How to confirm the suggested landmarks (do this, don't guess — applies to
N. crassa and H. sapiens; use it to sanity-check yeast too):**

1. **PTC adenine:** the PTC A2602 is universally conserved. Superpose the target
   large-subunit rRNA onto the *E. coli* 23S PTC region (or use a published
   alignment) and take the **structurally equivalent adenine's `N6`**. Sanity-check
   that it sits at the base of the tunnel, between the tRNA CCA ends.
2. **Tunnel-exit atom:** identify the **uL24 homolog** (universal name `uL24`;
   beware the *different* protein `eL24` — do not confuse them). Take the residue
   whose backbone `N` sits at the **cytosolic tunnel vestibule** (the O'Brien
   choice is E. coli uL24:51). Verify the `PTC → exit` vector actually threads the
   tunnel by rendering it against the rRNA.
3. **tRNA tails:** residue **76** (the 3′ CCA `A76`) of the A-/P-site tRNAs; the
   ribose-ring centroid. If the deposition lacks one tRNA, graft as in stage 2.

Record the confirmed IDs back into that organism's config file **and** its
provenance note. If you cannot confidently place a landmark for *N. crassa* or
*H. sapiens*, **stop and report** what's blocking you rather than shipping a
mis-oriented ribosome — a wrong orientation silently corrupts every downstream
synthesis run.

---

## 7. Acceptance tests (definition of done)

1. **E. coli regression (hard gate).** Running the full pipeline on *E. coli*
   (4V9D + 5JTE) must reproduce the shipped structures:
   - `out/ecoli/*_oriented.pdb` matches
     `topo/csp/structures/4v9d_50S_PtR_5jte_AtR_model.pdb` — same atoms, RMSD ≈ 0
     after best-fit (ideally the transform matches within numerical tolerance; a
     global sign/rotation ambiguity is acceptable **only** if the orientation
     asserts in §5 still hold, i.e. tunnel on +x, tails on +y).
   - After stages 4–5, the truncated CG bead set matches
     `topo/csp/structures/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb` in bead count
     and per-segID composition (use `truncate_ribosome.py`'s printed stats and the
     defaults `--r-cyl 30 --x-lo -8 --x-exit 58 --keep-segids PtR,AtR`).
   If it doesn't reproduce, fix stage 2/3 before touching the other organisms.
2. **Orientation asserts pass** for every organism (§5): PTC at origin, tunnel
   `v_x` along +x, tail `v_y` along +y, within tolerance.
3. **Truncation keeps the tunnel wall.** For each organism, the truncated CG PDB
   must retain the tunnel-lining proteins (e.g. uL4, uL22, uL23, uL24) and both
   tRNAs (`--keep-segids PtR,AtR`), and `load_ribosome` must parse it without a
   missing-bead-type error:
   ```python
   from topo.csp.ribosome import load_ribosome
   load_ribosome("out/<organism>/<tag>_cg_trunc.pdb")  # must not raise
   ```
4. **Provenance note** exists per organism (`out/<organism>/PROVENANCE.md`):
   PDB IDs + resolution, the chain→segID map, the four landmark atom IDs with how
   each was confirmed, the exact commands run, and any grafting/superposition done.

---

## 8. Deliverables & layout

```
sandbox/prep_rib/
  RUNBOOK.md             # this file
  gen_subunit.py         # stage 2 (you build)
  fix_orient.py          # stage 3 (you build)
  run.py  or  Makefile   # optional: drive all 5 stages per organism
  configs/
    ecoli.(ini|json)     # 4V9D + 5JTE, chain→segID map, landmarks
    yeast.(ini|json)     # 6Q8Y
    ncrassa.(ini|json)   # 7R81 (drop nascent chain; graft A-site tRNA)
    human.(ini|json)     # 8G61
  raw/<organism>/…       # downloaded cif (gitignore-able; never edit in place)
  out/<organism>/
    <tag>_model.pdb      # stage 2
    <tag>_oriented.pdb   # stage 3
    <tag>_cg.pdb         # stage 4
    <tag>_cg_trunc.pdb   # stage 5  ← the deliverable
    PROVENANCE.md
```

**Naming:** `<tag>` = `<pdbid>_<subunit>` style, e.g. `4v9d_50S_PtR_5jte_AtR` for
E. coli (match the shipped name), `6q8y_60S` for yeast, etc.

**Optional final step (only after all four organisms pass §7):** propose promoting
`gen_subunit.py` and `fix_orient.py` into `topo/csp/` as first-class modules
(`topo.csp.gen_subunit`, `topo.csp.fix_orient`) with CLIs mirroring
`cg_ribosome`/`truncate_ribosome`, and add a `structures/` entry + a `FILES.md`
row for each new organism. **Do not modify anything under `topo/csp/structures/`
or overwrite the shipped 4V9D files** — write only under `sandbox/prep_rib/out/`.

---

## 9. Environment & rules

- Run in the **`bioenv`** environment (the translation-dev env; `topo` is
  installed `-e`). Verify `python -c "import topo.csp; import MDAnalysis"` works
  before starting; if a library is missing, add it (prefer MDAnalysis/Biopython +
  numpy; parmed/pdbfixer optional).
- **Data safety:** never overwrite raw downloads or the shipped structures. All
  generated files go under `sandbox/prep_rib/`. Downloads are reproducible — you
  may gitignore `raw/`.
- **Reproducibility:** every stage is a script + a config; no manual PyMOL/VMD
  edits. If you must inspect visually, script the check or note it in PROVENANCE.
- Work **E. coli first** and pass the §7 regression **before** the eukaryotic
  organisms (yeast, *N. crassa*, human). That gate proves stages 2–3 are correct;
  only then generalize the landmarks.
- If a landmark for *N. crassa* or human is genuinely ambiguous, **ask / report**
  rather than shipping a guessed orientation.

Start by reading `topo/csp/FILES.md`, then scaffold `configs/ecoli.*` from the
known 4V9D+5JTE facts and drive stages 2→5 to reproduce the shipped E. coli files.
