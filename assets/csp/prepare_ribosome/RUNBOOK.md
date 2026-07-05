# Methodology & adding a new organism

How the reference ribosomes in this bundle are built, and how to add another. For
*using* the shipped structures and the quick-start commands see
[`README.md`](README.md); for the catalog of what ships see
[`MANIFEST.md`](MANIFEST.md); for per-organism provenance see
`structures/<organism>/PROVENANCE.md`.

Goal: turn a raw PDB deposition into a **rigid, oriented, truncated, coarse-grained
large subunit** whose exit tunnel lies on the **+x** axis (PTC at the origin, tRNA
tails on +y), ready to drop into `topo.csp` as the rigid ribosome.

## The 5-stage pipeline

| # | Stage | Script | Output |
|---|-------|--------|--------|
| 1 | **acquire** raw mmCIF from the PDB | `fetch_cifs.sh` (git-ignored `raw/`) | `raw/<org>/*.cif` |
| 2 | **gen** large subunit + A/P tRNAs → one all-atom PDB w/ TOPO segIDs | `gen_subunit.py` | `out/<org>/*_model.pdb` |
| 3 | **orient** PTC→origin, tunnel→+x, tRNA tails→+y | `fix_orient.py` | `out/<org>/*_oriented.pdb` |
| 4 | **coarse-grain** (protein→Cα; RNA→P/R/BR) | `cg_ribosome.py` | `out/<org>/*_cg.pdb` |
| 5 | **truncate** (tunnel cylinder + exit half-space) | `truncate_ribosome.py` | `out/<org>/*_cg_trunc.pdb` |

`run.py -c configs/<org>.ini` drives stages 2–5. The two shipped deliverables per
organism (`structures/<org>/*_cg.pdb` full CG + `*_cg_trunc.pdb` default crop) are
curated copies of the stage-4/5 outputs.

## Conventions

- **segIDs** (PDB cols 73–76): `23S`/`25S`/`26S`/`28S`, `5.8S`, `5S` for rRNAs;
  `L#` for large-subunit proteins; `PtR`/`AtR` for the P-/A-site tRNAs. Everything
  downstream keys off the segID.
- **CG mapping.** Protein → one Cα bead (resname kept). RNA → `P` (phosphate),
  `R` (ribose-ring centroid of **`C1' C2' C3' C4' C5'`**), `BR1`/`BR2` (base-ring
  centroids; `BR2` on purines only).
- **tRNA acceptor renumbered to residue 76.** `gen_subunit` shifts each tRNA so its
  3′-terminal acceptor (the CCA `A`) is residue 76 — the anchor the CSP engine looks
  up (`PtR`/`AtR` : 76 : `R`/`P`/`BR2`). The acceptor is always a purine A, so it has
  a `BR2` bead. (E.g. the human P-site tRNA is a 77-nt tRNA with its acceptor at 77 in
  the deposition → shifted to 76.)

## Orientation landmarks (stage 3)

`fix_orient` applies one rigid transform (proper rotation, `det(R)=+1`; asserts
verified) from four landmarks:

- `P_PTC`  — the PTC adenine `N6`  → the origin
- `P_exit` — the uL24-homolog backbone `N` → defines +x (the `PTC→exit` tunnel vector)
- `P_A` / `P_P` — centroid of the A-/P-site tRNA **3′-acceptor** ribose ring
  (`:last:ribose`, = residue 76 after renumbering) → the tail vector defines +y

`uL24` is the universal name; the **eukaryotic homolog is RPL26** (segID `L26`) — not
the different protein eL24.

### Per-organism landmarks (verified)

| Organism | Large subunit / main PDB | A-site tRNA | PTC (`N6`) | Exit (`N`, uL24) |
|----------|--------------------------|-------------|-----------|------------------|
| *E. coli* | 50S / **4V9D** | graft **5JTE** chain AW | `23S:2602` | `L24:51` |
| *S. cerevisiae* | 60S / **6Q8Y** | graft **8G61** chain At (via P-tRNA) | `25S:2971` | `L26:91` |
| *N. crassa* | 60S / **7R81** | **native** (auth `t1` / label BC) | `26S:2931` | `L26:91` |
| *H. sapiens* | 60S / **8G61** | **native** (chain At) | `28S:4548` | `L26:93` |

Notes: 6Q8Y has only P/E tRNAs → A-site grafted. **7R81 has a native A-site tRNA**
(it is *not* empty) — P-site = auth `u1` (label CC), A-site = auth `t1` (label BC),
verified by superposition; its bound **nascent peptide (chain `v1`) is dropped**.
8G61 has native A/P tRNAs and needs no graft.

## Adding a new organism

1. **Fetch** the cif(s): add to `fetch_cifs.sh` and run it (populates `raw/`).
2. **Config.** For a eukaryote, `python make_euk_config.py <org>` auto-selects the
   60S rRNAs + proteins (by proximity to the large rRNA) and the tRNAs; edit the
   per-organism spec at the top of that script (main cif, big/small rRNA chains,
   P-/A-site tRNA chains or graft donor, and the landmark residues). For a bacterium,
   hand-write the `[mol:*]` chain→segID map like `configs/ecoli.ini`.
3. **Verify the landmarks — do not guess.** rRNA/protein numbering differs between
   organisms, so confirm the PTC adenine and the uL24-homolog exit residue by
   structural homology:
   - `python helpers/scan_landmarks.py <cif> <p_tRNA_chain> <big_rRNA_chain> <uL24_chain>`
     superposes the P-site tRNA acceptor arm onto the validated E. coli oriented frame
     and reports the nearest-origin adenine (PTC) and the uL24-homolog residue at the
     tunnel exit (x ≈ +100 Å).
   - `python helpers/verify_landmarks.py …` confirms a specific candidate.
   - For a structure with two tRNAs of unknown site, superpose each acceptor arm onto
     the E. coli P- and A-site references — the better fit assigns the site.
4. **Run** `python run.py -c configs/<org>.ini` and check the stage-3 orientation
   asserts pass.
5. **Publish** the two shipped files (`*_cg.pdb`, `*_cg_trunc.pdb`) into
   `structures/<org>/` and write `structures/<org>/PROVENANCE.md` (PDB IDs +
   resolution, chain→segID map, the four landmark atoms and how each was confirmed,
   commands run, any graft/superposition).

## Acceptance criteria

1. **E. coli regression (hard gate).** The full pipeline on 4V9D + 5JTE reproduces
   the reference: `*_cg_trunc.pdb` matches the historical shipped bead set
   (4575 beads, exact per-segID composition) and the oriented model best-fits it at
   RMSD ≈ 0 (the 50S to ~0.03 Å; the grafted AtR within ~1.3 Å).
2. **Orientation asserts pass** for every organism: PTC at origin, tunnel `v_x` on +x,
   tail `v_y` on +y (within tolerance).
3. **Truncation keeps the tunnel wall** (uL4/uL22/uL23/uL24 homologs) and both tRNA
   res-76 acceptor anchors, and `topo.csp.ribosome.load_ribosome(<trunc>)` parses it
   without a missing-bead-type error.
4. **Provenance** exists per organism.
