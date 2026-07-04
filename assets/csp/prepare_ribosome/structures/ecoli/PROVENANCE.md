# Provenance — *E. coli* 50S reference ribosome

**Tag:** `4v9d_50S_PtR_5jte_AtR` · **Deliverable:** `4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb`

This is the **regression / hard-gate** organism: the pipeline reproduces the
hand-built structure shipped in `topo/csp/structures/4v9d_50S_PtR_5jte_AtR_model.pdb`.

## Source structures

| Role | PDB | Resolution | Notes |
|------|-----|-----------|-------|
| 50S large subunit + P-site tRNA | **4V9D** | 3.0 Å | *E. coli* 70S, **two copies** per asymmetric unit; the **D-copy** is the one the shipped file uses. |
| A-site tRNA donor | **5JTE** | 2.6 Å (cryo-EM) | A-site tRNA grafted (4V9D lacks a well-resolved A-site tRNA). |

Raw cif: `raw/ecoli/4v9d.cif`, `raw/ecoli/5jte.cif` (downloaded from RCSB;
`raw/` is reproducible and git-ignorable).

## Chain → segID map (4V9D auth chain, D-copy)

Determined by (a) exact sequence match and (b) **coordinate superposition** onto
the shipped file — the D-copy gives RMSD ≈ 0, the C-copy is 0.4–1.5 Å off, so the
copy is unambiguous.

| segID | 4V9D | segID | 4V9D | segID | 4V9D |
|-------|------|-------|------|-------|------|
| 23S | DA | L11 | DI | L27 | DW |
| 5S  | DB | L13 | DJ | L28 | DX |
| PtR | BV | L14 | DK | L29 | DY |
| L2  | DC | L15 | DL | L30 | DZ |
| L3  | DD | L16 | DM | L32 | D0 |
| L4  | DE | L17 | DN | L33 | D1 |
| L5  | DF | L18 | DO | L34 | D2 |
| L6  | DG | L19 | DP | L35 | D3 |
| L9  | DH | L20 | DQ | L36 | D4 |
|     |    | L21 | DR |     |    |
|     |    | L22 | DS |     |    |
|     |    | L23 | DT |     |    |
|     |    | L24 | DU |     |    |
|     |    | L25 | DV |     |    |

**A-site tRNA (AtR):** 5JTE chain **AW**, residues **2–76** (the aminoacyl LYS at
res 101 and the disordered 5′ res 1 are dropped). Grafted by superposing 5JTE's
**large subunit** onto 4V9D's (29 chains, 3316 backbone atoms, fit RMSD 3.1 Å)
and lifting AW into the target A-site frame.

## Landmarks (stage 3)

| Landmark | Atom | Confirmed |
|----------|------|-----------|
| `P_PTC`  | `23S:2602:N6`  | Given/authoritative; A2602 present with N6. |
| `P_exit` | `L24:51:N`     | Given/authoritative; uL24 res 51 backbone N present. |
| `P_A`    | `AtR:last:ribose` (res 76) (C1′–C5′ centroid) | tRNA res-76 present. |
| `P_P`    | `PtR:last:ribose` (res 76) (C1′–C5′ centroid) | tRNA res-76 present. |

> Note: the runbook §5 text calls the ribose set `C1' C2' C3' C4' O4'`, but the
> actual `cg_ribosome.py` R-bead uses `C1'..C5'`. We use **C1′–C5′** so the
> landmark centroid and the CG `R` bead agree (as the runbook itself requires).

## Commands

```bash
# stage 1: acquire
curl -sfL https://files.rcsb.org/download/4V9D.cif.gz | gunzip > raw/ecoli/4v9d.cif
curl -sfL https://files.rcsb.org/download/5JTE.cif.gz | gunzip > raw/ecoli/5jte.cif
# stages 2-5
python run.py -c configs/ecoli.ini      # gen -> orient -> cg -> truncate
```

Truncation defaults: `--r-cyl 30 --x-lo -8 --x-exit 58`, **no** `--keep-segids`
(this is what produced the shipped 4575-bead file; both tRNA acceptor ends —
res 76 of PtR and AtR, the elongation anchors — survive the cylinder crop).

## Verification (acceptance test 1)

| Check | Result |
|-------|--------|
| `*_oriented.pdb` vs shipped, best-fit RMSD (93 982 common atoms) | **0.164 Å** (50S+PtR 0.006–0.06 Å = PDB rounding; AtR graft 1.25 Å) |
| Orientation asserts (PTC@0, tunnel +x, tails +y) | **all pass** |
| `*_cg.pdb` total beads | **14 662** — exact match to shipped, all 33 segIDs identical |
| `*_cg_trunc.pdb` | **4575 beads, 0 per-segID differences** vs shipped |
| `load_ribosome(...)` | parses OK — `Ribosome`, 4575 beads |

The AtR graft lands 1.25 Å from the shipped A-site tRNA because the shipped
file's exact graft reference is not recoverable (likely a local/hand-refined
fit); a full-large-subunit superposition is the runbook-sanctioned method and
1.25 Å is negligible for a rigid excluded-volume tRNA whose only functional role
downstream is the res-76 anchor. The 50S tunnel geometry — everything the
truncation and synthesis depend on — reproduces exactly.
