# Provenance — *S. cerevisiae* (yeast) 60S reference ribosome

**Tag:** `6q8y_60S` · **Deliverable:** `6q8y_60S_model_cg_trunc.pdb`

## Source structures

| Role | PDB | Resolution | Notes |
|------|-----|-----------|-------|
| 60S large subunit + P-site tRNA | **6Q8Y** | 3.1 Å | Yeast 80S–Xrn1 complex; has **P-site (chain n)** and **E-site (chain m)** tRNAs but **no A-site tRNA**. Small subunit (18S + S-proteins), mRNA and Xrn1 are dropped. |
| A-site tRNA donor (graft) | **8G61** | 2.94 Å | Human 80S decoding complex; its A-site tRNA (chain `At`) is grafted (see below). |

## Chain → segID map (6Q8Y)

- **25S** rRNA = chain `BQ`, **5.8S** = `BS`, **5S** = `BR`
- **PtR** (P-site tRNA) = chain `n`
- **AtR** (A-site tRNA) = **grafted** from 8G61 chain `At`, residues 1–76
- **40 large-subunit (60S) proteins** — auto-selected by structural proximity to
  the 25S rRNA (closer to 25S than to 18S), non-ribosomal chains excluded. Each is
  labelled by its yeast L-name (segID). The uL24 homolog **L26 = chain `AK`**.

Full section list: `configs/yeast.ini`.

## Landmarks (stage 3) — given + verified

| Landmark | Atom | Confirmation |
|----------|------|-------------|
| `P_PTC`  | `25S:2971:N6` | Given/authoritative. **Verified:** residue is an adenine with N6; after superposing the P-site tRNA acceptor arm onto the validated E. coli oriented frame, its N6 lands **3.8 Å from the origin** (next-nearest adenine is 14.6 Å). |
| `P_exit` | `L26:91:N` (uL24 homolog, chain AK) | Given/authoritative. **Verified:** in the E. coli frame it lands at **x ≈ 98 Å** on the tunnel axis (the exit vestibule), so `PTC→exit` threads the tunnel. |
| `P_A`    | `AtR:last:ribose` (res 76) | Grafted tRNA res-76 present. |
| `P_P`    | `PtR:last:ribose` (res 76) | P-site tRNA res-76 present. |

uL24 homolog = ribosomal protein **L26** (not eL24) — as specified by the runbook
and consistent with the exit residue landing at the tunnel vestibule.

## A-site tRNA graft

6Q8Y lacks an A-site tRNA, so one is grafted from the human decoding complex 8G61
(chain `At`). Because donor and target are different species (no identical
sequence), the graft superposes the **conserved P-site-tRNA acceptor arm** of the
donor (8G61 `Pt`) onto the target P-site tRNA (`n`) — 54 backbone atoms, fit RMSD
**3.23 Å** — and lifts `At` with that transform. The A-/P-tRNA arrangement at the
ribosome is universally conserved, so this places the A-site tRNA in a
biologically sane A-site pose; only its res-76 `R` bead (the elongation anchor)
matters downstream, and it lands at x ≈ 7 Å, next to the P-site anchor.

## Commands

```bash
curl -sfL https://files.rcsb.org/download/6Q8Y.cif.gz | gunzip > raw/yeast/6q8y.cif
python make_euk_config.py yeast          # regenerate configs/yeast.ini
python run.py -c configs/yeast.ini       # gen -> orient -> cg -> truncate
```

## Verification

| Check | Result |
|-------|--------|
| Orientation asserts (PTC@0, tunnel +x, tails +y) | **all pass** |
| CG beads | 18 970 |
| Truncated deliverable | 5 199 beads |
| Tunnel wall retained | uL4 (`L4`), uL22 (`L17`), uL23 (`L25`), uL24 (`L26`) all present |
| Both tRNA anchors retained | `PtR:76:R` (x=4.9), `AtR:76:R` (x=7.1) |
| `load_ribosome(...)` | parses OK — 5 199 beads |
