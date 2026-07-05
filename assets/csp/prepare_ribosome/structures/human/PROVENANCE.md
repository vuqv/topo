# Provenance — *H. sapiens* (human) 60S reference ribosome

**Tag:** `8g61_60S` · **Deliverable:** `8g61_60S_model_cg_trunc.pdb`

## Source structure

| Role | PDB | Resolution | Notes |
|------|-----|-----------|-------|
| 60S large subunit + **both** tRNAs | **8G61** | 2.94 Å | Human 80S decoding complex; has native **A-site (chain `At`)** and **P-site (chain `Pt`, "RNA (77-MER)")** tRNAs + mRNA. Full 80S — the 40S (18S + S-proteins) and mRNA are dropped. **No graft needed.** |

Chosen because it is a complete, high-resolution 80S with both A- and P-site
tRNAs already bound (the runbook's preferred entry), so the tRNA-tail landmark and
the A-site anchor come directly from the deposition.

## Chain → segID map (8G61)

- **28S** rRNA = `L5`, **5.8S** = `L8`, **5S** = `L7`
- **PtR** = chain `Pt` (native), **AtR** = chain `At` (native, residues 1–76)
- **42 large-subunit (60S) proteins** — auto-selected by proximity to 28S (closer
  to 28S than 18S). uL24 homolog **L26 = chain `LY`** (human RPL26; the *different*
  protein eL24 = chain `LW`, "60S ribosomal protein L24", is correctly **not**
  used as the landmark).

Full section list: `configs/human.ini`.

## Landmarks (stage 3) — verified by structural homology

| Landmark | Atom | Confirmation |
|----------|------|-------------|
| `P_PTC`  | `28S:4548:N6` | **Verified:** adenine with N6; after P-site-tRNA superposition into the E. coli frame, its N6 lands **7.3 Å from the origin** (next-nearest adenine 14.3 Å), confirming the runbook's suggested residue. |
| `P_exit` | `L26:93:N` (RPL26 = uL24, chain LY) | **Verified:** lands at **x ≈ 96 Å** — the tunnel exit vestibule; residues 90–94 of RPL26 all cluster there (the homology-optimal residue is 93). Orientation is insensitive to ±2 residues on this loop (<3° over the ~100 Å tunnel vector). |
| `P_A`    | `AtR:last:ribose` → model res **76** | native A-site tRNA, standard 76-nt, 3′ acceptor A76. |
| `P_P`    | `PtR:last:ribose` → model res **76** | native P-site tRNA is a **77-nt tRNA** (Met-tRNA, entity "RNA (77-MER)"): 3′ acceptor is **A77 in the deposition**, **renumbered to 76 in the assembled model** (see box below). `:last:` picks the acceptor either way. |

> **Why the P-site tRNA has 77 residues.** Chain `Pt` is a 77-nucleotide tRNA
> (positions 1–77, all modelled, including modified G7M at 47), i.e. one
> nucleotide longer than the canonical 76-mer, so **in the deposition** its 3′
> acceptor is **A77**, not A76 (an aminoacyl Met is esterified to it = res 78,
> dropped). Two consequences, both handled:
>
> 1. **Orientation.** The tail landmark uses `:last:` (highest-numbered ribose
>    nucleotide = the true acceptor). A hardcoded `76` would have put the P
>    landmark on the penultimate C76 (5.2 Å away), rolling the ribosome **10.2°
>    about the tunnel axis**. Corrected A↔A tail spacing is 7.2 Å, matching the
>    ~6 Å P/A spacing in E. coli / N. crassa.
> 2. **Downstream anchors.** `gen_subunit` **renumbers the tRNA so its 3′
>    acceptor is residue 76** (Pt shifted −1: 1–77 → 0–76). The whole CSP engine
>    (`read_anchor`, `add_trna_tether`, `optimal_ptc_targets` /
>    `_ptc_bead_index`) looks up `PtR`/`AtR` `:76:{R,P,BR2}`. Without the shift,
>    human `PtR:76` was a **cytidine (no BR2 bead)**, which silently mis-placed
>    the P-anchor by ~5 Å and would have **crashed the (always-on) PTC-geometry
>    optimization `optimal_ptc_targets`** (missing `BR2`). After the shift,
>    `PtR:76` is the purine acceptor A with
>    R/P/BR2 all present; verified A–P anchor spacing 0.722 nm.

**uL24 vs eL24.** The exit landmark uses **RPL26 = uL24** (chain LY). The
similarly-named eL24 (chain LW, deposited as "60S ribosomal protein L24") is a
*different* protein and was deliberately avoided, per the runbook's warning.

## Commands

```bash
curl -sfL https://files.rcsb.org/download/8G61.cif.gz | gunzip > raw/human/8g61.cif
python make_euk_config.py human          # regenerate configs/human.ini
python run.py -c configs/human.ini       # gen -> orient -> cg -> truncate
```

## Verification

| Check | Result |
|-------|--------|
| Orientation asserts | **all pass** |
| CG beads | 21 244 |
| Truncated deliverable | 5 795 beads |
| Tunnel wall retained | uL4 (`L4`), uL22 (`L17`), uL23 (`L23a`), uL24 (`L26`) all present |
| Both tRNA anchors retained | `PtR:76:R` (x=2.1), `AtR:76:R` (x=3.5) |
| `load_ribosome(...)` | parses OK — 5 795 beads |
