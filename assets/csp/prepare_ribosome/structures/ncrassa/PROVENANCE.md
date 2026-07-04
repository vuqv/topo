# Provenance — *N. crassa* 60S reference ribosome

**Tag:** `7r81_60S` · **Deliverable:** `7r81_60S_model_cg_trunc.pdb`

## Source structure

| Role | PDB | Resolution | Notes |
|------|-----|-----------|-------|
| 60S large subunit + **both** tRNAs | **7R81** | 2.7 Å | *N. crassa* 80S translating ribosome (cycloheximide-arrested). Full 80S — the 40S (18S + S-proteins), mRNA (chain `w1`), RACK1 (`h2`) and the **bound nascent peptide (chain `v1`) are dropped**. |

**Runbook correction (verified).** The runbook expected 7R81's A-site to be empty
and told us to graft an A-site tRNA. It is **not** empty: the two tRNAs (auth
chains `t1`/`u1`, deposited under label_asym IDs `BC`/`CC`) both sit at the PTC —
their acceptor ends are 2.7–2.8 Å from the 4-residue nascent peptide, 6.0 Å apart,
and ~8 Å from the PTC adenine. So 7R81 carries a **native A-site and P-site tRNA**
and **no graft is needed**.

### Which tRNA is which (verified by homology, per the user's request to check)

Superposing each tRNA's acceptor arm onto the validated E. coli oriented P- and
A-site references:

| tRNA (auth / label) | fit to E. coli PtR | fit to E. coli AtR | PTC lands at | Assignment |
|---------------------|-------------------:|-------------------:|-------------:|------------|
| `u1` / `CC` | **1.28 Å** | 2.28 Å | 3.3 Å from origin | **P-site → PtR** |
| `t1` / `BC` | 1.99 Å | **1.27 Å** | 3.4 Å from origin | **A-site → AtR** |

## Chain → segID map (7R81)

- **26S** rRNA = `A1`, **5.8S** = `C1`, **5S** = `B1`
- **PtR** = chain `u1` (label CC), **AtR** = chain `t1` (label BC) — both native
- **41 large-subunit (60S) proteins** — auto-selected by proximity to 26S (closer
  to 26S than 18S); nascent peptide, RACK1 and mRNA excluded. uL24 homolog
  **L26 = chain `a1`** (deposition name "Ribosomal protein L26"; distinct from
  eL24 = chain `Y1`, "60S ribosomal protein L24", which is *not* the landmark).

Full section list: `configs/ncrassa.ini`.

## Landmarks (stage 3) — verified by structural homology

| Landmark | Atom | Confirmation |
|----------|------|-------------|
| `P_PTC`  | `26S:2931:N6` | **Verified** homolog of the yeast/E. coli PTC adenine: after P-site-tRNA superposition into the E. coli frame, this adenine's N6 lands **3.3 Å from the origin** (next-nearest 14.8 Å). |
| `P_exit` | `L26:91:N` (chain a1) | **Verified:** lands at **x ≈ 97 Å** on the tunnel axis; the residue is on the same L26 exit loop as yeast res 91 (residues 90–94 all cluster at the vestibule). |
| `P_A`    | `AtR:last:ribose` = t1 res 76 | native A-site tRNA res-76 present. |
| `P_P`    | `PtR:last:ribose` = u1 res 76 | native P-site tRNA res-76 present. |

## Commands

```bash
curl -sfL https://files.rcsb.org/download/7R81.cif.gz | gunzip > raw/ncrassa/7r81.cif
python make_euk_config.py ncrassa        # regenerate configs/ncrassa.ini
python run.py -c configs/ncrassa.ini     # gen -> orient -> cg -> truncate
```

## Verification

| Check | Result |
|-------|--------|
| Orientation asserts | **all pass** |
| CG beads | 18 978 |
| Truncated deliverable | 5 578 beads |
| Tunnel wall retained | uL4 (`L4`), uL22 (`L17`), uL23 (`L25`), uL24 (`L26`) all present |
| Both tRNA anchors retained | `PtR:76:R` (x=5.2), `AtR:76:R` (x=6.0) |
| `load_ribosome(...)` | parses OK — 5 578 beads |
| Nascent peptide dropped | chain `v1` not whitelisted (would corrupt the tunnel) |
