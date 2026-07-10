# Reference ribosome structures — manifest

Rigid, oriented, coarse-grained large-subunit ribosomes for `topo.csp`. Each is
**oriented in the TOPO tunnel frame** (PTC at the origin, exit tunnel on **+x**,
tRNA tails on **+y**), so `d = sqrt(y²+z²)` is the radial distance to the tunnel
axis. Two files ship per organism:

- `*_model_cg.pdb` — **full** coarse-grained large subunit (the re-truncation
  master; re-crop it with your own criteria, see `README.md`).
- `*_model_cg_trunc.pdb` — **default** truncation (`--r-cyl 30 --x-lo -8
  --x-exit 58`), ready to feed straight into `topo.csp` as the rigid ribosome.

All carry the TOPO segID convention (`23S/25S/26S/28S`, `5.8S`, `5S`, `L#`
proteins, `PtR`/`AtR` tRNAs) and the **tRNA 3′ acceptor renumbered to residue 76**
(the anchor the CSP engine looks up). Per-organism details, chain maps, landmark
confirmations and grafts are in each `structures/<org>/PROVENANCE.md`.

| Organism | Large subunit | Main PDB (res.) | A-site tRNA | PTC landmark | Exit landmark (uL24) | CG beads | Trunc beads |
|----------|---------------|-----------------|-------------|--------------|----------------------|---------:|------------:|
| *E. coli* | 50S | **4V9D** (3.0 Å) | graft **5JTE** | `23S:2602:N6` | `L24:51:N` | 14 662 | **4575** |
| *S. cerevisiae* | 60S | **6Q8Y** (3.1 Å) | graft **8G61** (cross-species) | `25S:2971:N6` | `L26:91:N` | 18 970 | **5199** |
| *N. crassa* | 60S | **7R81** (2.7 Å) | **native** (chain t1/BC) | `26S:2931:N6` | `L26:91:N` | 18 978 | **5578** |
| *H. sapiens* | 60S | **8G61** (2.94 Å) | **native** (chain At) | `28S:4548:N6` | `L26:93:N` (RPL26) | 21 244 | **5795** |

**E. coli** reproduces the hand-built reference shipped in `topo/csp/structures/`
(bit-exact truncated bead set; oriented best-fit RMSD 0.16 Å).

## Caveats to cite/note when using these

- **Yeast A-site tRNA is a cross-species graft** from the human 8G61 decoding
  complex (6Q8Y has only P/E tRNAs), superposed via the conserved P-site-tRNA
  acceptor arm (fit ~3.2 Å). Its acceptor `R` bead is a sane A-site anchor, but it
  is not a native yeast A-site tRNA.
- **Eukaryote exit landmarks are homology-assigned** (uL24 = RPL26 = segID `L26`,
  *not* eL24). The exact exit residue on the RPL26 vestibule loop is ±2 residues
  ambiguous; orientation is insensitive to this (<3° over the ~100 Å tunnel vector).
- **N. crassa** uses its **native** A- and P-site tRNAs (the 7R81 A-site is *not*
  empty, contrary to some expectations); the bound **nascent peptide is dropped**.
- **Human P-site tRNA is a 77-nt tRNA** (acceptor A77 in the deposition),
  renumbered so the acceptor is residue 76 like the others.

## Flexible exit-tunnel loop (`ribo_free_mask`)

By default the whole ribosome is rigid (mass-0). The `ribo_free_mask` CSP key frees
a portion of the exit-tunnel protein (uL24 family: E. coli **L24**, eukaryotic
**L26**) as a topo-style structure-based Go loop, reproducing O'Brien's mobile-L24
setup. Native contacts for the loop are built by `build_nonbonded_interaction` from
an **all-atom** chain carved out of the deposited CIF (the CG truncation is Cα-only
and cannot supply STRIDE H-bonds / heavy-atom contacts). One carved chain ships per
organism under `structures/<org>/`:

| Organism | Carved chain | Source CIF : auth chain | Resid range | Default `ribo_free_mask` |
|----------|--------------|-------------------------|:-----------:|--------------------------|
| *E. coli* | `structures/ecoli/L24_atomistic.pdb` | `4v9d.cif : DU` | 1–102 | `L24 : 42 - 59` |
| *S. cerevisiae* | `structures/yeast/L26_atomistic.pdb` | `6q8y.cif : AK` | 2–127 | `L26 : 82 - 99` |
| *N. crassa* | `structures/ncrassa/L26_atomistic.pdb` | `7r81.cif : a1` | 2–122 | `L26 : 82 - 99` |
| *H. sapiens* | `structures/human/L26_atomistic.pdb` | `8g61.cif : LY` | 1–134 | `L26 : 84 - 101` |

Each carved chain's residue numbering matches its truncated ribosome exactly, so
contacts align by resid. The mask ranges are the **same-window offset** of the
E. coli loop about the exit landmark (E. coli `L24:51` ≡ `L26:91` yeast/N. crassa,
`L26:93` human): an ~18-residue window centred on the apex. Note the eukaryotic
L26 β-hairpin is structurally **shorter** than E. coli L24's (a 5-residue turn vs a
12-residue loop), so this window also frees flanking β-strand residues — a
deliberate choice to keep the same physical extent across organisms.

Regenerate a carved chain with:

```
python helpers/carve_flexible_protein.py raw/ecoli/4v9d.cif DU \
    structures/ecoli/L24_atomistic.pdb
```

Use it in a CSP control file:

```ini
ribo_free_mask = L26 : 82 - 99
ribo_free_pdb  = assets/csp/prepare_ribosome/structures/yeast/L26_atomistic.pdb
```
