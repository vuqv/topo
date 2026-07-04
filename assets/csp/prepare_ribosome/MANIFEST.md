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
