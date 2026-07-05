# The ribosome structure (get one, or build your own)

The explicit-ribosome runner ({doc}`continuous_synthesis`) needs a **rigid,
oriented, truncated, coarse-grained large subunit** as its `ribosome` input — the
`ribosome_trunc.pdb` you see in the tutorials. TOPO **ships ready-made ones** for the
most-studied organisms, and includes the scripted pipeline that builds them, so you
can also prepare your own from any PDB deposition. (The cylinder runner needs none of
this.)

Every structure is oriented in the TOPO **tunnel frame**: PTC at the origin, exit
tunnel on **+x**, tRNA tails on **+y** — so the tunnel axis *is* the X-axis and the
radial distance to it is `sqrt(y²+z²)`.

## Ready-made reference ribosomes

They live in the repo at **`assets/csp/prepare_ribosome/structures/<organism>/`**, two
files each: `*_model_cg.pdb` (the full CG large subunit — a re-truncation master) and
`*_model_cg_trunc.pdb` (the default exit-tunnel crop, ready to run).

| Organism | Large subunit / PDB | PTC landmark | Exit landmark (uL24) | Trunc beads |
|----------|---------------------|--------------|----------------------|------------:|
| *E. coli* | 50S / **4V9D** (+5JTE A-tRNA) | `23S:2602:N6` | `L24:51:N` | 4 575 |
| *S. cerevisiae* | 60S / **6Q8Y** | `25S:2971:N6` | `L26:91:N` | 5 199 |
| *N. crassa* | 60S / **7R81** | `26S:2931:N6` | `L26:91:N` | 5 578 |
| *H. sapiens* | 60S / **8G61** | `28S:4548:N6` | `L26:93:N` (RPL26) | 5 795 |

Per-organism provenance (chain→segID maps, landmark confirmation, grafts) is in each
`structures/<organism>/PROVENANCE.md`; the full catalog is `MANIFEST.md`.

```{note}
**Caveats to cite when you use these.** The **yeast A-site tRNA is a cross-species
graft** from the human 8G61 decoding complex (6Q8Y has only P/E tRNAs); its acceptor
is a sane A-site anchor but not a native yeast tRNA. The **eukaryote exit landmark**
(uL24 = RPL26 = segID `L26`, *not* eL24) is homology-assigned (±2 residues on the
vestibule loop; orientation is insensitive to this). *N. crassa* uses its **native**
A/P tRNAs (its A-site is *not* empty); its bound nascent peptide is dropped.
```

## Use one in a run

The `ribosome` key of a `csp.ini` is a plain path — point it at a truncated
structure:

```ini
pdb_file  = my_protein.pdb
ribosome  = assets/csp/prepare_ribosome/structures/human/8g61_60S_model_cg_trunc.pdb
domain_def = my_protein.yaml
```

## Re-truncate with your own criteria

The default crop keeps the exit-tunnel shell (`r_cyl 30 Å`, `x_lo −8`, `x_exit 58`).
To keep more/less of the wall, a longer exit cap, or the whole tRNAs, run the truncator
on the shipped **full CG** (already oriented, so the X-axis is the tunnel):

```bash
cd assets/csp/prepare_ribosome
python truncate_ribosome.py \
    -i structures/human/8g61_60S_model_cg.pdb -o my_human_trunc.pdb \
    --r-cyl 35 --x-lo -8 --x-exit 70 --keep-segids PtR,AtR
```

## Build a new organism

The pipeline is five stages — **acquire → gen → orient → coarse-grain → truncate** —
driven by `run.py`:

```bash
cd assets/csp/prepare_ribosome
bash fetch_cifs.sh                     # download raw cifs (needs internet; gemmi)
python run.py -c configs/<organism>.ini
```

Adding a new organism means writing its config (landmark residues + chain→segID map)
and **verifying the PTC and exit landmarks by structural homology** — do not guess,
since rRNA/protein numbering differs between organisms. The bundle's
`RUNBOOK.md` is the methodology + "how to add an organism" guide, and
`helpers/scan_landmarks.py` / `verify_landmarks.py` confirm the landmarks against the
validated *E. coli* frame. See `assets/csp/prepare_ribosome/README.md` for the full
workflow.
