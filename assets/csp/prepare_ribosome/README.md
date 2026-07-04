# `prepare_ribosome` — reference ribosomes for `topo.csp`

Ready-to-use **rigid, oriented, coarse-grained large-subunit ribosomes** for the
continuous-synthesis engine (`topo.csp`), plus the scripted pipeline that builds
them. If you don't have a ribosome structure of your own, use one of the bundled
ones directly; if you do, run the pipeline on it.

Everything here is oriented in the TOPO **tunnel frame**: PTC at the origin, exit
tunnel on **+x**, tRNA tails on **+y**, so the tunnel axis *is* the X-axis.

## What ships (in `structures/<organism>/`)

For **E. coli, yeast, N. crassa, human** (see [`MANIFEST.md`](MANIFEST.md)):

| File | What | Use |
|------|------|-----|
| `*_model_cg.pdb` | **full** CG large subunit (oriented) | re-truncate with your own criteria |
| `*_model_cg_trunc.pdb` | **default** truncation (`r_cyl 30, x_lo −8, x_exit 58`) | drop straight into `topo.csp` |

The truncated file is what `topo.csp` consumes; the full CG is shipped so you can
re-crop without re-running the whole pipeline. Both carry the TOPO segID
convention and the **tRNA 3′ acceptor renumbered to residue 76** (the anchor the
CSP engine looks up). Provenance (PDB IDs, chain maps, landmark confirmations,
grafts, caveats) is in each `structures/<organism>/PROVENANCE.md`.

**Not shipped** (large / regenerable): the raw mmCIF depositions (18–35 MB each —
fetch with `fetch_cifs.sh`) and the all-atom / oriented / intermediate files
(recreate with `run.py`).

## Use a bundled ribosome in a CSP run

These are plain PDB files, not a package resource — point the `ribosome` key of
your CSP `.ini` at one by path (relative to your repo checkout, or an absolute
path):

```ini
# csp.ini
pdb_file  = my_protein.pdb
ribosome  = assets/csp/prepare_ribosome/structures/human/8g61_60S_model_cg_trunc.pdb
domain_def = my_protein.yaml
# ... MD keys ...
```

## Re-truncate with your own criteria (why the full CG ships)

The default crop keeps the exit-tunnel shell (`r_cyl 30 Å`, `x_lo −8`, `x_exit 58`).
To keep more/less of the wall, a longer exit cap, or the whole tRNAs, run the
truncator on the **full CG** — it's already oriented, so the X-axis is the tunnel:

```bash
python truncate_ribosome.py \
    -i structures/human/8g61_60S_model_cg.pdb \
    -o my_human_trunc.pdb \
    --r-cyl 35 --x-lo -8 --x-exit 70 --keep-segids PtR,AtR
```

## Regenerate, or build a new organism

Needs `gemmi` (`pip install gemmi`) in addition to `topo`.

```bash
bash fetch_cifs.sh                     # download raw cifs into ./raw/ (git-ignored)
python run.py -c configs/human.ini     # stages: gen -> orient -> cg -> truncate  (writes ./out/)
```

Pipeline stages (all four scripts live in this bundle):

| # | Stage | Script |
|---|-------|--------|
| 2 | **gen** large subunit + tRNAs → all-atom PDB w/ segIDs | `gen_subunit.py` |
| 3 | **orient** PTC→origin, tunnel→+x, tails→+y | `fix_orient.py` |
| 4 | **coarse-grain** | `python cg_ribosome.py` |
| 5 | **truncate** | `python truncate_ribosome.py` |

`run.py -c configs/<organism>.ini` drives all of it. To add an organism: fetch its
cif, generate a config (`make_euk_config.py` for a eukaryote, then verify landmarks
with `helpers/scan_landmarks.py` / `verify_landmarks.py`), and run. See
[`RUNBOOK.md`](RUNBOOK.md) for the full spec and [`helpers/view_trunc.tcl`](helpers/view_trunc.tcl)
to eyeball the X-axis alignment in VMD.

> **Note.** `run.py` writes intermediates and the regenerated CG/trunc files into
> `./out/` (git-ignored). The curated copies under `structures/` are the shipped
> assets; refresh them by copying the two `*_cg.pdb` / `*_cg_trunc.pdb` from
> `out/<organism>/` after a run.
