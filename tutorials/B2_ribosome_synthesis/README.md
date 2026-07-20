# Tutorial B.2 — Protein synthesis on a coarse-grained ribosome

**Goal:** synthesize a protein **residue by residue on TOPO's own coarse-grained
ribosome**, using codon-resolved kinetics with `topo-csp`. This is the
ribosome-based counterpart to [Tutorial B.1](B1_translation_cylinder.md),
which replaces the physical ribosome with an analytic cylindrical tunnel.

Everything here is **standalone TOPO**: the ribosome is a truncated CG bead model
(`ribosome_trunc.pdb`, loaded by `topo.csp.ribosome.load_ribosome`) — no external
`.cor`/`.psf`/`.prm` files. The nascent chain grows N→C through the exit region
and is ejected once complete.

![residue-by-residue synthesis of 4c5c](4c5c/img/process.gif)

***Synthesis on the ribosome** — the chain (coloured beads) grows through TOPO's coarse-grained ribosome (translucent grey); already-emerged N-terminal residues fold while the C-terminus is still being added inside.*

> The GIF is for the `4c5c/` system. After a run, view your own synthesis in VMD:
> stitch the per-length frames into one movie (adding the ribosome as scenery),
> then launch VMD *from inside* `synth_out/` (the generated `movie.tcl` loads
> `movie.{psf,dcd}` by basename, so it only resolves from that folder):
> ```bash
> topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb   # -> synth_out/movie.{psf,dcd,tcl}
> cd synth_out && vmd -e movie.tcl
> ```
> For the full movie reference (both synthesis models, all `topo-csp-movie` options),
> see [Visualizing the synthesis process](../usage/synthesis_visualization.md).

**Prerequisite:** the coarse-grained model of
[Tutorial A.1](A1_single_domain.md) and
the domain-scaling idea of
[Tutorial A.2](A2_multidomain.md).

## Two worked systems

| System | Residues | Notes |
|--------|----------|-------|
| `4c5c/` | 306, multidomain | Full-length synthesis `L = 1 → 306`, then ejection. The main validation case. |
| `P0CX28/` | 106, single-domain | The single-domain protein from Tutorial A.1, `L = 1 → 106`, `nscale = 2.5044`. |

## Files (per system)

| File | Role |
|------|------|
| `*_clean.pdb` | Target folded structure (defines the native contacts). |
| `*_stride.dat` | Precomputed STRIDE backbone hydrogen bonds. |
| `domain.yaml` | Domain definition / contact `nscale` (see [Domain definition file](../usage/domain_definition.rst)). |
| `*_mrna.txt` | The mRNA codon sequence that times each residue. |
| `trans_times.txt` | Per-codon translation times (kinetics input). |
| `ribosome_trunc.pdb` | TOPO's coarse-grained truncated ribosome. |
| `csp_val.ini` | Full-length run config (`L0=1 → L_max`, `AllBonds`, equilibrium-PTC, ejection). |
| `csp_debug.ini` | Short clamped debug profile (4c5c only). |

## Run it

From a system folder (e.g. `4c5c/`):

```bash
topo-csp -f csp_val.ini        # -> synth_out/
```

The runner uses **codon-resolved kinetics** (each elongation cycle split into
three kinetic sub-stages) with a per-stage **dt-halving stability guard**, rigid
`AllBonds` (the default), and **always-on equilibrium-PTC seeding** — each new
residue is placed one peptide bond from the previous C-terminus, clear of the
ribosome, so rigid bonds seed cleanly. A successful full-length run shows **no
`[stability]` lines** in the log.

## What it produces

- `synth_out/L_001/ … L_<L_max>/` — one folder per chain length, each with the
  MD trajectory and the nascent structure at that length.
- `synth_out/ejection/` — the post-synthesis release: `traj_final.pdb`, the
  built-once native-contact PDB, and run metadata.

## Validating the run

A successful run leaves **no `[stability]` lines** in the log (per-stage potential
energy stayed finite throughout), and `synth_out/ejection/traj_final.pdb` shows the
released chain clear of the tunnel — not penetrating the ribosome.
