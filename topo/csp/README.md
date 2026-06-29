# `topo.csp` — Continuous Synthesis Protocol (O'Brien)

The per-codon, three-stage co-translational synthesis protocol of
`continuous_synthesis_v6.py` (Yang Jiang, Dan Nissley, Ed O'Brien), ported to topo.
It is the **kinetic upgrade** of [`topo.translation.elongate`](../translation/README.md):
that driver grows the chain at a *fixed* `n_steps`; `topo.csp` times every residue
from its mRNA codon and splits it into the three sub-stages of the elongation cycle.

```bash
topo-csp -f csp.ini
python -m topo.csp -f csp.ini
```

See **[Tutorial 10](../../tutorials/10_csp_obrien/)** for a runnable example, and its
[`PLAN.md`](../../tutorials/10_csp_obrien/PLAN.md) for the porting design + the
done/remaining matrix.

## What it does

For nascent length `L = L0 .. L_max`, each residue is added through three sub-stages,
each run for an **exponentially-sampled** dwell time (first-passage-time sampling):

| stage | biology | dwell-time mean (s) | C-terminus restraint |
|-------|---------|---------------------|----------------------|
| 1 | peptidyl transfer / A-site delivery | `time_stage_1` | → A-anchor |
| 2 | translocation begins | `time_stage_2` (+ traffic) | → A-anchor |
| 3 | tRNA binding / waiting | `codon_total − stage1 − stage2` | → **P-anchor** |

The A→P restraint switch between stages 2 and 3 reproduces translocation. Stage 3's
final structure seeds the next residue. After `L_max`, optional **ejection** (release
the restraint) and **dissociation** (free drift-off) phases run.

Dwell times map to integration steps the O'Brien way:

```
t_sim_ns = dwell_s * 1e9 / scale_factor          # in-vivo → in-silico
steps    = int(t_sim_ns / (dt_ps * 1e-3))        # → MD steps
```

## Design

`topo.csp` adds **only** the kinetics and the 3-stage loop; **all MD is delegated**
to `topo.translation.elongate.run_length` (build-once-subset contacts, coordinate
seeding, rigid-ribosome scenery, tunnel wall, minimize/run/finalize). CSP drives the
position-restraint path (target switchable A↔P), so it forces `trna_tether = no`;
`rigid_ribosome` / `tunnel_wall` / excluded volume + electrostatics work as in
`topo.translation` build step v2.

| file | role |
|------|------|
| `kinetics.py` | pure timing core (no OpenMM): codon tables, FPT sampling, `scale_factor`→steps, 3-stage split, ribosome-traffic hook |
| `csp.py` | `CSPParams` / `CSPConfig`, `read_csp_config` (INI), `run_continuous_synthesis`, `csp()` CLI |
| `__init__.py`, `__main__.py` | public API; `python -m topo.csp` |

## Public API

```python
from topo.csp import run_continuous_synthesis, read_csp_config, CSPParams, kinetics

cfg = read_csp_config("csp.ini")
run_continuous_synthesis(cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max,
                         out_root=cfg.outdir, mrna=cfg.mrna, trans_times=cfg.trans_times,
                         domain_def=cfg.domain_def,
                         stride_output_file=cfg.stride_output_file, params=cfg.params)
```

## Control file (`csp.ini`)

A single `[OPTIONS]` section — a **superset of `elongate.ini`**. Required:
`pdb_file`, `ribosome`, `L0`. Non-uniform timing also requires `mrna` + `trans_times`.

**Kinetic keys**

| key | meaning |
|-----|---------|
| `mrna` | mRNA sequence file (raw nucleotides; one codon/residue + 1 stop) |
| `trans_times` | codon → mean in-vivo time (s) table |
| `scale_factor` | in-vivo seconds → in-silico ns compressor |
| `time_stage_1` / `time_stage_2` | mean peptidyl-transfer / translocation dwell (s) |
| `uniform_ta` / `uniform_mfpt` | ignore the mRNA; use a uniform mean codon time (s) |
| `ribosome_traffic` / `initiation_rate` | external traffic correction (falls back if binary absent) |
| `random_seed` | seed for the FPT sampler (reproducible schedules) |
| `max_steps_per_stage` / `min_steps_per_stage` | clamp each stage's step count (testing) |
| `ejection_steps` / `dissociation_steps` | post-synthesis free runs (0 = skip) |

**MD / ribosome keys** (same as `elongate.ini`; `n_steps` is **not** used — step
counts come from the kinetics): `dt`, `ref_t`, `tau_t`, `nstout`, `device`, `ppn`,
`constraints`, `restraint_k`, `buffer`, `minimize`, `rigid_ribosome`, `tunnel_wall`,
`tunnel_wall_x0`, `tunnel_wall_k`, `ptc_offset`, `nascent_only_output`.

**Units:** OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm². In-vivo dwell times are
in **seconds**.

## Output

`<outdir>/L_<L>/stage_<1,2,3>/` per residue (each a standalone topo run: `traj.dcd`,
`traj.log`, `traj.psf`, `traj.chk`, `traj_final.pdb`, `traj_runinfo.log`,
`native_1_<L>.pdb`), plus `ejection/` and optional `dissociation/`. Stage 3's
`traj_final.pdb` seeds the next residue.

## Movie

Stitch the per-stage trajectories into one VMD-playable movie (chain grows stage by
stage; static ribosome overlaid):

```bash
topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
vmd -e synth_out/movie.tcl
```

`topo.csp.make_movie` discovers the nested `L_<L>/stage_<s>/` segments and reuses the
shared stitching core `topo.translation.make_movie.stitch_segments`.

## Not yet ported

CHARMM PSF/TOP/PRM/COR ingestion (O'Brien's exact systems); multi-trajectory
multiprocessing; literal mid-residue peptide-bond toggling / explicit A/P tRNA bonded
geometry; `restart = 1`; quantitative validation. See
[Tutorial 10 `PLAN.md`](../../tutorials/10_csp_obrien/PLAN.md).
