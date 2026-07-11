# `topo.csp` — Continuous Synthesis Protocol (O'Brien)

The per-codon, three-stage protein synthesis protocol of
`continuous_synthesis_v6.py` (Yang Jiang, Dan Nissley, Ed O'Brien), ported to topo.
`topo.csp` (the CSP runner, `topo-csp`) times every residue from its mRNA codon and
splits it into the three sub-stages of the elongation cycle. It is a thin outer loop
over the shared low-level MD engine `topo.csp.core` (`run_length`, `RunParams`),
the rigid-ribosome scenery `topo.csp.ribosome`, and the timing core `topo.csp.kinetics`.

```bash
topo-csp -f csp.ini
python -m topo.csp -f csp.ini
```

See **[Tutorial 8](../../tutorials/08_ribosome_synthesis/)** for runnable, validated
examples — 4c5c (smoke `csp_debug.ini`, full-length `csp_val.ini`) and P0CX28 — and
[`review/TODO.md`](../../review/TODO.md) for the done/remaining matrix.

## What it does

For nascent length `L = L0 .. L_max`, each residue is added through three sub-stages,
each run for an **exponentially-sampled** dwell time (first-passage-time sampling):

| stage | biology | dwell-time mean (s) | C-terminus restraint |
|-------|---------|---------------------|----------------------|
| 1 | peptidyl transfer / A-site delivery | `time_stage_1` | → A-anchor |
| 2 | translocation begins | `time_stage_2` | → A-anchor |
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
to `topo.csp.core.run_length` (build-once-subset contacts, coordinate
seeding, rigid-ribosome scenery, tunnel wall, minimize/run/finalize). CSP drives the
position-restraint path (target switchable A↔P), so it forces `trna_tether = no`. The
supplied ribosome PDB is **always** loaded as rigid (mass-0) scenery (excluded volume +
electrostatics on) — there is no `rigid_ribosome` flag — and the tunnel wall is on with
its plane auto-derived from the structure.

| file | role |
|------|------|
| `kinetics.py` | pure timing core (no OpenMM): codon tables, FPT sampling, `scale_factor`→steps, 3-stage split |
| `protocol.py` | `RunParams` / `CSPConfig`, `read_csp_config` (INI), `run_continuous_synthesis`, `csp()` CLI |
| `cylinder.py` | parallel runner for the **cylinder** ribosome model (analytic exit tunnel, no beads): `CylinderParams` / `CylinderConfig`, `read_cylinder_config`, `run_cylinder_synthesis`, `cylinder()` CLI (`topo-cylinder`); same kinetics as CSP, one MD segment per residue |
| `__init__.py`, `__main__.py` | public API; `python -m topo.csp` |

## Public API

```python
from topo.csp import run_continuous_synthesis, read_csp_config, RunParams, kinetics

cfg = read_csp_config("csp.ini")
run_continuous_synthesis(cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max,
                         out_root=cfg.outdir, mrna=cfg.mrna, codon_time_table_path=cfg.codon_time_table_path,
                         domain_def=cfg.domain_def,
                         stride_output_file=cfg.stride_output_file, params=cfg.params)
```

## Control file (`csp.ini`)

A single `[OPTIONS]` section. Required: `pdb_file`, `ribosome`, `domain_def` (the
protein's contact-nscale definition). `L0` (default `1`) and `L_max` (default = full
length) are optional. Per-codon timing requires both `mrna` and a `codon_times` table
path (there is no bundled default -- pick one under `assets/csp/codon_dwell_times/`).
Setting `mrna = fastest`, `slowest` or `median` auto-builds a synonymous-codon mRNA from
the protein + table (see the `mrna` key below).

**Kinetic keys**

| key | meaning |
|-----|---------|
| `mrna` | mRNA sequence file (raw nucleotides; one codon/residue + 1 stop), **or** `fastest`/`slowest`/`median` to auto-generate a synonymous-codon mRNA — every residue encoded by its fastest/slowest/median-dwell-time codon per the `codon_times` table, written next to the PDB as `mrna_<mode>.txt`. A real filename must not be `fastest`/`slowest`/`median`. |
| `codon_times` | codon-time table path (per-codon; required, no bundled default -- pick one under `assets/csp/codon_dwell_times/`) **or** a positive number of s (uniform codon time). A table filename must not be a bare number. |
| `scale_factor` | in-vivo seconds → in-silico ns compressor |
| `time_stage_1` / `time_stage_2` | mean peptidyl-transfer / translocation dwell (s) |
| `random_seed` | seed for the FPT sampler (reproducible schedules) |
| `max_steps_per_stage` / `min_steps_per_stage` | clamp each stage's step count (testing) |
| `ejection_steps` / `dissociation_steps` | post-synthesis free runs (0 = skip) |
| `resume` | `auto` (default; resume iff an interrupted run is present under `outdir`), `yes` (require a resumable run, else error), `no` (always fresh). See *Resume* below. |

**MD / ribosome keys** (configure the shared `RunParams`; `n_steps` is **not** used — step
counts come from the kinetics): `dt`, `ref_t`, `tau_t`, `nstout`, `device`, `ppn`,
`constraints`, `restraint_k`, `minimize`, `tunnel_wall`, `nascent_ev_radii`,
`ribo_free_mask`, `ribo_free_pdb`. (Output is always nascent-only; `tunnel_wall_x0` and `tunnel_wall_k` are not keys -- the plane is auto-derived from the structure and the stiffness is a fixed model constant.)

**Flexible exit-tunnel loop.** `ribo_free_mask = SEG : lo - hi` frees a portion of one
ribosomal protein (E. coli `L24` / eukaryotic `L26`) as a topo-style Go loop while the
rest of the ribosome stays mass-0 (`append_flexible_l24_loop`); `ribo_free_pdb` gives the
all-atom chain its native contacts are built from. See the usage guide's *Flexible
exit-tunnel loop* section.

**Units:** OpenMM defaults — nm, ps, kJ/mol, K, kJ/mol/nm². In-vivo dwell times are
in **seconds**.

## Output

`<outdir>/L_<L>/stage_<1,2,3>/` per residue (each a standalone topo run: `traj.dcd`,
`traj.log`, `traj.psf`, `traj.chk`, `traj_final.pdb`, `traj_runinfo.log`,
`native_1_<L>.pdb`), plus `ejection/` and optional `dissociation/`. Stage 3's
`traj_final.pdb` seeds the next residue.

At the output root: `dwell_times.dat` (the **immutable plan** — the full per-residue
3-stage schedule, drawn once up front, with the PTC restraint geometry in its `#PTC`
header) and `progress.log` (the append-only `L_XXX RUNNING`/`DONE` **status** used by
resume). A launched run prints its **exact total planned step count** before the first
MD step.

## Resume

A production run is hours to days and today survives no interruption. Re-invoking
`topo-csp` on an interrupted `outdir` **continues from the last completed residue** with
a kinetic schedule *identical* to the uninterrupted run — the schedule and PTC geometry
are re-read from `dwell_times.dat` (no RNG redraw, no SLSQP re-solve), and the seed
conformation is reloaded from the previous residue's `traj_final.pdb`. The resume unit is
the **residue**: the one in flight at the crash is dropped and re-run from its persisted
schedule row (≤ 3 stages of redone work).

Controlled by the `resume` key (`auto`/`yes`/`no`) or the `--no-resume`/`--fresh` CLI
flag. `auto` (default) resumes iff a `progress.log` is present. Before resuming, every
prior length's `traj_final.pdb` must be present on disk — a missing one aborts with an
actionable message rather than silently leaving a hole. **Extending a run** (a larger
`L_max`) is a fresh run, not a resume: the schedule is fixed at first launch. The MD
micro-trajectories are *not* bit-reproducible across a restart (the thermostat RNG is
unseeded); the guarantee is the kinetic schedule + conformational continuity. See
`topo.csp.resume` and [`review/F_RESUME.md`](../../review/F_RESUME.md).

## Movie

Stitch the per-stage trajectories into one VMD-playable movie (chain grows stage by
stage; static ribosome overlaid):

```bash
topo-csp-movie -o synth_out --ribosome ribosome_trunc.pdb
vmd -e synth_out/movie.tcl
```

`topo.csp.movie` discovers the nested `L_<L>/stage_<s>/` segments (`stitch_movie`) and
stitches them with the shared `stitch_segments` core.

## Not yet ported

CHARMM PSF/TOP/PRM/COR ingestion (O'Brien's exact systems); multi-trajectory
multiprocessing; literal mid-residue peptide-bond toggling / explicit A/P tRNA bonded
geometry; quantitative validation. (Cross-length **resume** is now supported — see
*Resume* above.) See [`review/TODO.md`](../../review/TODO.md).
