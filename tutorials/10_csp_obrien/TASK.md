# TASK — O'Brien CSP port checklist

Status of each task in porting `continuous_synthesis_v6.py` → `topo.csp` + Tutorial 10.
`[x]` done · `[~]` partial/mapped · `[ ]` remaining.

## Module: `topo/csp/`

- [x] `kinetics.py` — codon timing core (pure, no OpenMM)
  - [x] `read_trans_times` — codon → seconds table (T→U normalised)
  - [x] `read_mrna` — split to codons, truncate at stop codon
  - [x] `codon_mfpt_list` / `uniform_mfpt_list` — per-residue intrinsic mFPT
  - [x] `sample_fpt` — exponential FPT sampling (v6 `sample_fpt_dist`)
  - [x] `seconds_to_steps` — `t·1e9/scale_factor / dt` (v6 mapping)
  - [x] `stage_dwell_times` — 3-stage split + stage-2 traffic correction + fast-codon floor
  - [x] `stage_steps` — dwell→steps with test clamps (`max/min_steps_per_stage`)
  - [x] `ribosome_traffic_times` — external binary hook, graceful fallback
  - [x] `build_mfpt_lists` — assemble intrinsic/real lists (uniform or per-codon)
- [x] `csp.py` — runner + INI + CLI
  - [x] `CSPParams` (composes `ElongationParams` + kinetic fields)
  - [x] `run_continuous_synthesis` — precompute once; 3-stage per-residue loop
  - [x] per-residue dwell-time table logged (mirrors v6 `*.out`)
  - [x] ejection + dissociation post-synthesis phases
  - [x] `CSPConfig` + `read_csp_config` (INI superset of `elongate.ini`)
  - [x] `csp()` CLI (`topo-csp -f csp.ini`, `-o`/`--device` overrides)
- [x] `__init__.py` / `__main__.py` (public API, `python -m topo.csp`)
- [x] `topo-csp` console script in `pyproject.toml`
- [x] `make_movie.py` — VMD movie of the per-stage trajectories (`topo-csp-movie`),
      reusing `topo.translation.make_movie.stitch_segments`
- [x] `README.md` module reference

## Reuse (not re-implemented)

- [x] `run_length`, `precompute_contacts`, `read_anchor`, `ElongationParams`,
      `TUNNEL_AXIS`, `TRNA_TETHER_BOND_NM` (from `topo.translation.elongate`)
- [x] `load_ribosome`, tunnel wall, excluded volume + electrostatics
      (from `topo.translation.ribosome`)
- [x] build / setup / finalize (from `topo.engine`); `strtobool` (`topo.utils.config`)

## Fixes

- [x] `add_cterm_restraint`: per-particle `k` (was global) — lets the position
      restraint coexist with the tunnel wall's global `k` in v2 (CSP's mode)

## Tutorial: `tutorials/10_csp_obrien/`

- [x] Inputs copied from Tutorial 7 (`P0CX28_clean.pdb`, `_stride.dat`,
      `domain.yaml`, `ribosome_trunc.pdb`)
- [x] `trans_times.txt` — Fluitt *E. coli* codon-time table (copied)
- [x] `P0CX28_mrna.txt` — synthetic mRNA (back-translated PDB sequence, +1 stop)
- [x] `csp.ini` — demo config (`L0=5..L_max=12`, `max_steps_per_stage=667` ≈ 2000/AA,
      v2 ribosome + wall, ejection)
- [x] `README.md` / `PLAN.md` / `TASK.md`
- [x] End-to-end run verified (~30 s CPU; chain 5→12 correct)

## Remaining

- [ ] CHARMM PSF/TOP/PRM/COR ingestion (run O'Brien's exact 4c5c + 50S systems)
- [ ] Multi-trajectory multiprocessing (`num_traj`/`tpn`/`ppn`, GPU device fan-out)
- [~] Literal 3-stage mechanics (peptide-bond toggling, explicit A/P tRNA bonded
      geometry) — currently mapped via the A→P restraint switch
- [ ] `restart = 1` resume of a partial trajectory
- [ ] Working external `ribosome_traffic` binary (exits 127 here; degrades to no traffic)
- [ ] Quantitative validation vs. v6 reference outputs
- [x] Tutorial entry in `tutorials/README.md` index
- [ ] Docs site page
- [ ] **Stage-1 blow-up at some residue additions** — 5/306 stages explode (PotE ~1e13)
      when the new bead is seeded on top of an existing bead; self-recovers next stage
      but corrupts those frames. See [`OBSERVATIONS.md`](OBSERVATIONS.md) #1. Fix:
      robustify seed placement / minimization (soft-core or large-force re-seed).
