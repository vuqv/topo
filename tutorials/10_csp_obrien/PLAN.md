# PLAN — Porting the O'Brien Continuous Synthesis Protocol to topo

This document records **what was ported, how, and what remains**. The runnable result
is the `topo.csp` package and this Tutorial 10.

## Source

O'Brien *Continuous Synthesis Protocol*, script **v6**
(`continuous_synthesis_v6.py`, 1697 lines; Yang Jiang, Dan Nissley, Ed O'Brien),
with the `example/continuous_synthesis/` inputs. The reference is OpenMM + parmed +
CHARMM (PSF/TOP/PRM/COR) with multiprocessing over trajectories.

## Goal

Reproduce the **scientific content** of v6 — codon-resolved translation kinetics, the
three-stage elongation cycle, the in-vivo→in-silico time mapping, ribosome traffic,
ejection/dissociation — in **topo style** (INI `[OPTIONS]` config, dataclasses,
`topo-*` console script, per-length output folders, OpenMM nm/ps/kJ units), while
**reusing** the existing `topo.translation` machinery rather than duplicating it.

## Design

A new package **`topo/csp/`** that adds only the O'Brien kinetics + 3-stage loop and
delegates all MD to `topo.translation.elongate.run_length`.

```
topo/csp/
├── kinetics.py   # PURE timing core (no OpenMM): codon times, FPT sampling,
│                 #   scale_factor → steps, 3-stage split, ribosome-traffic hook
├── csp.py        # CSPParams/CSPConfig, read_csp_config (INI), run_continuous_synthesis,
│                 #   csp() CLI — the outer loop calling run_length 3×/residue
├── __init__.py   # public API re-exports
├── __main__.py   # python -m topo.csp
└── README.md     # module reference
```

Console script `topo-csp = "topo.csp.csp:csp"` (registered in `pyproject.toml`).

### The 3 stages mapped onto `run_length`

`run_length(L, p_anchor=, a_anchor=, prev_final=, seed_override=, restrain=,
out_subdir=, n_steps_override=)` builds the length-`L` model, seeds, restrains
residue `L-1` to `p_anchor`, minimizes, runs, returns the final `(L,3)` coords. CSP
calls it three times per residue, chaining via `seed_override` and switching the
restraint target to reproduce A→P translocation:

| stage | dwell mean (s) | `run_length` call |
|-------|----------------|-------------------|
| 1 peptidyl transfer | `time_stage_1` | new residue placed at A-anchor; restrain → A; `n_steps = step1` |
| 2 translocation | `time_stage_2` (+ traffic) | `seed_override = stage1 final`; restrain → A; `n_steps = step2` |
| 3 tRNA binding | codon − stage 1 − stage 2 | `seed_override = stage2 final`; restrain → **P**; `n_steps = step3` |

Stage 3's final coords seed the next residue. The cold-start segment (`L == L0`) is
laid down the tunnel from the P-anchor. CSP uses the **position-restraint** path
(target switchable A↔P), so it forces `trna_tether = no`; `rigid_ribosome`,
`tunnel_wall`, excluded volume + electrostatics work exactly as in Tutorial 7's v2.

### Time → steps (v6 formulas, in `kinetics.py`)

```
dwell_s   = expovariate(1 / mean_s)              # first-passage-time sampling
t_sim_ns  = dwell_s * 1e9 / scale_factor         # in-vivo → in-silico
steps     = int(t_sim_ns / (dt_ps * 1e-3))       # → integration steps
```

with the per-residue means from the codon table (`intrinsic[L]`) and the stage-2
traffic correction `real[L-1] − intrinsic[L-1]` (guarded non-negative, v6 lines 76-80).

## Status — DONE

- [x] **Kinetics core** (`topo/csp/kinetics.py`): `read_trans_times`, `read_mrna`
      (stop-codon aware), `codon_mfpt_list`, `uniform_mfpt_list`, `sample_fpt`
      (exponential), `seconds_to_steps`, `stage_dwell_times` (3-stage split +
      traffic correction + fast-codon floor), `stage_steps` (+ test clamps),
      `ribosome_traffic_times` (external binary hook with graceful fallback),
      `build_mfpt_lists`. Unit-checked.
- [x] **Runner** (`topo/csp/csp.py`): `CSPParams` (composes `ElongationParams` +
      kinetic fields), `run_continuous_synthesis` (precompute + ribosome loaded
      once; per-residue 3-stage loop; per-residue dwell-time table logged like v6's
      `*.out`; ejection + dissociation phases).
- [x] **INI + CLI**: `CSPConfig`, `read_csp_config` (superset of `elongate.ini`,
      with v6-style validation), `csp()` console entry; `topo-csp` registered.
- [x] **topo-native inputs** (user decision): a PDB + topo's CG builder (no CHARMM
      ingestion). Tutorial reuses `P0CX28` + `ribosome_trunc.pdb` from Tutorial 7,
      plus the Fluitt `trans_times.txt` and a generated `P0CX28_mrna.txt`.
- [x] **Bug fix** — `add_cterm_restraint` (in `topo/translation/elongate.py`) now
      uses a **per-particle** `k` instead of a global, so it can coexist with the
      tunnel wall's global `k` (the position-restraint + wall combo CSP needs; a
      latent collision Tutorial 7 never hit because it uses the tRNA tether).
- [x] **Tutorial 10** (`csp.ini` + inputs + README/PLAN/TASK), runs end-to-end in
      ~30 s on CPU; chain verified to grow one residue per residue-step (5→12).

## Status — REMAINING (out of scope this pass)

- [ ] **CHARMM ingestion** — load O'Brien's exact PSF/TOP/PRM/COR (4c5c + 50S) so the
      *identical* reference systems can be run. topo currently builds its own CG
      model from a PDB; a CHARMM loader is a separate, larger effort.
- [ ] **Multi-trajectory parallelism** — v6's `num_traj`/`tpn`/`ppn` multiprocessing
      and GPU device fan-out. topo runs one trajectory serially; launch replicas
      externally (vary `random_seed`) for now.
- [ ] **Literal 3-stage mechanics** — mid-residue peptide-bond toggling and explicit
      A/P tRNA bonded geometry (bond + angles + improper to the tRNA R bead). We map
      translocation via the A→P restraint switch; the peptide bond is present from
      stage 1.
- [ ] **`restart = 1`** — resume a partially-finished trajectory by scanning prior
      stage outputs (v6 lines 1423-1457).
- [ ] **Working `ribosome_traffic`** — the external binary is wired up but exits 127
      (missing libs) in this environment; the run degrades gracefully to no traffic.
- [ ] **Quantitative validation** against v6 reference outputs (dwell-time
      distributions, co-translational folding observables).
- [ ] **Stage-1 seeding blow-up** (robustness bug, found via the movie): in a full
      run, ~5/306 stages — almost all **stage 1** (new-residue delivery) — explode to
      PotE ~10^13 kJ/mol when the new bead is seeded nearly on top of an existing bead;
      the next stage's minimization recovers it (final structures OK) but the affected
      frames are garbage. Fix: overlap-aware seed placement, soft-core/capped
      minimization, or a post-minimize large-force re-seed in `run_length`. Details in
      [`OBSERVATIONS.md`](OBSERVATIONS.md) #1.

## Verification performed

1. Unit checks of `kinetics.py` (table round-trip, `seconds_to_steps` hand value,
   capped `stage_steps`, uniform + per-codon `build_mfpt_lists`, traffic fallback).
2. `topo-csp -f csp.ini` end-to-end: config echo, ribosome (4576 beads), one-time
   contact precompute, per-residue codon-resolved 3-stage schedule, `L_005..L_012`
   × `stage_{1,2,3}` + `ejection/` outputs, finished ~30 s on CPU.
3. Chain growth: `L`'s `stage_3/traj_final.pdb` has exactly `L` CA atoms (5, 8, 12).
4. `pip install -e .` then `topo-csp -h` works.
