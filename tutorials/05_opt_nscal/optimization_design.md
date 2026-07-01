# Nscale optimization — design / tracking note

**Status:** production. The per-domain / per-interface *Q* scorer
([`topo.analysis.native_contacts`](../../topo/analysis/native_contacts.py)) and
the optimizer ([`topo.optimize`](../../topo/optimize/optimize.py), exposed as the
`topo-optimize` console command; the in-folder [`optimization.py`](optimization.py)
is now a thin shim) are both **implemented and tested** end-to-end (real
multicopy MD → split → score → decide). Each round writes one Q time series per
trajectory next to its DCD (`round_N/traj/Q_<k>.csv`, paired with `traj_<k>.dcd`;
columns `frame, Q_<unit>...`) for inspection. The only remaining open item is
optimization-level restart/resume (below). This note remains the design reference.

The goal: automatically choose the per-domain and per-interface `nscale`
(*n*<sub>scale</sub>) in `domain.yaml` — the smallest value on a discrete ladder
that keeps each domain/interface folded across many independent trajectories —
instead of hard-coding `nscale` by hand.

The *Q* method (how the fraction of native contacts is computed) and its
implementation both live in
[`topo.analysis.native_contacts`](../../topo/analysis/native_contacts.py) — see
its module docstring for the definition. It exposes importable functions plus a
`python -m topo.analysis.native_contacts` CLI, and emits one *Q* column per
domain and per interface. It was generalized from an earlier AK-specific
`calc_OP.py` template.

---

## Architecture

A dedicated **`optimization.py`** is the top-level driver. It owns the search
logic and *calls* the existing tools as sub-steps — it does **not** extend
`run_simulation.py`. Separation of concerns:

```
optimization.py  (the orchestrator — this is what the user runs)
   ├─ owns the nscale LADDER and the per-unit level state
   ├─ each round: writes domain.yaml, then INVOKES run_simulation.py (multi-copy MD)
   ├─ after each round: scores Q per domain / per interface (topo.analysis.native_contacts)
   ├─ makes the decision: stable? climb ladder? converged? fallback?
   └─ writes the final calibrated domain.yaml
```

`run_simulation.py` stays exactly as-is (the MD engine); `optimization.py`
invokes it once per round (as a subprocess, or by importing the same `topo`
entry points it uses).

## Components used (current setup)

| Component | Role in the loop |
|-----------|------------------|
| **`optimization.py`** (implemented) | **Top-level driver:** reads `optimize.ini`; owns ladder, per-unit levels, per-round Q, stability decision, convergence/fallback, final `domain.yaml`. |
| **`optimize.ini`** | Minimal user config (essentials + optimizer controls); the driver fills implicit defaults and expands it into each `round_N/md.ini`. |
| `domain.yaml` | Initial domains: `residues` + `class`. The driver writes a fresh `round_N/domain.yaml` with the current `nscale` each round. |
| `run_simulation.py` | The MD engine, **invoked by `optimization.py`** each round with the generated `round_N/md.ini`. **`n_copies = ntraj`** gives the N independent trajectories in **one** run (Tutorial 4). |
| **`topo.split_chains`** (package) | Splits the combined multi-copy DCD into N per-copy DCDs (`traj_<k>.dcd`). Memory-bounded chunked streaming (handles DCDs too large for memory); `optimization.py` calls it in-process. Also a `python -m topo.utils.multichain` CLI. |
| **`topo.analysis.native_contacts`** (package module) | Per-domain / per-interface *Q* from `domain.yaml` + reference PDB + CG PSF/DCD. `optimization.py` imports it in-process (builds native contacts once, reuses each round); also a `python -m topo.analysis.native_contacts` CLI. Generalized from an earlier AK-specific `calc_OP.py` template (since removed). |

## Inputs

The optimizer is driven by a single **minimal `optimize.ini`** (not a full
`md.ini`) with one flat `[OPTIONS]` section. The optimizer consumes its own
control keys (`ntraj`, `q_threshold`, `frame_fraction`, `max_rounds`,
`min_contacts`); every other key is an essential simulation parameter
(`pdb_file`, `domain_def`, `md_steps`, sampling, `ref_t`) passed through to each
round's `md.ini`.
Everything else (`dt`, `model`, thermostat, `device`, `ppn`, `minimize`, output
naming, ...) comes from the optimizer's **implicit protocol defaults**
(`IMPLICIT_DEFAULTS` in `optimization.py`). Notably `dt = 0.015` ps (the model's
15 fs timestep), *not* the bare `SimulationConfig` default of 0.01.

Each round the driver expands `optimize.ini` into a complete `round_N/md.ini`
(implicit defaults ← `optimize.ini` ← per-round values) and writes
`round_N/domain.yaml` with the current nscales.

```
optimize.ini             # MINIMAL config (see above); the only user input
  pdb_file               #   all-atom reference (native contacts + geometry)
  domain_def             #   initial domains: residues + class (nscale optimized)
LADDER                   # class -> list of nscale levels (Table 1, below)
ntraj          = 10      # independent trajectories per round (= n_copies)
Q_THRESHOLD    = 0.6688  # a frame is "folded" for a unit if Q > this
FRAME_FRACTION = 0.98    # a traj is "stable" if >= 98% of frames are folded
```

The nscale ladder, keyed by structural class (last value = median/level-3
fallback, used only when level 5 still fails):

```
LADDER = {
  alpha      : [1.1954, 1.4704, 1.7453, 2.0322, 2.5044,  fallback=1.7453],
  beta       : [1.4732, 1.8120, 2.1508, 2.5044, 2.5044,  fallback=2.1508],
  alpha_beta : [1.1556, 1.4213, 1.6871, 1.9644, 2.5044,  fallback=1.6871],
  interface  : [1.2747, 1.5679, 1.8611, 2.1670, 2.5044,  fallback=1.8611],
}
```

## Main loop (pseudo-code) — body of `optimization.py`

```
# ---- setup ----
domains    = read domains from domain.yaml          # each has: residues, class
interfaces = all pairs (D1, D2) of domains          # class = "interface"
units      = domains + interfaces                    # everything we must stabilize

level[u] = 0  for every unit u                       # start at level 1 (index 0)
MAX_LEVEL = 5

# ---- iterate rounds ----
for round in 1 .. MAX_LEVEL+1:

    # 1. assign current nscale to every unit from its class ladder
    for u in units:
        if round <= MAX_LEVEL:
            nscale[u] = LADDER[class(u)][ level[u] ]
        else:
            nscale[u] = LADDER[class(u)].fallback   # could-not-stabilize case

    # 2. write a domain.yaml carrying these nscales
    write_domain_yaml(round_dir/domain.yaml,
                      intra = {D: nscale[D]      for D in domains},
                      inter = {pair: nscale[pair] for pair in interfaces})

    # 3. run ntraj INDEPENDENT trajectories in ONE simulation (multi-copy).
    #    Set n_copies = ntraj in md.ini -> run_simulation.py replicates the chain
    #    into ntraj non-interacting copies (Tutorial 4), then split per chain.
    run_simulation.py -f round_dir/md.ini          # md.ini: n_copies = ntraj,
                                                    #         domain_def = round_dir/domain.yaml
    split_chains(round_dir/traj/traj.dcd,  # -> traj_0.dcd .. traj_{ntraj-1}.dcd
                         [traj_k.dcd for k in range(ntraj)])   # in-process, streaming

    # 4. score every trajectory -> per-unit Q time series.
    #    optimization.py imports topo.analysis.native_contacts and scores in-process
    #    (native contacts built once from the reference, reused every round), writing
    #    one CSV per trajectory next to its DCD: round_dir/traj/Q_<k>.csv with columns
    #    frame, Q_<dom>..., Q_<d1>-<d2>.... (Standalone equivalent of the scorer:
    #    python -m topo.analysis.native_contacts -d domain.yaml -r protein.pdb ...)
    for k in 0 .. ntraj-1:
        Q[k] = score(round_dir/traj/traj_k.dcd) -> round_dir/traj/Q_k.csv

    # 5. decide stability for each unit
    if round == MAX_LEVEL+1:
        break                                       # fallback already applied; stop

    all_stable = true
    for u in units:
        if is_stable(u):
            keep level[u] (and its nscale)        # stable: frozen at this level
        else:
            all_stable = false
            level[u] = level[u] + 1                 # unstable: climb ladder next round

    # 6. converged?
    if all_stable:
        break

# ---- result ----
write_domain_yaml(final/domain.yaml, nscale)      # the calibrated model
log "## Final nscales:", nscale
```

## Stability sub-routine

```
function is_stable(unit u):
    # u is stable only if ALL ntraj trajectories keep it folded >=98% of the time
    for i in 1 .. ntraj:
        Q_series        = read column "Q_<u>" from round_dir/Q_i.csv
        folded_fraction = fraction of frames where Q_series > Q_THRESHOLD
        if folded_fraction < FRAME_FRACTION:
            return false        # one bad trajectory => unit not stable
    return true
```

## Key behaviors

- **Per-unit, not global.** Each domain and interface climbs its ladder
  *independently*: if round N finds domain A unstable but B stable, only A's level
  increments; B stays put. Hence `level[]` is per-unit.
- **Monotonic climb.** Nscales only ever increase — we search for the
  *smallest* nscale that stabilizes each unit.
- **Two stop conditions.** (a) all units stable → done; or (b) a unit still
  unstable after level 5 → pinned to the median (level-3) nscale and accepted.
- **`class` drives the choice; input `nscale` is a placeholder** overwritten
  every round.

## Table 1 — nscale levels

| Structural Class | Level 1 | Level 2 | Level 3 | Level 4 | Level 5 |
| ------ | ------ | ------ | ------ | ------ | ------ |
| α          | 1.1954 | 1.4704 | 1.7453 | 2.0322 | 2.5044 |
| β          | 1.4732 | 1.8120 | 2.1508 | 2.5044 | 2.5044 |
| α/β        | 1.1556 | 1.4213 | 1.6871 | 1.9644 | 2.5044 |
| Interface  | 1.2747 | 1.5679 | 1.8611 | 2.1670 | 2.5044 |

(Calibrated on a training set of 18 small single-domain proteins.)

## Open items for the implementation

1. ~~**Generalize `calc_OP.py`** to read `domain.yaml` and emit one *Q* column
   per domain and per interface.~~ **Done** — implemented and tested as the
   package module [`topo.analysis.native_contacts`](../../topo/analysis/native_contacts.py)
   (interfaces = all unordered domain pairs; non-touching pairs report
   `Q = NaN`). The orchestrator imports it directly.
2. ~~**Per-round run config.**~~ **Done** — `optimization.py` templates a
   `round_N/md.ini` per round (`n_copies = ntraj`, `domain_def`, `output_dir`,
   `ref_t`), invokes `run_simulation.py`, then splits in-process with
   `topo.split_chains`. Copy independence is guaranteed by the
   multi-copy primitive (Tutorial 4).
3. ~~**`nscale` is required by the parser**~~ **Not an issue.** The parser reads
   `nscale` with a hard lookup ([nonbonded.py:649](../../topo/utils/nonbonded.py#L649)),
   but step 2 always writes a concrete value for every domain and interface from
   the current ladder decision, so the key is never missing.
4. **Restart / resume — still open.** Long scans should be resumable
   round-by-round (the legacy `opt_nscal_v3.py` supports `-r`). The prototype
   starts fresh each invocation; decide how to checkpoint `level[]` and completed
   rounds.
5. ~~**STRIDE caching.**~~ **Done.** STRIDE depends only on the fixed structure,
   so it is identical every round. If `stride_output_file` is set in `[OPTIONS]`
   it is resolved to an absolute path and passed through every round (STRIDE never
   runs). Otherwise round 1 runs `cwd=round_1/`, and the optimizer reuses the
   `{pdb_stem}_stride.dat` the model build cached there for rounds 2+ — STRIDE runs
   once instead of `max_rounds` times.
