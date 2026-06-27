# Tutorial 6 — Temperature annealing & quenching

**Goal:** drive a simulation through a **temperature protocol** instead of a
single constant temperature — hold the protein hot enough to unfold, then bring
it back down to `ref_t` to watch it **refold**. You will learn both supported
protocols (a delta **T-jump quench** and a **linear cooling ramp**), every config
key that controls them, and how the schedule interacts with restarts.

**Time:** two short runs, ~2–3 seconds each.

**Prerequisite:** do [Tutorial 1](https://vuqv.github.io/topo/tutorials/01_single_domain.html)
first — this reuses the same single-domain protein (`P0CX28`) and its calibrated
`domain.yaml`.

---

## The idea: a temperature *protocol*

An ordinary run (Tutorials 1–4) is **equilibrium**: one Langevin thermostat held
at `ref_t` for all `md_steps`. Annealing generalizes that to a **schedule** of
`(temperature, n_steps)` stages whose step counts always sum to `md_steps`:

- **Equilibrium** (`anneal = no`, the default) → one stage: `[(ref_t, md_steps)]`.
- **Annealing** (`anneal = yes`) → hold at a high temperature `t_high`, then
  bring the thermostat down to `ref_t` and finish the run there.

The same runner (`topo-mdrun`) handles both — only the temperature schedule
differs. `ref_t` is **always** the low / refold temperature; there is **no
separate `t_low` key**, so you never have to keep two "low temperature" settings
in sync.

When the run starts, the runner prints the exact schedule it will execute, e.g.

```
Temperature protocol: 600 K x 3000 -> 300 K x 3000
```

## Two ways down: jump vs. linear

`anneal_ramp` selects how the temperature goes from `t_high` to `ref_t`:

| `anneal_ramp` | What happens | Use it for |
|---------------|--------------|------------|
| `jump` (default) | Hold `t_high`, then **instantaneously** set the thermostat to `ref_t` (a delta T-jump). Folding then happens at a single, well-defined temperature. | **Folding kinetics / mechanism.** Clean folding times and pathways because rates aren't convolved with a changing T. Mirrors an experimental T-jump. |
| `linear` | Hold `t_high`, then **cool gradually** to `ref_t` over `anneal_ramp_steps`, in `anneal_ramp_increments` discrete steps, then hold `ref_t`. | **Refolding yield / structure recovery.** Slow cooling avoids kinetic traps (classic simulated annealing toward the native minimum). *Not* for clean kinetics — cooling and folding overlap. |

**Recommendation:** for studying how (and how fast) the protein refolds, use the
**delta `jump`**. Reach for `linear` only when the goal is "recover the native
state at all," not "measure the folding rate."

## All annealing options

All keys live in the single `[OPTIONS]` section, alongside the usual settings.
They are read **only when `anneal = yes`**.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `anneal` | bool | `no` | Turn the annealing protocol on. `no` → ordinary constant-`ref_t` equilibrium. |
| `t_high` | float [K] | — (**required** when `anneal = yes`) | High / unfolding temperature to hold before cooling. |
| `anneal_steps` | int | `0` | Steps held at `t_high` before cooling begins. |
| `anneal_ramp` | str | `jump` | `jump` (instantaneous drop to `ref_t`) or `linear` (gradual cool-down). |
| `anneal_ramp_steps` | int | `0` | **`linear` only.** Steps spent ramping `t_high → ref_t`. |
| `anneal_ramp_increments` | int | `20` | **`linear` only.** Number of discrete temperature steps in the ramp; the last increment lands exactly on `ref_t`. |

Reused, not new: `ref_t` (the low temperature), `tau_t` (thermostat coupling),
`md_steps` (the **total** across the whole protocol).

**The arithmetic is on you, lightly:** the steps must fit inside `md_steps`. The
runner spends `anneal_steps` at `t_high`, then (for `linear`) `anneal_ramp_steps`
cooling, and the **remainder** at `ref_t`:

```
remainder_at_ref_t = md_steps - anneal_steps - (anneal_ramp_steps if linear else 0)
```

If `anneal_steps (+ anneal_ramp_steps)` exceeds `md_steps`, the run stops
immediately with an error telling you to lengthen `md_steps` or shorten the
hold/ramp. If the remainder is zero, that's fine — the protocol just ends right
as it reaches `ref_t` (you'd typically leave some steps for the system to settle
and refold at `ref_t`).

## ⚠️ The most important physical knob: hold length vs. thermostat coupling

A Langevin thermostat doesn't change the temperature instantly even on a `jump`:
the kinetic energy relaxes toward the new setpoint over roughly **`1 / tau_t`**.
So the *setpoint* jumps, but the *system* takes ~`1/tau_t` to follow.

That has two consequences you must respect:

1. **`anneal_steps` must be ≫ the relaxation time** (`1/tau_t`, converted to
   steps via `dt`), or the system never actually reaches `t_high` and won't
   fully unfold.
2. **The hold must be long enough to lose native memory** — verify the protein
   genuinely unfolds during the hold (fraction of native contacts *Q* → 0, `Rg`
   plateaus). Otherwise refolding is biased by leftover native structure.

> **These tutorial configs cheat for speed:** they use `tau_t = 1.0` ps⁻¹
> (relaxation ≈ 1 ps ≈ 67 steps at `dt = 0.015`), so a 3000-step hold easily
> reaches 600 K. **Production runs typically use `tau_t ≈ 0.01` ps⁻¹**
> (relaxation ≈ 100 ps), so you would need an `anneal_steps` of *many* ×
> 100 ps — e.g. a 15 ns hold — for the same effect. Scale the hold with your
> friction, not with the demo numbers.

## Files in this folder

| File | Role |
|------|------|
| `P0CX28_clean.pdb` | Input structure (same single domain as Tutorial 1). |
| `domain.yaml` | Calibrated contact `strength` — folded at 300 K, unfolds at 600 K. |
| `md.ini` | **Delta T-jump quench** (`anneal_ramp = jump`). |
| `md_linear.ini` | **Linear cooling ramp** (`anneal_ramp = linear`). |
| `run_simulation.py` | The runner shim (`= python -m topo.mdrun`). |

## Step-by-step

### 1. Run the delta T-jump quench
```bash
python run_simulation.py -f md.ini
```
The console echoes the protocol:
```
Temperature annealing on: hold 600.0 K for 3000 steps, then jump to ref_t = 300.0 K
...
Temperature protocol: 600 K x 3000 -> 300 K x 3000
```
Outputs land in `traj_jump/` (so they don't collide with the linear run).

Inspect the temperature column of the log — you should see it hover near 600 K
for the first 3000 steps, then drop to ~300 K:
```python
import numpy as np
from topo.reporter.topo_reporter import readOpenMMReporterFile
d = readOpenMMReporterFile("traj_jump/traj.log")
S, T = np.array(d["Step"]), np.array(d["Temperature (K)"])
print("mean T during hold  :", round(T[S <= 3000].mean()))   # ~600 K
print("mean T after quench :", round(T[S >  3000].mean()))   # ~300 K
```
(Instantaneous CG temperatures are noisy for a small chain — judge the *mean* of
each block, not single frames.)

### 2. Run the linear cooling ramp
```bash
python run_simulation.py -f md_linear.ini
```
Now the schedule has a staircase of cooling stages (outputs in `traj_linear/`):
```
Temperature protocol: 600 K x 1500 -> 570 K x 300 -> 540 K x 300 -> ... -> 300 K x 300 -> 300 K x 1500
```
Plot `T` vs `Step` from `traj_linear/traj.log` and you'll see the staircase
descend from 600 K to 300 K and then flatten.

### 3. (Optional) Confirm unfolding/refolding
For a real study you'd track the fraction of native contacts *Q* with
`topo.analysis` (see Tutorial 5's machinery) over the trajectory: *Q* should fall
toward 0 during the 600 K hold and climb back toward 1 after the quench. With the
tiny demo `md_steps` here this is only illustrative — increase `md_steps` (and
`anneal_steps`) for a meaningful folding curve.

## Annealing + restart

Restarts (Tutorial 3) compose with annealing: the schedule is defined over
**absolute** step counts, so a restart **resumes mid-schedule**. If a 6000-step
jump protocol (`hold 0–3000 @ 600 K`, `quench 3000–6000 @ 300 K`) is interrupted
at step 4000, restarting with `restart = yes` skips the completed hold and the
already-done part of the quench, and runs only the remaining 2000 steps at
`ref_t` — exactly the right temperature for where it left off. As always,
`md_steps` is the **total**, and `output_dir` / `outname` / the protocol keys must
match between stages.

## Key takeaways

- Annealing is just a **temperature schedule**; equilibrium is its one-stage
  special case. Same runner, selected by `anneal`.
- **`ref_t` is the low temperature** — reused directly, no `t_high`/`t_low`
  bookkeeping.
- **`jump` for kinetics, `linear` for yield.**
- **`md_steps` is the total**; the hold (+ ramp) must fit inside it, and the
  remainder runs at `ref_t`.
- **Match `anneal_steps` to your `tau_t`** — the hold must be many thermal
  relaxation times, and long enough to actually unfold.

## Try next

- Sweep `t_high` (e.g. 500/600/700 K) and check at which temperature the hold
  actually unfolds the protein (*Q* → 0).
- Switch the demo to production-like coupling (`tau_t = 0.01`) and lengthen
  `anneal_steps` accordingly to see the realistic relaxation timescale.
- Generate many independent refolding trajectories at once by combining this
  protocol with **multi-copy** runs (Tutorial 4): set `n_copies > 1` to launch a
  batch of quenches in a single GPU job.
