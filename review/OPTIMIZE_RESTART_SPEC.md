# Spec — Restart / resume for the nscale optimizer (`topo.optimize`)

Status: **draft** · Author: (fill in) · Target module: [`topo/optimize/optimize.py`](../topo/optimize/optimize.py)

## 1. Motivation

The nscale optimizer runs a bounded round loop (up to `max_rounds`, default 6),
and **each round runs a full `ntraj`-copy MD** (GPU, production-length) followed
by per-trajectory Q scoring. A 6-round run of a multi-domain protein at
production `md_steps` is hours of GPU time.

Today the optimizer is **not resumable**: every invocation starts fresh at
round 1, level 0. Its own docstring notes this as the sole remaining limitation
([`optimize.py:40`](../topo/optimize/optimize.py#L40)), and the
`strength-optimization-status` memory lists restart/resume as the only open
item. A crash — GPU eviction, walltime kill on a scheduler, node failure —
throws away every completed round, **and** any progress inside the round that
was running. On HPC batch queues (routine preemption / walltime) this is the
difference between a usable tool and one that must finish inside a single
allocation.

**Goal:** an interrupted `topo-optimize` run can be relaunched and continue from
where it stopped — skipping fully completed rounds and *resuming the MD of a
round that was mid-flight* — reproducing the final calibrated model it would
have produced in one shot.

## 2. Two artifacts already carry the resume state

The design reuses state the run **already writes to disk** rather than inventing
a parallel state store. Two granularities, two artifacts:

### 2.1 Round-level state → `optimization.log`

`optimization.log` is already the authoritative narrative of the run
([`optimize.py:450-465`](../topo/optimize/optimize.py#L450)). Per round it emits:

- `## Round N:` header,
- per unit, the ladder tag + nscale **entering** round N —
  `level N` / `fallback` / `frozen` ([`optimize.py:548-554`](../topo/optimize/optimize.py#L548)),
- per unit, the folded-fraction list + `STABLE`/`unstable` — the **decision**
  ([`optimize.py:617-619`](../topo/optimize/optimize.py#L617)),
- terminal lines (`# Converged at round N`, `# Stopped ...`).

That is exactly the evolving state of the loop: `level[]` (the per-unit ladder
index — nscales are a pure function of it via
[`nscale_for`](../topo/optimize/optimize.py#L96)) plus the decision that advances
it. So the log already tells us **which round completed** and **the `level[]`
to resume with**. `frozen` and the `Scorer` native-contact lists rebuild
deterministically from the unchanged inputs.

### 2.2 Within-round MD state → `round_N/traj/traj.chk`

Each round runs one multi-copy `mdrun` with `output_dir = round_N/traj`,
`outname = "traj"` ([`optimize.py:557-573`](../topo/optimize/optimize.py#L557)),
so its OpenMM checkpoint is `round_N/traj/traj.chk`, written every `nstchk`
steps by the `CheckpointReporter` ([`engine.py:250`](../topo/engine.py#L250)).
`mdrun` already supports resuming from it: with `restart = yes` and the
checkpoint present, `setup_simulation` loads positions+velocities, skips
`done_steps`, and **appends** to `traj.dcd`/`traj.log`
([`engine.py:128-148`](../topo/engine.py#L128), `attach_reporters` append mode).
A missing checkpoint under `restart = yes` warns and cleanly fresh-starts that
round ([`engine.py:132-136`](../topo/engine.py#L132)).

The optimizer currently *forces* `restart = "no"` for every round
([`optimize.py:132`](../topo/optimize/optimize.py#L132)); resuming a mid-flight
round is just a matter of flipping that to `yes` for that one round and pointing
`mdrun` at the existing round dir — **no new checkpoint code is needed in
`mdrun`/`engine`; the file is already there.**

## 3. Design

### 3.1 Make the log the resume record

Two required changes so the log survives and is parsed deterministically:

1. **Append, don't truncate.** `run_optimizer` opens `optimization.log` with
   mode `"w"` ([`optimize.py:453`](../topo/optimize/optimize.py#L453)), which
   erases the very history we resume from. On `--resume`, open with `"a"` and
   emit a `# --- resumed from round N ---` banner.
2. **Emit a compact, machine-readable state line** at the end of each round's
   decision, so resume never has to parse free-text human lines:

   ```
   #STATE round=3 completed=yes converged=no level={"domain:D1":2,"domain:D2":0,"interface:D1-D2":1}
   ```

   One `grep`-able line per round, alongside the existing human-readable log
   narrative (kept as-is). The human lines remain for the user; the `#STATE`
   line is the parse anchor. (Alternative: parse the existing `level N`/`STABLE`
   lines — rejected as brittle; a purpose-built marker is a handful of lines and
   removes all format-guessing.)

   Also write a `#CONFIG fingerprint=<sha256>` line once at the top (§3.3).

This keeps the user's model — *the log has the status* — and adds no parallel
JSON file. The MD state lives in the round's `traj.chk`; the round/level state
lives in the log. Resume **combines the two**.

### 3.2 Resume algorithm

`--resume` (alias `--restart`), plumbed
[`parse_args`](../topo/optimize/optimize.py#L673) →
[`run_optimizer`](../topo/optimize/optimize.py#L407) → `_optimize_loop` as
`resume: bool`.

1. **No log / no `#STATE` line in `out_root`** → warn, start fresh (so a
   `--resume` of a never-started run is not an error — convenient for
   "always pass `--resume`" batch scripts).
2. **Parse the log.** Take the last `#STATE ... completed=yes` line → gives
   `R = last_completed_round` and the `level[]` entering round `R+1`. Verify the
   `#CONFIG` fingerprint matches the current inputs (§3.3); mismatch → hard error.
3. **Rebuild** `Scorer`, `unit_class`, `frozen` from the unchanged inputs; load
   `level[]` from the parsed `#STATE`.
4. **Skip** rounds `1..R` — their `round_N/` trajectories and `Q_*.csv` are on
   disk and their effect is captured in `level[]`. No MD, no scoring.
5. **In-progress round `R+1`** (a `## Round R+1:` header with no
   `completed=yes #STATE`):
   - if `round_{R+1}/traj/traj.chk` exists → re-invoke that round's `mdrun` with
     `restart = yes` (same `output_dir`/`outname`); engine resumes the MD from
     `done_steps` and appends to `traj.dcd`. Then split + score + decide as
     normal, writing that round's `#STATE`.
   - if the round dir is absent or has no checkpoint → run `R+1` fresh
     (`restart = no`), exactly as today.
6. **Continue** `R+2 …` normally (`restart = no`, fresh each round).

The only new per-round decision is *"is this the resumed round, and does it have
a checkpoint?"* — everything downstream is the existing loop body.

### 3.3 Run-identity fingerprint

Reusing on-disk rounds/checkpoints is valid only if the protocol is unchanged.
`#CONFIG fingerprint` = sha256 over:

- reference PDB content (not path — the run may move between nodes);
- seed `domain.yaml` content (residues + class; nscales in the seed are
  overwritten each round, but hash the whole file for a simple, conservative v1);
- optimizer controls: `ntraj`, `q_threshold`, `frame_fraction`, `min_contacts`
  (**exclude `max_rounds`** so an unconverged run can be resumed with a larger
  cap — validate `max_rounds` separately: it may only grow);
- the resolved per-round `sim_options` **except** driver-overridden keys
  (`pdb_file`, `domain_def`, `n_copies`, `output_dir`, `outname`, and CLI
  `--device`). Include `md_steps` (it defines trajectory comparability).

Mismatch → refuse with a message naming what differs; never mix incomparable
rounds. Platform note: an OpenMM `.chk` is **not portable** across GPU /
OpenMM-build / precision ([per the XML-restart TODO in §F](TODO.md)); a mid-round
resume must run on a compatible platform. If `traj.chk` fails to load, engine
already warns and fresh-starts that round — correct, just not free.

### 3.4 STRIDE cache across a resume

The STRIDE-cache optimization ([`optimize.py:578`](../topo/optimize/optimize.py#L578))
mutates `sim_options["stride_output_file"]` in memory after round 1; that is lost
on resume. Before the loop, if unset, re-point it at a completed round's cached
`round_1/<pdb_stem>_stride.dat`. Otherwise the first resumed round re-runs STRIDE
— correct but wasteful.

## 4. Success criteria

1. **Skip is real.** A resumed run does **not** launch `mdrun` for rounds
   `1..R` (assert via round-log timestamps or a mocked `run_md`).
2. **Mid-round MD resumes, not restarts.** Kill a run *during* round `R+1`'s MD
   (after ≥1 `nstchk`); on `--resume`, that round's `mdrun` is invoked with
   `restart = yes`, `traj.dcd` is **appended** (not rewritten), and the final
   per-copy frame count equals an uninterrupted round's.
3. **Same decision path.** Interrupt-then-resume yields the same `level[]`
   trajectory and the same `domain_optimized.yaml` as an uninterrupted reference
   run on a fixed-seed / compatible-platform config.
4. **Fingerprint guard.** Resuming with a changed `q_threshold` / `md_steps` /
   PDB errors out; changing only `max_rounds` (upward) is allowed.
5. **Fresh-start `--resume`.** `--resume` with no prior log starts at round 1
   without error; the log is **appended**, never truncated, across resumes.

## 5. Testing notes

- MD is stochastic (random initial velocities). Criterion 3's "same decision
  path" holds because resume **reuses** the interrupted round's trajectory
  (checkpoint-continued) rather than regenerating it; make the reference
  deterministic with a fixed MD seed if the runner supports it, else assert the
  reuse property directly (resumed `level[]` path == log-implied path).
- **Frame-boundary check** (criterion 2): verify the checkpoint-continued DCD
  has no duplicated/dropped frame at the restart boundary vs. the uninterrupted
  DCD — the `done_steps` skip + DCD-append path is the place a boundary
  off-by-one would hide, and Q scoring reads every frame.
- Unit-test the `#STATE`/`#CONFIG` line **emit ↔ parse** round-trip (unit-key ↔
  string encoding `"domain:D1"` / `"interface:D1-D2"`), and fingerprint
  stability (same inputs → same hash; changed control → different hash).
- Fast CI config: small protein, tiny `md_steps`, low `ntraj`, CPU device, kill
  after round K.

## 6. Open questions

- **Log as sole state store vs. a tiny sidecar.** This spec keeps state in the
  log (`#STATE`/`#CONFIG` markers) per the "the log has the status" intent. If
  robustness against manual log edits/rotation becomes a concern, mirror the
  last `#STATE` into a one-line `out_root/opt_state.json` — but that reintroduces
  a second source of truth to keep consistent. Recommend log-only for v1.
- **Partial `#STATE` write.** If killed mid-write of the `#STATE` line, that
  round parses as not-completed → its MD resumes from `traj.chk` (or reruns).
  Acceptable (bounded to one round); no atomic-write machinery needed for a
  single appended line.
- **Clobber guard.** Without `--resume`, should a non-empty `out_root` that
  already contains an `optimization.log` refuse to overwrite unless `--force`?
  Prevents silently destroying a resumable run, but may be too aggressive for
  existing workflows — decide.
- **`max_rounds` growth on resume** — excluded from the fingerprint here; confirm
  the intended semantics (allow raising the cap to extend an unconverged run;
  disallow lowering below `last_completed_round`).
