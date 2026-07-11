# Proposal: Resume for CSP Continuous Synthesis

*Date: 2026-07-11 · Scope: `topo/csp/protocol.py` (driver), `topo/csp/kinetics.py` (schedule), new `topo/csp/resume.py` · Status: draft for review*

## Summary

A production `topo-csp` run is hours to days of wall time and today survives no
interruption. This proposal adds resume with three deliberately small moving parts, plus
one separable refactor:

1. **Precompute the schedule up front (§3.1).** Draw every residue's stage step counts
   from the seeded generator in one pre-loop pass, persist the table (it already exists
   as `dwell_times.dat`), and have the main loop *read* step counts rather than draw
   them. This removes the hard part of resume — there is no live RNG cursor to
   serialize — and yields an exact **total-step count reported before any MD runs**, so
   a scheduler wall-clock request stops being guesswork.
2. **Track progress in a human-readable log (§3.2).** A `DONE`/`RUNNING` `progress.log`
   records how far the run got; resume drops the one `RUNNING` unit and continues from
   the last `DONE`.
3. **Reload the seed conformation from disk (§3.3).** The only per-residue state besides
   the schedule is the previous residue's final coordinates, already written every stage.

Separately, an **output-layout consolidation (§3.5)** collapses the three per-stage
subdirectories into one per-residue directory — a file-count/cleanliness win that is
independent of resume and ships as its own commit.

Terminology note: the single `random_seed` seeds one generator; we pre-draw all
per-residue *step counts* from it (not one seed per residue).

## 1. Problem and goals

A production run synthesizes `L = L0 .. L_max` residues, each as three MD sub-stages,
then post-synthesis ejection/dissociation. At production step counts
(~10⁵–10⁶ steps/stage) a full protein is **hours to days** of wall time. Two gaps:

- **No resume.** Any interruption — scheduler wall-clock limit, node failure, `Ctrl-C` —
  loses the whole run: re-invoking `topo-csp` restarts from `L0` and overwrites the
  output directory. This is the last open item on the CSP feature list
  (`review/TODO.md` §B/§D, memory `csp-obrien-status`).
- **No up-front cost.** The user cannot know the total step count / expected runtime
  before launching, making scheduler wall-clock requests guesswork.

**Goals.** Re-invoking `topo-csp` on an interrupted output directory (a) continues from
the last completed residue with a kinetic schedule **identical** to the uninterrupted
run, and (b) a launched run reports its **total step count** (and a nominal runtime)
before the first MD step.

## 2. Key insight: what state crosses a residue boundary

Reading the driver (`run_continuous_synthesis`, `protocol.py:71`), the main loop
(`for L in range(L0, L_max+1)`) carries **exactly two** pieces of mutable state between
iterations:

| State | Source | Recovery on resume |
|---|---|---|
| `prev_final` — the `(L-1, 3)` nm nascent coords seeding residue `L` | the previous residue's final structure, already on disk | **reload from disk** (§3.3) |
| the per-residue step schedule `(s1,s2,s3)` | drawn from `rng` **once per residue** in `stage_steps` → `stage_dwell_times` (`kinetics.py:485`) | **read from the persisted schedule** (§3.1) — no RNG at run time |

Everything else the loop consumes is derived **deterministically and cheaply** (relative
to MD) from the config at startup and is simply recomputed on resume:

- ribosome load, P/A anchors (`protocol.py:148–185`);
- `precompute_contacts` (contacts + STRIDE, cached via `stride_output_file`) (`:188`);
- `build_codon_time_lists` — the `intrinsic`/`real`/`codons` kinetic lists (`:204`).

**One exception is reloaded, not recomputed: the PTC restraint targets.**
`optimal_ptc_targets` (`core.py:129`) and the tunnel-wall plane derived from them
(`protocol.py:148–185`) are deterministic and cheap (a one-time SLSQP multistart over the
ribosome, no RNG, no per-residue cost), so recompute *would* work. But they are
persisted at first launch and **re-read on resume** (§3.1) for the same reason the
schedule is: a scipy/numpy upgrade between launch and resume could shift `a_target` /
`p_target` slightly, and the resumed residues would then be built against different
restraint targets than the residues already on disk. Reloading pins the geometry
identical across the interruption.

**Consequence:** no heavyweight checkpoint of the simulation state is needed. The natural
resume unit is **the residue**; the only cross-residue state is `prev_final` (on disk)
and the schedule (persisted once). What remains to record is just a pointer at how far
the run got — the `progress.log` of §3.2.

### 2.1 Reproducibility ceiling (why disk reload is enough)

The per-stage MD is already stochastic and *not* bit-reproducible across a process
boundary: the Langevin thermostat and `setVelocitiesToTemperature` draw from OpenMM's
own (unseeded) platform RNG, and stage 1 re-seeds + re-minimizes + redraws velocities
every residue regardless. So the achievable and *intended* guarantee is:

> **the kinetic schedule (per-residue dwell times → step counts) is identical to the
> uninterrupted run**, and the resumed conformation continues from the last completed
> residue.

The persisted schedule delivers the first half exactly (same file, re-read). The final
structure reload carries 0.001 nm (PDB) precision for `prev_final`, far below the model's
existing thermal noise — adequate and consistent with the current design. We do **not**
need the non-portable OpenMM `.chk` for this (TODO §D concern moot here).

## 3. Design

### 3.1 Precompute the schedule

Before the main loop, draw the whole schedule and persist it. **`dwell_times.dat` already
*is* this table** — today it is written *inside* the loop (`protocol.py:262–268`), one row
per residue, with `t1/t2/t3` (dwell seconds) and `steps1/2/3` (integration steps). The
change is to compute it in a pre-loop pass and make the file the source of truth the loop
reads back:

```python
# --- pre-loop: draw the full schedule once, persist it, report cost --------
schedule = []                       # schedule[L-L0] = (s1,s2,s3),(t1,t2,t3)
rng = random.Random(params.random_seed)
for L in range(L0, L_max + 1):
    schedule.append(kinetics.stage_steps(L, intrinsic, real, ..., rng=rng))
# persist the immutable plan: PTC targets (§2) as a header block + the schedule rows
write_schedule(dwell_log, schedule,
               header=ptc_header(a_target, p_target, wall_plane))  # a/p targets + wall x
total_steps = sum(s1 + s2 + s3 for (s1, s2, s3), _ in schedule)
print(f"[schedule] {L_max - L0 + 1} residues, {total_steps:,} planned MD steps"
      f"{est_walltime(total_steps, params)}")        # nominal; see caveat
```

The `a_target` / `p_target` (and the derived tunnel-wall plane) go in the **header** of
`dwell_times.dat` — full float precision, machine-readable comment lines — so the one
immutable-plan file carries both the schedule *and* the restraint geometry. On resume
`read_schedule` returns both, and neither `stage_steps` (RNG) nor `optimal_ptc_targets`
(SLSQP) runs again. (Not the human `run.log`, which prints the same numbers for the reader
but is not the machine source — see §3.2.)

The loop then reads step counts instead of drawing them, and records progress (§3.2):

```python
for L in range(resume_L, L_max + 1):
    (s1, s2, s3), (t1, t2, t3) = schedule[L - L0]    # read, no RNG
    append_progress(out_path, f"L_{L:03d}", "RUNNING")
    ... run 3 stages with s1/s2/s3 ...               # prev_final = f3
    append_progress(out_path, f"L_{L:03d}", "DONE")
```

This is a **behavior-preserving refactor of an uninterrupted run**: `rng` is consumed
**only** in `stage_steps` → `stage_dwell_times`, once per `L` in ascending order
(verified), so pre-drawing `L0..L_max` sequentially produces the identical stream the
current lazy loop does.

**Cost-report caveat** (state it in the console output too): `total_steps` is the
**planned** count at the base timestep and is **exact** as a step total. Wall-time is
only a **nominal estimate**, because (a) the dt-halving stability guard (`core.py:1162`)
re-runs a diverged stage at `dt/2 × 2` steps, so *actual* integration steps ≥ planned,
and (b) per-step cost grows as the nascent chain lengthens. Report total steps as exact;
label any time figure as an estimate with a calibration factor.

### 3.2 Track progress (`DONE` / `RUNNING`)

With the schedule persisted, "how far did we get" is the only mutable state left to
record. Rather than an opaque JSON blob, use a **human-readable append-only progress log**
at the output root — a single-user research run wants to `tail`/`cat` its progress, and an
explicit `RUNNING` marker *names* the partial unit to drop rather than inferring it:

```
<out_root>/progress.log
```

```
# csp progress log — schema 1
L_001 RUNNING
L_001 DONE
L_002 RUNNING
L_002 DONE
...
L_041 RUNNING
L_041 DONE
L_042 RUNNING        <- crash here: L_042 is the partial unit to drop on resume
```

**Protocol (crash-safe by ordering):**

- Append `L_XXX RUNNING` when a residue starts; append `L_XXX DONE` when **all three
  stages** finish. Post-synthesis phases use the same scheme
  (`ejection RUNNING`/`ejection DONE`, then `dissociation …`).
- The **`DONE` line is the commit point.** A single short append is effectively atomic
  on POSIX, and the ordering makes every crash recoverable: die after writing frames but
  before `DONE` → that unit is `RUNNING` → dropped and redone; die after `DONE` → resume
  at the next unit. No partial trajectory is ever mistaken for complete.
- **On resume:** take the *last* status per unit, drop every `RUNNING` entry **and its
  on-disk directory** (`L_XXX/`), and continue from `max(DONE) + 1`. There is at most one
  `RUNNING` to drop (the unit in flight when the run died); dropping one residue's
  directory is cheap and clean. The dropped residue is re-run from its *persisted*
  schedule row, so the redo matches the original plan.

**Config sanity (no fingerprint needed).** Because the schedule is **loaded, not
redrawn**, resume correctness does not depend on the config being bit-identical to the
launch config — a changed `random_seed` or a retuned kinetic knob cannot poison a resume,
since the RNG is never consulted again (the schedule file *is* the materialized draws).
So no SHA1 fingerprint guard is carried. The one hard constraint is structural, not
config-hash based: the persisted schedule covers exactly the original `L0..L_max`, so a
resume simply runs whatever residues that table still contains. **Extending a run
(a larger `L_max`) is a fresh run** — it would need draws past the materialized schedule.
If a launch config disagrees with the persisted schedule's `L_max`, that is caught
directly by reading the table, not by a config hash.

**Two files, clean separation.** `dwell_times.dat` is the **immutable plan** (all
residues, written once — §3.1); `progress.log` is the **mutable status**. Keeping them
apart avoids any rewrite hazard on the plan, and there is **no serialized RNG state**
anywhere — the schedule file *is* the materialized RNG output.

### 3.3 Startup / resume logic

Add a config knob (`RunParams` + INI `[OPTIONS]`):

```
resume = auto      # auto (default) | yes | no
```

- `auto` — resume iff `progress.log` exists (and the plan file it points at is readable);
  else fresh.
- `yes` — resume; error if no valid `progress.log`.
- `no` — always fresh (current behavior); refuse to clobber a non-empty out dir unless
  `--force`/`overwrite` (optional hardening, see §7).

At the top of `run_continuous_synthesis`, after the deterministic precompute (§2):

```python
resume_L = L0
if resume_enabled and progress_exists(out_path):
    prog = read_progress(out_path)                   # parse progress.log
    schedule, a_target, p_target, wall_plane = read_schedule(dwell_log)  # plan + PTC geom, re-read (not recomputed)
    verify_completed_units(out_path, prog)           # every prior length present
    drop_running_units(out_path, prog)               # rm the ≤1 RUNNING dir(s)
    resume_L = prog.last_done_residue + 1
    prev_final = load_final_pdb(final_path(out_path, prog.last_done_residue))
    print(f"[resume] {prog.last_done_residue} lengths present; continuing from L={resume_L}.")
else:
    a_target, p_target = optimal_ptc_targets(ribo)   # SLSQP, fresh start only
    wall_plane = derive_wall_plane(a_target)         # protocol.py:148–185
    schedule = precompute_and_write_schedule(..., a_target, p_target, wall_plane)  # §3.1 header + rows
    write_progress_header(out_path)
```

On a fresh start `optimal_ptc_targets` runs once and its result is written into the
schedule-file header; on resume that header is re-read and the SLSQP solve is skipped —
so the restraint geometry is guaranteed identical to the residues already on disk (§2).

(`final_path` resolves to `L_<L>/traj_final.pdb` under the consolidated layout of
§3.5, or `L_<L>/stage_3/traj_final.pdb` under the current per-stage layout — the resume
logic is identical either way.)

**Presence check before resuming (`verify_completed_units`).** `progress.log` records
intent; the disk records reality, and the two can diverge (a `DONE` residue's directory
deleted, moved, or lost to a scratch-purge). Because resume only reloads the *last* `DONE`
residue as `prev_final`, a hole anywhere earlier would go unnoticed and leave a permanently
broken assembled trajectory/movie. So before continuing, verify **every** length
`L0..last_done` is present:

```python
def verify_completed_units(out_path, prog):
    for L in range(L0, prog.last_done_residue + 1):
        fp = final_path(out_path, L)
        if not fp.is_file():
            raise SystemExit(f"[resume] L_{L:03d} is marked DONE but its final "
                             f"structure {fp} is missing — the output tree is "
                             f"incomplete. Refusing to resume (re-run fresh, or "
                             f"restore the missing length).")
```

A simple **presence** check is all a correct resume needs: every previous length's final
structure exists, or the run aborts with the offending length named. We deliberately do
**not** re-parse each final to count residues or otherwise validate its contents — that
byte-level consistency probe is more machinery than the guarantee requires. This runs
**whenever a resume is attempted**: it is the concrete meaning of `resume = yes`. For
`resume = auto` the same check applies once auto has decided to resume; the abort is
deliberate — silently falling back to a fresh start would overwrite the surviving good
data.

Two current-code fixes this requires:

1. **`dwell_times.dat` is opened `"w"`** (`protocol.py:230`) — it truncates the table on
   any re-run. In the adopted design it is written **once** in the pre-loop pass (fresh
   start) and **re-read, not rewritten** on resume. *(The current `"w"`-in-loop pattern
   is a latent data-loss bug on any re-invocation today.)*
2. **Partial residue `resume_L`:** the `RUNNING` residue's directory is removed wholesale
   by `drop_running_units`, then residue `resume_L` is **re-run from scratch** off its
   persisted `(s1,s2,s3)` row. (Even without the explicit drop the re-run is
   self-healing: every stage does `out_dir.mkdir(exist_ok=True)` and each stage's
   `setup_simulation`/`attach_reporters` truncates its own files. The drop just keeps the
   tree clean.) Worst-case redone work is one residue (≤ 3 stages).

### 3.4 Post-synthesis phases

`ejection` and `dissociation` (`protocol.py:322–344`) are single `run_length` calls with
**fixed** step counts (`params.ejection_steps` / `dissociation_steps`) — already
deterministic, no RNG. They appear in `progress.log` as their own units
(`ejection RUNNING`/`ejection DONE`, then `dissociation …`). On resume, if all residues
are `DONE`, skip elongation, reload the seed from the last completed unit's final
(`L_max` final, or `ejection/` if ejection already finished), and run only the remaining
phases. Their fixed step counts are added to the up-front `total_steps` report for
completeness.

### 3.5 Output-layout consolidation (separable refactor)

*Independent of the resume mechanism above — it changes only where bytes land, not what
state is tracked. Recommended as its own commit: it makes the `RUNNING`-drop of §3.2
cleaner (one directory per residue instead of three) and harmonizes CSP with the flat
`cylinder` loop, which already writes one directory per length.*

Collapse the three per-stage subdirectories into a single per-residue directory:

```
   before                          after
   L_042/stage_1/{traj.psf,        L_042/traj.psf            (one, shared)
                 native_1_42.pdb,        native_1_42.pdb     (one, shared)
                 traj.dcd,               traj_s1.dcd
                 traj_final.pdb,         traj_s2.dcd
                 traj_runinfo.log}       traj_s3.dcd
   L_042/stage_2/{... identical psf,     traj_final.pdb       (seed for L_043)
                  identical native,      traj_runinfo.log     (one)
                  traj.dcd, ...}
   L_042/stage_3/{...}
```

**Only one final per residue is kept — `traj_final.pdb` (the stage-3 final).** The
stage-1/2 `_final.pdb` files are dropped: nothing reads them (fresh runs chain stages in
memory via the returned `f1`/`f2` coords; resume and the next residue read only this final
as `prev_final`), and their last conformation already lives as the **last frame of
`traj_s1.dcd` / `traj_s2.dcd`**. Nothing on resume needs `f1`/`f2` from disk either
(resume is per-residue — §4 — and re-runs the whole residue from its schedule row).
`traj_final.pdb` is the one functional artifact — the seed for the next residue. Because
only stage 3 writes it, there is no `s3` in the name: one residue, one final.

**Why sharing `traj.psf` and `native_1_L.pdb` across the three stages is correct:** both
depend only on `L`. The nascent topology (atoms + bonds + contacts sliced to `L×L`) is
identical across a residue's three stages — the A/P-site differences are *forces*
(tether/restraint), not PSF atoms — and `write_subset_structure(full_pdb, L, …)` is a
pure function of `L`. Today they are regenerated **3× per residue with byte-identical
content** (`core.py:1023`, `core.py:1061`); sharing computes them once. **DCDs stay split
per stage** (`traj_s{1,2,3}.dcd`) to preserve the stage/kinetic-dwell boundaries, and one
shared PSF legitimately serves all three since the topology is identical.

**Benefits:** ~3× fewer files (≈5–6k → ≈2k for a 300-mer — real metadata/inode relief on
Lustre/GPFS); two redundant `write_subset_structure` + `dumpTopology` calls saved per
residue; one-directory-per-residue matches `cylinder`.

**Blast radius (update in lockstep — this is why it is its own commit):**

- `run_length` (`core.py`) output naming — decouple per-residue shared files (`traj.psf`,
  `native_1_L.pdb`) from per-stage files via an outname suffix (`cfg.outname = "traj_s{n}"`);
  guard so the shared files are written only once.
- **`run_length` return path (`core.py:1246–1249`)** — small behavioral change, not just
  renaming. Today step 6 reads `_final.pdb` back with `PDBFile(...)` to build the returned
  `(L,3)` coord array; refactor it to take the coords **directly from the context state**
  (`ctx.simulation.context.getState(getPositions=True)...[:L]`), removing the PDB
  round-trip. This decouples the return value from whether the file is written, so the
  `_final.pdb` write (currently unconditional in `_finalize_nascent` / `finalize_simulation`,
  `core.py:800–803`) can be gated on a new `persist_final` flag: **True** for stage 3 and
  the post-synthesis phases (their finals seed the next phase and are the resume-reload
  targets, §3.4), **False** for stages 1/2. Stages 1/2 then write no `_final.pdb` while
  still returning their `f1`/`f2` coords in memory as before.
- `protocol.py:288/298/311` — the `out_subdir=f"{ldir}/stage_N"` arguments (pass
  `persist_final=False` for the stage-1/2 calls, `True` for stage 3).
- **`movie.py:363–370`** — the one real downstream consumer; it globs
  `L_*/stage_N/{outname}.psf` and picks the traj. Must move to `L_*/traj_s{n}.dcd` + the
  shared `L_*/traj.psf`, or movie generation silently yields empty segments.
- The resume seed path (§3.3) — `final_path` resolves to `L_<L>/traj_final.pdb`.
- A `grep` of `tutorials/` and `docs/` for `stage_` path references before committing.
- `traj_runinfo.log` becomes one-per-residue (fold the three stage sections); per-stage
  wall-time/PE already also prints to stdout (the concise summary line), so no provenance
  is lost.

## 4. Resume granularity: per-residue (fixed decision)

**The resume unit is the residue, full stop.** Per-stage resume is explicitly **out of
scope and will not be implemented** — not "deferred," not "v2." The residue is the right
and only granularity because:

- it matches the natural state boundary — `prev_final` updates once per residue, and the
  schedule is indexed per residue;
- redone work on a crash is bounded to a single residue (its ≤ 3 stages re-run from the
  persisted schedule row), which is negligible against a run measured in hours to days;
- it is minimal — no per-stage index, no intra-residue bookkeeping, no reloading `f1`/`f2`
  from DCD last-frames. Stages 1 and 2 remain purely in-memory hand-offs, exactly as the
  current uninterrupted loop already treats them.

Committing to per-residue keeps the entire design centered on one recovery target
(`traj_final.pdb`) and one status per unit (`L_XXX DONE`). Nothing in `progress.log`,
`resume.py`, or the output layout carries a stage-level notion of resumability.

## 5. Implementation plan

The **resume mechanism** (steps 1–4) is localized to the driver and a new helper module;
`run_length`/`core.py` are untouched. The **layout consolidation** (step 5, §3.5) is a
separable refactor that *does* touch `core.py`/`movie.py` and ships as its own commit.

1. **`topo/csp/resume.py`** (new, ~80 lines) → `write_schedule`/`read_schedule` (the
   `dwell_times.dat` reader/writer — the header carries `a_target`/`p_target`/`wall_plane`,
   the rows carry the per-residue schedule), `write_progress_header`/`append_progress`/
   `read_progress`, `verify_completed_units`, `drop_running_units`, `load_final_pdb`,
   `est_walltime`. *Verify:* unit tests for schedule round-trip (write→read reproduces
   `(s1,s2,s3)` **and** the PTC targets to full precision), progress parse
   (last-status-per-unit), presence check (missing final → abort), and `RUNNING`-drop.
2. **`run_continuous_synthesis`** → (a) pre-loop schedule precompute + persist + cost
   report (§3.1); (b) loop reads `schedule[L-L0]` instead of drawing RNG; (c) resume
   branch (§3.3); (d) `RUNNING`/`DONE` append around each residue + post-synthesis phase.
   *Verify:* integration test in §6.
3. **`RunParams` + `read_csp_config`** → add `resume` key (`auto`/`yes`/`no`).
   *Verify:* INI parse test.
4. **CLI** (`csp()`) → optional `--no-resume` / `--fresh` flag mirroring the config knob
   (convenience; config is the source of truth).
5. **(Separate commit) Output-layout consolidation (§3.5)** → decouple per-residue shared
   files from per-stage files in `run_length`; update `protocol.py` `out_subdir` args;
   **update `movie.py:363–370` in lockstep**; point `final_path` at `traj_final.pdb`.
   *Verify:* a movie regenerated from a consolidated run plays back identically to one
   from the old layout; `grep tutorials/ docs/` for `stage_` clean.
6. **Docs** → `topo/csp/README.md` + Tutorial 12 note; `review/TODO.md` check off §B/§D.

## 6. Test plan (the success criterion)

Deterministic end-to-end resume equivalence, using the fast test clamp
(`max_steps_per_stage` small, as Tutorial 12_auto already does):

1. **Golden run:** synthesize `L0..L_max` (e.g. the 4c5c reproduction, small clamp),
   uninterrupted → record `dwell_times.dat` (the schedule) and each residue's final.
2. **Interrupted run:** same config, kill after residue `k` (mid-run, so `L_{k+1}` is
   left `RUNNING`). Re-invoke. *Assert:*
   - the resumed run's `dwell_times.dat` is **byte-identical** to the golden run — the
     schedule is re-read, not redrawn, so this holds by construction (the test pins it);
   - the partial `L_{k+1}/` directory was dropped and cleanly rebuilt;
   - the resumed run reaches `L_max` and produces all stage outputs;
   - `prev_final` continuity: residue `k+1`'s seed equals golden residue `k`'s final.
3. **Presence guard:** delete a mid-run length's final structure (e.g. `L_{k//2}`) whose
   `progress.log` still says `DONE`, then resume → asserts an actionable abort naming that
   length, **not** a silent continue-from-last-DONE that would leave the hole.
4. **Idempotent completed run:** resume a finished run (all units `DONE`) → no-op, exits
   cleanly.
5. **Schedule ≡ legacy:** the pre-loop schedule for a given seed reproduces the step
   counts the current in-loop code produces (guards the refactor in §3.1).

Add as `tests/csp/test_resume.py`. Current CSP has *no covering tests* (codegraph
blast-radius), so this also lays the first driver-level integration test.

## 7. Risks / open questions

- **No config guard on resume** — deliberately, resume does not fingerprint/validate the
  launch config against the original (the schedule is loaded, not redrawn, so a changed
  seed/knob cannot corrupt the tail). The cost is that a user who resumes with a genuinely
  different *force* config gets a run whose early lengths used the old config and later
  ones the new — a self-inflicted footgun, not a correctness bug the tool must police. If
  this bites in practice, a lightweight config-mismatch **warning** (not an abort) can be
  added later.
- **Wall-time estimate is nominal, not exact** — see §3.1. Report total steps as exact;
  label any time figure as an estimate (dt-halving retries and growing per-step cost both
  push actual time above the naive `steps × const`).
- **Extending a run** — a larger `L_max` on resume is intentionally rejected (§3.2); the
  schedule is fixed at first launch. If "extend an existing run" becomes a real need, it
  is a separate feature (regenerate-and-append the schedule tail).
- **Non-elongation determinism** — we explicitly do *not* claim bit-reproducible MD
  micro-trajectories (thermostat RNG is unseeded); the guarantee is the kinetic schedule
  + conformational continuity (§2.1). Bit-exact MD would be a separate, larger change
  (seed OpenMM platform RNG + save/restore `.chk` velocities).
- **Concurrent/duplicate invocation** — two `topo-csp` on the same out dir would race
  `progress.log`. Out of scope; a lockfile is a cheap optional add.
- **Layout consolidation is orthogonal** — §3.5 can ship before, after, or independently
  of the resume mechanism. Resume works with the current per-stage layout too (drop three
  stage dirs instead of one); the consolidation just makes the drop and the tree cleaner.

## Appendix — alternative considered (not adopted)

**Per-residue RNG-state checkpoint.** Keep drawing the schedule lazily in the loop and
serialize `rng.getstate()` into a checkpoint each residue, restoring it on resume.
Functionally equivalent for resume, but strictly more complex (RNG (de)serialization + a
test proving the restore reproduces the draw sequence) and gives **no** up-front cost
report. The precomputed-schedule approach dominates it: same resume guarantee, less
machinery, and the total-steps estimate for free.
