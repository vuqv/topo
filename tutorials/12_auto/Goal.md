# Goal — Reproduce O'Brien continuous synthesis in topo-style (autonomous)

> **Agent directive.** This is an autonomous task. Work the loop in
> [§4](#4-execution-loop) **until every box in [§3 Definition of Done](#3-definition-of-done-success-criteria)
> is checked or a [§6 stop condition](#6-stop-conditions) is hit.** Do not end your
> turn just because one step succeeded — re-read §3, find the first unmet criterion,
> and continue. Record progress in [§5](#5-milestones) after every milestone.

---

## 1. Objective

1. **Copy** the contents of `tutorials/10_csp_obrien/` into this folder
   (`tutorials/12_auto/`) as the working scaffold.
2. **Reproduce** the O'Brien `continuous_synthesis_v6.py` protocol **in topo-style**
   (an `csp.ini` + a `topo-csp` run, not by editing the legacy script) on the
   **4c5c** inputs in this folder — same per-codon, three-stage elongation kinetics
   and ribosome geometry — and get it to run to a **physically sane**, **quantitatively
   consistent** result vs. the reference run in `continuous_synthesis/output/`. A key
   physical target: the nascent chain must **not collapse onto the truncated ribosome**
   during synthesis, and on tether release must **diffuse out through the exit tunnel
   (along +x), without penetrating the tunnel wall** (see D5b).
3. **On success**, write a report explaining **why the original `tutorials/10_csp_obrien`
   does not work** (see D8).

## 2. Context & inputs (this folder)

All paths below are relative to this directory
(`tutorials/12_auto/`, i.e. `/storage/work/qzv5006/src/topo/tutorials/12_auto/`).

| Input | What it is |
|-------|-----------|
| `4c5c_model_clean.pdb` | all-atom structure of the target protein |
| `protein_cg_model/cg.cntrl` | CG model control (potential `bt`, `casm=0`, domain file) |
| `protein_cg_model/domain_def.dat` | domain definitions + per-domain Go scale factors |
| `continuous_synthesis/` | a **reference** continuous-synthesis run (input + output) to match |
| `continuous_synthesis/input/cont_synth_ecoli.cntrl` | the reference protocol's parameters (kinetics, ribosome, traffic) |
| `continuous_synthesis/input/setup/` | reference CG `.psf/.top/.prm/.cor`, mRNA, codon trans-times |
| `continuous_synthesis/output/` | reference trajectory + `info.log` to validate against |

**Reference source code** (read-only, for semantics — do **not** modify or run as the deliverable):
- `/storage/home/qzv5006/programs/cg_simtk_protein_folding/` — model-building reference
- `…/Continuous_synthesis_protocol/continuous_synthesis_v6.py` — the protocol being reproduced
- `…/Continuous_synthesis_protocol/` ribosome_traffic — upstream-queue delay binary

**Reuse first — do not re-implement.** The protocol is already ported:
- `topo/csp/` (`csp.py`, `kinetics.py`) + console script `topo-csp` — the runner to use.
- `tutorials/10_csp_obrien/` — the demo to copy here; has `csp.ini`, `PLAN.md`, `TASK.md`,
  `STAGES.md`, `OBSERVATIONS.md`. **Read these first**; `TASK.md` lists done vs. remaining,
  and `OBSERVATIONS.md` #1 records the stage-1 blow-up that 12_auto must avoid.
- `tutorials/11_reproduce_csp/` — prior attempt at the real CHARMM 4c5c + 50S system;
  mine it for the CHARMM-ingest path and any `cont_synth_ecoli.cntrl` → `csp.ini` mapping.
- **Ribosome structure for `topo-csp`.** The reference `setup/50S_tRNA_cg_truncated.*`
  is CHARMM-format; a topo-ready CG PDB of an equivalent truncated 50S+tRNA ribosome is
  `topo/translation/structures/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb` (≈ the same
  geometry). Use it for the topo run unless CHARMM ingestion of the exact `setup/` system
  is working.

The known gap (per `10_csp_obrien/TASK.md`) is **CHARMM `.psf/.top/.prm/.cor` ingestion**
so `topo.csp` can run O'Brien's exact system rather than a topo-rebuilt CG model. Closing
that gap (or building an equivalent topo CG model from `protein_cg_model/`) is the crux.

## 3. Definition of Done (success criteria)

The goal is achieved when **all** of the following hold. Each is verifiable — verify it,
don't assume it.

- [ ] **D0 — Scaffold copied.** Tutorial 10's contents are copied into `tutorials/12_auto/`
      and adapted to the 4c5c inputs. **Never overwrite this folder's existing raw inputs**
      (`4c5c_model_clean.pdb`, `protein_cg_model/`, `continuous_synthesis/`).
- [ ] **D1 — Config exists.** `tutorials/12_auto/csp.ini` reproduces the reference
      `cont_synth_ecoli.cntrl` parameters (kinetics `time_stage_1/2`, `scale_factor`,
      `mrna_seq`, `trans_times`, ribosome geometry, `ribosome_traffic`, `initiation_rate`)
      mapped to topo names. Any deliberate deviation is noted in `NOTES.md` with a reason.
- [ ] **D2 — Inputs prepared.** The protein CG model + ribosome inputs `topo-csp` needs
      are present (either via CHARMM ingestion of `continuous_synthesis/input/setup/*`,
      or a topo CG build from `protein_cg_model/` that matches the same beads/domains).
- [ ] **D3 — Run completes.** `topo-csp -f csp.ini` runs end-to-end with exit code 0,
      synthesizing the full nascent chain (1 → full length, here 306 residues) — or a
      documented short demo length if a full run is infeasible (state the length + why).
- [ ] **D4 — Outputs produced.** A trajectory + a per-residue dwell-time log are written
      under `tutorials/12_auto/synth_out/` (mirroring the reference `output/` layout).
- [ ] **D5 — Physically sane.** No stage blow-ups: per-stage potential energy stays
      finite (no PotE ≳ 1e12 spikes — cf. `10_csp_obrien/OBSERVATIONS.md` #1). The chain
      threads the exit tunnel (monotonic-ish +x egress), never falls back into the ribosome.
- [ ] **D5b — Clean ejection through the tunnel.** During synthesis the nascent chain
      must **not collapse onto / into the truncated ribosome**. After the last residue,
      when the C-terminus tether is released (ejection phase), the chain must **diffuse out
      of the exit tunnel along the tunnel axis (+x)** and fully clear the ribosome — it must
      **never penetrate or pass through the tunnel wall** or any ribosome bead. Verify in the
      trajectory: (a) no nascent bead crosses the `tunnel_wall` radius / overlaps ribosome
      beads at any frame (no steric clashes, energy stays finite through ejection); (b) the
      chain's centre of mass moves monotonically away from the PTC down +x until free; (c)
      the released chain ends up outside the tunnel, not lodged in or beside the ribosome.
      Record the exit trajectory (e.g. min nascent–ribosome distance over time, CoM-x vs.
      frame) in `NOTES.md`.
- [ ] **D6 — Quantitatively consistent with reference.** Compare against
      `continuous_synthesis/output/`: (a) total synthesized length matches; (b) per-codon
      dwell-time distribution / total synthesis time agree within a stated tolerance
      (e.g. within ~2× on summed in-vivo time, given stochastic FPT sampling); (c) final
      nascent-chain geometry is comparable (e.g. radius of gyration in range). Record the
      numbers in `NOTES.md`.
- [ ] **D7 — Reproducible & documented.** `tutorials/12_auto/README.md` gives the exact
      commands to reproduce the run; `NOTES.md` logs decisions, deviations, and the
      validation table. The tutorial is added to `tutorials/README.md`.
- [ ] **D8 — Post-mortem on Tutorial 10.** *On success of D0–D7*, write
      `tutorials/12_auto/WHY_10_FAILS.md` explaining **why the original `tutorials/10_csp_obrien`
      does not work**: contrast what 12_auto did to succeed against 10's setup, name the concrete
      root cause(s) (e.g. CG model/strength mismatch, missing CHARMM ingest, seed-placement
      blow-up in `OBSERVATIONS.md` #1, config-mapping error), and give the minimal fix that
      would make Tutorial 10 work. Back claims with evidence (logs, energies, diffs), not guesses.

## 4. Execution loop

Repeat until §3 is fully checked or §6 fires:

1. **Orient.** Re-read §3 and §5; pick the first unmet criterion.
2. **Plan the step.** Smallest action that advances it. For the hard parts (CHARMM
   ingest, CG build), check `tutorials/11_reproduce_csp/` and `topo/csp/` for existing code.
3. **Act.** Copy/adapt the scaffold, write/adjust `csp.ini`, prepare inputs, or run `topo-csp`.
4. **Verify.** Run the command; inspect logs/energies/trajectory. Compare to the
   reference output. Treat a finished run as a hypothesis to check, not a success.
5. **Record.** Tick the milestone in §5; append findings/decisions to `NOTES.md`.
6. **Iterate.** If it failed, diagnose (read the traceback / energy log), fix, retry.
   Prefer a short fast run (small `L_max`, `max_steps_per_stage`) to debug the pipeline,
   then scale up for the validation run. Do not silently lower success criteria.

**Debug-then-scale:** get an end-to-end short run green (D1–D5 on a demo length) before
attempting the full-length validation run (D3/D6).

> **Fast-test `scale_factor`.** `steps = dwell_s · 1e9 / scale_factor / dt`, so a *larger*
> `scale_factor` means *fewer* MD steps per stage = quicker runs. For the debug/short runs
> use **`scale_factor = 4331293 * 50 = 216564650`** (50× the production value `4331293`).
> Restore the original `4331293` for the D6 validation run — note which value produced
> each result in `NOTES.md`.

## 5. Milestones (update as you go)

- [ ] M0 — Tutorial 10 contents copied here (D0); read `{TASK,PLAN,STAGES,OBSERVATIONS}.md`
      + `topo/csp/README.md`; map `cont_synth_ecoli.cntrl` fields → `topo.csp` config names.
- [ ] M1 — Inputs prepared (D2): CG protein + ribosome usable by `topo-csp`.
- [ ] M2 — `csp.ini` drafted (D1).
- [ ] M3 — Short demo run green: completes + sane energies (D3/D5 at small `L_max`).
- [ ] M4 — Full / target-length run completes with outputs (D3/D4); ejection phase clears
      the chain through the tunnel without wall penetration or collapse (D5b).
- [ ] M5 — Validation vs. reference done; numbers in `NOTES.md` (D6).
- [ ] M6 — `README.md` + `tutorials/README.md` updated (D7).
- [ ] M7 — `WHY_10_FAILS.md` written (D8), only after M0–M6 pass.

## 6. Stop conditions (pause and ask the user)

Stop the loop and report concisely if any occurs:

- A required input is missing/corrupt and cannot be derived from this folder or the
  reference dirs.
- An external dependency is unavailable (e.g. the `ribosome_traffic` binary won't run, a
  GPU/OpenMM platform is missing) **and** no documented fallback in `topo.csp` suffices.
- The same run fails the same way after **3** distinct fix attempts — report the
  traceback, what you tried, and your best hypothesis.
- A validation criterion (D5/D6) cannot be met and the discrepancy looks like a genuine
  science question (model mismatch), not a bug — surface it with the numbers.
- Any step would overwrite reference data (`continuous_synthesis/output/`, `setup/`) or
  this folder's source inputs. **Never overwrite inputs — write only under `synth_out/`.**

When stopping, state: which milestone you reached, what's blocking, and the one decision
or input you need.
