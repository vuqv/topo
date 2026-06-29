# Observations / known artifacts (CSP movie)

Running notes from inspecting `topo-csp-movie` output. These are things to
investigate later, not settled conclusions.

## 1. Stage-1 blow-up at some residue additions (real, self-recovering)

**Seen (full-length GPU run, `synth_out/movie.dcd`, 907 frames):** around **frame
787** the protein looks collapsed; the next few frames "explode" (beads fly off); by
**frame 791** the chain is back in the tunnel, extended.

**Frame → residue/stage mapping** (`topo.csp.make_movie.find_stage_segments`):

| frame | segment | meaning |
|-------|---------|---------|
| 785–787 | L=93 **stage 3** | tRNA-binding wait; chain relaxed/collapsed at L=93 |
| 788–790 | L=94 **stage 1** | **new residue 94 delivered at the A-site** ← the "explosion" |
| 791–793 | L=94 stage 2 | recovered; chain extended in the tunnel again |

It sits exactly on the **L=93 → L=94 boundary** — the per-residue rebuild/seed step.

**This is a real numerical blow-up, not a rendering artifact.** The energy log
confirms it (`synth_out/L_094/stage_1/traj.log`):

```
L_093/stage_3 last frame:  PotE = -6.7 kJ/mol      Temp = 227 K     (fine)
L_094/stage_1 all frames:  PotE = 3.4e13 kJ/mol    Temp = 3.0e13 K  (EXPLODED)
L_094/stage_2 last frame:  PotE = 617 kJ/mol       Temp = 166 K     (recovered)
```

**It self-recovers.** The *next* stage's `minimizeEnergy` pulls the bad bead back, so
the chain returns to normal and the final structures are fine (`L_106/stage_3`:
PotE ≈ 306 kJ/mol, Temp ≈ 221 K).

**Where it happens.** Scanning all stage logs for `PotE > 1e6`, **5 of 306 stages**
blew up — and almost all are **stage 1** (new-residue delivery):

```
L_014/stage_1, L_014/stage_2, L_094/stage_1, L_095/stage_1, L_101/stage_1
```

**Likely cause.** Stage 1 places the **new bead at the A-anchor** (`buffer` nm beyond
it) while its bond partner (residue L-1) is wherever stage 3 left it. For most
residues `minimizeEnergy` relaxes the stretched bond / any excluded-volume overlap.
But when stage 3 happens to leave a bead **near the A-anchor**, the new bead is placed
nearly coincident with it → the excluded-volume `(sigma/r)^12` term is astronomically
large, the minimizer's line search can't escape the near-singularity, and the
subsequent dynamics (after `setVelocitiesToTemperature`) explode. The clamped demo
step counts (≤667/stage) don't cause it but do make it visible (no time to wash out
before the stage ends).

**Severity.** Cosmetic for the **final per-residue structures** (they recover via the
next stage's minimization) and for the **L=106 product**. But the **trajectory frames**
of the affected stage-1's are garbage, and there is a latent risk that an explosion
corrupts the carried-over coordinates of residues `1..L-1` before recovery — to be
confirmed (compare `L_094/stage_1` final residues 1..93 vs `L_093/stage_3` final).

**To investigate / fix:**
- Robustify stage-1 seeding so the new bead can't land on an existing bead:
  check the A-anchor neighbourhood and offset the placement, or grow the `buffer`
  adaptively until the min pair distance clears `sigma`.
- Harden minimization at residue addition: more iterations / tighter tolerance, or a
  **soft-core / capped** nonbonded for the minimization only, then switch to the real
  potential (standard "push apart overlaps first" trick).
- Add a post-minimization **energy/force sanity check** in `run_length` (the engine
  has `check_large_forces`; it is currently off for the seeded build) that re-minimizes
  or re-seeds if PotE is absurd, instead of integrating an exploded state.
- For the movie only: skip frames whose PotE exceeds a threshold (drop blown-up
  stage-1 frames) so the movie stays continuous.

**Status:** open. Add to `TASK.md` / `PLAN.md` "Remaining" — this is the most concrete
robustness bug found so far in the CSP port.
