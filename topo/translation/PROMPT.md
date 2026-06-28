# Implementation prompt — co-translational synthesis (`topo.translation`)

Paste this (or point me at it) to start a fresh session. It is the build plan;
the **full agreed design is in [`DESIGN.md`](DESIGN.md)** and the **structures /
conventions / done work are in [`FILES.md`](FILES.md)** — read both first.

> First instruction for the new session: **read `topo/translation/DESIGN.md` and
> `topo/translation/FILES.md` in full before writing any code.** Everything below
> is already agreed there; do not re-litigate the design.

Branch: `translation`. New module: `topo.translation` (sibling of `topo.mdrun` /
`topo.optimize`), reusing `topo.engine` (build → setup → finalize) where possible.

There are **two build steps: v1 then v2.** Do v1 fully (and test) before v2.

---

## Non-negotiable invariants (from DESIGN.md — respect in both versions)

1. **Particle ordering:** nascent chain = indices `0..L-1`; ribosome (v2) appended
   at `L..N-1`. Required, because the contact table is particle-indexed.
2. **Build-once-subset contacts (§3.5):** compute the native/non-native matrices
   **once** on the *full* native PDB → `R_full`, `eps_full` (`N_full×N_full`);
   each length uses the **top-left `L×L` block**. Never re-run STRIDE per length.
3. **Coordinates ≠ force field:** native PDB sets contact distances + bonded terms
   only; initial coordinates come from the previous step + new-residue placement.
4. **Restrain only the current C-terminus** to the P-anchor (hand-off is automatic
   because each step rebuilds and restrains only residue L).
5. **σ sources:** intra-chain non-native σ = TOPO structure-derived on the full
   structure, fixed across L; ribosome–NC σ = `model_parameters['radii']`
   (= collision diameter), `σ_ij=(σ_i+σ_j)/2`. (Validated against O'Brien — §3.5.)
6. **Placement buffer ≥ excluded-volume contact distance** so the new bead does
   not get a huge `(σ/r)¹²` kick (would crash). A short minimization each step
   also helps.

Test settings (§2.3/§2.4): **`n_steps_per_residue = 1000`** (constant), test
protein **P0CX28** (`tutorials/01_single_domain_quickstart/P0CX28_clean.pdb`, 106
res), **no tunnel-confinement term, no re-equilibration discard.** Pick a small
starting length `L0` (e.g. 5–10).

---

## STEP v1 — elongation loop, mechanics only (NO ribosome forces)

System = **nascent chain only**. The tRNA anchors are used purely as fixed
*coordinates* (for placement + the restraint target); no ribosome particles in
the System yet.

Implement `topo/translation/elongate.py` (module + CLI) that does:

1. **Precompute once:** run TOPO's contact builder on the full P0CX28 PDB →
   `R_full`, `eps_full`. Cache (and the STRIDE output).
2. **Read anchors** from the truncated ribosome
   (`structures/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb`): P-anchor = `PtR`
   resid-76 `R` bead coord; A-anchor = `AtR` resid-76 `R` bead coord. (These are
   just fixed points in v1.)
3. **Loop `L = L0 .. N_full`:**
   a. Build the length-`L` TOPO model on native residues 1..L — reuse the bonded
      terms / per-residue mass-charge-radii TOPO already builds, but **inject the
      `L×L` subset matrices** (`R_full[:L,:L]`, `eps_full[:L,:L]`) into the
      contact `CustomNonbondedForce` instead of recomputing. (You'll likely add a
      code path / option so the model can take precomputed matrices + a residue
      range, rather than calling `build_nonbonded_interaction` per step.)
   b. **Initial coordinates:**
      - `L == L0`: lay residues 1..L0 **extended along the tunnel axis (X)** from
        the P-anchor — residue L0 (C-term) at P-anchor, residue 1 (N-term) toward
        +x — at one CG bond length (0.381 nm) spacing.
      - `L > L0`: residues 1..L-1 from step (L-1)'s final structure; residue L at
        **A-anchor + buffer** (buffer per invariant 6).
   c. Add a **harmonic restraint** on residue L (C-terminus) → P-anchor
      (`CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")`, k ≈ 200
      kcal/mol/Å²). Restrain only residue L.
   d. Langevin integrator at `ref_t`; set Boltzmann velocities; **minimize**
      (relax clashes); **run `n_steps_per_residue` (=1000)**; save the final
      structure (per-length output folder) + trajectory/log.
   e. Seed L+1 with this final structure.

**v1 acceptance test:** run on P0CX28 from `L0` for several lengths; confirm it
runs without crashing, the chain grows N→C by one residue per step, the new
residue migrates A→P under the restraint, and per-length files are written. (No
ribosome physics expected yet — this validates the loop machinery.)

---

## STEP v2 — add the rigid ribosome forces

Add the truncated ribosome as fixed scenery and wire the two cross-interactions.

1. **Append the ribosome** (`*_cg_trunc.pdb`) at indices `L..N-1`, **mass = 0**
   (rigid; not integrated), coordinates as-is (do not shift). Now the anchors are
   real beads in the System.
2. **Contact force** (nascent native/non-native, the `L×L` table): restrict to an
   **interaction group `{nascent}×{nascent}`**; give ribosome beads an
   `addParticle` entry with **dummy `id=0`** (never evaluated). Table stays `L×L`.
3. **New ribosome–NC force:** a separate `CustomNonbondedForce`, pure
   `ε·(σ_ij/r)¹²`, `ε = 0.000132 kcal/mol`, per-particle σ = `model_parameters['radii']`,
   `σ_ij=(σ_i+σ_j)/2`, cutoff 2.0 nm / switch 1.8 nm, interaction group
   **`{nascent}×{ribosome}`**.
4. **Electrostatics:** extend the existing **Yukawa** force over all charged
   particles (nascent + ribosome phosphates −1e and charged residues).
5. Keep everything else from v1 (loop, placement, restraint, build-once-subset).

**v2 acceptance test:** same P0CX28 run, now with the rigid ribosome; confirm
stable runs (no clashes/explosions), the nascent chain stays in/around the tunnel
and extrudes toward +x as it grows.

---

## Notes / likely gotchas

- TOPO's `buildCoarseGrainModel` rebuilds contacts from the PDB and sets positions
  from it — for this work you need to (a) supply precomputed `L×L` matrices and
  (b) override positions. Plan a small extension or a dedicated build path in
  `topo.translation` rather than calling the stock builder per step.
- `CustomNonbondedForce` requires one `addParticle` per System particle — that's
  why ribosome beads need dummy entries in the contact force (v2).
- Rigid ribosome = mass 0 + `rm_cons_0_mass`-style handling so constraints on
  fixed atoms don't error.
- Open choice to make at run time: the numeric `L0`.
- Defer (not needed now): variable elongation schedule, in-vivo timescale
  mapping, tunnel-confinement term.
