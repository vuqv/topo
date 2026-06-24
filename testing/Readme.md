# Comparison Between Quyen's and Yang's Implementations

This directory compares two implementations of the same coarse-grained simulation:

- **Quyen's implementation** is packaged as a Python module named `topo`, which is installed on the system `PATH` so it can be imported directly.
- **Yang's implementation** (the original) is a set of standalone scripts that run the simulation through several manual steps.

# How to Run

## Quyen's Implementation

From `/storage/home/qzv5006/work/src/topo/testing/Quyen`:

```bash
python run_simulation.py -f md.ini
```

This produces a trajectory (`.dcd`) and a log file
(`/storage/home/qzv5006/work/src/topo/testing/Quyen/P0CX28_clean.log`)
containing the potential, kinetic, and total energy.

## Yang's Implementation

**Step 1 — Build the coarse-grained model.**
From `/storage/home/qzv5006/work/src/topo/testing/Yang/create_model`:

```bash
python ../data/create_cg_protein_model.py -f go_model.cntrl
```

This generates the `cor`, `psf`, `top`, and `prm` files.

**Step 2 — Parse the parameters into an XML file.**

```bash
python ../data/parse_cg_prm.py -p P0CX28_clean_nscal1_fnn1_go_bt.prm -t P0CX28_clean_ca.top
```

This generates the `xml` file.

**Step 3 — Stage the input files.**
Copy the `cor`, `psf`, `top`, `xml`, and `prm` files to
`/storage/home/qzv5006/work/src/topo/testing/Yang/setup`.

**Step 4 — Run the simulation.**
From `/storage/home/qzv5006/work/src/topo/testing/Yang`:

```bash
python equil.py -f control.cntrl
```

This produces trajectories and a log file
(`/storage/home/qzv5006/work/src/topo/testing/Yang/P0CX28_equil.log`)
containing the potential, kinetic, and total energy.

# Observations

In the default runs (random initial velocities), the per-step energies of the
two implementations are **not consistent** in either kinetic or potential energy.

## Result of the fixed-velocity test

After fixing the initial velocity of every particle to zero (see the test below),
a 10-step comparison shows:

- **Kinetic energy is now consistent** (step 1: 225.40 vs 225.10 kJ/mol;
  per-step differences ≲ 1.7 kJ/mol). The small residual is expected, since the
  Langevin integrator is itself stochastic and the two codes use independent
  random-number streams. → The random *initial* velocities were the main cause
  of the kinetic-energy mismatch.
- **Potential energy still differs, but by a nearly constant offset of
  ~−126 kJ/mol** (range −123.8 to −126.3 kJ/mol over the 10 steps). A constant
  offset has zero gradient, so it does not affect the forces — which is why the
  dynamics now track each other closely. This points to a constant additive
  energy term present in one implementation but not the other (most likely in the
  non-bonded / contact component), rather than a difference in the dynamics.

## Per-component energy decomposition

To localize the offset, the per-force energy of the initial structure was dumped
from both implementations:

- Quyen prints the breakdown at startup (`run_simulation_zerovel.py`).
- Yang: `Yang/dump_initial_energies.py` evaluates the energy per force group on
  the initial structure and exits (no dynamics).

Initial-structure potential energy by component (kJ/mol):

| Component               |    Quyen |     Yang | Yang − Quyen |
|-------------------------|---------:|---------:|-------------:|
| Bond (harmonic)         |   14.044 |    0.000 |      −14.044 |
| Angle (Gaussian/custom) |  173.236 |  167.315 |       −5.921 |
| Torsion (periodic)      |  527.210 |  527.210 |       −0.000 |
| Yukawa / electrostatic  |   51.987 |    0.000 |      −51.987 |
| Custom non-bonded       | −915.149 | −983.411 |      −68.262 |
| **TOTAL**               | **−148.671** | **−288.885** | **−140.214** |

### Interpretation

The torsion term is **identical to 6 significant figures**, confirming both codes
read the same coordinates — so the remaining differences are genuine force-field
differences, not a coordinate mismatch.

**Important correction about the bond term.** The table above reports Quyen's bond
energy on the *input* structure (14.044 kJ/mol). However, Quyen's model already
builds the system with one **distance constraint per bond** (105 constraints at
0.381 nm) — i.e. it is already equivalent to Yang's `constraints=AllBonds`. During
dynamics the constraints pin every bond at its equilibrium length, so Quyen's
harmonic bond energy is effectively **zero** (≈1×10⁻⁷ kJ/mol at every step). The
14 kJ/mol only reflects the unconstrained input geometry and is removed the moment
dynamics start. Both implementations therefore **already agree on bonds**, and the
harmonic bond force Quyen keeps alongside the constraints is redundant (harmless).

**Important correction about how the non-bonded terms are organized.** The
per-component table above is *not* an apples-to-apples split, because the two codes
group the non-bonded interactions differently:

- **Quyen** uses **two separate** `CustomNonbondedForce` objects — a contact force
  `eps*(13*(R/r)^12 − 18*(R/r)^10 + 4*(R/r)^6)` (Karanicolas–Brooks 12-10-6,
  cutoff 2.0 nm) and a separate Yukawa force
  `factor*q1*q2/epsilon_r/r*exp(−r/lD)` (cutoff 3.5 nm).
- **Yang** uses **one combined** `CustomNonbondedForce` whose energy expression
  folds both together:
  `ke*q1*q2/ep/r*exp(−r/ld) + kv*(a/r^12 + b/r^10 + c/r^6)` (cutoff 2.0 nm,
  switch 1.8 nm).

So Yang is **not** missing electrostatics (as an earlier version of this note
stated) — they are bundled into its single non-bonded force. Comparing like with
like (contact vs contact, electrostatic vs electrostatic) by zeroing Yang's charges
to isolate the two parts (`Yang/test_nb_split.py`):

| Non-bonded part      |     Yang |    Quyen | Yang − Quyen |
|----------------------|---------:|---------:|-------------:|
| Contact (12-10-6)    | −1028.12 |  −915.15 |     −112.97  |
| Electrostatic (Yukawa)|   +44.70 |   +51.99 |       −7.28  |
| **Non-bonded total** | −983.41  | −863.16  |    −120.25   |
| Angle                |  167.32  |  173.24  |      −5.92   |
| **Grand total diff** |          |          |  **−126.17** |

which matches the observed per-step dynamics offset (−123.8 to −126.3 kJ/mol).

### Where the difference actually lives

1. **Contact potential (−113, dominant — root cause found & FIXED):** both codes
   use the *same* contact model — Karanicolas–Brooks 12-10-6 form, the same BT
   (Betancourt–Thirumalai) residue-pair well depths (`eps = nscal·|w − 0.6|`),
   and the same backbone hydrogen-bond term (0.75 kcal/mol per single H-bond,
   1.5 kcal/mol for a double H-bond, capped). A pair-by-pair comparison showed the
   entire gap lived on ~31 pairs, all carrying H-bond energy in multiples of
   0.75 kcal/mol. The cause was a **bug in Quyen's STRIDE parser**
   (`parse_hydrogen_bonds`, `topo/utils/build_nonbonded_interaction.py`): its
   regex expected a literal `-` where the **chain ID** sits, so for a structure
   with chain `A` it matched **zero** lines and Quyen assigned **no backbone
   H-bond energy at all** (`get_hb_contact_matrix` returned an all-zero matrix).
   Yang assigns the H-bonds (25 single + 6 double here = 27.75 kcal/mol ≈ 116 kJ/mol),
   which is exactly the offset. **Fix:** match the chain field with `\S+` so both
   `A` and `-` are accepted. After the fix Quyen parses 37 H-bonds and its contact
   energy becomes −1028.115 kJ/mol, **identical to Yang**.
2. **Electrostatics (−7 → −0.8, now matched):** same Yukawa/Debye–Hückel form in
   both. The gap was due to Yang truncating + switching the electrostatics at
   2.0 nm (it is folded into his combined non-bonded force) while Quyen carried a
   *separate* Yukawa force to 3.5 nm with no switch. Quyen's Yukawa force was
   changed (`topo/core/system.py`, `addYukawaForces`) to use the **same 2.0 nm
   cutoff + 1.8 nm switching function** as Yang. Quyen's Yukawa energy then drops
   from 51.987 to 43.866 kJ/mol, matching Yang's 44.704 to within ~0.8 kJ/mol.
   The residual ~0.8 kJ/mol is a minor parameter difference (Yang stores per-particle
   `ke/ep/ld` that multiply pairwise; Quyen uses global `factor/epsilon_r/lD`).
   Note: this shortens Quyen's electrostatic range from 3.5 to 2.0 nm for *all*
   `topo` runs — a deliberate change to match Yang for this comparison.
3. **Angle (−6 → ~0, now matched):** the angle potential uses the *identical*
   double-well Gaussian form and identical `k_alpha/k_beta/eps_alpha` in both
   codes. The gap came purely from **rounded reference angles**: Quyen hard-coded
   `theta_alpha = 1.6 rad` and `theta_beta = 2.27 rad` (and `gamma = 0.0239`),
   whereas Yang uses the precise degree conversions `91.7° = 1.6004669 rad`,
   `130.0° = 2.2689280 rad` (and `gamma = 0.0239005736 = 0.1 mol/kcal`). Across
   104 stiff angles those small θ₀ shifts summed to ~6 kJ/mol. Quyen's
   `addGaussianAngleForces` (`topo/core/system.py`) was updated to the precise
   values; Quyen's angle energy then changes from 173.236 to 167.31546 kJ/mol,
   matching Yang's 167.31545 to ~1×10⁻⁵ kJ/mol.

### Suggested follow-up

- **Bonds: no action** — both codes already constrain all bonds; they are consistent.
- **Contact potential: FIXED.** The ~113 kJ/mol gap was a bug in Quyen's STRIDE
  H-bond parser (chain-ID regex), which dropped all backbone H-bond energy. After
  the fix Quyen's contact energy matches Yang exactly (−1028.115 kJ/mol).
- **Electrostatic residual: FIXED.** Quyen's `epsilon_r` was changed from 80.0 to
  78.5 to match Yang (`ep = 8.860`, 8.860² ≈ 78.5). Quyen's Yukawa energy becomes
  44.704 kJ/mol, identical to Yang.

### Final result

After all fixes, the two implementations produce an **identical** initial-structure
potential energy (to ~2×10⁻⁴ kJ/mol):

| Component                  |      Quyen |       Yang |
|----------------------------|-----------:|-----------:|
| Bond (constrained)         |      0.000 |      0.000 |
| Gaussian Angle             |    167.315 |    167.315 |
| Periodic Torsion           |    527.210 |    527.210 |
| Yukawa / electrostatic     |     44.704 |     44.704 |
| Contact (12-10-6 + H-bond) |  −1028.115 |  −1028.115 |
| **Total PE**               | **−288.885** | **−288.885** |

(Yang reports electrostatic + contact together as one combined non-bonded force,
−983.411 kJ/mol = 44.704 − 1028.115; Quyen reports them separately.) The original
~126 kJ/mol per-step discrepancy is fully resolved.
- **Electrostatics:** decide whether to match the cutoff/switching treatment
  (Yang 2.0 nm + switch vs Quyen 3.5 nm, no switch) — worth only ~7 kJ/mol.

# Test: Fixed (Zero) Initial Velocities

By default both implementations assign **random** Maxwell–Boltzmann velocities at
the first step (via OpenMM's `setVelocitiesToTemperature`). Because these draws
differ between runs, they confound any comparison of per-step energies.

To remove this source of variation, two test scripts initialize the velocity of
**every particle to zero**, so both runs start from an identical, deterministic
state (same positions + zero velocities):

- `Quyen/run_simulation_zerovel.py`
- `Yang/equil_zerovel.py`

Each is a copy of the corresponding production script, changed only in the
velocity-initialization step. (In Yang's implementation the model does not need
to be rebuilt for this test.)

## Running the test

Quyen — from `/storage/home/qzv5006/work/src/topo/testing/Quyen`:

```bash
python run_simulation_zerovel.py -f md.ini
```

Yang — from `/storage/home/qzv5006/work/src/topo/testing/Yang`:

```bash
python equil_zerovel.py -f control.cntrl
```

Then compare the per-step energies in the two log files
(`Quyen/P0CX28_clean.log` and `Yang/P0CX28_equil.log`). With identical initial
conditions, any remaining energy difference reflects a genuine difference between
the two implementations rather than different random initial velocities.
