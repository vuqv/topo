# Isolated Protein Simulation Specification

## Purpose

Run a molecular-dynamics simulation of a single isolated protein (optionally
replicated into several non-interacting copies in one box) with the TOPO
Cα coarse-grained model, driven entirely by a control file (`md.ini`). This is
the canonical runner `topo.mdrun` (console `topo-mdrun`), distinct from the
ribosome-coupled / continuous-synthesis (CSP) runners.

The runner is a thin orchestration layer over a shared engine
(`topo.engine.build_system` / `setup_simulation` / `finalize_simulation`) and a
temperature protocol (`topo.mdrun.protocol`). It builds the coarse-grained
system, brings an OpenMM `Simulation` to a ready-to-step state (platform,
integrator, starting state, restart), runs one of two temperature protocols
(constant-temperature equilibrium or annealing/quenching), and finalizes.

Scope boundaries:

- The structure-based (contact) non-bonded force is built by
  `build_nonbonded_interaction`; its internal behavior (native contacts, domain
  `nscale` scaling, IDR/disorder masking) is specified elsewhere
  (`nscale-optimization` and `intrinsically-disordered-regions` specs). This spec
  treats that build as a dependency: it must succeed, and its failure is fatal.
- "Isolated" means no external scenery (no ribosome, no wall); multiple copies,
  when requested, are mutually non-interacting.

## Requirements

### Requirement: Control file is a flat key/value list

The control file SHALL be parsed as a flat list of `key = value` lines. Section
headers are optional and carry no meaning — every setting is read into one flat
namespace regardless of which (if any) header it sits under. Inline comments
starting with `#` or `;` SHALL be ignored. A key set more than once SHALL be a
hard error (no last-one-wins). A file containing no settings SHALL be an error.
Underscores in integer step counts (e.g. `500_000`) SHALL be accepted.

#### Scenario: Duplicate key rejected

- **WHEN** the control file sets the same key twice
- **THEN** parsing fails with an error identifying the repeated key
- **AND** no simulation is run

#### Scenario: Header placement is irrelevant

- **WHEN** a setting appears under any section header, several, or none
- **THEN** it is applied identically

### Requirement: Configuration validation and coupled options

The parser SHALL enforce the physical dependencies between options:

- Pressure coupling (`pcoupl = yes`) SHALL require periodic boundary conditions
  (`pbc = yes` with a valid box).
- `ref_t` and `tau_t` SHALL carry temperature units only when temperature
  coupling is on (`tcoupl = yes`).
- Enabling `restart` SHALL force energy minimization off, regardless of the
  `minimize` setting.
- `pbc = yes` with a missing or unparseable `box_dimension` SHALL disable PBC
  (keeping the flag consistent with the absent box) rather than proceed with an
  undefined box.

#### Scenario: Pressure coupling without a box

- **WHEN** `pcoupl = yes` but `pbc` is off
- **THEN** configuration fails with an assertion error

#### Scenario: Restart disables minimization

- **WHEN** `restart = yes` and `minimize = yes`
- **THEN** minimization is forced off (the loaded state, not the input geometry,
  is simulated)

#### Scenario: Invalid box turns PBC off

- **WHEN** `pbc = yes` but `box_dimension` is empty or unparseable
- **THEN** PBC is disabled and the run proceeds without a periodic box

### Requirement: Periodic box specification

When PBC is enabled with a valid `box_dimension`, a scalar `L` SHALL define a
cubic box `[L, L, L]` and an explicit `[x, y, z]` list SHALL define a
rectangular box. The box vectors SHALL be set on both the topology and the
system. When PBC is off, `box_dimension` SHALL be `None`.

#### Scenario: Scalar box is cubic

- **WHEN** `pbc = yes` and `box_dimension = 20`
- **THEN** the periodic box is cubic 20 × 20 × 20 nm

#### Scenario: List box is rectangular

- **WHEN** `pbc = yes` and `box_dimension = [10, 12, 15]`
- **THEN** the periodic box is rectangular with those edge lengths

### Requirement: Coarse-grained model build

The system SHALL be built as an alpha-carbon-only model from the input structure
(`pdb_file`), comprising bonds, angles, periodic torsions, Yukawa
electrostatics, and the structure-based (contact) non-bonded force. Per-residue
mass and charge SHALL be assigned by residue type; the excluded-volume radius
SHALL be the per-residue structure-derived Karanicolas–Brooks `Rmin/2` from the
contact build.

#### Scenario: Model terms present

- **WHEN** the model is built
- **THEN** it contains CA bonds, angles, periodic torsions, Yukawa
  electrostatics, and the contact non-bonded force

### Requirement: Bond treatment is rigid or flexible, never both

Bond treatment SHALL be selected by `constraints`: `AllBonds` (default) makes
every bond a rigid distance constraint at its equilibrium length and adds no
harmonic bond force; `None` (or `'None'`/`'none'`/empty) makes bonds flexible
via a harmonic bond force and adds no constraints. The two modes are mutually
exclusive. Any other value SHALL raise an error.

#### Scenario: Rigid bonds

- **WHEN** `constraints = AllBonds`
- **THEN** each bond is added as a distance constraint and no harmonic bond
  force is created

#### Scenario: Flexible bonds

- **WHEN** `constraints = None`
- **THEN** a harmonic bond force is created and no constraints are added

#### Scenario: Unknown constraint mode

- **WHEN** `constraints` is set to an unrecognized value
- **THEN** the build fails with an error

### Requirement: Structure-based non-bonded build is mandatory

The contact (structure-based) non-bonded interaction SHALL be built during model
construction. If `build_nonbonded_interaction` raises, the build SHALL fail with
a fatal error naming the domain/STRIDE inputs; it SHALL NOT be silently swallowed
(which would run an incomplete force field with no native contacts).

#### Scenario: Contact build failure is fatal

- **WHEN** the structure-based non-bonded build raises an exception
- **THEN** model construction raises a `RuntimeError` and no simulation runs

### Requirement: Build-time structure sanity checks

On a fresh (non-restart) build the system object SHALL be validated: bond
distances SHALL always be checked, and an initial-energy / large-force check
SHALL run on the input structure. On restart the large-force / initial-energy
check SHALL be skipped, because the loaded checkpoint state — not the input
geometry — is what is simulated.

#### Scenario: Fresh build checks input structure

- **WHEN** a fresh run builds the model
- **THEN** bond distances and initial large forces are checked

#### Scenario: Restart skips the initial-energy check

- **WHEN** `restart = yes`
- **THEN** the build skips the large-force / initial-energy check on the input
  structure

### Requirement: Multi-copy non-interacting replication

When `n_copies > 1` the single-chain model SHALL be replicated into that many
mutually non-interacting copies within one simulation, each translated along x
by `copy_shift` nm relative to the previous. A multi-copy run SHALL additionally
write a `<outname>_multi.psf` topology for the combined trajectory. The default
`n_copies = 1` is a single chain.

#### Scenario: Replicate into copies

- **WHEN** `n_copies = 10`
- **THEN** the system holds 10 non-interacting copies, x-shifted by `copy_shift`
- **AND** `<outname>_multi.psf` is written alongside `<outname>.psf`

### Requirement: Center-of-mass motion removal is opt-in

A `CMMotionRemover` SHALL be added only when `nstcomm` is set. When `nstcomm` is
unset (default) no COM-motion removal is added, because COM removal suits a
single chain but couples the drift of multiple independent chains.

#### Scenario: COM removal off by default

- **WHEN** `nstcomm` is not set
- **THEN** no `CMMotionRemover` is added

#### Scenario: COM removal when requested

- **WHEN** `nstcomm = 1000`
- **THEN** a `CMMotionRemover` is added at that frequency

### Requirement: Topology and native-reference dumps

The build SHALL write the single-chain topology as `<outname>.psf` and the input
structure coarse-grained to Cα resolution as `<pdbstem>_CA.pdb` (named from the
input structure, as an alignment reference for the trajectory).

#### Scenario: Reference artifacts written

- **WHEN** the model is built from `P0CX28_clean.pdb`
- **THEN** `<outname>.psf` and `P0CX28_clean_CA.pdb` are written to the output
  directory

### Requirement: Integrator and platform selection

The simulation SHALL use a Langevin integrator created at `ref_t` with time
constant `tau_t`, timestep `dt`, and the configured constraint tolerance. The
protocol driver sets the per-stage temperature before stepping, so the
integrator's initial temperature does not affect the run. The platform SHALL be
CUDA (mixed precision, device 0) when `device = GPU`, otherwise the CPU platform
using `ppn` threads.

#### Scenario: GPU platform

- **WHEN** `device = GPU`
- **THEN** the CUDA platform is used with mixed precision on device 0

#### Scenario: CPU platform

- **WHEN** `device = CPU` and `ppn = 4`
- **THEN** the CPU platform is used with 4 threads

### Requirement: Fresh-run starting coordinates

For a fresh run, starting coordinates SHALL come from `init_position` when set,
otherwise from the structure used to build the system (`pdb_file`). Input
coordinates SHALL be used as-is with no automatic translation (so multi-copy and
geometry-dependent layouts are preserved). When `init_position` is given, its
atom count SHALL match the system's particle count, else the run SHALL exit with
an error.

#### Scenario: Default coordinate source

- **WHEN** `init_position` is not set
- **THEN** starting coordinates are taken from `pdb_file` unchanged

#### Scenario: init_position atom-count mismatch

- **WHEN** `init_position` has a different atom count than the built system
- **THEN** the run exits with an error

### Requirement: Fresh-run velocities and minimization

For a fresh run, velocities SHALL be drawn from the Boltzmann distribution at
`ref_t`. Energy minimization of the input structure SHALL run when `minimize` is
on (the default for a fresh run) and SHALL be skipped on restart.

#### Scenario: Boltzmann velocities

- **WHEN** a fresh run starts
- **THEN** initial velocities are assigned from the Boltzmann distribution at
  `ref_t`

### Requirement: Restart from checkpoint

When `restart = yes` and the checkpoint file exists, positions AND velocities
SHALL be loaded from the checkpoint and the step count resumed from the saved
state; the loaded-state potential energy SHALL be reported. When `restart = yes`
but the checkpoint is missing, the runner SHALL warn and start fresh (coordinates
from `init_position` or `pdb_file`, Boltzmann velocities at `ref_t`).

#### Scenario: Restart with checkpoint present

- **WHEN** `restart = yes` and the checkpoint exists
- **THEN** positions, velocities, and step count are restored from it

#### Scenario: Restart with checkpoint missing

- **WHEN** `restart = yes` but the checkpoint file is absent
- **THEN** the runner warns and starts a fresh run

### Requirement: Constant-temperature equilibrium protocol

When `anneal = no` (default) the run SHALL be a single constant-temperature
production stage of `md_steps` at `ref_t`, written to `<outname>.dcd` and
`<outname>.log`.

#### Scenario: Equilibrium run

- **WHEN** `anneal = no` and `md_steps = 1_000_000`
- **THEN** the run integrates 1,000,000 steps at `ref_t` into `<outname>.dcd`/
  `.log`

### Requirement: Annealing requires its parameters

When `anneal = yes` the configuration SHALL require temperature coupling
(`tcoupl = yes`), a `t_high` (high/unfolding temperature), and `anneal_steps`
greater than 0. `anneal_ramp` SHALL be either `jump` or `linear`. `anneal_steps`
is separate from `md_steps`; the grand total is `quench_steps() + md_steps`, and
`ref_t` is the low/refold temperature (there is no separate `t_low`).

#### Scenario: Annealing without t_high

- **WHEN** `anneal = yes` but `t_high` is not set
- **THEN** configuration fails with an error

#### Scenario: Invalid ramp mode

- **WHEN** `anneal_ramp` is neither `jump` nor `linear`
- **THEN** configuration fails with an error

### Requirement: Annealing protocol — quench then production

When `anneal = yes` the run SHALL have two phases. The quench phase writes
`<outname>_quench.dcd`/`.log`; the production phase then runs `md_steps` at
`ref_t` and writes `<outname>.dcd`/`.log`. For `anneal_ramp = jump`, the quench
phase holds at `t_high` for `anneal_steps` and the drop to `ref_t` occurs at the
phase boundary (production starts thermostatted at `ref_t`). For
`anneal_ramp = linear`, the hold at `t_high` is followed by discrete cooling
increments (`anneal_ramp_increments` stages over `anneal_ramp_steps`) ending
exactly at `ref_t`. Between phases the step/time clock SHALL be reset so
production is a clean run from step 0 while positions and velocities carry over.

#### Scenario: Jump quench

- **WHEN** `anneal = yes` and `anneal_ramp = jump`
- **THEN** the quench holds at `t_high` for `anneal_steps` (into `_quench.*`)
- **AND** production runs `md_steps` at `ref_t` (into `<outname>.*`)

#### Scenario: Linear ramp quench

- **WHEN** `anneal = yes` and `anneal_ramp = linear`
- **THEN** the quench holds at `t_high` then cools in
  `anneal_ramp_increments` steps to `ref_t` over `anneal_ramp_steps`
- **AND** the final increment lands exactly on `ref_t`

### Requirement: Restart under annealing resumes production only

A restart under an annealing configuration SHALL skip the quench phase entirely
and resume the production phase from the checkpoint's step count, because the
checkpoint only ever holds production state (the quench writes no checkpoint) and
the quench is a one-time preparation.

#### Scenario: Restart skips the quench

- **WHEN** `restart = yes` under an annealing config with a production checkpoint
- **THEN** the quench phase is skipped and production resumes from the saved step

### Requirement: Output reporters

Each run phase SHALL attach a trajectory reporter (`<outname>[suffix].dcd` every
`nstxout` steps, honoring the periodic box when PBC is on), a fixed-width log
reporter (`<outname>[suffix].log` every `nstlog` steps), and — for phases that
own the checkpoint — a checkpoint reporter (`<outname>.chk` every `nstchk`, which
defaults to `nstxout` when unset). The quench phase SHALL use the `_quench`
suffix and SHALL NOT write a checkpoint, so the checkpoint only ever holds
production state.

#### Scenario: Production reporters

- **WHEN** the production phase runs
- **THEN** it writes `<outname>.dcd`, `<outname>.log`, and `<outname>.chk` at
  their configured frequencies

#### Scenario: Quench writes no checkpoint

- **WHEN** the quench phase runs
- **THEN** it writes `<outname>_quench.dcd`/`.log` but no checkpoint

### Requirement: Finalization artifacts

On completion the runner SHALL save the final checkpoint, write the last
conformation as `<outname>_final.pdb` (usable to seed a later run via
`init_position`, using the possibly-multi-copy topology and honoring the
periodic box when PBC is on), and close out run metadata (`<outname>_runinfo.log`).

#### Scenario: Final outputs

- **WHEN** the run finishes
- **THEN** `<outname>.chk`, `<outname>_final.pdb`, and the run-metadata log are
  written
