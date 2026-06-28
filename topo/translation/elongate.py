"""Protein synthesis: the elongation loop (``topo.translation``).

This is the nascent-chain elongation driver described in ``DESIGN.md`` (the
consecutive rebuild-and-continue protocol, §2.1/§2.5) and ``PROMPT.md``.

**Build step v1 — mechanics only, no ribosome forces.** The simulated System is
the *nascent chain only*. The truncated ribosome is used purely as a source of
two fixed points — the **P-anchor** (P-site tRNA residue-76 ``R`` bead) and the
**A-anchor** (A-site tRNA residue-76 ``R`` bead) — for new-residue placement and
the C-terminus restraint target.

**Build step v2 — add the rigid ribosome forces** (enable with
``rigid_ribosome = yes``). The truncated ribosome is appended to the System as
fixed (mass-0) scenery and the ribosome<->nascent excluded-volume +
electrostatics are wired up (:mod:`topo.translation.ribosome`). Because the
P-anchor is itself a ribosome bead, the C-terminus is held a short distance into
the tunnel from it (``ptc_offset_nm``) so it does not clash with that bead.

What the loop does, for ``L = L0 .. N_full`` (``L`` = current nascent length):

1. **Precompute once** (before the loop): run TOPO's contact builder on the
   *full* native PDB -> ``R_full``, ``eps_full`` (``N_full x N_full``). STRIDE is
   run at most once and cached. Each length uses the top-left ``L x L`` block;
   STRIDE / heavy-atom analysis are never re-run per length (DESIGN §3.5).
2. **Build** the length-``L`` TOPO model on native residues ``1..L`` (bonds,
   angles, torsions, Yukawa, contacts), injecting the ``L x L`` contact subset
   instead of recomputing (``buildCoarseGrainModel(precomputed_contacts=...)``).
3. **Seed coordinates.** ``L == L0``: lay residues ``1..L0`` extended along the
   tunnel axis (+x) from the P-anchor (C-terminus at the P-anchor, N-terminus
   toward +x), one CG bond length apart. ``L > L0``: residues ``1..L-1`` from the
   previous step's final structure; the new residue ``L`` at the A-anchor + buffer.
4. **Restrain only residue ``L``** (the current C-terminus) to the P-anchor with
   a harmonic ``CustomExternalForce`` (``k = 83680 kJ/mol/nm^2`` = 200 kcal/mol/A^2).
   The hand-off is automatic: each rebuilt step restrains only its own C-terminus.
5. **Minimize** (relax the placement / the stretched new bond), draw Boltzmann
   velocities at ``ref_t``, **run ``n_steps_per_residue`` steps**, write the
   per-length outputs, and seed ``L+1`` from this step's final structure.

Use it as a CLI, driven by an INI control file (see :func:`read_elongate_config`
and ``topo/translation/README.md``)::

    topo-elongate -f elongate.ini
    python -m topo.translation -f elongate.ini

or call :func:`run_elongation` from your own script. Build / setup / finalize are
reused from :mod:`topo.engine`; this module is the outer loop over length ``L``.
"""
from __future__ import annotations

import argparse
import configparser
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit

import topo
from topo import engine
from topo.core.system import system as TopoSystem
from topo.translation.ribosome import (Ribosome, load_ribosome, append_ribosome,
                                       add_trna_tether, add_tunnel_wall,
                                       TRNA_TETHER_BOND_NM, TUNNEL_WALL_X0_NM,
                                       TUNNEL_WALL_K)
from topo.utils.config import strtobool
from topo.utils.nonbonded import build_nonbonded_interaction

# --- physical constants / conversions -------------------------------------
# CG protein bond length (nm); cold-start beads are laid one bond length apart so
# no bond is pre-stretched (matches model_parameters['topo']['bond_length_protein']).
CG_BOND_LENGTH_NM = 0.381
# Tunnel / exit axis: the working ribosome is X-aligned (FILES.md), so the chain
# extrudes toward +x.
TUNNEL_AXIS = np.array([1.0, 0.0, 0.0])
# Default C-terminus position-restraint force constant (kJ/mol/nm^2 = OpenMM units;
# 83680 = 200 kcal/mol/A^2). Used only in the non-tether (position-restraint) mode.
RESTRAINT_K_KJ = 83680.0


# --------------------------------------------------------------------------
# Anchors
# --------------------------------------------------------------------------
def read_anchor(pdb_file: str, segid: str, resid: int = 76,
                bead: str = "R") -> np.ndarray:
    """Return the coordinate (nm) of one ribosome bead, e.g. a tRNA anchor.

    Parses the (CG) ribosome PDB for the single atom whose **segID** (columns
    73-76), **residue number** and **atom name** match ``segid`` / ``resid`` /
    ``bead``. Used to pick the P-anchor (``segid='PtR'``) and A-anchor
    (``segid='AtR'``) from the truncated ribosome (both keep residue 76 — see
    FILES.md). PDB coordinates are in angstrom; the returned vector is in nm.

    Raises
    ------
    ValueError
        If no matching atom (or more than one) is found.
    """
    matches = []
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            res_seq = line[22:26].strip()
            seg = line[72:76].strip()
            if name == bead and seg == segid and res_seq == str(resid):
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                matches.append(np.array([x, y, z]) / 10.0)  # angstrom -> nm
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one atom (segid={segid!r}, resid={resid}, "
            f"name={bead!r}) in {pdb_file}, found {len(matches)}.")
    return matches[0]


# --------------------------------------------------------------------------
# Build-once-subset contacts
# --------------------------------------------------------------------------
def precompute_contacts(full_pdb: str,
                        domain_def: Optional[str] = None,
                        stride_output_file: Optional[str] = None
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Run TOPO's contact builder once on the full native PDB (DESIGN §3.5).

    Returns ``(R_full, eps_full)`` -- the ``N_full x N_full`` well-position
    (nm) and well-depth (kJ/mol) matrices. STRIDE is run at most once here (and
    cached by :func:`build_nonbonded_interaction`); each length later reuses the
    top-left ``L x L`` block, so neither STRIDE nor the heavy-atom analysis is
    ever re-run per length.
    """
    print("=" * 66)
    print("[ Precompute contacts (build-once-subset) ]")
    print("=" * 66)
    print(f"Running TOPO contact builder on full native structure: {full_pdb}")
    R_full, eps_full = build_nonbonded_interaction(full_pdb, domain_def,
                                                   stride_output_file)
    print(f"  full contact matrices: {R_full.shape}")
    return R_full, eps_full


# --------------------------------------------------------------------------
# Subset native structure (residues 1..L)
# --------------------------------------------------------------------------
def write_subset_structure(full_pdb: str, L: int, out_pdb: str) -> None:
    """Write a CA-only PDB of the first ``L`` residues of ``full_pdb``.

    This length-``L`` native structure supplies the bonded terms and per-residue
    mass/charge/radii for the length-``L`` build; its coordinates set the native
    contact distances and CG bond lengths only (the *simulated* coordinates come
    from the seeding scheme, not from here). Residues are taken in file order so
    particle ``i`` (``0..L-1``) corresponds to native residue ``i+1``, matching
    the ``L x L`` contact subset.
    """
    import MDAnalysis as mda  # local import: heavy, only needed for slicing

    u = mda.Universe(full_pdb)
    ca = u.select_atoms("protein and name CA")
    if len(ca) < L:
        raise ValueError(f"requested L={L} but {full_pdb} has only {len(ca)} CA atoms.")
    subset = ca[:L]
    subset.write(out_pdb)


# --------------------------------------------------------------------------
# Length-L model build (dedicated build-once-subset path)
# --------------------------------------------------------------------------
def build_length_model(sub_pdb: str, R_L: np.ndarray, eps_L: np.ndarray,
                       constraints="AllBonds", model: str = "topo"):
    """Build a length-``L`` TOPO model, injecting the ``L x L`` contact subset.

    This is the dedicated "build-once-subset" path the design anticipates
    (DESIGN §3.5 / PROMPT notes): it mirrors
    :meth:`topo.models.buildCoarseGrainModel` *exactly* — reusing the same
    :class:`topo.core.system.system` methods for bonds, angles, torsions, Yukawa
    electrostatics and per-residue mass/charge/radii — but **supplies the
    precomputed contact matrices directly** instead of recomputing them (so
    STRIDE / heavy-atom analysis are not re-run per length). PBC is off and the
    build-time large-force check is skipped (the simulated coordinates come from
    the seeding scheme and are minimized explicitly, so the native-structure
    energy is irrelevant here).

    Parameters
    ----------
    sub_pdb : str
        Length-``L`` native CA structure (bonded terms + per-residue properties).
    R_L, eps_L : numpy.ndarray
        The ``L x L`` contact well-position (nm) and well-depth (kJ/mol) blocks.
    constraints : str or None
        ``'AllBonds'`` (rigid, default) or ``None`` (flexible harmonic bonds) --
        same semantics as :meth:`buildCoarseGrainModel`.

    Returns
    -------
    topo.core.system.system
        The built model (``.system``, ``.topology``, ``.positions`` ready to use).
    """
    print("=" * 66)
    print(f"[ Length-L build ] {sub_pdb}  (build-once-subset contacts)")
    print("=" * 66)

    topo_model = TopoSystem(sub_pdb, model)

    # CA topology (atoms, bonds, angles) -- identical to buildCoarseGrainModel.
    topo_model.getCAlphaOnly()
    topo_model.getAtoms()
    topo_model.getBonds()
    topo_model.getAngles()

    n_ca = topo_model.n_atoms
    for nm, m in (("distance", R_L), ("energy", eps_L)):
        if m.shape != (n_ca, n_ca):
            raise ValueError(
                f"{nm} subset matrix has shape {m.shape} but the length-L model has "
                f"{n_ca} CA atoms; expected ({n_ca}, {n_ca}).")

    # Resolve the bond-constraint mode (same vocabulary as buildCoarseGrainModel).
    if constraints is None or str(constraints).strip().lower() == "none":
        use_constraints = False
    elif str(constraints).strip().lower() == "allbonds":
        use_constraints = True
    else:
        raise ValueError(f"Invalid constraints option: {constraints!r}. "
                         f"Expected 'AllBonds' or None.")
    print(f"  chains={topo_model.n_chains}  CA atoms={n_ca}  "
          f"bonds={topo_model.n_bonds}  angles={topo_model.n_angles}  "
          f"bonds: {'rigid (AllBonds)' if use_constraints else 'flexible (harmonic)'}")

    if use_constraints:
        for bond in topo_model.bonds:
            topo_model.system.addConstraint(bond[0].index, bond[1].index,
                                            topo_model.bonds[bond][0])

    # Per-residue particle properties.
    topo_model.setCAMassPerResidueType()
    topo_model.setCAChargePerResidueType()
    topo_model.setCARadiusPerResidueType()

    # Bonded terms.
    topo_model.setBondForceConstants()
    if not use_constraints:
        topo_model.addHarmonicBondForces()
    topo_model.addGaussianAngleForces()
    topo_model.getTorsions()
    topo_model.addPeriodicTorsionForce()

    # Electrostatics (no PBC for the nascent-chain runs).
    topo_model.addYukawaForces(use_pbc=False)

    # Structure-based contacts: the injected L x L subset (not recomputed).
    topo_model.distance_matrix = R_L
    topo_model.energy_matrix = eps_L
    topo_model.addCustomNonBondedForce(R_L, eps_L, use_pbc=False)

    # Assemble the OpenMM System. Skip the large-force check (native-structure
    # energy is irrelevant; the seeded structure is minimized explicitly later).
    topo_model.createSystemObject(minimize=False, check_bond_distances=True,
                                  check_large_forces=False)
    return topo_model


# --------------------------------------------------------------------------
# Coordinate seeding
# --------------------------------------------------------------------------
# Small transverse zig-zag amplitude (nm) applied to the cold-start layout. A
# perfectly collinear chain makes every bond angle exactly 180 degrees, where the
# angle-force gradient is singular (-> NaN on the first minimization step), so the
# beads are offset by +/- this amount perpendicular to the tunnel axis to break
# collinearity. Tiny relative to the 0.381 nm bond, so the chain stays "extended
# along the axis"; the bond-length constraints absorb the negligible change.
COLD_START_ZIGZAG_NM = 0.03


def cold_start_positions(L0: int, p_anchor: np.ndarray) -> np.ndarray:
    """Extended cold-start layout for the first length ``L0`` (DESIGN §2.5).

    Residues ``1..L0`` are laid along the tunnel axis (+x) from the P-anchor:
    the C-terminus (residue ``L0``) sits *at* the P-anchor and the N-terminus
    (residue 1) points toward the exit (+x), one CG bond length apart. A small
    alternating transverse zig-zag (``COLD_START_ZIGZAG_NM``) breaks the exact
    collinearity that would make the angle force singular. Returns an ``(L0, 3)``
    array in nm (row ``i`` = residue ``i+1``).
    """
    positions = np.empty((L0, 3))
    for i in range(L0):  # i = 0..L0-1  -> native residue i+1
        # residue L0 (i = L0-1) at offset 0 (the P-anchor); residue 1 furthest +x.
        offset = (L0 - 1 - i) * CG_BOND_LENGTH_NM
        pos = p_anchor + offset * TUNNEL_AXIS
        # Alternate a small displacement along +y / -y to avoid a straight chain.
        pos = pos + ((-1) ** i) * COLD_START_ZIGZAG_NM * np.array([0.0, 1.0, 0.0])
        positions[i] = pos
    return positions


def seed_positions(prev_final: np.ndarray, a_anchor: np.ndarray,
                   buffer_nm: float) -> np.ndarray:
    """Seed length ``L`` from the previous final structure + the new residue.

    Residues ``1..L-1`` keep their coordinates from step ``L-1``'s final
    structure (``prev_final``, shape ``(L-1, 3)`` nm); the new C-terminal residue
    ``L`` is placed at the A-anchor offset by ``buffer_nm`` along +x (the buffer
    clears excluded volume so the new bead does not get a huge ``(sigma/r)^12``
    kick — DESIGN §2.5 / invariant 6). Returns an ``(L, 3)`` array in nm.
    """
    new_residue = a_anchor + buffer_nm * TUNNEL_AXIS
    return np.vstack([prev_final, new_residue[None, :]])


# --------------------------------------------------------------------------
# C-terminus restraint
# --------------------------------------------------------------------------
def add_cterm_restraint(system: mm.System, particle_index: int,
                        anchor_nm: np.ndarray, k: float) -> mm.Force:
    """Add a harmonic restraint pulling one particle toward a fixed anchor.

    ``U = k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)`` via a ``CustomExternalForce``.
    Only the current C-terminus (``particle_index = L-1``) is restrained -- to the
    P-anchor -- so the tether hand-off between lengths is automatic (each rebuilt
    step restrains only its own C-terminus). ``k`` is the force constant in OpenMM
    units (kJ/mol/nm^2).
    """
    restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", k)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    restraint.addParticle(int(particle_index),
                          [float(anchor_nm[0]), float(anchor_nm[1]), float(anchor_nm[2])])
    system.addForce(restraint)
    return restraint


# --------------------------------------------------------------------------
# PDB writing helper
# --------------------------------------------------------------------------
def _write_pdb(topology, positions_nm: np.ndarray, path: str) -> None:
    """Write a PDB from a topology and an ``(N, 3)`` nm coordinate array."""
    coords = [mm.Vec3(float(x), float(y), float(z)) for x, y, z in positions_nm] * unit.nanometer
    with open(path, "w") as fh:
        mm.app.PDBFile.writeFile(topology, coords, fh)


class NascentDCDReporter:
    """A DCD reporter that writes only the first ``n_keep`` atoms each frame.

    Used in build step v2 so the (large, static) ribosome beads are **not** written
    to the trajectory every frame -- only the nascent chain (indices ``0..n_keep-1``)
    is saved. Mirrors :class:`openmm.app.DCDReporter` but slices the positions and
    uses a fixed ``n_keep``-atom topology, so the DCD header records ``n_keep`` atoms
    and pairs with the nascent-only PSF.
    """

    def __init__(self, file, reportInterval, nascent_topology, n_keep, append=False):
        self._reportInterval = reportInterval
        self._topology = nascent_topology
        self._n = n_keep
        self._append = append
        self._out = open(file, "r+b" if append else "wb")
        self._dcd = None

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        # No PBC in these runs; positions only.
        return {"steps": steps, "periodic": False, "include": ["positions"]}

    def report(self, simulation, state):
        if self._dcd is None:
            self._dcd = mm.app.DCDFile(
                self._out, self._topology, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, self._append)
        positions = state.getPositions(asNumpy=True)[:self._n]
        self._dcd.writeModel(positions)

    def __del__(self):
        try:
            self._out.close()
        except Exception:
            pass


def _dump_topology_psf(cgModel, path: str) -> None:
    """Write a PSF for the (nascent + ribosome) system via parmed.

    The model's own ``dumpTopology`` keys per-atom mass/charge off its nascent-only
    lists, so it cannot describe the v2 system. parmed reads masses/bonds straight
    from the OpenMM topology + System instead (charges default to 0 in the PSF --
    cosmetic; the real electrostatics live in the Yukawa force).
    """
    import parmed as pmd
    pmd.openmm.load_topology(cgModel.topology, system=cgModel.system).save(
        path, overwrite=True)


def _finalize_nascent(cfg, ctx, nascent_topology, n_keep: int,
                      start_epoch: float) -> None:
    """Finalize a v2 length writing a **nascent-only** final structure.

    Like :func:`topo.engine.finalize_simulation` but the written ``_final.pdb`` is
    only the first ``n_keep`` (nascent) atoms -- the rigid ribosome is dropped. The
    saved checkpoint still holds the **full** system (needed for a correct restart).
    """
    sim = ctx.simulation
    sim.saveCheckpoint(ctx.checkpoint)
    final_pdb = cfg.output_path("_final.pdb")
    pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
    pos = pos[:n_keep].value_in_unit(unit.nanometer)
    _write_pdb(nascent_topology, pos, final_pdb)
    print(f"Wrote last nascent conformation to {final_pdb}")
    topo.runinfo.write_run_end(ctx.runinfo_path, simulation=sim,
                               start_epoch=start_epoch, final_structure=final_pdb)
    print("--- Finished in %s seconds ---" % (time.time() - start_epoch))


# --------------------------------------------------------------------------
# Per-length configuration
# --------------------------------------------------------------------------
@dataclass
class ElongationParams:
    """Run parameters shared by every length (set once from the CLI)."""
    n_steps: int = 1000
    dt_ps: float = 0.015
    ref_t: float = 300.0
    tau_t: float = 0.01
    nstout: int = 50
    device: str = "CPU"
    ppn: int = 1
    # Flexible (harmonic) bonds by default -- NOT the package default of rigid
    # 'AllBonds'. The elongation step places the new residue at the A-anchor while
    # restraining it to the P-anchor (~0.9-1.1 nm away), so the new bond starts far
    # from its equilibrium length. A rigid distance constraint cannot be seeded that
    # far off (the constraint solver / minimizer diverges); a harmonic bond absorbs
    # the stretch, lets the minimization relax it, and the restraint then translocates
    # the residue A->P -- exactly the DESIGN §2.5 mechanism.
    constraints: object = None
    restraint_k: float = RESTRAINT_K_KJ   # kJ/mol/nm^2 (= 200 kcal/mol/A^2)
    buffer_nm: float = 0.4
    minimize: bool = True
    # Build step v2: append the truncated ribosome as rigid (mass-0) scenery and
    # wire the ribosome<->nascent excluded-volume + electrostatics. Default off
    # (v1 = nascent chain only); enable with `rigid_ribosome = yes` in the INI.
    rigid_ribosome: bool = False
    # How far into the tunnel (+x) from the P-anchor *bead* to hold the
    # C-terminus (nm). The P-anchor is the PtR residue-76 R bead, so in v2 a
    # zero offset would restrain the C-terminus on top of that ribosome bead --
    # a near-coincident excluded-volume clash that explodes. A small offset holds
    # the C-terminus just inside the tunnel, clearing the P-tRNA bead. None ->
    # auto: 0.0 in v1 (no ribosome), 0.4 nm in v2 (or the tether bond length when
    # the tRNA tether is on).
    ptc_offset_nm: Optional[float] = None
    # v2: tether the C-terminus to the P-site tRNA R bead the O'Brien way -- a bond
    # plus a CA(L-1)-CA(L)-tRNA orienting angle -- instead of a generic position
    # restraint. The angle aims the nascent chain down the tunnel (toward the exit)
    # so it extrudes N-first and folds outside, rather than balling up at the PTC.
    # Only used in v2 (needs the real tRNA bead); ignored in v1.
    trna_tether: bool = True
    # v2: O'Brien's one-sided planar tunnel wall on the nascent chain -- keeps beads
    # at x >= tunnel_wall_x0, so the chain can only extrude forward (+x, toward the
    # exit) and cannot fold back past the synthesis point into the ribosome interior.
    # x0_nm (nm) is the C-terminal-AA addition plane (the PTC / P-site where each new
    # residue is placed and tethered) ~ P-anchor x + tether bond length ~ 1.05 nm;
    # k in kJ/mol/nm^2. Applied throughout synthesis and the post-elongation phase.
    tunnel_wall: bool = True
    tunnel_wall_x0_nm: float = TUNNEL_WALL_X0_NM
    tunnel_wall_k: float = TUNNEL_WALL_K
    # v2 output: write only the nascent chain to the trajectory / PSF / final
    # structure (the rigid ribosome is static, so writing it every frame wastes
    # storage). The checkpoint still holds the full system. No effect in v1 (which
    # is already nascent-only). Set False to write the full system in v2.
    nascent_only_output: bool = True
    # Post-elongation phase (runs after the chain reaches its final length, for
    # post_elongation_steps steps; 0 = skip). 'ejection' releases the C-terminus
    # tether and lets the protein move (-> 'ejection/'); 'stallation' keeps the
    # restraint so the chain stays stalled on the ribosome (-> 'stallation/').
    post_elongation: str = "stallation"
    post_elongation_steps: int = 0


def _make_cfg(out_dir: Path, sub_pdb: str, seed_pdb: str,
              params: ElongationParams) -> topo.SimulationConfig:
    """Build a per-length :class:`SimulationConfig` for the engine helpers.

    Each length is a self-contained standalone run (its own output folder), so
    this mirrors a single ``topo-mdrun`` invocation: a constant-temperature
    production run of ``n_steps`` at ``ref_t``. The contact matrices are injected
    at build time (not via this config), so ``domain_def`` is left unset here.
    """
    cfg = topo.SimulationConfig()
    cfg.md_steps = params.n_steps
    cfg.dt = params.dt_ps * unit.picoseconds
    cfg.nstxout = params.nstout
    cfg.nstlog = params.nstout
    cfg.nstchk = params.nstout
    cfg.tcoupl = True
    cfg.ref_t = params.ref_t * unit.kelvin
    cfg.tau_t = params.tau_t / unit.picoseconds
    cfg.pbc = False
    cfg.box_dimension = None
    cfg.pdb_file = sub_pdb
    cfg.init_position = seed_pdb
    cfg.constraints = params.constraints
    cfg.output_dir = str(out_dir)
    cfg.outname = "traj"
    cfg.device = params.device
    cfg.ppn = params.ppn
    cfg.restart = False
    cfg.minimize = False  # we minimize the seeded structure explicitly below
    return cfg


# --------------------------------------------------------------------------
# Single length step
# --------------------------------------------------------------------------
def run_length(L: int, *, full_pdb: str, R_full: np.ndarray, eps_full: np.ndarray,
               p_anchor: np.ndarray, a_anchor: np.ndarray,
               prev_final: Optional[np.ndarray], out_root: Path,
               params: ElongationParams,
               ribo: Optional[Ribosome] = None,
               seed_override: Optional[np.ndarray] = None,
               restrain: bool = True,
               out_subdir: Optional[str] = None,
               n_steps_override: Optional[int] = None,
               label: Optional[str] = None) -> np.ndarray:
    """Build, seed, (restrain,) minimize and run one length-``L`` system.

    Used both for an elongation step and for the post-synthesis phase (§post-
    synthesis below). When ``ribo`` is given (build step v2), the rigid ribosome is
    appended as fixed (mass-0) scenery with the ribosome<->nascent cross-interactions
    (:func:`topo.translation.ribosome.append_ribosome`).

    Parameters that tailor the standard elongation step into the post-synthesis
    phase:

    - ``seed_override`` : use these ``(L, 3)`` nm nascent coordinates directly
      (no cold-start / new-residue placement) -- e.g. the fully-synthesized
      structure for an ejection / stall run.
    - ``restrain`` : if False, do **not** add the C-terminus restraint (ejection
      = the tether is released and the chain is free to move).
    - ``out_subdir`` : output folder name under ``out_root`` (default
      ``L_<L>``); e.g. ``ejection`` / ``stall``.
    - ``n_steps_override`` : run this many steps instead of ``params.n_steps``.
    - ``label`` : header text for the console banner.

    Returns the final **nascent** ``(L, 3)`` nm coordinate array.
    """
    print()
    print("#" * 66)
    print("# " + (label or (f"Nascent length L = {L}"
                            + ("  (+ rigid ribosome)" if ribo else ""))))
    print("#" * 66)

    out_dir = out_root / (out_subdir or f"L_{L:03d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. length-L native structure (bonded terms + per-residue properties) ---
    sub_pdb = str(out_dir / f"native_1_{L}.pdb")
    write_subset_structure(full_pdb, L, sub_pdb)

    # 2. inject the L x L contact subset (build-once-subset) -----------------
    R_L = np.ascontiguousarray(R_full[:L, :L])
    eps_L = np.ascontiguousarray(eps_full[:L, :L])
    cgModel = build_length_model(sub_pdb, R_L, eps_L, constraints=params.constraints)

    # 3. seed nascent coordinates --------------------------------------------
    if seed_override is not None:   # post-synthesis: use the given structure as-is
        if seed_override.shape[0] != L:
            raise ValueError(
                f"seed_override has {seed_override.shape[0]} residues but L = {L}.")
        nascent_pos = seed_override
    elif prev_final is None:        # cold start (L == L0)
        nascent_pos = cold_start_positions(L, p_anchor)
    else:                           # continue from previous length + new residue
        if prev_final.shape[0] != L - 1:
            raise ValueError(
                f"prev_final has {prev_final.shape[0]} residues but L-1 = {L - 1}.")
        nascent_pos = seed_positions(prev_final, a_anchor, params.buffer_nm)

    # v2 nascent-only output writes only the nascent chain to the trajectory /
    # PSF / final structure. Capture a nascent (L-atom) topology BEFORE appending
    # the ribosome (append mutates cgModel.topology), and write the nascent PSF now
    # (the model's dumpTopology keys off its nascent-only per-atom lists).
    nascent_only_v2 = ribo is not None and params.nascent_only_output
    nascent_topology = None
    if nascent_only_v2:
        nascent_topology = mm.app.PDBFile(sub_pdb).topology
        cgModel.dumpTopology(str(out_dir / "traj.psf"))

    # 3b. v2: append the rigid ribosome (mass-0 scenery + cross-interactions).
    if ribo is not None:
        append_ribosome(cgModel, ribo)
        positions = np.vstack([nascent_pos, ribo.coords_nm])
    else:
        positions = nascent_pos

    # 4. tether the current C-terminus (residue L) ---------------------------
    # (skipped for an ejection run: restrain=False -> the tether is released).
    # v2 + trna_tether: O'Brien peptidyl-tRNA linkage (bond + CA-CA-tRNA orienting
    # angle to the P-site R bead) -- aims the chain down the tunnel. Otherwise a
    # generic harmonic position restraint of residue L to the P-target point.
    if restrain:
        if ribo is not None and params.trna_tether:
            prev_index = (L - 2) if L >= 2 else None
            add_trna_tether(cgModel, L - 1, prev_index, ribo, L)
        else:
            add_cterm_restraint(cgModel.system, L - 1, p_anchor, params.restraint_k)

    # 4b. v2 tunnel wall: keep nascent beads at x >= x0 (no leaking through the
    # truncated-ribosome cutout). Applied in elongation and post-elongation alike.
    if ribo is not None and params.tunnel_wall:
        add_tunnel_wall(cgModel.system, range(L),
                        x0_nm=params.tunnel_wall_x0_nm, k=params.tunnel_wall_k)

    # 5. set up, minimize, run, finalize (reuse topo.engine) -----------------
    # v1: write a small seed PDB and feed it via init_position. v2: the full
    # system includes the (large, static) ribosome, so seed coordinates are set
    # directly on built.positions instead of writing an N-atom seed PDB.
    if ribo is None:
        seed_pdb = str(out_dir / "seed.pdb")
        _write_pdb(cgModel.topology, positions, seed_pdb)
        cfg = _make_cfg(out_dir, sub_pdb, seed_pdb, params)
        cgModel.dumpTopology(cfg.output_path(".psf"))
        built_positions = cgModel.positions
    else:
        cfg = _make_cfg(out_dir, sub_pdb, None, params)
        if not nascent_only_v2:
            # Full-system PSF (nascent-only PSF was already written above).
            _dump_topology_psf(cgModel, cfg.output_path(".psf"))
        built_positions = ([mm.Vec3(*r) for r in positions]) * unit.nanometer

    # Post-synthesis phases run for their own step count.
    if n_steps_override is not None:
        cfg.md_steps = n_steps_override

    built = engine.BuiltSystem(cgModel=cgModel, system=cgModel.system,
                               topology=cgModel.topology, positions=built_positions)

    start = time.time()
    ctx = engine.setup_simulation(cfg, built)
    if params.minimize:
        print("Minimizing seeded structure (relax placement / new bond)...")
        ctx.simulation.minimizeEnergy()
        # Re-draw velocities for the relaxed structure.
        ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)

    engine.attach_reporters(cfg, ctx.simulation, suffix="", total_steps=cfg.md_steps)
    if nascent_only_v2:
        # Swap the full-system DCD reporter for a nascent-only one (the rigid
        # ribosome is static -- no need to write its ~thousands of beads/frame).
        ctx.simulation.reporters[1] = NascentDCDReporter(
            cfg.output_path(".dcd"), cfg.nstxout, nascent_topology, L)

    ctx.simulation.step(cfg.md_steps)

    if nascent_only_v2:
        _finalize_nascent(cfg, ctx, nascent_topology, L, start)
    else:
        engine.finalize_simulation(cfg, ctx, built.topology, start)

    # 6. final NASCENT coordinates seed the next length ----------------------
    final_pdb = cfg.output_path("_final.pdb")
    final = mm.app.PDBFile(final_pdb).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]


# --------------------------------------------------------------------------
# Elongation loop
# --------------------------------------------------------------------------
def run_elongation(full_pdb: str, ribosome_pdb: str, *,
                   L0: int, L_max: Optional[int] = None,
                   out_root: str = "synth_out",
                   domain_def: Optional[str] = None,
                   stride_output_file: Optional[str] = None,
                   params: Optional[ElongationParams] = None) -> None:
    """Run the full nascent-chain elongation loop ``L = L0 .. L_max`` (v1).

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein (the nascent chain at full length).
    ribosome_pdb : str
        Truncated CG ribosome PDB. Always the source of the P-/A-anchor
        coordinates; in v2 (``params.rigid_ribosome``) it is also appended as the
        rigid (mass-0) ribosome scenery.
    L0 : int
        Starting nascent-chain length (cold-start layout).
    L_max : int, optional
        Final length; defaults to the full residue count ``N_full``.
    out_root : str
        Root output directory; each length writes to ``<out_root>/L_<L>/``.
    domain_def, stride_output_file : str, optional
        Passed to the one-time contact precompute (n_scale / STRIDE).
    params : ElongationParams, optional
        Per-length run parameters (defaults to the test settings).
    """
    if params is None:
        params = ElongationParams()
    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # Anchors (fixed points from the rigid truncated ribosome).
    p_anchor = read_anchor(ribosome_pdb, "PtR", resid=76, bead="R")
    a_anchor = read_anchor(ribosome_pdb, "AtR", resid=76, bead="R")
    print(f"P-anchor (PtR 76 R): {p_anchor} nm")
    print(f"A-anchor (AtR 76 R): {a_anchor} nm")

    # Where the C-terminus is held / seeded, measured into the tunnel (+x) from the
    # P-anchor bead. The P-anchor is the PtR-76 R bead, so it must not sit on top of
    # it (clash). Auto: 0 in v1; in v2, the tether bond length when the tRNA tether
    # is on (the bond sets the distance), else 0.4 nm for the position restraint.
    tether_on = params.rigid_ribosome and params.trna_tether
    offset = params.ptc_offset_nm
    if offset is None:
        if tether_on:
            offset = TRNA_TETHER_BOND_NM
        elif params.rigid_ribosome:
            offset = 0.4
        else:
            offset = 0.0
    p_target = p_anchor + offset * TUNNEL_AXIS
    if offset:
        print(f"C-terminus restraint/cold-start target: P-anchor + {offset} nm "
              f"into tunnel (+x) = {p_target} nm")

    # Build step v2: load the rigid ribosome once (identical in every length).
    ribo = None
    if params.rigid_ribosome:
        ribo = load_ribosome(ribosome_pdb, model="topo")
        print(f"Rigid ribosome: {ribo.n} beads from {ribosome_pdb} "
              f"(appended as mass-0 scenery; ribosome<->nascent forces on).")

    # Build-once-subset contacts on the full native structure.
    R_full, eps_full = precompute_contacts(full_pdb, domain_def, stride_output_file)
    N_full = R_full.shape[0]
    if L_max is None:
        L_max = N_full
    if not (1 <= L0 <= L_max <= N_full):
        raise ValueError(f"require 1 <= L0 <= L_max <= N_full; got L0={L0}, "
                         f"L_max={L_max}, N_full={N_full}.")

    print()
    print(f"Elongating {full_pdb}: L = {L0} .. {L_max} (N_full = {N_full}), "
          f"{params.n_steps} steps/residue.")

    prev_final: Optional[np.ndarray] = None
    for L in range(L0, L_max + 1):
        prev_final = run_length(
            L, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
            p_anchor=p_target, a_anchor=a_anchor, prev_final=prev_final,
            out_root=out_path, params=params, ribo=ribo)

    print()
    print(f"Done. Elongated {L0} -> {L_max}. Per-length outputs under {out_path}/")

    # Post-elongation phase: once the chain reaches its final length, either release
    # the C-terminus tether and let the protein move (ejection) or keep it tethered
    # (stallation). Continues the same length-L_max system from the final structure.
    if params.post_elongation_steps > 0:
        phase = params.post_elongation.strip().lower()
        if phase not in ("ejection", "stallation"):
            raise ValueError(f"post_elongation must be 'ejection' or 'stallation', "
                             f"got {params.post_elongation!r}.")
        restrain = phase == "stallation"
        print()
        print(f"=== Post-elongation: {phase} (L = {L_max}, {params.post_elongation_steps} "
              f"steps, C-terminus restraint {'ON' if restrain else 'OFF'}) "
              f"-> {out_path / phase}/ ===")
        run_length(
            L_max, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
            p_anchor=p_target, a_anchor=a_anchor, prev_final=None,
            out_root=out_path, params=params, ribo=ribo,
            seed_override=prev_final, restrain=restrain, out_subdir=phase,
            n_steps_override=params.post_elongation_steps,
            label=f"Post-elongation: {phase} (L = {L_max})")
        print(f"Done. {phase.capitalize()} written to {out_path / phase}/")


# --------------------------------------------------------------------------
# INI control file
# --------------------------------------------------------------------------
@dataclass
class ElongateConfig:
    """Parsed contents of an elongation control file (``elongate.ini``).

    Bundles the run inputs (structures, the ``L0..L_max`` schedule, output
    directory, one-time-precompute options) with the per-length
    :class:`ElongationParams`.
    """
    pdb_file: str
    ribosome: str
    L0: int
    L_max: Optional[int] = None
    outdir: str = "synth_out"
    domain_def: Optional[str] = None
    stride_output_file: Optional[str] = None
    params: ElongationParams = field(default_factory=ElongationParams)
    config_file: Optional[str] = None


def read_elongate_config(config_file: str, verbose: bool = True) -> ElongateConfig:
    """Parse an elongation control file (INI) into an :class:`ElongateConfig`.

    The file has a single ``[OPTIONS]`` section. Required keys: ``pdb_file``,
    ``ribosome``, ``L0``. All other keys are optional and fall back to the
    defaults in :class:`ElongationParams` / :class:`ElongateConfig`:

    - ``pdb_file`` -- full native PDB of the target protein (the nascent chain).
    - ``ribosome`` -- truncated CG ribosome PDB (source of the P-/A-anchors; rigid
      scenery in v2).
    - ``L0`` -- starting nascent-chain length (cold-start layout).
    - ``L_max`` -- final length (blank -> full residue count).
    - ``outdir`` -- root output directory (per-length subfolders ``L_<L>/``).
    - ``domain_def`` -- domain YAML for contact ``n_scale`` (one-time precompute).
    - ``stride_output_file`` -- precomputed STRIDE (else STRIDE runs once if on PATH).
    - ``n_steps`` -- integration steps per residue (constant schedule).
    - ``dt`` -- time step (ps); ``ref_t`` -- temperature (K); ``tau_t`` -- Langevin
      friction (1/ps); ``nstout`` -- trajectory/log/checkpoint write frequency.
    - ``device`` -- 'CPU' or 'GPU'; ``ppn`` -- CPU threads (device = CPU).
    - ``constraints`` -- 'None' (flexible, default) or 'AllBonds' (rigid).
    - ``restraint_k`` -- C-terminus position-restraint constant (kJ/mol/nm^2).
    - ``buffer`` -- new-residue placement buffer beyond the A-anchor (nm).
    - ``minimize`` -- yes/no, per-step energy minimization.
    - ``rigid_ribosome`` -- yes/no, v2: append the rigid ribosome + its forces.
    - ``trna_tether`` -- v2: tether the C-terminus to the P-site tRNA R bead (bond +
      CA-CA-tRNA orienting angle, O'Brien-style) vs. a position restraint (default yes).
    - ``tunnel_wall`` -- v2: one-sided planar wall keeping the chain at
      ``x >= tunnel_wall_x0`` so it only extrudes forward (default yes).
    - ``tunnel_wall_x0`` -- wall plane (nm; default ~1.05 = the C-terminal-AA
      addition plane / PTC ~ P-anchor x + tether bond length).
    - ``tunnel_wall_k`` -- wall force constant (kJ/mol/nm^2; default 8368 = 20 kcal/mol/A^2).
    - ``ptc_offset`` -- hold/seed the C-terminus this far (+x) from the P-anchor bead
      (nm); blank -> auto (0 in v1; tether bond length or 0.4 nm in v2).
    - ``nascent_only_output`` -- v2: write only the nascent chain to the
      trajectory/PSF/final (default yes; the rigid ribosome is static).
    - ``post_elongation`` -- 'ejection' (release the tether) or 'stallation' (keep it);
      runs only if ``post_elongation_steps > 0``.
    - ``post_elongation_steps`` -- steps for the post-elongation phase (0 = skip).

    Inline comments starting with ``#`` or ``;`` are ignored; underscores in
    ``n_steps`` (e.g. ``1_000``) are allowed. Paths are resolved relative to the
    current working directory (as for ``md.ini``).

    **Units:** OpenMM defaults throughout -- length nm, time ps, energy kJ/mol,
    temperature K, angle rad, force constants kJ/mol/nm^2.
    """
    def log(msg: str) -> None:
        if verbose:
            print(msg)

    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read elongation config file: {config_file!r}")
    if "OPTIONS" not in cp:
        raise ValueError(f"{config_file}: missing required [OPTIONS] section.")
    o = cp["OPTIONS"]

    def opt(key: str) -> Optional[str]:
        """Return a stripped option value, or None if absent/blank."""
        v = o.get(key, None)
        if v is None:
            return None
        v = v.strip()
        return v if v != "" else None

    def req(key: str) -> str:
        v = opt(key)
        if v is None:
            raise ValueError(f"{config_file}: required option '{key}' is missing or blank.")
        return v

    log(f"Reading elongation parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    ribosome = req("ribosome")
    L0 = int(req("L0"))
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None
    outdir = opt("outdir") or "synth_out"
    domain_def = opt("domain_def")
    stride_output_file = opt("stride_output_file")

    # Per-length run parameters: start from defaults, override what is present.
    p = ElongationParams()
    if opt("n_steps") is not None:
        p.n_steps = int(str(opt("n_steps")).replace("_", ""))
    if opt("dt") is not None:
        p.dt_ps = float(opt("dt"))
    if opt("ref_t") is not None:
        p.ref_t = float(opt("ref_t"))
    if opt("tau_t") is not None:
        p.tau_t = float(opt("tau_t"))
    if opt("nstout") is not None:
        p.nstout = int(opt("nstout"))
    if opt("device") is not None:
        p.device = opt("device")
    if opt("ppn") is not None:
        p.ppn = int(opt("ppn"))
    # 'None' (case-insensitive) / blank -> flexible bonds (the runner default).
    cons = opt("constraints")
    p.constraints = None if (cons is None or cons.lower() == "none") else cons
    if opt("restraint_k") is not None:
        p.restraint_k = float(opt("restraint_k"))
    if opt("buffer") is not None:
        p.buffer_nm = float(opt("buffer"))
    if opt("minimize") is not None:
        p.minimize = bool(strtobool(opt("minimize")))
    if opt("rigid_ribosome") is not None:
        p.rigid_ribosome = bool(strtobool(opt("rigid_ribosome")))
    if opt("trna_tether") is not None:
        p.trna_tether = bool(strtobool(opt("trna_tether")))
    if opt("tunnel_wall") is not None:
        p.tunnel_wall = bool(strtobool(opt("tunnel_wall")))
    if opt("tunnel_wall_x0") is not None:
        p.tunnel_wall_x0_nm = float(opt("tunnel_wall_x0"))
    if opt("tunnel_wall_k") is not None:
        p.tunnel_wall_k = float(opt("tunnel_wall_k"))
    if opt("ptc_offset") is not None:
        p.ptc_offset_nm = float(opt("ptc_offset"))
    if opt("nascent_only_output") is not None:
        p.nascent_only_output = bool(strtobool(opt("nascent_only_output")))
    if opt("post_elongation") is not None:
        p.post_elongation = opt("post_elongation")
    if opt("post_elongation_steps") is not None:
        p.post_elongation_steps = int(str(opt("post_elongation_steps")).replace("_", ""))

    log(f"  inputs: pdb_file={pdb_file}, ribosome={ribosome}")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}, "
        f"n_steps={p.n_steps}")
    log(f"  ribosome forces (v2): {'on (rigid)' if p.rigid_ribosome else 'off (v1)'}"
        + (f"; tether: {'tRNA bond+angle' if p.trna_tether else 'position restraint'}"
           f"; tunnel wall: {('x>=%.2f nm' % p.tunnel_wall_x0_nm) if p.tunnel_wall else 'off'}"
           f"; output: {'nascent-only' if p.nascent_only_output else 'full system'}"
           if p.rigid_ribosome else ""))
    log(f"  mechanics: constraints={'None (flexible)' if p.constraints is None else p.constraints}, "
        f"restraint_k={p.restraint_k} kJ/mol/nm^2, buffer={p.buffer_nm} nm, minimize={p.minimize}")
    if p.post_elongation_steps > 0:
        _pe = p.post_elongation.strip().lower()
        log(f"  post-elongation: {_pe} "
            f"({'release tether' if _pe == 'ejection' else 'keep tether'}) "
            f"for {p.post_elongation_steps} steps -> {_pe}/")
    log(f"  integrator: dt={p.dt_ps} ps, ref_t={p.ref_t} K, tau_t={p.tau_t} /ps, nstout={p.nstout}")
    log(f"  hardware/output: device={p.device}, ppn={p.ppn}, outdir={outdir}")

    return ElongateConfig(pdb_file=pdb_file, ribosome=ribosome, L0=L0, L_max=L_max,
                          outdir=outdir, domain_def=domain_def,
                          stride_output_file=stride_output_file, params=p,
                          config_file=config_file)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def elongate(argv: Optional[List[str]] = None) -> None:
    """Console entry point: ``topo-elongate -f elongate.ini``.

    The simulation is controlled by an INI file (see :func:`read_elongate_config`).
    ``-o`` / ``--device`` are optional overrides handy for sweeps; everything else
    lives in the control file.
    """
    parser = argparse.ArgumentParser(
        prog="topo-elongate",
        description="Protein synthesis elongation loop (build step v1: "
                    "nascent chain only, no ribosome forces). Grows the nascent "
                    "chain N->C one residue per step, restraining the current "
                    "C-terminus to the ribosome P-anchor. Controlled by an INI "
                    "file: topo-elongate -f elongate.ini",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input", "-f", dest="config", type=str,
                        help="elongation control file (INI, [OPTIONS] section).")
    parser.add_argument("-o", "--outdir", default=None,
                        help="override the output directory from the config file.")
    parser.add_argument("--device", default=None, choices=["CPU", "GPU"],
                        help="override the compute device from the config file.")

    # A bare `topo-elongate` (no arguments) prints help.
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(argv)
    if not args.config:
        parser.error("an elongation control file is required: -f elongate.ini")

    print(f"OpenMM version: {mm.__version__}")

    cfg = read_elongate_config(args.config)
    if args.outdir:
        cfg.outdir = args.outdir
    if args.device:
        cfg.params.device = args.device

    run_elongation(
        cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
        domain_def=cfg.domain_def, stride_output_file=cfg.stride_output_file,
        params=cfg.params)


if __name__ == "__main__":
    elongate()
