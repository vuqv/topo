"""Core nascent-chain MD machinery (``topo.csp.core``).

This module holds the **per-length / per-stage molecular-dynamics building blocks**
for co-translational synthesis — the foundation that the higher-level runners build
on. It is a *library*, not a runner: it exposes no CLI and no outer loop.

Its consumers are:

* :mod:`topo.csp.protocol` — the O'Brien Continuous Synthesis Protocol (CSP), which
  calls :func:`run_length` three times per residue (one per kinetic sub-stage), and
* the Tutorial-9 cylinder runner (`tutorials/09_translation_cylinder/cylinder.py`),
  which reuses these blocks with the ribosome modelled as an analytic cylinder.

(The standalone fixed-rate elongation runner — ``run_elongation`` / the
``topo-elongate`` CLI — and its Tutorial 7 were removed: a fixed per-residue step
count is not a physically meaningful synthesis model. Use CSP (``topo-csp``) for
codon-resolved kinetics.)

The central primitive is :func:`run_length`, which builds, seeds, restrains, runs
and finalizes **one length-``L`` segment** of nascent-chain MD:

1. **Build** the length-``L`` TOPO model on native residues ``1..L`` (bonds,
   angles, torsions, Yukawa, contacts), injecting the precomputed ``L x L`` contact
   subset instead of recomputing (``buildCoarseGrainModel(precomputed_contacts=...)``).
2. **Seed coordinates.** ``L == L0``: lay residues ``1..L0`` extended along the
   tunnel axis (+x) from the P-anchor, one CG bond length apart. ``L > L0``:
   residues ``1..L-1`` from the previous segment's final structure; the new residue
   ``L`` at the A-anchor + buffer.
3. **Restrain only residue ``L``** (the current C-terminus) to a chosen anchor with
   a harmonic ``CustomExternalForce`` (``k = restraint_k``).
4. **Minimize**, draw Boltzmann velocities at ``ref_t``, **run the requested step
   count** under the per-stage *stability guard* (see :data:`STABILITY_POTE_LIMIT_KJ`),
   write the per-segment outputs, and return the final NASCENT coordinates to seed
   the next segment.

Contacts are precomputed **once** on the full native PDB (:func:`precompute_contacts`);
STRIDE / heavy-atom analysis are never re-run per length. The optional rigid ribosome
scenery and tunnel wall are wired up via :mod:`topo.csp.ribosome`. Build / setup /
finalize are reused from :mod:`topo.engine`.
"""
from __future__ import annotations

import math
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
from topo.csp.ribosome import (Ribosome, load_ribosome, append_ribosome,
                                       add_trna_tether, add_tunnel_wall,
                                       TRNA_TETHER_BOND_NM, TUNNEL_WALL_X0_NM,
                                       TUNNEL_WALL_K)
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

# --- per-stage stability guard --------------------------------------------
# The seeded structure is integrated at a 15 fs timestep with flexible (un-
# constrained) bonds. For a few configurations -- when a newly formed native
# contact introduces a stiff Go well -- 15 fs is too large and the dynamics diverge
# (potential energy -> ~1e13 kJ/mol), corrupting that stage's trajectory frames
# (see tutorials/12_auto/WHY_10_FAILS.md). O'Brien's reference v6 avoids
# this with rigid AllBonds constraints; topo keeps flexible bonds (needed to seed
# the far-placed new residue at the A-site) and instead detects a diverging stage
# and re-runs it with a halved timestep and proportionally more steps -- which
# preserves the physical in-vivo dwell time (dwell = n_steps * dt) exactly while
# stabilising the integration. Up to STABILITY_MAX_ATTEMPTS halvings are tried
# (dt -> dt/2^k). A finite, physically sane CG stage has |PotE| of order 1e1-1e4
# kJ/mol, so the 1e9 limit cleanly separates "diverged" from "fine" with margin.
STABILITY_POTE_LIMIT_KJ = 1.0e9
STABILITY_MAX_ATTEMPTS = 6


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

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein.
    L : int
        Number of N-terminal residues to keep (the current nascent length).
    out_pdb : str
        Path of the CA-only length-``L`` PDB to write.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If ``full_pdb`` has fewer than ``L`` CA atoms.
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

    ``k`` is a **per-particle** parameter (not a global) so this force can coexist
    with other ``CustomExternalForce``\\s that also call their global constant ``k``
    -- e.g. the tunnel wall (:func:`topo.csp.ribosome.add_tunnel_wall`). Two
    forces sharing a global parameter name with different default values is an
    OpenMM error; per-particle ``k`` avoids the clash (this combination -- position
    restraint + tunnel wall in v2 -- is what :mod:`topo.csp` uses).

    Parameters
    ----------
    system : openmm.System
        The System to add the restraint force to (modified in place).
    particle_index : int
        Index of the restrained particle (the current C-terminus, ``L-1``).
    anchor_nm : numpy.ndarray
        The fixed anchor position ``(x0, y0, z0)`` in nm (the P-anchor target).
    k : float
        Harmonic force constant in OpenMM units (kJ/mol/nm^2).

    Returns
    -------
    openmm.Force
        The added ``CustomExternalForce`` restraint.
    """
    restraint = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    for p in ("k", "x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    restraint.addParticle(int(particle_index),
                          [float(k), float(anchor_nm[0]), float(anchor_nm[1]),
                           float(anchor_nm[2])])
    system.addForce(restraint)
    return restraint


# --------------------------------------------------------------------------
# PDB writing helper
# --------------------------------------------------------------------------
def _write_pdb(topology, positions_nm: np.ndarray, path: str) -> None:
    """Write a PDB from a topology and an ``(N, 3)`` nm coordinate array.

    Parameters
    ----------
    topology : openmm.app.Topology
        The topology describing the atoms to write.
    positions_nm : numpy.ndarray
        An ``(N, 3)`` array of coordinates in nm.
    path : str
        Output PDB path.

    Returns
    -------
    None
    """
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
        """Open the output DCD and store the reporting settings.

        Parameters
        ----------
        file : str
            Path of the DCD trajectory to write.
        reportInterval : int
            Number of integration steps between written frames.
        nascent_topology : openmm.app.Topology
            The ``n_keep``-atom (nascent-only) topology; sets the DCD atom count.
        n_keep : int
            Number of leading atoms (the nascent chain) to write each frame.
        append : bool, optional
            If True, open the file for appending instead of truncating.
        """
        self._reportInterval = reportInterval
        self._topology = nascent_topology
        self._n = n_keep
        self._append = append
        self._out = open(file, "r+b" if append else "wb")
        self._dcd = None

    def describeNextReport(self, simulation):
        """Tell the simulation when and what this reporter needs next.

        Parameters
        ----------
        simulation : openmm.app.Simulation
            The running simulation (queried for the current step).

        Returns
        -------
        dict
            ``{'steps', 'periodic', 'include'}`` -- steps until the next report,
            no PBC, positions only.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        # No PBC in these runs; positions only.
        return {"steps": steps, "periodic": False, "include": ["positions"]}

    def report(self, simulation, state):
        """Write one frame of the first ``n_keep`` atoms to the DCD.

        Lazily creates the underlying :class:`openmm.app.DCDFile` on the first
        call, then writes the sliced (nascent-only) positions.

        Parameters
        ----------
        simulation : openmm.app.Simulation
            The running simulation (supplies the integrator / step info).
        state : openmm.State
            The current state (source of the positions).

        Returns
        -------
        None
        """
        if self._dcd is None:
            self._dcd = mm.app.DCDFile(
                self._out, self._topology, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, self._append)
        positions = state.getPositions(asNumpy=True)[:self._n]
        self._dcd.writeModel(positions)

    def __del__(self):
        """Close the output file handle when the reporter is garbage-collected."""
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

    Parameters
    ----------
    cgModel : topo.core.system.system
        The built (nascent + ribosome) model providing ``.topology`` and
        ``.system``.
    path : str
        Output PSF path.

    Returns
    -------
    None
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

    Parameters
    ----------
    cfg : topo.SimulationConfig
        The per-length config (supplies output paths and the runinfo path).
    ctx : topo.engine.SimulationContext
        The active simulation context (simulation, checkpoint, runinfo paths).
    nascent_topology : openmm.app.Topology
        The ``n_keep``-atom nascent-only topology for the written ``_final.pdb``.
    n_keep : int
        Number of leading (nascent) atoms to keep in the final structure.
    start_epoch : float
        Wall-clock start time (``time.time()``) for the elapsed-time report.

    Returns
    -------
    None
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
    # NOTE: whether to append the truncated ribosome as rigid (mass-0) scenery is no
    # longer a flag here -- the CSP runner always loads the supplied ribosome PDB as
    # rigid scenery (passed to run_length via its `ribo` argument). run_length itself
    # keys off that `ribo` argument, not a field.
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
    # exit) and cannot fold back past the synthesis point into the truncated-ribosome
    # void below the PTC. Applied throughout synthesis + post-phase. Neither the plane
    # nor the stiffness is a user/INI knob: x0_nm is **auto-derived from the ribosome
    # structure** by the CSP runner (the lower / P-site C-terminus hold plane =
    # min(P,A anchor x) + ptc_offset) so it can never go stale when the structure
    # changes (left None here; the runner fills it in before run_length), and
    # tunnel_wall_k is a **fixed model constant** (O'Brien's 20 kcal/mol/A^2). Only the
    # `tunnel_wall` on/off toggle is exposed.
    tunnel_wall: bool = True
    tunnel_wall_x0_nm: Optional[float] = None
    tunnel_wall_k: float = TUNNEL_WALL_K
    # Note: when a ribosome is present the output is **always nascent-only** -- only the
    # nascent chain is written to the trajectory / PSF / final structure (the rigid
    # ribosome is static, so writing it every frame would waste storage; the checkpoint
    # still holds the full system, and the movie overlays the ribosome separately). This
    # is the default CSP behavior and is no longer a flag (keyed off `ribo is not None`).
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

    Parameters
    ----------
    out_dir : pathlib.Path
        Output directory for this length's run.
    sub_pdb : str
        Length-``L`` native CA PDB (the model's ``pdb_file``).
    seed_pdb : str or None
        Seed-coordinate PDB fed via ``init_position`` (v1); ``None`` in v2 where
        seed coordinates are set directly on ``built.positions``.
    params : ElongationParams
        Shared per-length run parameters (steps, dt, temperature, device, ...).

    Returns
    -------
    topo.SimulationConfig
        The configured per-length simulation config.
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
    (:func:`topo.csp.ribosome.append_ribosome`).

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

    # With a ribosome present, output is always nascent-only: write only the nascent
    # chain to the trajectory / PSF / final structure. Capture a nascent (L-atom)
    # topology BEFORE appending the ribosome (append mutates cgModel.topology), and
    # write the nascent PSF now (dumpTopology keys off its nascent-only per-atom lists).
    nascent_only = ribo is not None
    nascent_topology = None
    if nascent_only:
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
        if params.tunnel_wall_x0_nm is None:
            raise ValueError(
                "tunnel_wall is on but tunnel_wall_x0_nm is unset -- it is auto-derived "
                "from the ribosome structure by the CSP runner (run_continuous_synthesis); "
                "set params.tunnel_wall_x0_nm before calling run_length directly.")
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
        if not nascent_only:
            # Full-system PSF (nascent-only PSF was already written above).
            _dump_topology_psf(cgModel, cfg.output_path(".psf"))
        built_positions = ([mm.Vec3(*r) for r in positions]) * unit.nanometer

    # Post-synthesis phases run for their own step count.
    if n_steps_override is not None:
        cfg.md_steps = n_steps_override

    built = engine.BuiltSystem(cgModel=cgModel, system=cgModel.system,
                               topology=cgModel.topology, positions=built_positions)

    start = time.time()
    # Stability-guarded stage run. A few configurations are unstable at the
    # configured 15 fs timestep with flexible bonds: when a new native contact
    # forms a stiff Go well, the dynamics diverge (PotE -> ~1e13) and corrupt the
    # stage's frames (OBSERVATIONS.md #1). O'Brien's reference v6 sidesteps this
    # with rigid AllBonds constraints; topo keeps flexible bonds (needed to seed
    # the far-placed new residue), so instead -- *only if a stage diverges* -- we
    # re-run it with a halved timestep and proportionally more steps. That keeps
    # the physical in-vivo dwell time identical (dwell = n_steps * dt) while making
    # the integration stable. The common case runs once at the configured dt.
    # Divergence is judged from the **maximum** |PotE| over the whole stage (a
    # mid-run blow-up can cool back below the limit by the final frame yet still
    # have ruined those frames), and the chunked stepping aborts a diverging stage
    # early. Each retry's fresh setup_simulation + attach_reporters truncates the
    # per-stage output, so a successful attempt cleanly overwrites the aborted one.
    base_dt = cfg.dt
    base_steps = cfg.md_steps
    ctx = None
    max_pe = float("nan")
    for attempt in range(STABILITY_MAX_ATTEMPTS):
        cfg.dt = base_dt / (2 ** attempt)
        cfg.md_steps = base_steps * (2 ** attempt)
        if attempt > 0:
            print(f"[stability] stage diverged (max|PotE| = {max_pe:.3g} kJ/mol > "
                  f"{STABILITY_POTE_LIMIT_KJ:g}); re-running with dt = "
                  f"{cfg.dt.value_in_unit(unit.picoseconds):g} ps x {cfg.md_steps} "
                  f"steps (identical dwell time; attempt {attempt + 1}/"
                  f"{STABILITY_MAX_ATTEMPTS}).")
        ctx = engine.setup_simulation(cfg, built)
        if params.minimize:
            if attempt == 0:
                print("Minimizing seeded structure (relax placement / new bond)...")
            ctx.simulation.minimizeEnergy()
            # Re-draw velocities for the relaxed structure.
            ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)

        engine.attach_reporters(cfg, ctx.simulation, suffix="", total_steps=cfg.md_steps)
        if nascent_only:
            # Swap the full-system DCD reporter for a nascent-only one (the rigid
            # ribosome is static -- no need to write its ~thousands of beads/frame).
            ctx.simulation.reporters[1] = NascentDCDReporter(
                cfg.output_path(".dcd"), cfg.nstxout, nascent_topology, L)

        # Step in chunks so a divergence is caught (and the stage aborted) mid-run.
        max_pe = 0.0
        chunk = max(cfg.nstxout, cfg.md_steps // 20, 1)
        done = 0
        diverged = False
        while done < cfg.md_steps:
            n = min(chunk, cfg.md_steps - done)
            ctx.simulation.step(n)
            done += n
            pe = abs(ctx.simulation.context.getState(getEnergy=True).getPotentialEnergy(
                ).value_in_unit(unit.kilojoule_per_mole))
            max_pe = max(max_pe, pe)
            if not math.isfinite(pe) or pe > STABILITY_POTE_LIMIT_KJ:
                diverged = True
                break
        if not diverged:
            break
    else:
        print(f"[stability][warning] stage '{out_subdir or ('L_%03d' % L)}' still "
              f"diverges after {STABILITY_MAX_ATTEMPTS} attempts (max|PotE| = "
              f"{max_pe:.3g} kJ/mol). Continuing, but this stage's frames are suspect.")
    # Restore the configured values (cfg is per-call, but keep it tidy for finalize).
    cfg.dt = base_dt
    cfg.md_steps = base_steps

    if nascent_only:
        _finalize_nascent(cfg, ctx, nascent_topology, L, start)
    else:
        engine.finalize_simulation(cfg, ctx, built.topology, start)

    # 6. final NASCENT coordinates seed the next length ----------------------
    final_pdb = cfg.output_path("_final.pdb")
    final = mm.app.PDBFile(final_pdb).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]

