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
import os
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
                                       TUNNEL_WALL_K, RIBO_NC_EPS_KJ)
from topo.utils.nonbonded import build_nonbonded_interaction

# --- console verbosity -----------------------------------------------------
# Continuous synthesis runs ~3 short MD stages per residue, so the per-stage
# banners (build block, "chains=...", "Running simulation...", "[tracking]...",
# "Minimizing...", "Wrote last...", "Finished in...") add up to tens of lines per
# amino acid. By default we print one concise summary line per stage instead.
# Set TOPO_CSP_VERBOSE=1 to restore the full per-stage detail (e.g. when debugging
# a single length). Stability-guard warnings are always printed regardless.
VERBOSE = os.environ.get("TOPO_CSP_VERBOSE", "").strip().lower() in ("1", "true", "yes", "on")


def _vprint(*args, **kwargs) -> None:
    """print() only when TOPO_CSP_VERBOSE is set (see :data:`VERBOSE`)."""
    if VERBOSE:
        print(*args, **kwargs)


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
# Optimal PTC restraint-target geometry (step 2 -- equilibrium-bond seeding)
# --------------------------------------------------------------------------
# All quantities here are in OpenMM default units: length nm, energy kJ/mol, angle rad.
# O'Brien resting-geometry restraint parameters (continuous_synthesis_v6.py), in nm:
#   aminoacyl/peptidyl-tRNA bond lengths + the equilibrium peptide bond.
_PTC_D_A_NM, _PTC_D_P_NM = 0.427, 0.476        # tRNA bond lengths (nm; O'Brien 4.27/4.76 A)
_PTC_PEPTIDE_NM = 0.381                          # equilibrium peptide bond (nm; 3.81 A)
# Orienting-angle / improper restraint force constant: O'Brien 25 kcal/mol/rad^2 in kJ/mol/rad^2.
_PTC_ANGLE_K_KJ = 25.0 * 4.184                   # = 104.6 kJ/mol/rad^2
# (The excluded-volume well depth is topo.csp.ribosome.RIBO_NC_EPS_KJ, kJ/mol.)


def _ptc_bead_index(ribo, segid: str, resid: int, name: str) -> int:
    """Index of one named ribosome bead in a :class:`Ribosome` (or raise)."""
    for i in range(ribo.n):
        if ribo.segids[i] == segid and ribo.resids[i] == resid and ribo.names[i] == name:
            return i
    raise ValueError(f"ribosome bead {segid}:{resid}@{name} not found.")


def optimal_ptc_targets(ribo, *, aa_rmin_2_nm: float = 0.5,
                        n_starts: int = 60
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Optimal A-site / P-site C-terminus restraint **target points** (nm).

    Returns ``(a_target, p_target)`` -- two **fixed points in space** (nm; NOT bonds;
    the CSP path restrains the C-terminus to a fixed point via
    :func:`add_cterm_restraint`) that sit exactly one peptide bond
    (``_PTC_PEPTIDE_NM`` = 0.381 nm) apart and clear of the ribosome excluded volume.
    Seeding the new residue at ``a_target`` while the previous residue rests at
    ``p_target`` makes the always-present peptide bond start at its equilibrium
    length, so a rigid ``AllBonds`` constraint seeds/minimizes cleanly at 15 fs and
    the dt-halving stability guard is not triggered (tutorials/14 step 2).

    The points are found (all in OpenMM units: nm, kJ/mol, rad) by minimizing the soft
    O'Brien restraint energy -- the **A/P tRNA bonds** (``_PTC_D_A_NM`` / ``_PTC_D_P_NM``,
    harmonic ``kb = RESTRAINT_K_KJ``) and orienting **angles/improper**
    (``_PTC_ANGLE_K_KJ``; ``continuous_synthesis_v6.py``) -- plus the ribosome
    ``eps*(sigma/r)^12`` excluded volume (``eps = RIBO_NC_EPS_KJ``,
    ``sigma_ij = 1/2 (sigma_AA + sigma_bead)``), subject to:

    - the **peptide bond** as the only **hard** distance constraint
      (``_PTC_PEPTIDE_NM`` = 0.381 nm -- it is topo's AllBonds constraint length, so it
      must hold exactly);
    - an **exit-side inequality** keeping each point at ``x >= x`` of its tRNA R bead
      (between the tRNA and the +x exit port, not buried in the ribosome); and
    - O'Brien's two excluded-volume exclusions (new AA<->AtR:76@R, prev AA<->PtR:76@R),

    over ``n_starts`` deterministic full-sphere starts. The tRNA bond lengths are soft
    (finite-k, as in O'Brien's model), so ``|A-RA|`` / ``|P-RP|`` come out *near*
    0.427 / 0.476 nm rather than exactly -- and the solve stays feasible even on
    geometries where those distances cannot be met exactly together with the peptide bond.

    Parameters
    ----------
    ribo : Ribosome
        The parsed rigid CG ribosome (supplies the AtR/PtR 76 R/P/BR2 beads, all bead
        coordinates and the excluded-volume radii).
    aa_rmin_2_nm : float, optional
        Nascent-bead Rmin/2 (nm) used in the seed excluded-volume term, combined with the
        ribosome per-bead Rmin/2 by O'Brien's SUM rule (R = aa_rmin_2_nm + Rmin/2_ribo) --
        the SAME EV the simulation applies. Default 0.5 (a conservative value ~ the largest
        per-residue Karanicolas-Brooks radius, so the seed clears the wall for essentially
        every residue -> fewer dt-halving blow-ups).
    n_starts : int, optional
        Number of deterministic multistart initial orientations (default 60).

    Returns
    -------
    (numpy.ndarray, numpy.ndarray)
        ``a_target`` and ``p_target``, each a ``(3,)`` coordinate in nm.
    """
    from scipy.optimize import minimize  # lazy: only needed for this geometry mode

    RB = ribo.coords_nm                                     # all bead coords (nm)
    RBr = np.asarray(ribo.Rmin_2_nm)                        # bead Rmin/2 (nm)
    iRA = _ptc_bead_index(ribo, "AtR", 76, "R"); RA = RB[iRA]
    iRP = _ptc_bead_index(ribo, "PtR", 76, "R"); RP = RB[iRP]
    PA = RB[_ptc_bead_index(ribo, "AtR", 76, "P")]
    PP = RB[_ptc_bead_index(ribo, "PtR", 76, "P")]
    U2A = RB[_ptc_bead_index(ribo, "AtR", 76, "BR2")]       # topo BR2 == O'Brien PU2
    U2P = RB[_ptc_bead_index(ribo, "PtR", 76, "BR2")]

    # Pair contact distance for the seed EV, CONSISTENT with the simulation's NC<->ribosome
    # force (topo.csp.ribosome.append_ribosome): O'Brien's SUM rule R_ij = Rmin/2_nascent +
    # Rmin/2_ribo (NOT the average 0.5*(sig_i+sig_j)). aa_rmin_2_nm is the nascent bead's Rmin/2
    # (a conservative value ~ the max per-residue Karanicolas-Brooks radius, so the seed is
    # placed clear of the wall for essentially every residue -> fewer dt-halving blow-ups).
    R_pair = aa_rmin_2_nm + RBr                              # pair contact dist (nm), sum rule
    mA = np.ones(ribo.n, bool); mA[iRA] = False             # O'Brien exclusions
    mP = np.ones(ribo.n, bool); mP[iRP] = False
    ka, eps = _PTC_ANGLE_K_KJ, RIBO_NC_EPS_KJ               # kJ/mol/rad^2 , kJ/mol
    kb = RESTRAINT_K_KJ                                      # tRNA bond k (kJ/mol/nm^2; O'Brien 200 kcal/mol/A^2)

    def _ang(p, q, r):
        """Bond angle p-q-r (rad), the angle at vertex ``q``.

        Parameters
        ----------
        p, q, r : numpy.ndarray
            ``(3,)`` coordinates; ``q`` is the central (vertex) point.

        Returns
        -------
        float
            Angle in radians, in ``[0, pi]``.
        """
        u = p - q; v = r - q
        return np.arccos(np.clip(u.dot(v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))

    def _dih(p, q, r, s):
        """Dihedral (torsion) angle p-q-r-s (rad).

        Parameters
        ----------
        p, q, r, s : numpy.ndarray
            The four ``(3,)`` coordinates defining the torsion.

        Returns
        -------
        float
            Signed dihedral angle in radians, in ``(-pi, pi]``.
        """
        b1, b2, b3 = q - p, r - q, s - r
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        m = np.cross(n1, b2 / np.linalg.norm(b2))
        return np.arctan2(m.dot(n2), n1.dot(n2))

    def _ev(pt, mask):                                      # excluded volume (kJ/mol)
        """O'Brien 12-10-6 excluded-volume energy of one seed point vs. the ribosome.

        Parameters
        ----------
        pt : numpy.ndarray
            ``(3,)`` seed coordinate (nm) whose burial is scored.
        mask : numpy.ndarray
            Boolean mask over ribosome beads to sum over (the point's own tRNA R
            bead is excluded, per O'Brien).

        Returns
        -------
        float
            Summed excluded-volume energy (kJ/mol); minimizing it places the seed
            at the least-buried point of the actual wall.
        """
        # O'Brien 12-10-6 (same form as the simulation NC<->ribosome force), summed over
        # ribosome beads. Minimizing it places the seed at the least-buried point of the
        # actual wall (the -18/+4 terms give a negligible eps-deep attractive tail).
        d = np.linalg.norm(RB[mask] - pt, axis=1)
        x = R_pair[mask] / d
        return np.sum(eps * (13.0 * x ** 12 - 18.0 * x ** 10 + 4.0 * x ** 6))

    d2r = np.deg2rad
    # Only the peptide bond is a HARD equality constraint -- it is topo's AllBonds
    # constraint length, so a_target/p_target MUST be exactly 0.381 nm apart or the
    # rigid bond is seeded stretched. The two tRNA bond lengths are SOFT harmonic
    # penalties in the objective (O'Brien's bonds are finite-k harmonic, not rigid):
    # they pull the pair toward the PTC but give a little, so the result stays feasible
    # even on geometries where 0.427/0.476 cannot be met exactly. Inequality
    # constraints keep each point on the EXIT side of its tRNA R bead (x >= x_tRNA),
    # i.e. between the tRNA and the +x exit port -- never buried back in the ribosome
    # (the working ribosome is +x-aligned; FILES.md / tutorials/14). Without this a
    # feasible-but-buried minimum (A behind the A-tRNA) can win on excluded volume alone.
    cons = [{"type": "eq", "fun": lambda x: np.linalg.norm(x[:3] - x[3:]) - _PTC_PEPTIDE_NM},
            {"type": "ineq", "fun": lambda x: x[0] - RA[0]},   # A.x >= A-tRNA.x (exit side)
            {"type": "ineq", "fun": lambda x: x[3] - RP[0]}]   # P.x >= P-tRNA.x (exit side)

    def _obj(x):                                            # kJ/mol
        """Total restraint energy to minimize for the A/P seed geometry.

        Sums the soft A-tRNA and P-tRNA bond penalties, the orienting
        angle/dihedral penalties about each tRNA, and the ribosome excluded
        volume of both seed points.

        Parameters
        ----------
        x : numpy.ndarray
            ``(6,)`` state vector: the A-site point ``x[:3]`` and the P-site
            point ``x[3:]`` (nm).

        Returns
        -------
        float
            Total energy (kJ/mol).
        """
        A, P = x[:3], x[3:]
        E = kb * (np.linalg.norm(A - RA) - _PTC_D_A_NM) ** 2     # soft A-tRNA bond
        E += kb * (np.linalg.norm(P - RP) - _PTC_D_P_NM) ** 2    # soft P-tRNA bond
        E += ka * (_ang(A, RA, PA) - d2r(106)) ** 2 + ka * (_ang(A, RA, U2A) - d2r(127)) ** 2
        E += ka * (_ang(P, RP, PP) - d2r(117)) ** 2 + ka * (_ang(P, RP, U2P) - d2r(130)) ** 2
        E += ka * (_dih(A, RA, PA, U2A) - d2r(128)) ** 2 + ka * (_dih(P, RP, PP, U2P) - d2r(-161)) ** 2
        return E + _ev(A, mA) + _ev(P, mP)

    # Deterministic full-sphere (Fibonacci) multistart over the new-residue offset
    # direction: covers out-of-plane orientations so the global minimum is found
    # robustly and unit-invariantly (a planar ring of starts misses it and makes the
    # result depend on the variable scale). The previous residue is seeded toward the
    # P-tRNA. Keep the lowest-objective feasible solution.
    best = None
    golden = np.pi * (1.0 + 5.0 ** 0.5)
    for t in range(n_starts):
        z = 1.0 - 2.0 * (t + 0.5) / n_starts
        rho = np.sqrt(max(0.0, 1.0 - z * z))
        dirn = np.array([rho * np.cos(golden * t), rho * np.sin(golden * t), z])
        sA = RA + _PTC_D_A_NM * dirn
        sP = sA + _PTC_PEPTIDE_NM * (RP - sA) / np.linalg.norm(RP - sA)
        r = minimize(_obj, np.concatenate([sA, sP]), method="SLSQP",
                     constraints=cons, options={"maxiter": 300, "ftol": 1e-12})
        if r.success and (best is None or r.fun < best.fun):
            best = r
    if best is None:
        raise RuntimeError("optimal_ptc_targets: no feasible restraint geometry found.")
    return best.x[:3], best.x[3:]                           # nm (already)


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
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run TOPO's contact builder once on the full native PDB (DESIGN §3.5).

    Returns ``(R_full, eps_full, rmin_2_full)`` -- the ``N_full x N_full`` well-position
    (nm) and well-depth (kJ/mol) matrices, plus the per-residue Karanicolas-Brooks
    collision radius ``Rmin/2`` (nm, shape ``(N_full,)``) that O'Brien uses as the
    nascent excluded-volume radius (the structure-derived ``A_i`` values). STRIDE is
    run at most once here (and cached by :func:`build_nonbonded_interaction`); each
    length later reuses the top-left ``L x L`` block (and ``rmin_2_full[:L]``), so
    neither STRIDE nor the heavy-atom analysis is ever re-run per length.
    """
    print("=" * 66)
    print("[ Precompute contacts (build-once-subset) ]")
    print("=" * 66)
    print(f"Running TOPO contact builder on full native structure: {full_pdb}")
    R_full, eps_full, rmin_2_full = build_nonbonded_interaction(
        full_pdb, domain_def, stride_output_file, return_rmin_2=True)
    print(f"  full contact matrices: {R_full.shape}; "
          f"K-B Rmin/2 range {rmin_2_full.min():.3f}-{rmin_2_full.max():.3f} nm")
    return R_full, eps_full, rmin_2_full


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
                       constraints="AllBonds", model: str = "topo",
                       nascent_rmin_2: Optional[np.ndarray] = None):
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
    nascent_rmin_2 : numpy.ndarray or None
        Per-residue collision radius Rmin/2 (nm), length ``L``, for the nascent chain
        -- the structure-derived Karanicolas-Brooks values (Option A). Used as the
        particle excluded-volume radius (``rf_sigma``). This is deliberately **not**
        taken from ``model_parameters`` (whose fixed per-AA protein Rmin_2 is the rigid
        *ribosome* scenery value, not the mobile chain's). ``rf_sigma`` only feeds
        ``dumpForceFieldData``, so ``None`` is harmless (the radius is left unset).

    Returns
    -------
    topo.core.system.system
        The built model (``.system``, ``.topology``, ``.positions`` ready to use).
    """
    _vprint("=" * 66)
    _vprint(f"[ Length-L build ] {sub_pdb}  (build-once-subset contacts)")
    _vprint("=" * 66)

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
    _vprint(f"  chains={topo_model.n_chains}  CA atoms={n_ca}  "
            f"bonds={topo_model.n_bonds}  angles={topo_model.n_angles}  "
            f"bonds: {'rigid (AllBonds)' if use_constraints else 'flexible (harmonic)'}")

    if use_constraints:
        for bond in topo_model.bonds:
            topo_model.system.addConstraint(bond[0].index, bond[1].index,
                                            topo_model.bonds[bond][0])

    # Per-residue particle properties. Mass/charge are per-AA (from model_parameters).
    # The nascent excluded-volume radius is the per-residue K-B Rmin/2 (structure-derived,
    # Option A), NOT the fixed per-AA model_parameters value (that is the rigid ribosome
    # scenery radius). rf_sigma only feeds dumpForceFieldData, so None is harmless.
    topo_model.setCAMassPerResidueType()
    topo_model.setCAChargePerResidueType()
    if nascent_rmin_2 is not None:
        topo_model.setParticlesRadii(list(np.asarray(nascent_rmin_2, dtype=float)))

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
    topo_model.rmin_matrix = R_L
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
                   buffer_nm: float,
                   seed_point: Optional[np.ndarray] = None) -> np.ndarray:
    """Seed length ``L`` from the previous final structure + the new residue.

    Residues ``1..L-1`` keep their coordinates from step ``L-1``'s final
    structure (``prev_final``, shape ``(L-1, 3)`` nm). The new C-terminal residue
    ``L`` is placed either at ``seed_point`` directly (the equilibrium-bond A-site
    target from :func:`optimal_ptc_targets`, when given -- so the always-present
    peptide bond ``L-1<->L`` starts at its equilibrium length and a rigid ``AllBonds``
    constraint seeds cleanly; tutorials/14 step 2), or at the A-anchor offset by
    ``buffer_nm`` along +x (default; the buffer clears excluded volume so the new bead
    does not get a huge ``(sigma/r)^12`` kick -- DESIGN §2.5 / invariant 6). Returns an
    ``(L, 3)`` array in nm.
    """
    if seed_point is not None:
        new_residue = np.asarray(seed_point, dtype=float)
    else:
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
    _vprint(f"Wrote last nascent conformation to {final_pdb}")
    topo.runinfo.write_run_end(ctx.runinfo_path, simulation=sim,
                               start_epoch=start_epoch, final_structure=final_pdb)
    _vprint("--- Finished in %s seconds ---" % (time.time() - start_epoch))


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
    # Optimize the PTC (peptidyl-transferase-center) restraint/seed geometry
    # (tutorials/14 step 2; opt-in). When True, the CSP runner seeds each new residue at
    # the optimal A-site target point one peptide bond (0.381 nm) from the previous
    # C-terminus (:func:`optimal_ptc_targets`) and uses those optimized A/P points as the
    # C-terminus restraint targets and the tunnel-wall plane, so the always-present
    # peptide bond starts at its equilibrium length -- letting a rigid
    # `constraints='AllBonds'` build seed/minimize cleanly at 15 fs without the
    # dt-halving stability guard firing. Default False keeps the validated far-seed +
    # flexible-bond + dt-halving behavior (Tutorials 12 & 13 unchanged). Enable with
    # `optimize_ptc_geometry = yes` + `constraints = AllBonds` in the INI.
    #
    # NOTE: this depends on the ribosome carrying tRNA beads under the expected names
    # (segids ``AtR``/``PtR``, resid 76, beads ``R``/``P``/``BR2``); :func:`optimal_ptc_targets`
    # reads them directly. A ribosome PDB with no tRNA -- or with differently-named tRNA
    # segments -- makes this (and the plain anchor lookup) fail. See the tRNA-presence
    # TODO in ``review/TODO.md``.
    optimize_ptc_geometry: bool = False
    # Nascent per-residue excluded-volume radius (Rmin/2) for the NC<->ribosome 12-10-6 force:
    #   "kb"     -> per-residue Karanicolas-Brooks collision diameter from the native structure
    #               (Option A; O'Brien's actual nascent A_i radii). DEFAULT.
    #   "per_aa" -> per-amino-acid sidechain radii OBRIEN_SC_RMIN_2_NM (Option B; O'Brien's
    #               ribosomal-protein S<aa1> values). Optional fallback.
    # See tutorials/15_claude_fix/OPTION_A_vs_B.md. The ribosome side always uses O'Brien's
    # Rmin/2 from model_parameters (via load_ribosome), independent of this choice.
    nascent_ev_radii: str = "kb"
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
    # CSP runs one config per stage; silence the engine's per-run platform /
    # metadata banners unless TOPO_CSP_VERBOSE restores full detail.
    cfg.quiet = not VERBOSE
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
               seed_point: Optional[np.ndarray] = None,
               tether_segid: str = "PtR",
               tether_prev_segid: Optional[str] = None,
               nascent_rmin_2: Optional[np.ndarray] = None,
               minimize_override: Optional[bool] = None,
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
    - ``label`` : short stage tag (e.g. ``"stage 1 peptidyl-transfer"``) shown in
      the concise per-stage summary line and, under ``TOPO_CSP_VERBOSE``, in the
      verbose console banner.

    Returns the final **nascent** ``(L, 3)`` nm coordinate array.
    """
    _vprint()
    _vprint("#" * 66)
    _vprint("# " + (f"L={L}  {label}" if label else
                    (f"Nascent length L = {L}"
                     + ("  (+ rigid ribosome)" if ribo else ""))))
    _vprint("#" * 66)

    out_dir = out_root / (out_subdir or f"L_{L:03d}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. length-L native structure (bonded terms + per-residue properties) ---
    sub_pdb = str(out_dir / f"native_1_{L}.pdb")
    write_subset_structure(full_pdb, L, sub_pdb)

    # 2. inject the L x L contact subset (build-once-subset) -----------------
    R_L = np.ascontiguousarray(R_full[:L, :L])
    eps_L = np.ascontiguousarray(eps_full[:L, :L])
    # Nascent per-residue K-B Rmin/2 (Option A), sliced to length L. Used both for the
    # nascent particle excluded-volume radius (in build_length_model) and for the nascent
    # side of the NC<->ribosome EV (append_ribosome below).
    nasc_rm = None if nascent_rmin_2 is None else np.asarray(nascent_rmin_2)[:L]
    cgModel = build_length_model(sub_pdb, R_L, eps_L, constraints=params.constraints,
                                 nascent_rmin_2=nasc_rm)

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
        nascent_pos = seed_positions(prev_final, a_anchor, params.buffer_nm,
                                     seed_point=seed_point)

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
    # nascent_rmin_2 (per-residue Karanicolas-Brooks Rmin/2, full-structure array) gives the
    # nascent side of the NC<->ribosome excluded volume (Option A: structure-derived per-residue
    # radii); slice to the current length. If None, append_ribosome falls back to the per-AA
    # table (Option B).
    if ribo is not None:
        append_ribosome(cgModel, ribo, nascent_rmin_2=nasc_rm)
        positions = np.vstack([nascent_pos, ribo.coords_nm])
    else:
        positions = nascent_pos

    # 4. tether the current C-terminus (residue L) ---------------------------
    # (skipped for an ejection run: restrain=False -> the tether is released).
    # v2 + trna_tether: peptidyl-tRNA linkage (bond + CA-CA-tRNA orienting
    # angle to the P-site R bead) -- aims the chain down the tunnel. Otherwise a
    # generic harmonic position restraint of residue L to the P-target point.
    if restrain:
        if ribo is not None and params.trna_tether:
            # tRNA tether for the current C-terminus (residue L) to this
            # stage's site (A-site stages 1-2, P-site stage 3): bond + 2 orienting
            # angles + improper + a backbone angle aiming the chain down the tunnel.
            prev_index = (L - 2) if L >= 2 else None
            add_trna_tether(cgModel, L - 1, prev_index, ribo, L, segid=tether_segid)
            # Optionally also tether the previous residue L-1 to its site (the P-site
            # in stage 1, where L sits at A and L-1 rests at P) -- pins both
            # ends of the new peptide bond at the equilibrium-PTC geometry (feature 3).
            if tether_prev_segid is not None and L >= 2:
                pprev_index = (L - 3) if L >= 3 else None
                add_trna_tether(cgModel, L - 2, pprev_index, ribo, L,
                                segid=tether_prev_segid)
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
    # stage's frames (OBSERVATIONS.md #1). Rigid AllBonds constraints sidestep this
    # by removing the fast bond mode; topo keeps flexible bonds (needed to seed
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
    # Per-call minimize gate. Default = params.minimize; a caller may set minimize_override
    # to skip an unnecessary minimization -- e.g. CSP stage 2 continues from stage 1's already
    # relaxed final at the SAME restraint target, so re-minimizing is redundant (the seeded
    # structure setup_simulation loads is already at a local minimum). Retries still skip it
    # too; a diverging stage is instead stabilized by the dt-halving below.
    do_minimize = params.minimize if minimize_override is None else bool(minimize_override)
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
        diverged = False
        max_pe = 0.0
        if do_minimize:
            if attempt == 0:
                _vprint("Minimizing seeded structure (relax placement / new bond)...")
            try:
                ctx.simulation.minimizeEnergy()
                # Re-draw velocities for the relaxed structure.
                ctx.simulation.context.setVelocitiesToTemperature(cfg.ref_t)
            except Exception as exc:
                # A NaN during minimization (e.g. a bead seeded deep inside the stiff
                # ribosome 12-10-6 wall) -> treat as divergence and retry with a
                # halved timestep (the seeded geometry is the same; smaller dt + more steps
                # lets the wall push the bead out without overshooting).
                print(f"[stability] minimization failed ({type(exc).__name__}: "
                      f"{str(exc).splitlines()[0][:80]}); treating as divergence.")
                diverged = True

        if not diverged:
            engine.attach_reporters(cfg, ctx.simulation, suffix="", total_steps=cfg.md_steps)
            if nascent_only:
                # Swap the full-system DCD reporter for a nascent-only one (the rigid
                # ribosome is static -- no need to write its ~thousands of beads/frame).
                ctx.simulation.reporters[1] = NascentDCDReporter(
                    cfg.output_path(".dcd"), cfg.nstxout, nascent_topology, L)

            # Step in chunks so a divergence is caught (and the stage aborted) mid-run.
            chunk = max(cfg.nstxout, cfg.md_steps // 20, 1)
            done = 0
            while done < cfg.md_steps:
                n = min(chunk, cfg.md_steps - done)
                try:
                    ctx.simulation.step(n)
                except Exception as exc:
                    # OpenMM raises "Particle coordinate is NaN" when a stage blows up
                    # outright (stiff EV + 15 fs). Catch it so the dt-halving guard can
                    # retry instead of crashing the whole synthesis.
                    print(f"[stability] integration blew up ({type(exc).__name__}: "
                          f"{str(exc).splitlines()[0][:80]}).")
                    diverged = True
                    break
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

    # Potential energy of the last integrated step (informative per-stage health
    # check; queried before finalize while the context is still current).
    final_pe = ctx.simulation.context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)

    if nascent_only:
        _finalize_nascent(cfg, ctx, nascent_topology, L, start)
    else:
        engine.finalize_simulation(cfg, ctx, built.topology, start)

    # One concise, column-aligned line per stage (verbose banners above are gated
    # on VERBOSE). Fixed-width fields so stacked stage lines line up:
    #   L, stage tag, steps, wall-time, final-step potential energy.
    print(f"  L={L:>3d}  {(label or 'run'):<26s}  {base_steps:>5d} steps  "
          f"{time.time() - start:>6.2f} s  PE={final_pe:>+13.4e} kJ/mol")

    # 6. final NASCENT coordinates seed the next length ----------------------
    final_pdb = cfg.output_path("_final.pdb")
    final = mm.app.PDBFile(final_pdb).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(final)[:L]

