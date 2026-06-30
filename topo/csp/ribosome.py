"""Rigid ribosome scenery + cross-interactions for elongation build step v2.

Build step v1 (:mod:`topo.csp.core`) simulates the **nascent chain
only** and uses the truncated ribosome merely as two fixed anchor *coordinates*.
**Build step v2** adds the truncated ribosome to the System as **rigid (mass-0)
scenery** and wires the two ribosome <-> nascent-chain interactions, following
``DESIGN.md`` §2.2/§3.2 and ``PROMPT.md`` §v2:

1. **Append** the ribosome beads at indices ``L..N-1`` with **mass = 0** (frozen;
   not integrated), coordinates as-is. The P-/A-anchors are now real beads.
2. **Contact force** (the nascent ``L×L`` native/non-native table): give the
   ribosome beads a dummy in-range ``id = 0`` ``addParticle`` entry (never read)
   and restrict the force to the interaction group ``{nascent}×{nascent}`` -- so
   the table stays ``L×L`` and ribosome beads are never evaluated by it.
3. **Ribosome-NC excluded volume:** a separate ``CustomNonbondedForce`` with the
   pure repulsion ``ε·(σ_ij/r)¹²`` (``ε = 0.000132 kcal/mol``,
   ``σ_ij = ½(σ_i+σ_j)`` from the per-particle ``model_parameters['radii']``),
   cutoff 2.0 nm / switch 1.8 nm, interaction group ``{nascent}×{ribosome}``.
4. **Electrostatics:** extend the existing Yukawa force with the ribosome charges
   (rRNA phosphates −1e, charged residues ±1e) and restrict it to
   ``{nascent}×{nascent}`` + ``{nascent}×{ribosome}`` (the rigid ribosome's
   intra-interactions are constant and never computed -- DESIGN §2.2).

The ribosome is held rigid, so **no intra-ribosome forces are ever computed**.
Because the nascent chain uses flexible bonds (no constraints) and the ribosome
has none, no distance constraint ever involves a mass-0 particle, so the
``mass-0`` freezing needs no special constraint handling.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit
from openmm.app.element import carbon as _CARBON

from topo.parameters.model_parameters import parameters as MODEL_PARAMS

# Ribosome-NC excluded-volume well depth: 0.000132 kcal/mol -> kJ/mol.
RIBO_NC_EPS_KJ = 0.000132 * 4.184

# All force constants below are in OpenMM units (kJ/mol/nm^2); the kcal/mol/A^2 ->
# kJ/mol/nm^2 factor is 4.184 * 100 = 418.4, shown for provenance.
# O'Brien peptidyl-tRNA tether (P-site resting geometry), from
# Continuous_synthesis_protocol/continuous_synthesis_v6.py:
#   bond   C-term(CA) -- PtR:76 R      d = 4.76 A,  k = 200 kcal/mol/A^2 = 83680 kJ/mol/nm^2
#   angle  CA(L-1) -- CA(L) -- R       double-Gaussian Ep (TOPO backbone-angle form)
TRNA_TETHER_BOND_NM = 0.476
TRNA_TETHER_BOND_K = 83680.0   # kJ/mol/nm^2 (= 200 kcal/mol/A^2)
# O'Brien orienting-restraint force constants (continuous_synthesis_v6.py):
#   harmonic angle   k = 25 kcal/mol/rad^2,  improper k = 25 kcal/mol/rad^2.
TRNA_TETHER_ANGLE_K = 25.0 * 4.184      # kJ/mol/rad^2 (= 104.6)
TRNA_TETHER_IMPROPER_K = 25.0 * 4.184   # kJ/mol/rad^2
# Per-site resting geometry of the peptidyl-/aminoacyl-tRNA tether (O'Brien), tuple
#   (bond N--R nm, angle N-R-P deg, angle N-R-PU2 deg, improper N-R-P-PU2 deg).
# topo bead names: R == tRNA:76 R, P == O'Brien's P (R-1), BR2 == O'Brien's PU2 (R+2).
_TRNA_SITE_GEOM = {
    "AtR": (0.427, 106.0, 127.0, 128.0),    # A-site (aminoacyl-tRNA): stages 1-2
    "PtR": (0.476, 117.0, 130.0, -161.0),   # P-site (peptidyl-tRNA):  stage 3 / prev-AA
}
# O'Brien improper-dihedral energy form (parse_cg_prm.py): periodic harmonic on |theta-theta0|.
_IMPROPER_ENERGY = ("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); "
                    "pi = 3.1415926535")

# Planar tunnel wall (O'Brien): a one-sided restraint that keeps the nascent chain
# at x >= x0, so it can only extrude *forward* (+x, toward the exit) and cannot move
# backward past the synthesis point into the ribosome interior / the truncated
# underside. x0 is the C-terminal-AA addition plane -- the PTC / P-site where each new
# residue is placed and tethered = P-anchor x + tether bond length (0.5705 + 0.476 ~
# 1.05 nm for this ribosome). (O'Brien quote 58 A, but that is their coordinate frame.)
TUNNEL_WALL_X0_NM = 1.05
TUNNEL_WALL_K = 8368.0         # kJ/mol/nm^2 (= 20 kcal/mol/A^2)


@dataclass
class Ribosome:
    """Parsed rigid CG ribosome: per-bead coordinates (nm) and force parameters."""
    coords_nm: np.ndarray       # (M, 3)
    radii_nm: List[float]       # excluded-volume sigma per bead (model_parameters)
    charges: List[float]
    names: List[str]
    resnames: List[str]
    resids: List[int]
    segids: List[str]

    @property
    def n(self) -> int:
        """Number of ribosome beads.

        Returns
        -------
        int
            The bead count, taken as the length of ``radii_nm``.
        """
        return len(self.radii_nm)


def _bead_type(name: str, resname: str) -> str:
    """Parameter-lookup key for a CG bead (FILES.md mapping).

    Protein Cα beads (atom name ``CA``) look up by residue name; RNA beads look up
    by atom name with trailing digits stripped (``P``, ``R``, ``BR1``/``BR2`` → ``BR``).

    Parameters
    ----------
    name : str
        The PDB atom name of the bead (e.g. ``"CA"``, ``"P"``, ``"BR1"``).
    resname : str
        The PDB residue name of the bead (used for protein Cα beads).

    Returns
    -------
    str
        The key into ``model_parameters`` for this bead: ``resname`` for Cα
        beads, otherwise ``name`` with trailing digits stripped.
    """
    if name == "CA":
        return resname
    return name.rstrip("0123456789")


def load_ribosome(pdb_file: str, model: str = "topo") -> Ribosome:
    """Parse a (truncated) CG ribosome PDB into a :class:`Ribosome`.

    Reads each ATOM/HETATM record, derives its bead type (:func:`_bead_type`), and
    looks up the excluded-volume radius (σ) and charge from
    ``model_parameters[model]``. Coordinates are converted from angstrom to nm.

    Parameters
    ----------
    pdb_file : str
        Path to the (truncated) coarse-grained ribosome PDB file.
    model : str, optional
        Key selecting the parameter set in ``model_parameters`` used for radii
        and charge lookups. Default is ``"topo"``.

    Returns
    -------
    Ribosome
        The parsed ribosome with per-bead coordinates (nm), radii, charges,
        atom names, residue names, residue ids, and segment ids.

    Raises
    ------
    ValueError
        If a parsed bead type is missing from ``model_parameters[model]``, or if
        no ATOM/HETATM records are found in the file.
    """
    params = MODEL_PARAMS[model]
    coords, radii, charges = [], [], []
    names, resnames, resids, segids = [], [], [], []
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            resid = int(line[22:26])
            seg = line[72:76].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            btype = _bead_type(name, resname)
            if btype not in params or not isinstance(params[btype], dict):
                raise ValueError(
                    f"ribosome bead type {btype!r} (atom {name!r}, residue "
                    f"{resname!r}) is not defined in model_parameters[{model!r}].")
            coords.append((x / 10.0, y / 10.0, z / 10.0))
            radii.append(params[btype]["radii"])
            charges.append(params[btype]["charge"])
            names.append(name); resnames.append(resname)
            resids.append(resid); segids.append(seg)
    if not coords:
        raise ValueError(f"no ATOM records parsed from ribosome file {pdb_file!r}.")
    return Ribosome(np.asarray(coords, dtype=float), radii, charges,
                    names, resnames, resids, segids)


def append_ribosome(nascent_model, ribo: Ribosome) -> Tuple[List[int], List[int]]:
    """Append the rigid ribosome to a built nascent model (system + topology).

    Mutates ``nascent_model`` in place (its ``.system`` and ``.topology``):
    appends mass-0 ribosome particles; extends the contact and Yukawa forces with
    ribosome entries and the appropriate interaction groups; adds the ribosome-NC
    ``(σ/r)¹²`` force; and extends the topology with a ribosome chain per segID.

    Must be called **after** :func:`topo.csp.core.build_length_model`
    (the nascent forces must already exist). Returns ``(nascent_indices,
    ribosome_indices)``.
    """
    system = nascent_model.system
    topology = nascent_model.topology
    L = nascent_model.n_atoms
    M = ribo.n
    N = L + M
    nascent_idx = list(range(L))
    ribo_idx = list(range(L, N))

    # 1. Rigid (mass-0) ribosome particles. OpenMM does not integrate mass-0
    #    particles, so they stay fixed; no constraints involve them.
    for _ in range(M):
        system.addParticle(0.0)

    # 2. Contact force: dummy id=0 per ribosome bead (never read) + restrict to
    #    nascent x nascent so the L x L table stays nascent-only.
    cf = nascent_model.custom_non_bonded_force
    for _ in range(M):
        cf.addParticle((0,))
    cf.addInteractionGroup(nascent_idx, nascent_idx)

    # 3. Yukawa electrostatics: add ribosome charges; restrict to nascent-nascent
    #    + nascent-ribosome (no intra-ribosome electrostatics).
    yf = nascent_model.yukawaForce
    for q in ribo.charges:
        yf.addParticle((q,))
    yf.addInteractionGroup(nascent_idx, nascent_idx)
    yf.addInteractionGroup(nascent_idx, ribo_idx)

    # 4. Ribosome-NC excluded volume: pure (sigma/r)^12, nascent x ribosome only.
    nc = mm.CustomNonbondedForce("eps*(sigma/r)^12; sigma=0.5*(sigma1+sigma2)")
    nc.addGlobalParameter("eps", RIBO_NC_EPS_KJ)
    nc.addPerParticleParameter("sigma")
    for s in nascent_model.rf_sigma:          # nascent collision diameters (nm)
        nc.addParticle((s,))
    for s in ribo.radii_nm:                    # ribosome collision diameters (nm)
        nc.addParticle((s,))
    nc.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    nc.setUseSwitchingFunction(True)
    nc.setSwitchingDistance(1.8 * unit.nanometer)
    nc.setCutoffDistance(2.0 * unit.nanometer)
    nc.addInteractionGroup(nascent_idx, ribo_idx)
    # OpenMM (CPU platform) requires every CustomNonbondedForce to share an
    # identical exclusion list. The contact / Yukawa forces carry the nascent
    # bonded (1-2, 1-3) exclusions, so copy them here. They concern nascent-nascent
    # pairs only -- irrelevant to this force's {nascent}x{ribosome} group -- but the
    # lists must match.
    for k in range(cf.getNumExclusions()):
        i, j = cf.getExclusionParticles(k)
        nc.addExclusion(i, j)
    system.addForce(nc)

    # 5. Topology: append a ribosome chain per segID, grouping beads into residues.
    _append_topology(topology, ribo)

    return nascent_idx, ribo_idx


def _append_topology(topology, ribo: Ribosome) -> None:
    """Append the ribosome beads to an OpenMM topology (one chain per segID).

    Beads are grouped into one chain per segment id and one residue per
    ``(segID, resid)`` pair, then each bead is added as a carbon atom under its
    residue. Mutates ``topology`` in place.

    Parameters
    ----------
    topology : openmm.app.Topology
        The topology to extend with the ribosome chains, residues, and atoms.
    ribo : Ribosome
        The parsed ribosome supplying segment ids, residue ids/names, and atom
        names for each bead.

    Returns
    -------
    None
    """
    chains: dict = {}
    residues: dict = {}
    for i in range(ribo.n):
        seg = ribo.segids[i] or "RB"
        chain = chains.get(seg)
        if chain is None:
            chain = topology.addChain(id=seg)
            chains[seg] = chain
        rkey = (seg, ribo.resids[i])
        res = residues.get(rkey)
        if res is None:
            res = topology.addResidue(ribo.resnames[i], chain, id=str(ribo.resids[i]))
            residues[rkey] = res
        topology.addAtom(ribo.names[i], _CARBON, res)


def bead_system_index(ribo: Ribosome, L_nascent: int, segid: str,
                      resid: int, bead: str = "R") -> Optional[int]:
    """System index of a named ribosome bead (e.g. the P-site tRNA R anchor).

    Ribosome beads are appended after the ``L_nascent`` nascent particles in load
    order, so the system index is ``L_nascent + (position in the Ribosome arrays)``.
    Returns None if no such bead exists (e.g. a pyrimidine has no ``BR2``).
    """
    for i in range(ribo.n):
        if ribo.segids[i] == segid and ribo.resids[i] == resid and ribo.names[i] == bead:
            return L_nascent + i
    return None


def add_trna_tether(nascent_model, cterm_index: int, prev_index,
                    ribo: Ribosome, L_nascent: int,
                    segid: str = "PtR", resid: int = 76) -> None:
    """Tether a nascent residue to a tRNA, the full O'Brien orienting way.

    Reproduces the aminoacyl-/peptidyl-tRNA linkage of ``continuous_synthesis_v6.py``
    (``A_site_tRNA_binding`` / ``translocation_AtR``) for the resting geometry of the
    given site (``segid``: ``"AtR"`` A-site / ``"PtR"`` P-site). For the restrained
    residue N (= ``cterm_index``):

    - **bond** ``N -- tRNA:R``  (length/k from :data:`_TRNA_SITE_GEOM` / :data:`TRNA_TETHER_BOND_K`);
    - **orienting angles** ``N -- R -- P`` and ``N -- R -- PU2``  (harmonic,
      :data:`TRNA_TETHER_ANGLE_K`) -- fix the residue's bearing in the tRNA frame;
    - **improper** ``N -- R -- P -- PU2``  (periodic-harmonic on \\|θ−θ0\\|,
      :data:`TRNA_TETHER_IMPROPER_K`) -- fixes the out-of-plane sense;
    - **backbone orienting angle** ``prev -- N -- R`` (double-Gaussian backbone form) --
      aims the chain down the tunnel toward the exit.

    Together (vs. the single bond + 1 angle of the old version) these reproduce
    O'Brien's full orientation control: the chain extrudes N-first down the tunnel
    rather than balling up at the PTC. The peptide bond ``N-1<->N`` is left in place
    (topo keeps the always-bonded chain -- a deliberate deviation, DIFFERENCES.md).

    ``prev_index`` is the CA(N-1) particle index, or None (skips the backbone angle).
    ``P`` (O'Brien R-1) and ``PU2`` (O'Brien R+2) are resolved by topo bead name
    (``"P"`` / ``"BR2"``); a site missing either bead skips that angle/improper.

    Parameters
    ----------
    nascent_model : object
        The built nascent model whose ``.system`` receives the tether forces and
        whose ``.gaussianAngleForce`` receives the backbone orienting angle.
    cterm_index : int
        System index of the restrained nascent CA(N) bead.
    prev_index : int or None
        System index of CA(N-1), or None (skips the backbone angle).
    ribo : Ribosome
        The parsed ribosome (supplies the tRNA R/P/BR2 beads).
    L_nascent : int
        Number of nascent particles preceding the ribosome beads (index offset).
    segid : str, optional
        tRNA segment id: ``"AtR"`` (A-site) or ``"PtR"`` (P-site, default).
    resid : int, optional
        tRNA acceptor residue id. Default is ``76``.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If ``segid`` is not a known tRNA site, or its ``R`` bead is not found.
    """
    import math
    if segid not in _TRNA_SITE_GEOM:
        raise ValueError(f"tRNA tether: unknown site {segid!r} (expected one of "
                         f"{sorted(_TRNA_SITE_GEOM)}).")
    system = nascent_model.system
    R_idx = bead_system_index(ribo, L_nascent, segid, resid, "R")
    if R_idx is None:
        raise ValueError(f"tRNA tether: {segid} resid {resid} R bead not found.")
    P_idx = bead_system_index(ribo, L_nascent, segid, resid, "P")
    U2_idx = bead_system_index(ribo, L_nascent, segid, resid, "BR2")
    bond_nm, ang_P_deg, ang_U2_deg, imp_deg = _TRNA_SITE_GEOM[segid]
    d2r = math.radians

    # 1. bond N -- R (holds the residue at the PTC at its site's resting length).
    bond = mm.HarmonicBondForce()
    bond.addBond(int(cterm_index), R_idx, bond_nm, TRNA_TETHER_BOND_K)
    system.addForce(bond)

    # 2. orienting harmonic angles N -- R -- P and N -- R -- PU2 (tRNA frame).
    haf = mm.HarmonicAngleForce()
    if P_idx is not None:
        haf.addAngle(int(cterm_index), R_idx, P_idx, d2r(ang_P_deg), TRNA_TETHER_ANGLE_K)
    if U2_idx is not None:
        haf.addAngle(int(cterm_index), R_idx, U2_idx, d2r(ang_U2_deg), TRNA_TETHER_ANGLE_K)
    if haf.getNumAngles() > 0:
        system.addForce(haf)

    # 3. improper N -- R -- P -- PU2 (out-of-plane sense; O'Brien periodic-harmonic form).
    if P_idx is not None and U2_idx is not None:
        imf = mm.CustomTorsionForce(_IMPROPER_ENERGY)
        imf.addPerTorsionParameter("k")
        imf.addPerTorsionParameter("theta0")
        imf.addTorsion(int(cterm_index), R_idx, P_idx, U2_idx,
                       [TRNA_TETHER_IMPROPER_K, d2r(imp_deg)])
        system.addForce(imf)

    # 4. backbone orienting angle prev -- N -- R (aims the chain down the tunnel).
    if prev_index is not None and nascent_model.gaussianAngleForce is not None:
        nascent_model.gaussianAngleForce.addAngle(int(prev_index), int(cterm_index), R_idx)


def add_tunnel_wall(system, nascent_indices, x0_nm: float = TUNNEL_WALL_X0_NM,
                    k: float = TUNNEL_WALL_K) -> mm.Force:
    """Add O'Brien's one-sided planar tunnel wall on the nascent chain.

    ``U = k * min(x - x0, 0)^2`` per nascent bead -- penalizes ``x < x0`` only, so the
    chain is kept at ``x >= x0`` and can only extrude forward (+x). ``x0_nm`` is the
    plane position (nm; the C-terminal-AA addition plane / PTC) and ``k`` is the force
    constant in OpenMM units (kJ/mol/nm^2).

    Parameters
    ----------
    system : openmm.System
        The system to which the planar-wall force is added.
    nascent_indices : iterable of int
        System indices of the nascent-chain beads subjected to the wall.
    x0_nm : float, optional
        Plane position (nm); beads are penalized only for ``x < x0``. Default is
        :data:`TUNNEL_WALL_X0_NM`.
    k : float, optional
        Force constant (kJ/mol/nm^2) of the one-sided restraint. Default is
        :data:`TUNNEL_WALL_K`.

    Returns
    -------
    openmm.Force
        The created ``CustomExternalForce`` (already added to ``system``).
    """
    force = mm.CustomExternalForce("k*r^2; r=min(x-x0, 0)")
    force.addGlobalParameter("k", k)
    force.addGlobalParameter("x0", x0_nm)
    for i in nascent_indices:
        force.addParticle(int(i), [])
    system.addForce(force)
    return force
