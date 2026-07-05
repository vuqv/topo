"""Rigid ribosome scenery + cross-interactions for elongation build step v2.

Build step v1 (:mod:`topo.csp.core`) simulates the **nascent chain
only** and uses the truncated ribosome merely as two fixed anchor *coordinates*.
**Build step v2** adds the truncated ribosome to the System as **rigid (mass-0)
scenery** and wires the two ribosome <-> nascent-chain interactions, following
the RNC design (see ``topo/csp/README.md`` and ``docs/usage/continuous_synthesis.md``):

1. **Append** the ribosome beads at indices ``L..N-1`` with **mass = 0** (frozen;
   not integrated), coordinates as-is. The P-/A-anchors are now real beads.
2. **Contact force** (the nascent ``L×L`` native/non-native table): give the
   ribosome beads a dummy in-range ``id = 0`` ``addParticle`` entry (never read)
   and restrict the force to the interaction group ``{nascent}×{nascent}`` -- so
   the table stays ``L×L`` and ribosome beads are never evaluated by it.
3. **Ribosome-NC excluded volume:** a separate ``CustomNonbondedForce`` reproducing
   O'Brien's NC<->ribosome interaction -- the **12-10-6** form
   ``ε[13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶]`` (``ε = 0.000132 kcal/mol``) with the **sum**
   combination rule ``R_ij = Rmin/2_i + Rmin/2_j`` (O'Brien's convention; nascent Rmin/2 =
   per-residue K-B (Option A) or per-AA ``OBRIEN_SC_RMIN_2_NM``, ribosome Rmin/2 from
   ``model_parameters`` via ``load_ribosome``), cutoff
   2.0 nm / switch 1.8 nm, interaction group ``{nascent}×{ribosome}``. (Earlier topo used a
   pure ``ε(σ/r)¹²`` + average rule that was ~1000× too soft -- see
   ``tutorials/15_claude_fix/TOPO_OBrien_NCribosome_nonbonded_compare.md``.)
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
# kJ/mol/nm^2 factor is 4.184 * 100 = 418.4, shown for provenance. They are the
# O'Brien NOMINAL stiffnesses in the E = k (x-x0)^2 convention (no implicit 1/2).
# The CustomForce terms (improper, tunnel wall, C-terminus position restraint) use
# them as-is; the two built-in Harmonic{Bond,Angle}Force tether terms (add_trna_tether)
# evaluate E = 1/2 k (x-x0)^2, so they are passed 2*k there to realize the same
# nominal strength -- keep all four terms consistent if you change these.
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

# O'Brien per-amino-acid sidechain-bead Rmin/2 (nm), the S<aa1> NONBONDED types from his
# CG protein/ribosome .prm (combine_ribo_L24_Yang.prm). These are the fixed per-AA
# excluded-volume radii used for the ribosomal-protein beads. The authoritative values
# now live in model_parameters[...] (the protein Rmin_2 entries, == these S<aa> values),
# so this table is derived from there to keep a single source; the HSD/HSE/HSP histidine
# tautomers alias HIS. Retained as a named table for the Option-B nascent fallback in
# :func:`append_ribosome` (nascent per-AA Rmin/2 when the per-residue K-B array is absent).
_MP_TOPO = MODEL_PARAMS["topo"]
_AA20 = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
         "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
OBRIEN_SC_RMIN_2_NM = {aa: _MP_TOPO[aa]["Rmin_2"] for aa in _AA20}
OBRIEN_SC_RMIN_2_NM.update({t: OBRIEN_SC_RMIN_2_NM["HIS"] for t in ("HSD", "HSE", "HSP")})
# O'Brien's 12-10-6 Go/excluded-volume leading coefficients (U = eps[13(R/r)^12 -
# 18(R/r)^10 + 4(R/r)^6]); well minimum -eps at r=R. Same form parse_cg_prm.py emits.
_NC_126_ENERGY = ("eps*(13*(R/r)^12 - 18*(R/r)^10 + 4*(R/r)^6); R = rm1 + rm2")

# Planar tunnel wall (O'Brien): a one-sided restraint that keeps the nascent chain
# at x >= x0, so it can only extrude *forward* (+x, toward the exit) and cannot move
# backward past the synthesis point into the ribosome interior / the truncated
# underside. The plane x0 is always derived per structure by the CSP runner from the
# PTC-optimized A/P target points (the lower of their x) and passed in explicitly, so
# there is no module-level default for it -- only the wall stiffness k has one.
TUNNEL_WALL_K = 8368.0         # kJ/mol/nm^2 (= 20 kcal/mol/A^2)


@dataclass
class Ribosome:
    """Parsed rigid CG ribosome: per-bead coordinates (nm) and force parameters."""
    coords_nm: np.ndarray       # (M, 3)
    Rmin_2_nm: List[float]      # per-bead collision radius Rmin/2 (nm; sum-rule EV)
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
            The bead count, taken as the length of ``Rmin_2_nm``.
        """
        return len(self.Rmin_2_nm)


def _bead_type(name: str, resname: str) -> str:
    """Parameter-lookup key for a CG bead.

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
    looks up its collision radius Rmin/2 and charge from ``model_parameters[model]``
    -- the fixed per-type scenery table: **RNA beads** (P/R/BR) by O'Brien's per-type
    Rmin/2, **protein Cα beads** by O'Brien's per-AA sidechain Rmin/2 (his ribosome
    ``S<aa>`` values; the one structure-based Go protein L24 is treated the same
    per-AA way -- his per-residue B-type values are not reproduced). Coordinates are
    converted from angstrom to nm. This reproduces O'Brien's ribosome excluded volume
    from topo's own CG PDB -- no .cor/.psf/.prm needed. The mobile nascent chain uses
    its own per-residue K-B Rmin/2 (Option A), never this table (see model_parameters).

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
            # Rmin/2 (collision radius, nm) and charge from model_parameters -- the fixed
            # per-type scenery table: RNA P/R/BR per-type, ribosomal protein per-AA
            # (O'Brien's S<aa> values). Correct for the *rigid* ribosome; the mobile
            # nascent chain uses its own per-residue K-B Rmin/2 (Option A), never this.
            coords.append((x / 10.0, y / 10.0, z / 10.0))
            radii.append(params[btype]["Rmin_2"])
            charges.append(params[btype]["charge"])
            names.append(name); resnames.append(resname)
            resids.append(resid); segids.append(seg)
    if not coords:
        raise ValueError(f"no ATOM records parsed from ribosome file {pdb_file!r}.")
    return Ribosome(np.asarray(coords, dtype=float), radii, charges,
                    names, resnames, resids, segids)


def append_ribosome(nascent_model, ribo: Ribosome,
                    nascent_rmin_2: Optional[np.ndarray] = None
                    ) -> Tuple[List[int], List[int]]:
    """Append the rigid ribosome to a built nascent model (system + topology).

    Mutates ``nascent_model`` in place (its ``.system`` and ``.topology``):
    appends mass-0 (fixed) ribosome particles; extends the contact and Yukawa
    forces with ribosome entries and restricts their interaction groups; adds the
    O'Brien-consistent ribosome-NC excluded-volume force (the 12-10-6 form
    ``U = eps[13(R/r)¹² - 18(R/r)¹⁰ + 4(R/r)⁶]`` with the *sum* combination rule
    ``R = rm_i + rm_j``, acting on {nascent}×{ribosome} only); and extends the
    topology with one ribosome chain per segID.

    Must be called **after** :func:`topo.csp.core.build_length_model`
    (the nascent forces must already exist).

    Parameters
    ----------
    nascent_model
        The built nascent model whose ``.system`` and ``.topology`` are mutated.
    ribo : Ribosome
        The parsed rigid ribosome (coords, per-bead Rmin/2, charges, topology).
    nascent_rmin_2 : np.ndarray, optional
        Per-nascent-bead Rmin/2 (nm) for the NC-ribosome excluded volume. When
        given (length must equal the nascent atom count), these are O'Brien's
        structure-derived per-residue Karanicolas-Brooks collision radii
        ("Option A"). When ``None``, falls back to per-AA sidechain radii
        ``OBRIEN_SC_RMIN_2_NM`` ("Option B").

    Returns
    -------
    (nascent_indices, ribosome_indices) : Tuple[List[int], List[int]]
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

    # 4. Ribosome-NC excluded volume, O'Brien-consistent (nascent x ribosome only).
    #    Reproduces O'Brien's continuous_synthesis_v6 NC<->ribosome interaction (the rnc.xml
    #    CustomNonbondedForce): the 12-10-6 form U = eps[13(R/r)^12 - 18(R/r)^10 + 4(R/r)^6]
    #    with the SUM combination rule R_ij = Rmin/2_i + Rmin/2_j (NOT the average 0.5*(s1+s2)),
    #    same eps (RIBO_NC_EPS_KJ = 0.000132 kcal/mol). Per TOPO_OBrien_NCribosome_nonbonded_
    #    compare.md, the old pure (sigma/r)^12 + average rule was ~1000x too soft. Nascent
    #    per-bead Rmin/2 uses O'Brien's per-AA sidechain values (OBRIEN_SC_RMIN_2_NM, "Option B");
    #    ribosome per-bead Rmin/2 are O'Brien's values (ribo.Rmin_2_nm: RNA per-type,
    #    protein per-AA). Both sides are thus on O'Brien's Rmin/2 convention.
    nc = mm.CustomNonbondedForce(_NC_126_ENERGY)
    nc.addGlobalParameter("eps", RIBO_NC_EPS_KJ)
    nc.addPerParticleParameter("rm")    # per-bead Rmin/2 (nm); pair R = rm1 + rm2
    # Nascent per-bead Rmin/2:
    #  - Option A (default when nascent_rmin_2 given): O'Brien's structure-derived per-residue
    #    Karanicolas-Brooks collision radius (A1..An in his .prm) -- this is what O'Brien's
    #    nascent chain actually uses for the nascent<->ribosome excluded volume.
    #  - Option B (fallback): per-AA sidechain radii OBRIEN_SC_RMIN_2_NM (his ribosomal-protein
    #    S<aa1> values); kept for reference / when the K-B array is unavailable.
    nascent_atoms = list(topology.atoms())[:L]   # nascent CA beads (ribosome not appended yet)
    if nascent_rmin_2 is not None:
        if len(nascent_rmin_2) != L:
            raise ValueError(f"nascent_rmin_2 has {len(nascent_rmin_2)} entries but L={L}.")
        for rm in nascent_rmin_2:                  # Option A: per-residue K-B Rmin/2 (nm)
            nc.addParticle((float(rm),))
    else:
        for atom in nascent_atoms:                # Option B: per-AA sidechain Rmin/2 (nm)
            rn = atom.residue.name
            if rn not in OBRIEN_SC_RMIN_2_NM:
                raise ValueError(f"NC-ribosome EV: residue {rn!r} has no O'Brien Rmin/2 "
                                 f"(OBRIEN_SC_RMIN_2_NM); cannot set the nascent excluded-volume radius.")
            nc.addParticle((OBRIEN_SC_RMIN_2_NM[rn],))
    for s in ribo.Rmin_2_nm:                      # ribosome: O'Brien Rmin/2 (nm)
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


def anchor_coord(ribo: Ribosome, segid: str, resid: int = 76,
                 bead: str = "R") -> np.ndarray:
    """Coordinate (nm) of a named ribosome bead from a loaded :class:`Ribosome`.

    The :class:`Ribosome`-object analog of :func:`topo.csp.core.read_anchor` (which parses
    a PDB): used to pick the P-/A-anchors (``segid='PtR'/'AtR'``, resid 76, ``R`` bead)
    directly from the loaded :class:`Ribosome` arrays. Raises if the bead is absent or
    non-unique.
    """
    matches = [i for i in range(ribo.n)
               if ribo.segids[i] == segid and ribo.resids[i] == resid
               and ribo.names[i] == bead]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one bead (segid={segid!r}, resid={resid}, "
                         f"name={bead!r}) in the ribosome, found {len(matches)}.")
    return ribo.coords_nm[matches[0]]


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

    # OpenMM's built-in HarmonicBondForce/HarmonicAngleForce evaluate E = 1/2 k (x-x0)^2,
    # so we pass 2*k to realize the O'Brien nominal stiffness (E = k (x-x0)^2, the
    # CHARMM/no-1/2 convention in which TRNA_TETHER_*_K are defined). This keeps the
    # tether bond/angle at the SAME full strength as the improper, tunnel wall and
    # C-terminus position restraint, which use CustomForces (no implicit 1/2).

    # 1. bond N -- R (holds the residue at the PTC at its site's resting length).
    bond = mm.HarmonicBondForce()
    bond.addBond(int(cterm_index), R_idx, bond_nm, 2.0 * TRNA_TETHER_BOND_K)
    system.addForce(bond)

    # 2. orienting harmonic angles N -- R -- P and N -- R -- PU2 (tRNA frame).
    haf = mm.HarmonicAngleForce()
    if P_idx is not None:
        haf.addAngle(int(cterm_index), R_idx, P_idx, d2r(ang_P_deg), 2.0 * TRNA_TETHER_ANGLE_K)
    if U2_idx is not None:
        haf.addAngle(int(cterm_index), R_idx, U2_idx, d2r(ang_U2_deg), 2.0 * TRNA_TETHER_ANGLE_K)
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


def add_tunnel_wall(system, nascent_indices, x0_nm: float,
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
    x0_nm : float
        Plane position (nm); beads are penalized only for ``x < x0``. Always supplied
        by the CSP runner from the PTC-optimized A/P targets (structure-dependent);
        there is no default.
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
