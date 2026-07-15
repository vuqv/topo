"""
Build non-bonded interaction matrices for the TOPO coarse-grained model.

This module constructs distance and energy matrices used by OpenMM for
structure-based (native) contacts and repulsive non-native contacts.
Inputs: PDB structure, STRIDE hydrogen-bond output, and domain definition YAML.

Main workflow
-------------
1. **Hydrogen bonds** (STRIDE): parsed from stride output; 1 H-bond -> 0.75 kcal/mol,
   2+ H-bonds capped at 1.5 kcal/mol per pair.
2. **Backbone-sidechain (BS)** and **sidechain-sidechain (SS)** contacts: distance-based
   (default cutoff 4.5 Angstrom); BS energy 0.37 kcal/mol; SS energies from BT potential
   and domain scaling.
3. **Domain scaling**: YAML defines intra_domains (residue lists + nscale) and
   optional inter_domains (pair nscales). ``nscale`` is the native-contact scaling
   factor applied to the sidechain-sidechain well depths (not an absolute energy --
   the contacts still interact at nscale = 1.0). Single-domain proteins can omit
   inter_domains.
4. **Non-native contacts**: repulsive well with per-residue Rmin/2 from the
   nearest non-contact CA-CA distance (sum rule); energy ENERGY_PARAMS['non_native'].

Typical use
-----------
>>> from topo.utils.nonbonded import build_nonbonded_interaction
>>> rmin_matrix, energy_matrix = build_nonbonded_interaction(
...     "protein.pdb", "domain.yaml", "stride.dat"
... )
>>> # rmin_matrix (well positions), energy_matrix in nm and kJ/mol for OpenMM
"""
import os
import subprocess
import warnings
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd
import re
from pathlib import Path
from collections import defaultdict
import yaml
from typing import Dict, List, Tuple, Set, Optional

from topo.utils.external import find_executable
# import logging

# Configure logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
KCAL_TO_KJ = 4.184
DEFAULT_CUTOFF = 4.5  # Angstrom; distance cutoff for backbone-sidechain and sidechain-sidechain contacts
LOCAL_SEPARATION = 2   # Residues within this sequence separation are excluded from non-local contacts
RMIN_SCALE_FACTOR = 2**(1/6)  # 2^(1/6): scales the nearest-contact Ca distance to the collision diameter Rmin (LJ-like)
DISTANCE_TO_NM = 10.0  # Angstrom to nanometer (CA distances read in Angstrom)

# Energy parameters (kcal/mol values, stored in kJ/mol for OpenMM)
ENERGY_PARAMS = {
    'hydrogen_bond': 0.75 * KCAL_TO_KJ,       # per H-bond; 2+ H-bonds capped at 1.5 kcal/mol total
    'backbone_sidechain': 0.37 * KCAL_TO_KJ,
    'non_native': 0.000132 * KCAL_TO_KJ,
    'bt_shift': 0.6   # reference shift for the BT contact potential (kcal/mol)
}


# -----------------------------------------------------------------------------
# Chain-ID handling (multi-chain support)
# -----------------------------------------------------------------------------
def _norm_chain(chain) -> str:
    """
    Normalize a chain identifier to a canonical string.

    Residues are keyed by ``(chain, resid)`` so that structures with multiple
    chains (where the same residue number can appear in more than one chain) are
    handled correctly. Different sources spell a blank chain differently:
    MDAnalysis reports ``''`` or a space, while STRIDE prints ``'-'``. All of
    these are normalized to the empty string so the two sources agree.
    """
    chain = ('' if chain is None else str(chain)).strip()
    return '' if chain in ('', '-') else chain


def _residue_chain(res) -> str:
    """
    Return the normalized chain ID for an MDAnalysis residue.

    Prefers the PDB ``chainID`` (matching STRIDE's chain column); falls back to
    ``segid`` if ``chainID`` is unavailable.
    """
    chain = None
    try:
        chain_ids = res.atoms.chainIDs
        if len(chain_ids) > 0:
            chain = chain_ids[0]
    except (AttributeError, IndexError):
        chain = None
    if _norm_chain(chain) == '':
        # Fall back to segid when chainID is blank/missing.
        try:
            chain = res.segid
        except AttributeError:
            chain = None
    return _norm_chain(chain)


def get_residue_mapping(universe: mda.Universe) -> Tuple[Dict[Tuple[str, int], int], Dict[int, str], int]:
    """
    Build residue index and name mappings from an MDAnalysis protein universe.

    Uses the order of residues in ``universe.select_atoms("protein").residues``
    to define 0-based indices. Useful for filling contact/energy matrices by
    residue index instead of PDB resid.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis universe with a loaded protein structure.

    Returns
    -------
    key_to_index : dict
        Maps ``(chain, resid)`` -> 0-based index in residue list, where ``chain``
        is the normalized chain ID (see :func:`_norm_chain`) and ``resid`` is the
        PDB residue ID. Keying on ``(chain, resid)`` instead of ``resid`` alone is
        required for multi-chain structures, where the same residue number can
        appear in more than one chain.
    index_to_resname : dict
        Maps 0-based index -> three-letter residue name (e.g. 'ALA', 'GLY').
    n_residues : int
        Total number of protein residues.

    Example
    -------
    >>> u = mda.Universe("protein.pdb")
    >>> key_to_idx, idx_to_name, n = get_residue_mapping(u)
    >>> key_to_idx[('A', 1)]   # first residue's index (chain A, resid 1)
    0
    >>> idx_to_name[0]
    'MET'
    """
    residues = universe.select_atoms("protein").residues
    n_residues = len(residues)
    key_to_index = {(_residue_chain(res), res.resid): idx
                    for idx, res in enumerate(residues)}
    index_to_resname = {idx: res.resname for idx, res in enumerate(residues)}
    return key_to_index, index_to_resname, n_residues

def parse_hydrogen_bonds(stride_output_file: str) -> List[Tuple]:
    """
    Parse hydrogen bonds from a STRIDE output file (donor-acceptor pairs).

    Only DNR (donor) lines are read so that each physical H-bond is counted once.
    STRIDE writes each bond twice (DNR and ACC); using DNR only avoids double-counting
    while still allowing pairs with multiple H-bonds to appear multiple times in the list.

    Parameters
    ----------
    stride_output_file : str
        Path to the STRIDE output file (e.g. ``stride.dat``).

    Returns
    -------
    list of tuple
        Each element is a pair ``((chain1, resid1), (chain2, resid2))`` with the
        tuple sorted so the pair is canonical. ``chain`` is the normalized chain
        ID and ``resid`` the PDB residue ID. Repeated pairs indicate multiple
        H-bonds between the same residue pair. The ``(chain, resid)`` keys match
        those produced by :func:`get_residue_mapping`, so multi-chain structures
        (where a residue number repeats across chains) are resolved correctly.

    Raises
    ------
    FileNotFoundError
        If the STRIDE file cannot be opened.

    Notes
    -----
    STRIDE line format: ``DNR  RES chain  resid  seq ->  RES chain  resid  seq``,
    where ``chain`` is a chain letter (e.g. ``A``) or ``-`` for a blank chain.
    The regex captures resname, chain, and the first resid for donor and acceptor.

    Example
    -------
    >>> pairs = parse_hydrogen_bonds("stride.dat")
    >>> pairs[0]  # e.g. (('A', 3), ('A', 43))
    (('A', 3), ('A', 43))
    """
    try:
        with open(stride_output_file, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"STRIDE output file not found: {stride_output_file}")
        raise

    # Fields after the residue name are: chain ID, resid, seq. The chain may be a
    # chain letter (e.g. "A") or "-" for a blank chain; capturing it (\S+) lets us
    # key residues by (chain, resid) and supports multi-chain structures. The
    # original pattern matched only a literal "-" for the chain, so chained
    # structures parsed zero H-bonds (dropping all backbone hydrogen-bond energy).
    pattern = re.compile(
        r"(?:DNR|ACC)\s+\w+\s+(\S+)\s+(\d+)\s+\d+\s+->\s+\w+\s+(\S+)\s+(\d+)\s+\d+"
    )

    hb_pairs = []
    # STRIDE reports each H-bond twice (DNR and ACC). Count each bond once by using only DNR lines.
    # Then pairs with multiple physical H-bonds appear as multiple DNR lines and are counted correctly.
    for line in lines:
        if line.startswith("DNR"):
            match = pattern.search(line)
            if match:
                chain1, res1_pdb, chain2, res2_pdb = match.groups()
                donor = (_norm_chain(chain1), int(res1_pdb))
                acceptor = (_norm_chain(chain2), int(res2_pdb))
                pair = tuple(sorted([donor, acceptor]))
                hb_pairs.append(pair)
    return hb_pairs

def build_hb_contact_matrix(hb_pairs: List[Tuple],
                            key_to_index: Dict[Tuple[str, int], int],
                            n_residues: int) -> np.ndarray:
    """
    Build the hydrogen-bond contact matrix from a list of H-bond pairs.

    Counts how many H-bonds exist between each residue pair and stores 0, 1, or 2
    (values greater than 2 are capped at 2). The model uses 0.75 kcal/mol per single
    H-bond and 1.5 kcal/mol total for multiple H-bonds; capping at 2 enforces that.

    Parameters
    ----------
    hb_pairs : list of tuple
        List of canonical pairs from :func:`parse_hydrogen_bonds`. Each pair is
        ``((chain1, resid1), (chain2, resid2))``.
    key_to_index : dict
        Maps ``(chain, resid)`` -> 0-based residue index, from
        :func:`get_residue_mapping`. Used to place each H-bond at the correct
        matrix index regardless of chain (multi-chain safe).
    n_residues : int
        Number of residues (determines matrix shape).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues), dtype int
        Symmetric matrix. Entry [i, j] is 0, 1, or 2 (number of H-bonds between
        residue index i and j, capped at 2). Indices are the 0-based residue order
        from :func:`get_residue_mapping`.

    Notes
    -----
    H-bond ``(chain, resid)`` keys are resolved to matrix indices via
    ``key_to_index``. A pair whose residue is not found in the mapping (e.g. a
    hetero residue STRIDE reported but that is not in the protein selection) is
    skipped with a warning rather than mis-placed.
    """
    hb_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)

    # Count hydrogen bonds between residue pairs
    pair_counts = defaultdict(int)
    for donor, acceptor in hb_pairs:
        pair = tuple(sorted([donor, acceptor]))
        pair_counts[pair] += 1

    # Fill contact matrix; cap at 2 (model: 1 H-bond -> 0.75 kcal/mol, 2+ H-bonds -> 1.5 kcal/mol only)
    for (key1, key2), count in pair_counts.items():
        if key1 not in key_to_index or key2 not in key_to_index:
            print(f"Warning: H-bond residue {key1} or {key2} not found in residue "
                  f"mapping; skipping this H-bond.")
            continue
        i = key_to_index[key1]
        j = key_to_index[key2]
        val = min(count, 2)
        hb_contact_matrix[i, j] = val
        hb_contact_matrix[j, i] = val

    return hb_contact_matrix


def print_pairs_with_multiple_hb(hb_pairs: List[Tuple]) -> Dict[Tuple, int]:
    """
    Count H-bonds per residue pair, print pairs with more than one, and return counts.

    Parameters
    ----------
    hb_pairs : list of tuple
        List of canonical H-bond pairs from :func:`parse_hydrogen_bonds`.

    Returns
    -------
    dict
        Maps canonical pair ``((chain1, resid1), (chain2, resid2))`` to number
        of H-bonds (before any cap). Useful for logging or downstream analysis.

    Example
    -------
    >>> pairs = parse_hydrogen_bonds("stride.dat")
    >>> counts = print_pairs_with_multiple_hb(pairs)
    Residue pairs with more than 1 hydrogen bond:
      A3 -- A43:  2 H-bonds
      ...
      (total: 20 pairs)
    >>> counts[(('A', 3), ('A', 43))]
    2
    """
    pair_counts = defaultdict(int)
    for donor, acceptor in hb_pairs:
        pair = tuple(sorted([donor, acceptor]))
        pair_counts[pair] += 1
    multi = [(p, c) for p, c in pair_counts.items() if c > 1]
    multi.sort(key=lambda x: (-x[1], x[0]))
    if multi:
        print("Residue pairs with more than 1 hydrogen bond:")
        for (c1, i1), (c2, i2) in [p for p, _ in multi]:
            c = pair_counts[tuple(sorted([(c1, i1), (c2, i2)]))]
            print(f"  {c1}{i1} -- {c2}{i2}:  {c} H-bonds")
        print(f"  (total: {len(multi)} pairs)")
    return dict(pair_counts)


def get_hb_contact_matrix(stride_output_file: str,
                          key_to_index: Dict[Tuple[str, int], int],
                          n_residues: int) -> np.ndarray:
    """
    Build the hydrogen-bond contact matrix directly from a STRIDE output file.

    Convenience function that parses the file and builds the matrix. Values are
    0, 1, or 2 (multiple H-bonds capped at 2).

    Parameters
    ----------
    stride_output_file : str
        Path to the STRIDE output file.
    key_to_index : dict
        Maps ``(chain, resid)`` -> 0-based residue index, from
        :func:`get_residue_mapping` (multi-chain safe).
    n_residues : int
        Number of residues (must match the system used for STRIDE).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues), dtype int
        Symmetric H-bond contact matrix. See :func:`build_hb_contact_matrix`.

    Example
    -------
    >>> key_to_index, _, n_residues = get_residue_mapping(u)
    >>> hb_matrix = get_hb_contact_matrix("stride.dat", key_to_index, n_residues)
    >>> hb_matrix.shape
    (283, 283)
    >>> (hb_matrix > 1).sum() // 2   # number of pairs with 2 H-bonds (symmetric)
    20
    """
    hb_pairs = parse_hydrogen_bonds(stride_output_file)
    return build_hb_contact_matrix(hb_pairs, key_to_index, n_residues)

def get_bs_contact_matrix(u: mda.Universe, cutoff: float = DEFAULT_CUTOFF) -> np.ndarray:
    """
    Build backbone–sidechain contact matrix from structure (distance-based).

    Two residues are in contact if any backbone atom of one is within ``cutoff``
    of any sidechain atom of the other (and vice versa), and they are not
    within LOCAL_SEPARATION along the sequence. The matrix is symmetric and
    stores the total number of directional contacts (i → j and j → i).

    Parameters
    ----------
    u : mda.Universe
        MDAnalysis universe with protein structure (must have backbone/sidechain).
    cutoff : float, optional
        Distance cutoff in Angstrom (default :const:`DEFAULT_CUTOFF` = 4.5).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues), dtype int
        Symmetric matrix. Entry [i, j] is the number of directional backbone–sidechain
        contacts between residue index i and j (0, 1, or 2). Indices are 0-based
        residue order from :func:`get_residue_mapping`.

    Notes
    -----
    Backbone: protein backbone atoms excluding H. Sidechain: protein non-backbone
    excluding H. Pairs within abs(resid_i - resid_j) <= LOCAL_SEPARATION are excluded,
    but only when both residues are in the same chain; residues in different chains
    are never excluded by sequence separation (multi-chain safe).

    **Why the count is 0, 1, or 2 (and not just a 0/1 flag).** A BS contact is
    *directional*: the backbone belongs to one residue and the sidechain to the
    other, so the roles are not interchangeable. A pair (i, j) therefore admits two
    independent contacts, and either, both, or neither may be present:

    - **0** — neither backbone reaches the other's sidechain. (The pair may still be
      a native contact via an H-bond or an SS contact; see
      :func:`build_nonbonded_interaction`.)
    - **1** — exactly one direction: backbone of i is within `cutoff` of a sidechain
      atom of j, *or* backbone of j is within `cutoff` of a sidechain atom of i.
    - **2** — both directions at once: backbone of i touches sidechain of j *and*
      backbone of j touches sidechain of i (mutually interdigitated).

    That is exactly what the two-step construction below computes: `contacts_bs`
    holds directed (backbone-residue -> sidechain-residue) pairs, which are stored in
    an **asymmetric** matrix where `M[i, j] = 1` means "backbone of i within cutoff of
    sidechain of j". The returned value is the symmetrized count `M + M.T`, whose
    entries are the number of directions realized. The energy is linear in that count
    (0, 1, or 2 x ENERGY_PARAMS['backbone_sidechain']), so a count of 2 marks a
    tighter, mutually buried pair and earns twice the well depth.

    Two easy confusions worth flagging:

    - This is **not** the same rule as the H-bond count, which also tops out at 2
      (see :func:`get_hb_contact_matrix`) but for an unrelated reason: there the 2
      counts *two H-bonds* reported by STRIDE, not two directions. The BS 2 is a
      structural maximum, not a cap.
    - Glycine has no sidechain atoms, so it can only ever be the backbone partner —
      any BS contact involving a glycine caps at 1 by construction.
    """
    key_to_index, _, n_residues = get_residue_mapping(u)

    # The MDAnalysis 'backbone' keyword only matches N, CA, C, O, so the terminal
    # carboxyl oxygens (OXT/OT1/OT2) would otherwise fall into the "not backbone"
    # sidechain set. They are backbone (terminal) atoms, not sidechain; counting
    # them as sidechain creates spurious sidechain contacts for the C-terminal
    # residue (and gives glycine, which has no sidechain, a fake one).
    backbone = u.select_atoms('protein and (backbone or name OXT OT1 OT2) and not name H*')
    sidechain = u.select_atoms("protein and not backbone and not name H* OXT OT1 OT2")

    dists_bs = distance_array(backbone.positions, sidechain.positions)

    # Per-atom (chain, resid) keys (chain-aware so resid numbers can repeat across chains)
    bb_keys = [(_residue_chain(a.residue), a.resid) for a in backbone]
    sc_keys = [(_residue_chain(a.residue), a.resid) for a in sidechain]

    # Build directional residue-residue contact list, keyed by (chain, resid).
    # Sequence-local exclusion applies only within the same chain.
    contacts_bs = set()
    for i, k1 in enumerate(bb_keys):
        for j, k2 in enumerate(sc_keys):
            if dists_bs[i, j] <= cutoff:
                if k1[0] == k2[0] and abs(k1[1] - k2[1]) <= LOCAL_SEPARATION:
                    continue
                contacts_bs.add((k1, k2))  # preserve direction: bb → sc

    # Build asymmetric matrix first
    bs_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    for k1, k2 in contacts_bs:
        bs_contact_matrix[key_to_index[k1], key_to_index[k2]] = 1

    # Symmetric count = directional bb→sc plus sc→bb. The contact list is already
    # separation-filtered, so a plain symmetrization is correct (and chain-safe).
    bs_symmetric_count = bs_contact_matrix + bs_contact_matrix.T

    return bs_symmetric_count

def get_ss_contact_matrix(u: mda.Universe, cutoff: float = DEFAULT_CUTOFF) -> np.ndarray:
    """
    Build sidechain–sidechain contact matrix from structure (distance-based).

    Two residues are in contact if any sidechain atom of one is within ``cutoff``
    of any sidechain atom of the other, and they are not within LOCAL_SEPARATION
    along the sequence. The matrix is binary (0 or 1) and symmetric.

    Parameters
    ----------
    u : mda.Universe
        MDAnalysis universe with protein structure.
    cutoff : float, optional
        Distance cutoff in Angstrom (default :const:`DEFAULT_CUTOFF` = 4.5).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues), dtype int
        Symmetric binary matrix. Entry [i, j] is 1 if residue pair (i, j) has
        at least one sidechain–sidechain contact within cutoff, else 0. Indices
        are 0-based residue order from :func:`get_residue_mapping`.
    """
    key_to_index, _, n_residues = get_residue_mapping(u)
    # Exclude terminal carboxyl oxygens (OXT/OT1/OT2): they are backbone atoms,
    # not sidechain, so they must not generate sidechain-sidechain contacts. See
    # the note in get_bs_contact_matrix.
    sidechain = u.select_atoms("protein and not backbone and not name H* OXT OT1 OT2")

    dists_ss = distance_array(sidechain.positions, sidechain.positions)
    ss_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)

    # Per-atom (chain, resid) keys (chain-aware so resid numbers can repeat across chains)
    sc_keys = [(_residue_chain(a.residue), a.resid) for a in sidechain]

    # Build residue-residue contact list, keyed by (chain, resid). Sequence-local
    # exclusion applies only within the same chain.
    contacts_ss = set()
    for i, k1 in enumerate(sc_keys):
        for j, k2 in enumerate(sc_keys):
            if dists_ss[i, j] <= cutoff:
                if k1[0] == k2[0] and abs(k1[1] - k2[1]) <= LOCAL_SEPARATION:
                    continue
                contacts_ss.add(tuple(sorted([k1, k2])))

    # Fill contact matrix
    for k1, k2 in contacts_ss:
        ss_contact_matrix[key_to_index[k1], key_to_index[k2]] = 1
        ss_contact_matrix[key_to_index[k2], key_to_index[k1]] = 1

    return ss_contact_matrix

def load_bt_potential(bt_file: str = 'bt_potential.csv') -> pd.DataFrame:
    """
    Load BT (Betancourt–Thirumalai) potential matrix from CSV.

    If ``bt_file`` is a filename only (no path), the file is loaded from
    ``topo/parameters/data/`` so the package works regardless of current
    working directory. If ``bt_file`` is an absolute path, that path is used.

    Parameters
    ----------
    bt_file : str, optional
        Filename (e.g. ``'bt_potential.csv'``) or absolute path to the CSV.
        Default ``'bt_potential.csv'`` is looked up in ``topo/parameters/data/``.

    Returns
    -------
    pd.DataFrame
        Matrix indexed by residue name (rows and columns). Values are in kJ/mol,
        computed as ``KCAL_TO_KJ * |raw_value - bt_shift|`` from the CSV.

    Raises
    ------
    FileNotFoundError
        If the CSV file cannot be found.

    Example
    -------
    >>> df = load_bt_potential()  # uses topo/parameters/data/bt_potential.csv
    >>> df.loc['ALA', 'GLY']  # energy for ALA–GLY pair in kJ/mol
    1.234
    """
    bt_path = Path(bt_file)
    if not bt_path.is_absolute():
        data_dir = Path(__file__).resolve().parent.parent / 'parameters' / 'data'
        bt_path = data_dir / bt_path.name
    try:
        df = pd.read_csv(bt_path, index_col=0)
        return KCAL_TO_KJ * np.abs(df - ENERGY_PARAMS['bt_shift'])
    except FileNotFoundError:
        print(f"BT potential file not found: {bt_path}")
        raise

def get_ss_interaction_energy(u: mda.Universe, bt_file: str = 'bt_potential.csv') -> np.ndarray:
    """
    Build residue–residue sidechain interaction energy matrix from BT potential.

    Uses the BT potential CSV (residue name vs residue name) and the residue
    sequence in ``u`` to fill an (n_residues × n_residues) matrix of pairwise
    energies in kJ/mol.

    Parameters
    ----------
    u : mda.Universe
        MDAnalysis universe; residue order determines matrix indices.
    bt_file : str, optional
        Passed to :func:`load_bt_potential` (default ``'bt_potential.csv'``).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues)
        Matrix of pairwise sidechain interaction energies in kJ/mol. Entry
        [i, j] is the BT energy for the residue types at indices i and j.
    """
    eps_ss = load_bt_potential(bt_file)
    _, index_to_resname, n_residues = get_residue_mapping(u)
    
    sc_interaction_energy = np.zeros((n_residues, n_residues))
    for i in range(n_residues):
        for j in range(n_residues):
            sc_interaction_energy[i, j] = eps_ss.loc[index_to_resname[i], 
                                                   index_to_resname[j]]
    
    return sc_interaction_energy

def parse_residue_list(residue_items: List) -> List[int]:
    """
    Parse a list of residue specifiers into a flat list of residue numbers.

    Accepts integers and strings. Strings may be a single number (e.g. ``"5"``)
    or a range ``"start-end"`` (inclusive), which is expanded to all integers
    in [start, end].

    Parameters
    ----------
    residue_items : list
        Elements are int or str. Examples: ``[1, 2, "5-10", 15]``.

    Returns
    -------
    list of int
        Sorted residue numbers are not guaranteed; order follows the input
        and range expansion.

    Example
    -------
    >>> parse_residue_list([1, "3-5", 7])
    [1, 3, 4, 5, 7]
    """
    residues = []
    for item in residue_items:
        if isinstance(item, int):
            residues.append(item)
        elif isinstance(item, str):
            if '-' in item:
                start, end = map(int, item.split('-'))
                residues.extend(range(start, end + 1))
            else:
                residues.append(int(item))
    return residues

def read_yaml_config(filepath: str) -> Tuple[Dict, Dict, Dict]:
    """
    Read and parse domain definition YAML (intra/inter nscales, residue lists).

    ``nscale`` is the native-contact scaling factor applied to the sidechain-sidechain
    contact well depths -- a multiplier, not an absolute energy (the contacts still
    interact at ``nscale = 1.0``). Per-domain values scale intra-domain contacts;
    per-interface values scale contacts between two domains.

    Required keys: ``intra_domains``, ``n_residues``. Optional: ``inter_domains``
    (omit for single-domain proteins; then inter_nscales will be empty).
    Residues not listed in any domain are assigned to domain ``'X'`` with
    intra nscale 1.0 and inter 1.0 to all other domains.

    The per-domain scaling key is ``nscale``. The legacy key ``strength`` is still
    accepted as a deprecated alias (a one-time deprecation notice is printed).

    Parameters
    ----------
    filepath : str
        Path to the YAML file (e.g. ``domain.yaml``).

    Returns
    -------
    domain_to_residues : dict
        Domain name -> list of residue numbers (1-based).
    intra_nscales : dict
        Domain name -> float (intra-domain contact nscale).
    inter_nscales : dict
        (domain1, domain2) -> float (inter-domain nscale); symmetric keys
        (d1, d2) and (d2, d1) are both set.

    Raises
    ------
    FileNotFoundError
        If the YAML file cannot be opened.
    KeyError
        If a domain entry has neither an ``nscale`` nor a legacy ``strength`` value.

    Example
    -------
    YAML format::

        n_residues: 110
        intra_domains:
          A: { residues: [1-50], nscale: 1.0 }
          B: { residues: [51-110], nscale: 1.0 }
        inter_domains:
          A-B: 0.5

    >>> dom, intra, inter = read_yaml_config("domain.yaml")
    >>> dom["A"][:3]
    [1, 2, 3]
    >>> inter[("A", "B")]
    0.5
    """
    try:
        with open(filepath, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Domain configuration file not found: {filepath}")
        raise

    intra = config['intra_domains']
    # inter_domains optional: single-domain proteins have no inter-domain pairs
    inter = config.get('inter_domains', {})
    n_residues = int(config['n_residues'])

    domain_to_residues = {}
    intra_nscales = {}
    all_residues = set()
    used_legacy_key = False

    # Parse intra-domain configurations
    for domain, values in intra.items():
        raw_residues = values['residues']
        residues = parse_residue_list(raw_residues)
        domain_to_residues[domain] = residues
        # `nscale` is the current key; `strength` is the deprecated alias.
        if 'nscale' in values:
            intra_nscales[domain] = values['nscale']
        elif 'strength' in values:
            intra_nscales[domain] = values['strength']
            used_legacy_key = True
        else:
            raise KeyError(
                f"domain '{domain}' in {filepath} needs an 'nscale' value "
                f"(the native-contact scaling factor).")
        all_residues.update(residues)

    if used_legacy_key:
        warnings.warn(
            f"{filepath}: the domain.yaml key 'strength' is deprecated -- rename it to "
            f"'nscale' (the native-contact scaling factor). 'strength' is still accepted "
            f"for now but may be removed in a future release.",
            DeprecationWarning, stacklevel=2)
        print(f"[deprecation] {filepath}: 'strength' -> use 'nscale' (still works for now).")

    # Handle unassigned residues
    full_residues = set(range(1, n_residues + 1))
    unassigned_residues = sorted(full_residues - all_residues)
    if unassigned_residues:
        domain_to_residues['X'] = unassigned_residues
        intra_nscales['X'] = 1.0

    # Parse inter-domain configurations
    inter_nscales = {}
    for pair_str, nscale in inter.items():
        d1, d2 = pair_str.strip().split('-')
        inter_nscales[(d1, d2)] = nscale
        inter_nscales[(d2, d1)] = nscale  # ensure symmetry

    # Add inter-domain interactions for domain X
    if 'X' in domain_to_residues:
        for other in domain_to_residues:
            if other != 'X':
                inter_nscales[('X', other)] = 1.0
                inter_nscales[(other, 'X')] = 1.0

    return domain_to_residues, intra_nscales, inter_nscales

def get_scaling_ss_matrix(domain_def: str) -> np.ndarray:
    """
    Build scaling matrix for sidechain–sidechain energies by domain.

    Reads domain definitions from YAML and builds an (n × n) matrix where
    entry [i, j] is the intra-domain nscale if residues i and j are in
    the same domain, or the inter-domain nscale if they are in different
    domains (defaulting to 1.0 — no scaling — when no inter nscale is
    defined for that domain pair).

    Parameters
    ----------
    domain_def : str
        Path to domain definition YAML (see :func:`read_yaml_config`).

    Returns
    -------
    np.ndarray, shape (n, n)
        Symmetric matrix. Rows/columns are ordered by sorted residue list
        (all residues that appear in domain_to_residues). Values are floats
        (typically 0.0 to 1.0) used to scale SS contact energies.
    """
    domain_to_residues, intra_nscales, inter_nscales = read_yaml_config(domain_def)
    
    # Build residue to domain mapping
    residue_to_domain = {}
    residue_list = []
    for domain, residues in domain_to_residues.items():
        for res in residues:
            residue_to_domain[res] = domain
            residue_list.append(res)
    
    residue_list = sorted(set(residue_list))
    res_to_idx = {res: i for i, res in enumerate(residue_list)}
    n = len(residue_list)
    matrix = np.zeros((n, n))
    
    # Fill scaling matrix
    for i_res in residue_list:
        i_idx = res_to_idx[i_res]
        dom_i = residue_to_domain[i_res]
        for j_res in residue_list:
            j_idx = res_to_idx[j_res]
            dom_j = residue_to_domain[j_res]
            
            if dom_i == dom_j:
                matrix[i_idx, j_idx] = intra_nscales[dom_i]
            else:
                key = (dom_i, dom_j)
                if key in inter_nscales:
                    matrix[i_idx, j_idx] = inter_nscales[key]
                elif (dom_j, dom_i) in inter_nscales:
                    matrix[i_idx, j_idx] = inter_nscales[(dom_j, dom_i)]
                else:
                    # Default to 1.0 (identity = no scaling), not 0.0. The scaling
                    # matrix only modulates the nscale of contacts that already
                    # exist in the native contact map; an unspecified inter-domain
                    # pair should leave those native contacts unchanged rather than
                    # silently delete them. This matches the neutral default used
                    # for single-domain systems and for the auto 'X' domain. To
                    # decouple two domains, set their pair to 0.0 explicitly.
                    matrix[i_idx, j_idx] = 1.0
    
    return matrix

def calculate_rmin_2_values(binary_contact_matrix: np.ndarray, ca_distances: np.ndarray,
                           n_residues: int) -> List[float]:
    """
    Compute per-residue collision radius Rmin/2 for non-native contacts.

    For residue i, ``Rmin/2[i] = 0.5 * RMIN_SCALE_FACTOR * (minimum CA-CA distance to
    residues that are (1) not in contact with i and (2) not within LOCAL_SEPARATION in
    sequence)``. The minimum CA-CA distance is taken as the LJ sigma (closest approach);
    ``RMIN_SCALE_FACTOR = 2^(1/6)`` scales it to the full collision *diameter* Rmin, and
    the ``0.5`` gives **Rmin/2** -- the per-residue collision *radius* (structure-derived
    Karanicolas-Brooks). If no such residue j exists, ``Rmin/2[i] = 0.0``. Non-native pairs
    combine by the SUM rule ``Rmin_ij = Rmin/2[i] + Rmin/2[j]``.

    Parameters
    ----------
    binary_contact_matrix : np.ndarray, shape (n_residues, n_residues)
         Binary contact matrix (1 = in contact, 0 = not).
    ca_distances : np.ndarray, shape (n_residues, n_residues)
        CA–CA distances (same units as used later; typically Angstrom).
    n_residues : int
        Number of residues.

    Returns
    -------
    list of float
        ``Rmin/2[i]`` per residue i; non-native well position via the sum rule
        ``rmin_matrix[i, j] = Rmin/2[i] + Rmin/2[j]``.
    """
    rmin_2 = []
    for i in range(n_residues):
        not_in_contact_with_i = [
            j for j in range(n_residues) 
            if abs(i - j) > LOCAL_SEPARATION and binary_contact_matrix[i, j] == 0
        ]
        if not_in_contact_with_i:
            distance_to_i = ca_distances[i, not_in_contact_with_i]
            rmin_2.append(0.5 * RMIN_SCALE_FACTOR * np.min(distance_to_i))
        else:
            rmin_2.append(0.0)  # fallback value
    return rmin_2

def build_nonbonded_interaction(
    pdb_file: str,
    domain_def: Optional[str] = None,
    stride_output_file: Optional[str] = None,
    return_rmin_2: bool = False,
):
    """
    Build the well-position (Rmin) and energy matrices for TOPO non-bonded contacts.

    Combines hydrogen bonds (STRIDE), backbone–sidechain and sidechain–sidechain
    contacts (distance cutoffs), and domain-based scaling. Native contacts get their
    well at the measured CA–CA distance with energies from H-bond, BS, and scaled SS
    terms; non-native pairs get their well at the sum-rule collision distance
    ``Rmin/2_i + Rmin/2_j`` and a small repulsive energy (non_native). Output matrices
    are (n_residues × n_residues), symmetric, in nm and kJ/mol for use with OpenMM.
    With ``return_rmin_2=True`` also returns the per-residue Rmin/2 array (nm).

    Parameters
    ----------
    pdb_file : str
        Path to PDB structure (all-atom; used for residue list and CA distances).
    domain_def : str, optional
        Path to domain YAML (intra_domains, n_residues, optional inter_domains).
        If None, all residues are treated as a single domain (scaling 1.0 everywhere).
    stride_output_file : str, optional
        Path to STRIDE output for hydrogen bond list. If None, the function checks
        whether the ``stride`` program is available; if yes, it runs
        ``stride -h pdb_file``, caches the output to ``{prefix}_stride.dat`` (prefix
        taken from the PDB filename), and parses that file. STRIDE success is judged
        by output content (presence of DNR records), not the process exit code,
        since some STRIDE builds return non-zero even on success. If stride is not
        found, raises an error (provide a precomputed STRIDE file or install STRIDE).

    Returns
    -------
    rmin_matrix : np.ndarray, shape (n_residues, n_residues)
        Pairwise well position Rmin (nm). Native: CA–CA distance; non-native:
        Rmin/2_i + Rmin/2_j (sum rule).
    energy_matrix : np.ndarray, shape (n_residues, n_residues)
        Pairwise well depth in kJ/mol. Native: sum of H-bond (0.75/1.5 kcal/mol),
        backbone–sidechain (0.37 kcal/mol), and scaled SS; non-native: ENERGY_PARAMS['non_native'].

    Raises
    ------
    FileNotFoundError
        If pdb_file (or stride_output_file when given) cannot be opened.
    RuntimeError
        If stride_output_file is None and the ``stride`` program is not found or not executable.

    Example
    -------
    >>> dist, energy = build_nonbonded_interaction("2ww4.pdb", "domain.yaml", "stride.dat")
    >>> dist, energy = build_nonbonded_interaction("2ww4.pdb")  # single domain, run stride if available
    >>> dist.shape
    (283, 283)
    """
    print("Loading protein structure...")
    # A placeholder "1 A^3" CRYST1 record is common in CG/processed PDBs and makes
    # MDAnalysis warn that the unit cell is set to None. The box is irrelevant here
    # (only coordinates/topology are used), so the warning is safe to silence.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*CRYST1 record.*",
            category=UserWarning,
        )
        u = mda.Universe(pdb_file)
    resid_to_index, index_to_resname, n_residues = get_residue_mapping(u)

    # Resolve STRIDE output: use the file if provided; otherwise, if the `stride`
    # program is available, run it and cache the output to "{prefix}_stride.dat"
    # (prefix taken from the PDB filename, e.g. 1AKE_A.pdb -> 1AKE_A_stride.dat),
    # then parse that file.
    stride_path = stride_output_file
    if stride_path is None:
        try:
            stride_exe = find_executable("stride")
        except RuntimeError as exc:
            raise RuntimeError(
                "stride_output_file was not supplied and STRIDE could not be located. "
                "Either provide a precomputed STRIDE output file (stride_output_file=...) "
                "or make STRIDE available. " + str(exc)
            ) from None
        prefix = os.path.splitext(os.path.basename(pdb_file))[0]
        stride_path = "{}_stride.dat".format(prefix)
        print("Running STRIDE on structure (stride -h {} -> {}).".format(pdb_file, stride_path))
        try:
            result = subprocess.run(
                [stride_exe, "-h", pdb_file],
                capture_output=True,
                text=True,
                timeout=60,
            )
        except FileNotFoundError:
            raise RuntimeError(
                "stride executable was found but could not be run. "
                "Ensure stride is executable and on PATH."
            ) from None
        # NOTE: some STRIDE builds return a non-zero exit code even on success, so we
        # validate by output content (presence of DNR hydrogen-bond records, which
        # parse_hydrogen_bonds keys on) rather than trusting the process return code.
        if "DNR" not in result.stdout:
            raise RuntimeError(
                "STRIDE produced no hydrogen-bond (DNR) records for {} (exit code {}). "
                "stderr: {}".format(pdb_file, result.returncode, result.stderr)
            )
        with open(stride_path, "w") as f:
            f.write(result.stdout)

    hb_contact_matrix = get_hb_contact_matrix(stride_path, resid_to_index, n_residues)
    # Matrix values: 0, 1, or 2 only (2+ H-bonds capped; energy 0.75 kcal/mol per single, 1.5 for multiple)
    hb_interaction_energy = ENERGY_PARAMS['hydrogen_bond'] * hb_contact_matrix

    bs_contact_matrix = get_bs_contact_matrix(u, cutoff=DEFAULT_CUTOFF)
    # Matrix values: 0, 1, or 2 — the number of *directions* realized (bb_i–sc_j and
    # bb_j–sc_i are independent contacts), not an H-bond-style cap. Energy is linear
    # in that count: 0.37 kcal/mol per direction. See get_bs_contact_matrix.
    bs_interaction_energy = bs_contact_matrix * ENERGY_PARAMS['backbone_sidechain']

    if domain_def is None:
        # Single domain: all residues scaled by 1.0
        scaling_matrix = np.ones((n_residues, n_residues))
    else:
        scaling_matrix = get_scaling_ss_matrix(domain_def)
    ss_contact_matrix = get_ss_contact_matrix(u, cutoff=DEFAULT_CUTOFF)
    ss_interaction_energy = get_ss_interaction_energy(u)
    
    # Element-wise multiplication
    scaled_ss_interaction_energy = scaling_matrix * ss_contact_matrix * ss_interaction_energy
    
    # Total interaction energy for native contacts
    eps_ij = hb_interaction_energy + bs_interaction_energy + scaled_ss_interaction_energy
    
    # Build binary contact matrix
    contact_matrix = hb_contact_matrix + bs_contact_matrix + ss_contact_matrix
    binary_contact_matrix = (contact_matrix > 0).astype(int)
    
    # Build the per-pair well-position (Rmin) matrix: native distance for contacts,
    # sum-rule collision distance for non-native pairs.
    ca_atoms = u.select_atoms('protein and name CA')
    ca_distances = distance_array(ca_atoms, ca_atoms)
    rmin_matrix = np.zeros_like(ca_distances)

    # Native contacts: well position = the measured native CA-CA distance (Go).
    contact_mask = binary_contact_matrix == 1
    rmin_matrix[contact_mask] = ca_distances[contact_mask]

    # Non-native pairs: per-residue Rmin/2 (structure-derived), combined by the SUM rule.
    rmin_2 = calculate_rmin_2_values(binary_contact_matrix, ca_distances, n_residues)
    for i in range(n_residues):
        for j in range(n_residues):
            if binary_contact_matrix[i, j] == 0:
                rmin_matrix[i, j] = rmin_2[i] + rmin_2[j]        # sum rule
                eps_ij[i, j] = ENERGY_PARAMS['non_native']

    # Convert to nm for OpenMM compatibility
    rmin_matrix /= DISTANCE_TO_NM
    
    # Report contact statistics over unique residue pairs (upper triangle, i < j).
    # Native contacts are residue pairs flagged in the binary contact matrix;
    # all remaining pairs interact through the non-native (excluded-volume) term.
    n_native = int(np.triu(binary_contact_matrix, k=1).sum())
    n_pairs = n_residues * (n_residues - 1) // 2
    n_non_native = n_pairs - n_native
    print(f"  native contacts: {n_native}  |  non-native pairs: {n_non_native}  "
          f"(of {n_pairs} residue pairs)")
    if return_rmin_2:
        # Per-residue Karanicolas-Brooks collision radius Rmin/2 (nm) -- O'Brien's nascent
        # excluded-volume radius (the A1..An values in his .prm), used for the nascent side
        # of the NC<->ribosome excluded volume (see tutorials/15_claude_fix/
        # TOPO_OBrien_NCribosome_nonbonded_compare.md). calculate_rmin_2_values already
        # returns Rmin/2 (Angstrom); just convert to nm.
        rmin_2_nm = np.asarray(rmin_2, dtype=float) / DISTANCE_TO_NM   # Angstrom -> nm
        return rmin_matrix, eps_ij, rmin_2_nm
    return rmin_matrix, eps_ij

# Main execution
if __name__ == "__main__":
    try:
        # With all inputs:
        R_ij, eps_ij = build_nonbonded_interaction('2ww4.pdb', 'domain.yaml', 'stride.dat')
        # Optional: single domain, no domain YAML; STRIDE run automatically if stride in PATH:
        # R_ij, eps_ij = build_nonbonded_interaction('2ww4.pdb')
        print(f"Distance matrix shape: {R_ij.shape}")
        print(f"Energy matrix shape: {eps_ij.shape}")
    except Exception as e:
        print(f"Error building non-bonded interactions: {e}")
        raise