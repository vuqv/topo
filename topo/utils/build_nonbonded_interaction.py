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
3. **Domain scaling**: YAML defines intra_domains (residue lists + strength) and
   optional inter_domains (pair strengths). Single-domain proteins can omit
   inter_domains.
4. **Non-native contacts**: repulsive well with sigma from nearest non-contact
   CA-CA distance; energy ENERGY_PARAMS['non_native'].

Typical use
-----------
>>> from topo.utils.build_nonbonded_interaction import build_nonbonded_interaction
>>> distance_matrix, energy_matrix = build_nonbonded_interaction(
...     "protein.pdb", "domain.yaml", "stride.dat"
... )
>>> # distance_matrix, energy_matrix in nm and kJ/mol for OpenMM
"""
import os
import shutil
import subprocess
import tempfile
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd
import re
from pathlib import Path
from collections import defaultdict
import yaml
from typing import Dict, List, Tuple, Set, Optional
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
SIGMA_SCALE_FACTOR = 2**(1/6)  # Factor for repulsive sigma (LJ-like)
DISTANCE_TO_NM = 10.0  # Angstrom to nanometer (CA distances read in Angstrom)

# Energy parameters (kcal/mol values, stored in kJ/mol for OpenMM)
ENERGY_PARAMS = {
    'hydrogen_bond': 0.75 * KCAL_TO_KJ,       # per H-bond; 2+ H-bonds capped at 1.5 kcal/mol total
    'backbone_sidechain': 0.37 * KCAL_TO_KJ,
    'non_native': 0.000132 * KCAL_TO_KJ,
    'yang_shift': 0.6   # Used in BT potential (kcal/mol)
}


def get_residue_mapping(universe: mda.Universe) -> Tuple[Dict[int, int], Dict[int, str], int]:
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
    resid_to_index : dict
        Maps PDB residue ID (resid) -> 0-based index in residue list.
    index_to_resname : dict
        Maps 0-based index -> three-letter residue name (e.g. 'ALA', 'GLY').
    n_residues : int
        Total number of protein residues.

    Example
    -------
    >>> u = mda.Universe("protein.pdb")
    >>> resid_to_idx, idx_to_name, n = get_residue_mapping(u)
    >>> resid_to_idx[1]   # first residue's index
    0
    >>> idx_to_name[0]
    'MET'
    """
    residues = universe.select_atoms("protein").residues
    n_residues = len(residues)
    resid_to_index = {res.resid: idx for idx, res in enumerate(residues)}
    index_to_resname = {idx: res.resname for idx, res in enumerate(residues)}
    return resid_to_index, index_to_resname, n_residues

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
        Each element is a pair ``((resname1, resid1), (resname2, resid2))`` with the
        tuple sorted so the pair is canonical. Resid is the first number in the
        STRIDE line (PDB residue ID). Repeated pairs indicate multiple H-bonds
        between the same residue pair.

    Raises
    ------
    FileNotFoundError
        If the STRIDE file cannot be opened.

    Notes
    -----
    STRIDE line format: ``DNR  RES -  resid  seq ->  RES -  resid  seq  ...``
    The regex captures resname and the first resid for donor and acceptor.

    Example
    -------
    >>> pairs = parse_hydrogen_bonds("stride.dat")
    >>> pairs[0]  # e.g. ((u'THR', 3), (u'ILE', 43))
    (('THR', 3), ('ILE', 43))
    """
    try:
        with open(stride_output_file, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"STRIDE output file not found: {stride_output_file}")
        raise
    
    pattern = re.compile(
        r"(?:DNR|ACC)\s+(\w+)\s+-\s+(\d+)\s+\d+\s+->\s+(\w+)\s+-\s+(\d+)\s+\d+"
    )
    
    hb_pairs = []
    # STRIDE reports each H-bond twice (DNR and ACC). Count each bond once by using only DNR lines.
    # Then pairs with multiple physical H-bonds appear as multiple DNR lines and are counted correctly.
    for line in lines:
        if line.startswith("DNR"):
            match = pattern.search(line)
            if match:
                res1, res1_pdb, res2, res2_pdb = match.groups()
                donor = (res1, int(res1_pdb))
                acceptor = (res2, int(res2_pdb))
                pair = tuple(sorted([donor, acceptor]))
                hb_pairs.append(pair)
    return hb_pairs

def build_hb_contact_matrix(hb_pairs: List[Tuple], n_residues: int) -> np.ndarray:
    """
    Build the hydrogen-bond contact matrix from a list of H-bond pairs.

    Counts how many H-bonds exist between each residue pair and stores 0, 1, or 2
    (values greater than 2 are capped at 2). The model uses 0.75 kcal/mol per single
    H-bond and 1.5 kcal/mol total for multiple H-bonds; capping at 2 enforces that.

    Parameters
    ----------
    hb_pairs : list of tuple
        List of canonical pairs from :func:`parse_hydrogen_bonds`. Each pair is
        ``((resname1, resid1), (resname2, resid2))`` with resid as PDB residue ID.
    n_residues : int
        Number of residues (determines matrix shape).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues), dtype int
        Symmetric matrix. Entry [i, j] is 0, 1, or 2 (number of H-bonds between
        residue index i and j, capped at 2). Indices are 0-based and must match
        the residue order used elsewhere (e.g. from :func:`get_residue_mapping`);
        resids in hb_pairs are converted to 0-based index as ``resid - 1``.

    Notes
    -----
    Resid in hb_pairs is the PDB residue number (1-based). It is converted to
    matrix index as ``resid - 1``. Ensure your residue list order matches the
    PDB/STRIDE numbering when interpreting i, j.
    """
    hb_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    
    # Count hydrogen bonds between residue pairs
    pair_counts = defaultdict(int)
    for donor, acceptor in hb_pairs:
        pair = tuple(sorted([donor, acceptor]))
        pair_counts[pair] += 1
    
    # Fill contact matrix; cap at 2 (model: 1 H-bond -> 0.75 kcal/mol, 2+ H-bonds -> 1.5 kcal/mol only)
    for (res1, res2), count in pair_counts.items():
        i = res1[1] - 1  # 0-based index
        j = res2[1] - 1  # 0-based index
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
        Maps canonical pair ``((resname1, resid1), (resname2, resid2))`` to number
        of H-bonds (before any cap). Useful for logging or downstream analysis.

    Example
    -------
    >>> pairs = parse_hydrogen_bonds("stride.dat")
    >>> counts = print_pairs_with_multiple_hb(pairs)
    Residue pairs with more than 1 hydrogen bond:
      THR3 -- ILE43:  2 H-bonds
      ...
      (total: 20 pairs)
    >>> counts[(('THR', 3), ('ILE', 43))]
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
        for (r1, i1), (r2, i2) in [p for p, _ in multi]:
            c = pair_counts[tuple(sorted([(r1, i1), (r2, i2)]))]
            print(f"  {r1}{i1} -- {r2}{i2}:  {c} H-bonds")
        print(f"  (total: {len(multi)} pairs)")
    return dict(pair_counts)


def get_hb_contact_matrix(stride_output_file: str, n_residues: int) -> np.ndarray:
    """
    Build the hydrogen-bond contact matrix directly from a STRIDE output file.

    Convenience function that parses the file and builds the matrix. Values are
    0, 1, or 2 (multiple H-bonds capped at 2).

    Parameters
    ----------
    stride_output_file : str
        Path to the STRIDE output file.
    n_residues : int
        Number of residues (must match the system used for STRIDE).

    Returns
    -------
    np.ndarray, shape (n_residues, n_residues), dtype int
        Symmetric H-bond contact matrix. See :func:`build_hb_contact_matrix`.

    Example
    -------
    >>> n_residues = 283
    >>> hb_matrix = get_hb_contact_matrix("stride.dat", n_residues)
    >>> hb_matrix.shape
    (283, 283)
    >>> (hb_matrix > 1).sum() // 2   # number of pairs with 2 H-bonds (symmetric)
    20
    """
    hb_pairs = parse_hydrogen_bonds(stride_output_file)
    return build_hb_contact_matrix(hb_pairs, n_residues)

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
    excluding H. Pairs with |resid_i - resid_j| <= LOCAL_SEPARATION are excluded.
    """
    resid_to_index, _, n_residues = get_residue_mapping(u)
    
    backbone = u.select_atoms('protein and backbone and not name H*')
    sidechain = u.select_atoms("protein and not backbone and not name H*")
    
    dists_bs = distance_array(backbone.positions, sidechain.positions)
    
    # Build directional residue-residue contact list
    contacts_bs = set()
    for i, atom1 in enumerate(backbone):
        for j, atom2 in enumerate(sidechain):
            if (dists_bs[i, j] <= cutoff and 
                abs(atom1.resid - atom2.resid) > LOCAL_SEPARATION):
                pair = (atom1.resid, atom2.resid)  # preserve direction: bb → sc
                contacts_bs.add(pair)
    
    # Build asymmetric matrix first
    bs_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    for contact in sorted(contacts_bs, key=lambda x: x[0]):
        bs_contact_matrix[resid_to_index[contact[0]], 
                        resid_to_index[contact[1]]] = 1
    
    # Build symmetric count matrix
    bs_symmetric_count = np.zeros((n_residues, n_residues), dtype=int)
    for i in range(n_residues):
        for j in range(n_residues):
            if abs(i - j) > LOCAL_SEPARATION:
                count = bs_contact_matrix[i, j] + bs_contact_matrix[j, i]
                bs_symmetric_count[i, j] = count
                bs_symmetric_count[j, i] = count
    
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
    resid_to_index, _, n_residues = get_residue_mapping(u)
    sidechain = u.select_atoms("protein and not backbone and not name H*")
    
    dists_ss = distance_array(sidechain.positions, sidechain.positions)
    ss_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    
    # Build residue-residue contact list
    contacts_ss = set()
    for i, atom1 in enumerate(sidechain):
        for j, atom2 in enumerate(sidechain):
            if (dists_ss[i, j] <= cutoff and 
                abs(atom1.resid - atom2.resid) > LOCAL_SEPARATION):
                pair = tuple(sorted([atom1.resid, atom2.resid]))
                contacts_ss.add(pair)
    
    # Fill contact matrix
    for contact in sorted(contacts_ss, key=lambda x: x[0]):
        ss_contact_matrix[resid_to_index[contact[0]], 
                        resid_to_index[contact[1]]] = 1
        ss_contact_matrix[resid_to_index[contact[1]], 
                        resid_to_index[contact[0]]] = 1
    
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
        computed as ``KCAL_TO_KJ * |raw_value - yang_shift|`` from the CSV.

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
        return KCAL_TO_KJ * np.abs(df - ENERGY_PARAMS['yang_shift'])
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
    Read and parse domain definition YAML (intra/inter strengths, residue lists).

    Required keys: ``intra_domains``, ``n_residues``. Optional: ``inter_domains``
    (omit for single-domain proteins; then inter_strengths will be empty).
    Residues not listed in any domain are assigned to domain ``'X'`` with
    intra strength 1.0 and inter 1.0 to all other domains.

    Parameters
    ----------
    filepath : str
        Path to the YAML file (e.g. ``domain.yaml``).

    Returns
    -------
    domain_to_residues : dict
        Domain name -> list of residue numbers (1-based).
    intra_strengths : dict
        Domain name -> float (intra-domain contact strength).
    inter_strengths : dict
        (domain1, domain2) -> float (inter-domain strength); symmetric keys
        (d1, d2) and (d2, d1) are both set.

    Raises
    ------
    FileNotFoundError
        If the YAML file cannot be opened.

    Example
    -------
    YAML format::

        n_residues: 110
        intra_domains:
          A: { residues: [1-50], strength: 1.0 }
          B: { residues: [51-110], strength: 1.0 }
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
    intra_strengths = {}
    all_residues = set()
    
    # Parse intra-domain configurations
    for domain, values in intra.items():
        raw_residues = values['residues']
        residues = parse_residue_list(raw_residues)
        domain_to_residues[domain] = residues
        intra_strengths[domain] = values['strength']
        all_residues.update(residues)
    
    # Handle unassigned residues
    full_residues = set(range(1, n_residues + 1))
    unassigned_residues = sorted(full_residues - all_residues)
    if unassigned_residues:
        domain_to_residues['X'] = unassigned_residues
        intra_strengths['X'] = 1.0
    
    # Parse inter-domain configurations
    inter_strengths = {}
    for pair_str, strength in inter.items():
        d1, d2 = pair_str.strip().split('-')
        inter_strengths[(d1, d2)] = strength
        inter_strengths[(d2, d1)] = strength  # ensure symmetry
    
    # Add inter-domain interactions for domain X
    if 'X' in domain_to_residues:
        for other in domain_to_residues:
            if other != 'X':
                inter_strengths[('X', other)] = 1.0
                inter_strengths[(other, 'X')] = 1.0
    
    return domain_to_residues, intra_strengths, inter_strengths

def get_scaling_ss_matrix(domain_def: str) -> np.ndarray:
    """
    Build scaling matrix for sidechain–sidechain energies by domain.

    Reads domain definitions from YAML and builds an (n × n) matrix where
    entry [i, j] is the intra-domain strength if residues i and j are in
    the same domain, or the inter-domain strength if they are in different
    domains (or 0 if no inter strength is defined).

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
    domain_to_residues, intra_strengths, inter_strengths = read_yaml_config(domain_def)
    
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
                matrix[i_idx, j_idx] = intra_strengths[dom_i]
            else:
                key = (dom_i, dom_j)
                if key in inter_strengths:
                    matrix[i_idx, j_idx] = inter_strengths[key]
                elif (dom_j, dom_i) in inter_strengths:
                    matrix[i_idx, j_idx] = inter_strengths[(dom_j, dom_i)]
                else:
                    matrix[i_idx, j_idx] = 0.0
    
    return matrix

def calculate_sigma_values(binary_contact_matrix: np.ndarray, ca_distances: np.ndarray,
                          n_residues: int) -> List[float]:
    """
    Compute repulsive sigma (distance) for each residue for non-native contacts.

    For residue i, sigma[i] is set to SIGMA_SCALE_FACTOR times the minimum
    CA–CA distance to residues that are (1) not in contact with i (binary_contact_matrix[i, j] == 0)
    and (2) not within LOCAL_SEPARATION in sequence. This defines a soft repulsion
    distance for non-native pairs. If no such residue j exists, sigma[i] = 0.0.

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
        sigma[i] for each residue i. Used to set distance_matrix and repulsive
        well for non-native contacts (e.g. 0.5 * (sigma[i] + sigma[j])).
    """
    sigma = []
    for i in range(n_residues):
        not_in_contact_with_i = [
            j for j in range(n_residues) 
            if abs(i - j) > LOCAL_SEPARATION and binary_contact_matrix[i, j] == 0
        ]
        if not_in_contact_with_i:
            distance_to_i = ca_distances[i, not_in_contact_with_i]
            sigma.append(SIGMA_SCALE_FACTOR * np.min(distance_to_i))
        else:
            sigma.append(0.0)  # fallback value
    return sigma

def build_nonbonded_interaction(
    pdb_file: str,
    domain_def: Optional[str] = None,
    stride_output_file: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build distance and energy matrices for TOPO non-bonded (native + repulsive) contacts.

    Combines hydrogen bonds (STRIDE), backbone–sidechain and sidechain–sidechain
    contacts (distance cutoffs), and domain-based scaling. Native contacts get
    CA–CA distances and energies from H-bond, BS, and scaled SS terms; non-native
    pairs get sigma-based distances and a small repulsive energy (non_native).
    Output matrices are (n_residues × n_residues), symmetric, in nm and kJ/mol
    for use with OpenMM.

    Parameters
    ----------
    pdb_file : str
        Path to PDB structure (all-atom; used for residue list and CA distances).
    domain_def : str, optional
        Path to domain YAML (intra_domains, n_residues, optional inter_domains).
        If None, all residues are treated as a single domain (scaling 1.0 everywhere).
    stride_output_file : str, optional
        Path to STRIDE output for hydrogen bond list. If None, the function checks
        whether the ``stride`` program is available and executable; if yes, runs
        ``stride -h pdb_file`` and uses that output. If stride is not found, raises
        an error (provide a precomputed STRIDE file or install STRIDE).

    Returns
    -------
    distance_matrix : np.ndarray, shape (n_residues, n_residues)
        Pairwise distance in nm. Native: CA–CA distance; non-native: 0.5 * (sigma_i + sigma_j).
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
    u = mda.Universe(pdb_file)
    resid_to_index, index_to_resname, n_residues = get_residue_mapping(u)

    # Resolve STRIDE output: use file if given, else run stride -h pdb_file if executable
    stride_path = stride_output_file
    temp_stride_path = None
    if stride_path is None:
        stride_exe = shutil.which("stride")
        if stride_exe is None:
            raise RuntimeError(
                "stride_output_file was not supplied and the 'stride' program was not found. "
                "Either provide a precomputed STRIDE output file (stride_output_file=...) or "
                "install STRIDE and ensure it is on PATH."
            )
        print("Running STRIDE on structure (stride -h {}).".format(pdb_file))
        try:
            result = subprocess.run(
                [stride_exe, "-h", pdb_file],
                capture_output=True,
                text=True,
                check=True,
                timeout=60,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                "STRIDE failed: {} {}".format(e.stderr or "", e.stdout or "")
            ) from e
        except FileNotFoundError:
            raise RuntimeError(
                "stride executable was found but could not be run. "
                "Ensure stride is executable and on PATH."
            ) from None
        fd, temp_stride_path = tempfile.mkstemp(suffix=".stride.dat", text=True)
        try:
            with os.fdopen(fd, "w") as f:
                f.write(result.stdout)
            stride_path = temp_stride_path
        except Exception:
            if temp_stride_path is not None:
                os.unlink(temp_stride_path)
            raise

    print("Building hydrogen bond contact matrix...")
    try:
        hb_contact_matrix = get_hb_contact_matrix(stride_path, n_residues)
    finally:
        if temp_stride_path is not None:
            try:
                os.unlink(temp_stride_path)
            except OSError:
                pass
    # Matrix values: 0, 1, or 2 only (2+ H-bonds capped; energy 0.75 kcal/mol per single, 1.5 for multiple)
    hb_interaction_energy = ENERGY_PARAMS['hydrogen_bond'] * hb_contact_matrix

    print("Building backbone-sidechain contact matrix...")
    bs_contact_matrix = get_bs_contact_matrix(u, cutoff=DEFAULT_CUTOFF)
    bs_interaction_energy = bs_contact_matrix * ENERGY_PARAMS['backbone_sidechain']

    print("Building sidechain-sidechain contact matrix...")
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
    
    # Calculate distance matrix
    ca_atoms = u.select_atoms('protein and name CA')
    ca_distances = distance_array(ca_atoms, ca_atoms)
    distance_matrix = np.zeros_like(ca_distances)
    
    # Set distances for native contacts
    contact_mask = binary_contact_matrix == 1
    distance_matrix[contact_mask] = ca_distances[contact_mask]
    
    # Calculate sigma values for non-native contacts
    sigma = calculate_sigma_values(binary_contact_matrix, ca_distances, n_residues)
    
    # Set distances and energies for non-native contacts
    for i in range(n_residues):
        for j in range(n_residues):
            if binary_contact_matrix[i, j] == 0:
                distance_matrix[i, j] = 0.5 * (sigma[i] + sigma[j])
                eps_ij[i, j] = ENERGY_PARAMS['non_native']
    
    # Convert to nm for OpenMM compatibility
    distance_matrix /= DISTANCE_TO_NM
    
    print("Non-bonded interaction matrices built successfully")
    return distance_matrix, eps_ij

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