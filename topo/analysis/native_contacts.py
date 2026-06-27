#!/usr/bin/env python3
"""
Per-domain and per-interface fraction of native contacts (Q).

Computes the fraction of native contacts Q for every domain (intra-domain) and
every domain-domain interface (inter-domain) declared in a ``domain.yaml``, for
each frame of a coarse-grained (CG, Cα-only) trajectory.

Method
------
1. Native contacts are defined ONCE from the all-atom reference structure:
   residues i, j form a native contact when they are more than
   ``local_separation`` residues apart in sequence (abs(i - j) > 3 by default) AND
   at least one pair of heavy atoms (one per residue) is within ``cutoff``
   (4.5 Å by default). The reference Cα–Cα distance of each contact is stored.
2. Contacts are grouped by domain using the residue sets in ``domain.yaml``:
   - intra-domain  : both residues in the same domain
   - interface     : one residue in each of two domains
3. For every CG frame, a native contact (i, j) is counted as *formed* when the
   instantaneous Cα–Cα distance satisfies d_CG(i, j) <= tolerance * d_native(i, j)
   (tolerance = 1.2). Q is the fraction of that group's native contacts formed.

Use as a library (e.g. from the strength optimizer) by importing the functions,
or as a CLI:

    python -m topo.analysis.native_contacts -d domain.yaml -r ref.pdb \\
        -p cg.psf -f cg.dcd -o Q.csv

Residue numbering: ``domain.yaml`` uses 1-based residue numbers that run
1..n_residues in CG bead order; internally these map to 0-based indices (n - 1),
which is also the Cα/bead order of the reference and CG models.
"""

from pathlib import Path
import argparse
import itertools
import sys
import warnings

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_capped_distance

# Reuse the package's domain.yaml residue-list parser (no duplicate copy).
from ..utils.nonbonded import parse_residue_list


def load_universe(*args, **kwargs):
    """``mda.Universe`` with the harmless placeholder-CRYST1 warning suppressed.

    Some PDBs carry a dummy ``CRYST1  1 1 1`` record; MDAnalysis then warns
    "1 A^3 CRYST1 record, this is usually a placeholder. Unit cell dimensions
    will be set to None." That is expected for these non-periodic CG inputs, so
    we silence only that one UserWarning and let any others through.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*CRYST1 record",
                                category=UserWarning)
        return mda.Universe(*args, **kwargs)


# --------------------------------------------------------------------------- #
# domain.yaml parsing
# --------------------------------------------------------------------------- #
def load_domains(domain_yaml, n_res):
    """Read domains from domain.yaml as 0-based residue-index sets.

    Returns
    -------
    domains : dict[str, np.ndarray]
        Domain name -> sorted array of 0-based residue indices.
    n_residues : int
        ``n_residues`` declared in the file (validated against the structure).
    """
    import yaml  # local import: only needed here

    with open(domain_yaml) as fh:
        cfg = yaml.safe_load(fh)

    if "intra_domains" not in cfg:
        raise ValueError(f"{domain_yaml}: missing required 'intra_domains' section")

    n_residues = int(cfg.get("n_residues", n_res))

    domains = {}
    assigned = set()
    for name, values in cfg["intra_domains"].items():
        resnums = parse_residue_list(values["residues"])     # 1-based
        idx = np.array(sorted(r - 1 for r in resnums), dtype=int)  # -> 0-based
        if idx.min() < 0 or idx.max() >= n_residues:
            raise ValueError(
                f"Domain {name}: residue numbers must be within 1..{n_residues}"
            )
        domains[name] = idx
        assigned.update(idx.tolist())

    unassigned = set(range(n_residues)) - assigned
    if unassigned:
        print(
            f"[warning] {len(unassigned)} residue(s) are not assigned to any "
            f"domain; they contribute to Q_protein but to no per-domain Q.",
            file=sys.stderr,
        )

    return domains, n_residues


# --------------------------------------------------------------------------- #
# Native contacts from the all-atom reference
# --------------------------------------------------------------------------- #
def reference_residue_geometry(ref_universe):
    """Return reference Cα positions (0-based residue order) and a mapping from
    each heavy atom to its 0-based protein-residue index.

    Building an explicit ``universe-resindex -> 0-based position`` map makes the
    indexing robust to hetero residues or non-sequential numbering.
    """
    protein_residues = ref_universe.select_atoms("protein").residues
    resindex_to_pos = {ri: pos for pos, ri in enumerate(protein_residues.resindices)}

    ca = ref_universe.select_atoms("protein and name CA")
    if ca.n_atoms != len(protein_residues):
        raise ValueError(
            f"Reference has {len(protein_residues)} protein residues but "
            f"{ca.n_atoms} Cα atoms; need exactly one Cα per residue."
        )
    ca_positions = ca.positions.copy()  # row k == residue k (0-based)
    return ca_positions, resindex_to_pos


def build_native_contacts(group1, group2, heavy_positions, heavy_res,
                          ca_positions, cutoff, local_separation):
    """Build the native-contact list for one domain or interface.

    Parameters
    ----------
    group1, group2 : set[int] or None
        0-based residue-index sets:
        - group1 is None, group2 is None -> all residues (whole protein)
        - group1 set,    group2 is None  -> intra-domain (both in group1)
        - group1 set,    group2 set      -> interface (one in each)
    heavy_positions : (A, 3) array
        Coordinates of all reference heavy (non-H) atoms.
    heavy_res : (A,) int array
        0-based residue index of each heavy atom.
    ca_positions : (N, 3) array
        Reference Cα coordinates in 0-based residue order.
    cutoff : float
        Heavy-atom distance defining a contact (Å).
    local_separation : int
        Minimum sequence separation (abs(i - j) > local_separation).

    Returns
    -------
    pairs : (M, 2) int array
        Native-contact residue index pairs (i < j), 0-based.
    native_dist : (M,) float array
        Reference Cα–Cα distance for each contact.
    """
    # All heavy-atom pairs within the cutoff (efficient neighbour search).
    atom_pairs = self_capped_distance(
        heavy_positions, max_cutoff=cutoff, box=None, return_distances=False
    )

    ri = heavy_res[atom_pairs[:, 0]]
    rj = heavy_res[atom_pairs[:, 1]]

    # Sequence-separation filter (also removes intra-residue atom pairs).
    keep = np.abs(ri - rj) > local_separation
    ri, rj = ri[keep], rj[keep]

    # Group membership filter.
    if group1 is None and group2 is None:
        allowed = np.ones(ri.shape, dtype=bool)
    elif group2 is None:
        in1_i = np.isin(ri, list(group1))
        in1_j = np.isin(rj, list(group1))
        allowed = in1_i & in1_j
    else:
        g1, g2 = list(group1), list(group2)
        i_in1, i_in2 = np.isin(ri, g1), np.isin(ri, g2)
        j_in1, j_in2 = np.isin(rj, g1), np.isin(rj, g2)
        allowed = (i_in1 & j_in2) | (i_in2 & j_in1)

    ri, rj = ri[allowed], rj[allowed]

    # Collapse atom-level pairs to unique residue-level contacts (i < j).
    lo = np.minimum(ri, rj)
    hi = np.maximum(ri, rj)
    pairs = np.unique(np.stack([lo, hi], axis=1), axis=0)

    if pairs.shape[0] == 0:
        return pairs.reshape(0, 2), np.zeros(0, dtype=float)

    native_dist = np.linalg.norm(
        ca_positions[pairs[:, 0]] - ca_positions[pairs[:, 1]], axis=1
    )
    return pairs, native_dist


# --------------------------------------------------------------------------- #
# Q per frame
# --------------------------------------------------------------------------- #
def fraction_native_contacts(ca_positions, pairs, native_dist, tolerance):
    """Q for one CG frame: fraction of native contacts with
    d_CG <= tolerance * d_native. Returns NaN if the group has no contacts.
    """
    if pairs.shape[0] == 0:
        return np.nan
    d_cg = np.linalg.norm(
        ca_positions[pairs[:, 0]] - ca_positions[pairs[:, 1]], axis=1
    )
    formed = d_cg <= tolerance * native_dist
    return float(np.count_nonzero(formed)) / pairs.shape[0]


# --------------------------------------------------------------------------- #
# CLI / driver
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(
        prog="python -m topo.analysis.native_contacts",
        description="Per-domain and per-interface fraction of native contacts (Q) "
                    "for a CG trajectory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-d", "--domain", required=True, help="domain.yaml")
    p.add_argument("-r", "--reference", required=True,
                   help="All-atom reference PDB (defines native contacts).")
    p.add_argument("-p", "--psf", required=True, help="CG PSF topology.")
    p.add_argument("-f", "--dcd", required=True, help="CG DCD trajectory.")
    p.add_argument("-o", "--output", required=True, help="Output CSV path.")
    p.add_argument("--cutoff", type=float, default=4.5,
                   help="Heavy-atom distance defining a native contact (Å).")
    p.add_argument("--local-separation", type=int, default=3,
                   help="Minimum sequence separation (abs(i - j) > this).")
    p.add_argument("--tolerance", type=float, default=1.2,
                   help="Cα-distance stretch factor for a 'formed' contact.")
    p.add_argument("--start-frame", type=int, default=0,
                   help="First frame to score (inclusive).")
    p.add_argument("--end-frame", type=int, default=None,
                   help="Last frame to score (exclusive). Default: all frames.")
    return p.parse_args()


def main():
    args = parse_args()

    for label, path in (("reference", args.reference), ("psf", args.psf),
                        ("dcd", args.dcd), ("domain", args.domain)):
        if not Path(path).exists():
            raise FileNotFoundError(f"Missing {label} file: {path}")

    # --- reference structure & native-contact geometry ---
    u_ref = load_universe(args.reference)
    ca_positions, resindex_to_pos = reference_residue_geometry(u_ref)
    n_res = ca_positions.shape[0]

    heavy = u_ref.select_atoms("protein and not name H*")
    heavy_positions = heavy.positions.copy()
    # Map each heavy atom's universe resindex to a 0-based protein-residue index.
    heavy_res = np.fromiter(
        (resindex_to_pos[ri] for ri in heavy.resindices), dtype=int,
        count=heavy.n_atoms,
    )

    # --- domains & interfaces ---
    domains, n_residues = load_domains(args.domain, n_res)
    if n_residues != n_res:
        raise ValueError(
            f"domain.yaml n_residues={n_residues} but reference has {n_res} "
            f"protein residues."
        )
    interfaces = list(itertools.combinations(sorted(domains), 2))

    # --- build the native-contact list for every unit (once) ---
    print("=== Native contacts from all-atom reference ===")
    print(f"Definition: heavy-atom distance <= {args.cutoff} Å; "
          f"sequence separation > {args.local_separation}")

    units = {}  # column name -> (pairs, native_dist)

    pairs, dnat = build_native_contacts(None, None, heavy_positions, heavy_res,
                                  ca_positions, args.cutoff, args.local_separation)
    units["Q_protein"] = (pairs, dnat)
    print(f"  whole protein           : {pairs.shape[0]} native contacts")

    for name, idx in domains.items():
        pairs, dnat = build_native_contacts(set(idx.tolist()), None, heavy_positions,
                                      heavy_res, ca_positions, args.cutoff,
                                      args.local_separation)
        units[f"Q_{name}"] = (pairs, dnat)
        print(f"  domain {name:<16}: {pairs.shape[0]} native contacts")

    for a, b in interfaces:
        pairs, dnat = build_native_contacts(set(domains[a].tolist()),
                                      set(domains[b].tolist()), heavy_positions,
                                      heavy_res, ca_positions, args.cutoff,
                                      args.local_separation)
        units[f"Q_{a}-{b}"] = (pairs, dnat)
        note = "" if pairs.shape[0] else "  (no contacts -> Q = NaN)"
        print(f"  interface {a}-{b:<11}: {pairs.shape[0]} native contacts{note}")

    # --- CG trajectory ---
    # MDAnalysis emits a DeprecationWarning announcing a v3.0 change to how the
    # DCDReader handles timesteps. It does not affect the coordinates we read, so
    # silence it to keep the scoring output clean.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*DCDReader currently makes independent timesteps.*",
            category=DeprecationWarning,
        )
        u_cg = mda.Universe(args.psf, args.dcd)
    if u_cg.atoms.n_atoms != n_res:
        raise ValueError(
            f"CG model has {u_cg.atoms.n_atoms} beads but reference has {n_res} "
            f"residues; the CG model must be one Cα bead per residue."
        )

    n_frames_total = u_cg.trajectory.n_frames
    start = args.start_frame
    end = n_frames_total if args.end_frame is None else args.end_frame
    if not (0 <= start < end <= n_frames_total):
        raise ValueError(
            f"Invalid frame range [{start}, {end}) for {n_frames_total} frames."
        )

    print(f"\n=== Scoring {end - start} frames "
          f"(of {n_frames_total}) for {len(units)} Q series ===")

    cg_atoms = u_cg.atoms
    records = {name: np.empty(end - start, dtype=float) for name in units}
    frames = np.arange(start, end, dtype=int)

    for k, _ts in enumerate(u_cg.trajectory[start:end]):
        pos = cg_atoms.positions
        for name, (pairs, dnat) in units.items():
            records[name][k] = fraction_native_contacts(
                pos, pairs, dnat, args.tolerance
            )

    # --- write CSV (Frame first, then Q_protein, domains, interfaces) ---
    df = pd.DataFrame({"Frame": frames, **records})
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    print("\n=== Mean Q per unit ===")
    for name in units:
        print(f"  {name:<20}: {np.nanmean(records[name]):.4f}")
    print(f"\nSaved: {out_path.resolve()}")


if __name__ == "__main__":
    main()
