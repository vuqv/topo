#!/usr/bin/env python3
"""
Coarse-grain a ribosome (or any protein + RNA structure) to the TOPO convention.

Mapping
-------
* **Protein** -> one bead per residue at the **Cα** atom. The residue name is
  **kept unchanged** (ALA, GLY, ...); the bead atom name is ``CA``. This matches
  TOPO's Cα model and the per-residue parameter lookup in
  :mod:`topo.parameters.model_parameters` (keyed by residue name).

* **RNA** -> **3 beads for pyrimidines (C, U)** and **4 beads for purines
  (A, G)**, following the O'Brien-lab representation and the ``P`` / ``R`` / ``BR``
  bead types defined in :mod:`topo.parameters.model_parameters`:

    - ``P``  : the phosphate (placed at the phosphate ``P`` atom; q = -1e),
    - ``R``  : the ribose-ring centroid (C1',C2',C3',C4',O4'),
    - ``BR`` : the centroid of **each** conjugated base ring
               (pyrimidine: 1 ring -> ``BR1``; purine: 2 rings -> ``BR1`` + ``BR2``).

  Each nucleotide is kept as **one residue** (original chain + residue number and
  the nucleotide name A/U/G/C are preserved); its 3/4 beads are the residue's
  atoms (``P``, ``R``, ``BR1``[, ``BR2``]). The **bead type** used for parameter
  lookup is the atom name with trailing digits stripped (``P``, ``R``, ``BR``).

The PDB **segID** (cols 73-76) is **preserved** on every bead. In the ribosome
inputs it labels the biological molecule each bead belongs to: ``AtR`` (A-site
tRNA), ``PtR`` (P-site tRNA), ``23S`` / ``5S`` (rRNAs), and ``L2``..``L36``
(large-subunit ribosomal proteins). This lets you select molecules downstream
(e.g. the P-site tRNA the nascent chain attaches to, or the tunnel proteins).

Anything that is neither a standard amino acid nor a standard nucleotide
(ions, water, ligands) is skipped with a warning.

The ribosome is coarse-grained **before** truncation (see ``FILES.md``): every
bead centroid is computed from the complete all-atom residue.

Usage
-----
    python -m topo.csp.cg_ribosome -i ribosome.pdb -o ribosome_cg.pdb

or as a library::

    from topo.csp.cg_ribosome import coarse_grain
    coarse_grain("ribosome.pdb", "ribosome_cg.pdb")
"""
from __future__ import annotations

import argparse
import sys
from collections import OrderedDict

# Standard amino acids (protein) and ribonucleotides. Kept local so the script
# runs standalone; these mirror topo.parameters.protein_list / nucleic_list.
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# Map common RNA residue names to a canonical A/U/G/C.
PURINES = {"A", "G", "RA", "RG", "ADE", "GUA"}
PYRIMIDINES = {"C", "U", "RC", "RU", "CYT", "URA"}

# Ring atom sets (PDB v3 names).
RIBOSE_RING = ["C1'", "C2'", "C3'", "C4'", "O4'"]
BASE_RING_6 = ["N1", "C2", "N3", "C4", "C5", "C6"]   # pyrimidine ring; purine 6-ring
PURINE_RING_5 = ["C4", "C5", "N7", "C8", "N9"]       # purine 5-ring (fused)


def _centroid(coords):
    """Compute the geometric centroid of a set of 3D coordinates.

    Parameters
    ----------
    coords : list of tuple of float
        Sequence of ``(x, y, z)`` coordinates. May be empty.

    Returns
    -------
    tuple of float or None
        The component-wise mean ``(x, y, z)`` of the input coordinates, or
        ``None`` if ``coords`` is empty.
    """
    if not coords:
        return None
    n = len(coords)
    return (sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n)


def _parse_pdb(path):
    """Parse ATOM/HETATM records into ordered residues.

    Returns a list of residues, each:
        {"chain", "resseq", "icode", "resname", "atoms": OrderedDict(name -> (x,y,z))}
    Only the first altloc of each atom name in a residue is kept.
    """
    residues = OrderedDict()      # key -> residue dict
    order = []
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resseq = line[22:26].strip()
            icode = line[26]
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            segid = line[72:76].strip()      # cols 73-76: biological-molecule label
            key = (chain, resseq, icode, resname, segid)
            if key not in residues:
                residues[key] = {"chain": chain, "resseq": resseq, "icode": icode,
                                 "resname": resname, "segid": segid,
                                 "atoms": OrderedDict()}
                order.append(key)
            atoms = residues[key]["atoms"]
            if name not in atoms:          # keep first altloc only
                atoms[name] = (x, y, z)
    return [residues[k] for k in order]


def _beads_for_residue(res, warn):
    """Build the coarse-grained beads for a single residue.

    Proteins map to one ``CA`` bead (residue name unchanged). Nucleotides map
    to a phosphate bead (``P``, when present), a ribose-ring centroid bead
    (``R``), and one base-ring centroid bead for pyrimidines (``BR1``) or two
    for purines (``BR1`` + ``BR2``). Missing or insufficient atoms trigger a
    warning via ``warn`` and the affected bead is skipped.

    Parameters
    ----------
    res : dict
        A residue dict as produced by :func:`_parse_pdb`, with keys
        ``"resname"``, ``"chain"``, ``"resseq"`` and ``"atoms"`` (an
        ``OrderedDict`` mapping atom name to ``(x, y, z)``).
    warn : callable
        Callback invoked with a single message string to record a warning
        about missing atoms or an unrecognized residue.

    Returns
    -------
    list of tuple
        A list of ``(atom_name, element, (x, y, z))`` beads. Empty if the
        residue type is not recognized (the caller skips it).
    """
    rn = res["resname"]
    atoms = res["atoms"]

    # ---- protein: one Cα bead, residue name unchanged ----
    if rn in AMINO_ACIDS:
        if "CA" in atoms:
            return [("CA", "C", atoms["CA"])]
        warn(f"protein residue {rn} {res['chain']}{res['resseq']} has no CA atom; skipped")
        return []

    # ---- RNA: 3 (pyrimidine) or 4 (purine) beads ----
    is_purine = rn in PURINES or ("N9" in atoms)
    is_pyrimidine = rn in PYRIMIDINES
    if is_purine or is_pyrimidine:
        beads = []

        # P bead: at the phosphate P atom (absent on a 5'-terminal nucleotide).
        if "P" in atoms:
            beads.append(("P", "P", atoms["P"]))

        # R bead: ribose-ring centroid.
        ribose = [atoms[a] for a in RIBOSE_RING if a in atoms]
        c = _centroid(ribose)
        if c is None or len(ribose) < 3:
            warn(f"nucleotide {rn} {res['chain']}{res['resseq']} has too few ribose "
                 f"ring atoms ({len(ribose)}/5); R bead skipped")
        else:
            beads.append(("R", "C", c))

        # BR bead(s): centroid of each conjugated base ring.
        ring6 = _centroid([atoms[a] for a in BASE_RING_6 if a in atoms])
        if is_purine:
            ring5 = _centroid([atoms[a] for a in PURINE_RING_5 if a in atoms])
            if ring6 is not None:
                beads.append(("BR1", "C", ring6))
            if ring5 is not None:
                beads.append(("BR2", "C", ring5))
            if ring6 is None or ring5 is None:
                warn(f"purine {rn} {res['chain']}{res['resseq']} missing base-ring atoms")
        else:  # pyrimidine: single ring -> BR1
            if ring6 is not None:
                beads.append(("BR1", "C", ring6))
            else:
                warn(f"pyrimidine {rn} {res['chain']}{res['resseq']} missing base-ring atoms")
        return beads

    warn(f"unrecognized residue {rn} {res['chain']}{res['resseq']} (not protein/RNA); skipped")
    return []


def _pdb_atom_line(serial, name, resname, chain, resseq, icode, xyz, element, segid=""):
    """Format one standard PDB ATOM record (fixed columns).

    Columns (1-indexed): serial 7-11, name 13-16, altLoc 17, resName 18-20,
    chainID 22, resSeq 23-26, iCode 27, x/y/z 31-54, occ 55-60, temp 61-66,
    segID 73-76, element 77-78.

    Parameters
    ----------
    serial : int
        Atom serial number; written modulo 100000.
    name : str
        Atom (bead) name. Names shorter than four characters are right-padded
        with a leading space so the name occupies columns 14-16.
    resname : str
        Residue name; truncated to three characters.
    chain : str
        Chain identifier; the first character is used (a blank if falsy).
    resseq : str
        Residue sequence number; truncated to four characters.
    icode : str
        Insertion code; a blank is written if it is empty/whitespace.
    xyz : tuple of float
        Cartesian ``(x, y, z)`` coordinates in angstrom.
    element : str
        Element symbol; truncated to two characters.
    segid : str, optional
        Segment identifier preserved in columns 73-76; truncated to four
        characters. Defaults to an empty string.

    Returns
    -------
    str
        A single newline-terminated PDB ``ATOM`` record.
    """
    # Atom name: 1-3 char names get a leading space (cols 13 blank, 14-16 name).
    name4 = name[:4] if len(name) >= 4 else (" " + name).ljust(4)
    return ("ATOM  %5d %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n" % (
        serial % 100000, name4, " ", resname[:3], (chain or " ")[:1], resseq[:4],
        icode if icode.strip() else " ", xyz[0], xyz[1], xyz[2], 1.0, 0.0,
        segid[:4], element[:2]))


def coarse_grain(input_pdb: str, output_pdb: str, verbose: bool = True) -> dict:
    """Coarse-grain ``input_pdb`` to the TOPO convention and write ``output_pdb``.

    Returns a stats dict (bead/residue counts, warnings).
    """
    warnings_list = []

    def warn(msg):
        """Record a warning message for the current coarse-graining run.

        Parameters
        ----------
        msg : str
            The warning message to append to the enclosing ``warnings_list``.
        """
        warnings_list.append(msg)

    residues = _parse_pdb(input_pdb)

    stats = dict(protein_beads=0, rna_nucleotides=0, P=0, R=0, BR=0,
                 purines=0, pyrimidines=0, residues_skipped=0)
    segids = OrderedDict()        # preserved biological-molecule labels -> bead count

    lines = []
    serial = 0
    prev_chain = None
    for res in residues:
        beads = _beads_for_residue(res, warn)
        if not beads:
            stats["residues_skipped"] += 1
            continue

        # TER between chains for readability.
        if prev_chain is not None and res["chain"] != prev_chain:
            lines.append("TER\n")
        prev_chain = res["chain"]

        # tally
        if res["resname"] in AMINO_ACIDS:
            stats["protein_beads"] += 1
        else:
            stats["rna_nucleotides"] += 1
            if any(b[0] == "BR2" for b in beads):
                stats["purines"] += 1
            else:
                stats["pyrimidines"] += 1

        sid = res.get("segid", "")
        for name, element, xyz in beads:
            serial += 1
            segids[sid] = segids.get(sid, 0) + 1   # beads per biological molecule
            btype = name.rstrip("0123456789")   # CA/P/R/BR
            if btype in ("P", "R", "BR"):
                stats[btype] += 1
            lines.append(_pdb_atom_line(serial, name, res["resname"], res["chain"],
                                        res["resseq"], res["icode"], xyz, element,
                                        res.get("segid", "")))
    lines.append("END\n")

    with open(output_pdb, "w") as fh:
        fh.write(f"REMARK  TOPO coarse-grained model of {input_pdb}\n")
        fh.write("REMARK  protein: 1 bead/residue at CA (resname kept). "
                 "RNA: P (phosphate), R (ribose), BR1/BR2 (base ring centroids).\n")
        fh.writelines(lines)

    if verbose:
        print(f"Coarse-grained {input_pdb} -> {output_pdb}")
        print(f"  protein Cα beads : {stats['protein_beads']}")
        print(f"  RNA nucleotides  : {stats['rna_nucleotides']} "
              f"(purines {stats['purines']}, pyrimidines {stats['pyrimidines']})")
        print(f"  RNA beads        : P={stats['P']}  R={stats['R']}  BR={stats['BR']}")
        print(f"  total beads      : {serial}")
        print(f"  segIDs preserved : {len(segids)}  -> {dict(segids)}")
        if stats["residues_skipped"]:
            print(f"  residues skipped : {stats['residues_skipped']}")
        if warnings_list:
            print(f"  warnings         : {len(warnings_list)} "
                  f"(e.g. {warnings_list[0]})")
    stats["total_beads"] = serial
    stats["segids"] = segids
    stats["warnings"] = warnings_list
    return stats


def main(argv=None):
    """Command-line entry point for coarse-graining a structure.

    Parses ``--input``/``-i`` and ``--output``/``-o`` arguments (both
    required) and calls :func:`coarse_grain` to convert the all-atom input
    PDB to a TOPO coarse-grained PDB (protein ``CA`` beads; RNA ``P``/``R``/
    ``BR`` beads).

    Parameters
    ----------
    argv : list of str, optional
        Argument vector to parse. If ``None`` (the default), arguments are
        read from ``sys.argv``.
    """
    p = argparse.ArgumentParser(
        prog="python -m topo.csp.cg_ribosome",
        description="Coarse-grain a protein+RNA structure to the TOPO convention "
                    "(protein: CA; RNA: P/R/BR beads).")
    p.add_argument("-i", "--input", required=True, help="input all-atom PDB")
    p.add_argument("-o", "--output", required=True, help="output coarse-grained PDB")
    args = p.parse_args(argv)
    coarse_grain(args.input, args.output)


if __name__ == "__main__":
    main()
