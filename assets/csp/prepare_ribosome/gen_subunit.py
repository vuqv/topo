#!/usr/bin/env python3
"""
Stage 2 of the ribosome-prep pipeline: **gen** the large subunit + tRNAs.

Assemble a ribosome large subunit (rRNAs + large-subunit proteins) plus the
A-/P-site tRNAs into **one all-atom PDB** carrying the TOPO **segID** convention
in columns 73-76 (``23S`` / ``25S`` / ``28S`` / ``5S`` / ``5.8S`` for rRNAs,
``L2``..``L##`` for proteins, ``PtR`` / ``AtR`` for the P-/A-site tRNAs). The
segID column is what every downstream stage (``cg_ribosome``, ``truncate_ribosome``,
``ribosome.load_ribosome``) keys off, so it is populated for every atom.

What it does
------------
* Reads a per-organism **INI config** (``configs/<organism>.ini``): one
  ``[mol:<segid>]`` section per biological molecule to keep, each giving the
  **source chain** in the deposition and (optionally) an output chain letter.
* **Keeps only the ribosome scenery** — you whitelist exactly the molecules to
  keep by segID; everything else in the deposition (small subunit, mRNA,
  factors, a bound **nascent chain**, antibiotics, ions, water) is simply never
  extracted and therefore dropped.
* **Canonicalises** modified residues to their standard parent (modified
  nucleotides -> ``A``/``U``/``G``/``C``; modified amino acids -> the standard
  20) and keeps only standard protein/RNA residues, so the coarse-grainer
  recognises every residue.
* Within each molecule keeps only residues of the molecule's **dominant polymer
  type**, which cleanly drops an **aminoacyl amino acid** carried on a charged
  tRNA (e.g. the LYS on 5JTE's A-site tRNA) or any stray het.
* **Renumbers each tRNA so its 3'-terminal acceptor is residue 76** -- the
  convention the whole downstream CSP engine assumes (``read_anchor`` /
  ``add_trna_tether`` / ``optimal_ptc_targets`` look up ``PtR``/``AtR``
  ``:76:{R,P,BR2}``). A no-op for standard 76-nt tRNAs; it shifts e.g. the human
  77-nt P-site tRNA (deposition acceptor A77) down to 76 so its purine acceptor
  -- with an ``R``, ``P`` and ``BR2`` bead -- lands where the engine looks.
* **A-site tRNA graft:** a ``[mol:AtR]`` section may point at a *different* cif
  (``cif = ...``) and ask for a rigid **large-subunit superposition**
  (``superpose = lsu``) that aligns the donor's large subunit onto the target's
  before lifting the tRNA into the target's A-site frame.

Usage
-----
    python gen_subunit.py -c configs/ecoli.ini

Writes the assembled all-atom model to the ``out =`` path from the config and
prints a per-segID atom/residue count so a human can eyeball the result.
"""
from __future__ import annotations

import argparse
import configparser
import os
import sys
from collections import OrderedDict

import numpy as np
import gemmi

# ---------------------------------------------------------------------------
# Residue tables
# ---------------------------------------------------------------------------
STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
STD_NT = {"A", "U", "G", "C"}
ONE2THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU",
    "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def classify_residue(res):
    """Classify a gemmi residue and return ``(kind, canonical_name)``.

    ``kind`` is ``"aa"`` (standard/modified amino acid), ``"nt"`` (standard/
    modified ribonucleotide), or ``None`` (water/ion/ligand -> drop). The
    canonical name is a standard 3-letter amino acid or a single ``A/U/G/C``
    nucleotide; modified residues are mapped to their standard parent so the
    coarse-grainer recognises them.
    """
    n = res.name
    if n in STD_AA:
        return "aa", n
    if n in STD_NT:
        return "nt", n

    info = gemmi.find_tabulated_residue(n)
    if info is not None and info.is_water():
        return None, None

    atomset = {a.name for a in res}
    has_ribose = "C1'" in atomset and "C4'" in atomset and "O4'" in atomset

    # --- nucleotide (tabulated as nucleic acid, or has a ribose ring) ---
    if (info is not None and info.is_nucleic_acid()) or has_ribose:
        if info is not None and info.one_letter_code.upper() in STD_NT:
            return "nt", info.one_letter_code.upper()
        # fall back to base-atom sniffing for untabulated modifications
        if "N9" in atomset:                              # purine
            return "nt", ("G" if ("O6" in atomset or "N2" in atomset) else "A")
        return "nt", ("C" if "N4" in atomset else "U")   # pyrimidine

    # --- amino acid (tabulated, or a protein backbone with no ribose) ---
    backbone = {"N", "CA", "C", "O"} <= atomset
    if (info is not None and info.is_amino_acid()) or backbone:
        if info is not None and info.one_letter_code.upper() in ONE2THREE:
            return "aa", ONE2THREE[info.one_letter_code.upper()]
        return "aa", (n if n in STD_AA else "ALA")

    return None, None                                    # ion / ligand -> drop


def _parse_resrange(spec):
    """Parse a ``"lo-hi"`` residue range string into an ``(lo, hi)`` int tuple."""
    lo, hi = spec.split("-")
    return int(lo), int(hi)


def extract_molecule(model, chain_id, resrange=None):
    """Extract one chain as a list of canonicalised residues.

    Returns ``[(seqid, icode, canon_name, kind, [(atom_name, element, xyz)])]``
    in chain order, keeping only residues of the chain's **dominant** polymer
    type (so an aminoacyl amino acid on a tRNA is dropped) and, if ``resrange``
    is given, only residues within ``[lo, hi]``. Hydrogens are dropped.
    """
    chain = None
    for c in model:
        if c.name == chain_id:
            chain = c
            break
    if chain is None:
        raise ValueError(f"chain {chain_id!r} not found in structure")

    raw = []                              # (seqid, icode, kind, canon, atoms)
    counts = {"aa": 0, "nt": 0}
    for res in chain:
        kind, canon = classify_residue(res)
        if kind is None:
            continue
        if resrange is not None and not (resrange[0] <= res.seqid.num <= resrange[1]):
            continue
        atoms = []
        for a in res:
            if a.is_hydrogen():
                continue
            atoms.append((a.name, a.element.name,
                          np.array([a.pos.x, a.pos.y, a.pos.z])))
        raw.append((res.seqid.num, res.seqid.icode, kind, canon, atoms))
        counts[kind] += 1

    dominant = "aa" if counts["aa"] >= counts["nt"] else "nt"
    return [(sid, ic, canon, kind, atoms)
            for (sid, ic, kind, canon, atoms) in raw if kind == dominant]


def _kabsch(P, Q):
    """Return proper-rotation ``R`` and translation ``t`` best-fitting ``P`` to ``Q``."""
    P = np.asarray(P); Q = np.asarray(Q)
    pc = P.mean(0); qc = Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum() / len(P)))
    return R, t, rmsd


def superpose_via_trna(donor_model, donor_ref_chain, main_model, target_ref_chain):
    """Best-fit a donor tRNA onto the target's P-site tRNA; return ``(R, t)``.

    Used for **cross-species** A-site grafts where donor and target share no
    identical sequence: the P-/A-site tRNA arrangement at the ribosome is
    universally conserved, so aligning the donor's P-site tRNA onto the target's
    P-site tRNA (backbone atoms, matched by residue number) places the donor's
    A-site tRNA into the target's A-site frame when lifted with the same ``(R, t)``.
    """
    # Restrict to the acceptor arm (res 1-7 + 66-76): the rigid amino-acyl stem
    # that sits at the PTC and is the part whose pose we actually need to match.
    accept = set(range(1, 8)) | set(range(66, 77))

    def bb(model, cid):
        for c in model:
            if c.name == cid:
                d = {}
                for r in c:
                    if r.seqid.num not in accept:
                        continue
                    for a in r:
                        if a.name in ("P", "C1'", "C4'"):
                            d[(r.seqid.num, a.name)] = np.array(
                                [a.pos.x, a.pos.y, a.pos.z])
                return d
        raise ValueError(f"chain {cid} not found")

    d_bb = bb(donor_model, donor_ref_chain)
    t_bb = bb(main_model, target_ref_chain)
    P = [xyz for key, xyz in d_bb.items() if key in t_bb]
    Q = [t_bb[key] for key in d_bb if key in t_bb]
    if len(P) < 6:
        raise RuntimeError(f"via_trna graft: only {len(P)} shared backbone atoms")
    R, t, rmsd = _kabsch(P, Q)
    print(f"  [graft] P-site-tRNA superposition ({donor_ref_chain}->{target_ref_chain}): "
          f"{len(P)} atoms, fit RMSD = {rmsd:.2f} A")
    return R, t


def _chain_seq_and_bb(chain):
    """Return ``(sequence, {(seqid, bb_atom): xyz})`` for a gemmi chain.

    ``bb_atom`` is ``CA`` for amino acids and ``P`` for nucleotides; the
    sequence uses one-letter codes (upper = protein, lower = RNA).
    """
    three2one = {v: k for k, v in ONE2THREE.items()}
    seq = []
    bb = {}
    for res in chain:
        kind, canon = classify_residue(res)
        if kind == "aa":
            seq.append(three2one.get(canon, "X"))
            key = "CA"
        elif kind == "nt":
            seq.append(canon.lower())
            key = "P"
        else:
            continue
        for a in res:
            if a.name == key:
                bb[(res.seqid.num, key)] = np.array([a.pos.x, a.pos.y, a.pos.z])
    return "".join(seq), bb


def superpose_lsu(donor_model, target_lsu_chains, target_model):
    """Best-fit the donor's large subunit onto the target's; return ``(R, t)``.

    ``target_lsu_chains`` is the list of target chain IDs that make up the large
    subunit (rRNAs + proteins, tRNAs excluded). Donor chains are matched to
    target chains by exact one-letter sequence; all shared backbone atoms
    (CA/P, matched by residue number) drive one global Kabsch fit.
    """
    tgt_seq = {}
    tgt_bb = {}
    for cid in target_lsu_chains:
        for c in target_model:
            if c.name == cid:
                s, bb = _chain_seq_and_bb(c)
                tgt_seq[cid] = s
                tgt_bb[cid] = bb
    donor = {}
    for c in donor_model:
        s, bb = _chain_seq_and_bb(c)
        if s:
            donor[c.name] = (s, bb)

    P, Q = [], []
    used = 0
    for cid in target_lsu_chains:
        ts = tgt_seq.get(cid)
        if not ts:
            continue
        match = [dc for dc, (ds, _) in donor.items() if ds == ts]
        if not match:
            continue
        _, dbb = donor[match[0]]
        for key, dxyz in dbb.items():
            if key in tgt_bb[cid]:
                P.append(dxyz); Q.append(tgt_bb[cid][key])
        used += 1
    if len(P) < 3:
        raise RuntimeError("large-subunit superposition: too few matched atoms "
                           f"({len(P)}); check that donor and target share the LSU")
    R, t, rmsd = _kabsch(P, Q)
    print(f"  [graft] LSU superposition: {used} chains, {len(P)} backbone atoms, "
          f"fit RMSD = {rmsd:.2f} A")
    return R, t


# ---------------------------------------------------------------------------
# PDB writing
# ---------------------------------------------------------------------------
def _pdb_line(serial, name, resname, chain, resseq, icode, xyz, element, segid):
    """Format one fixed-column PDB ATOM record (segID in cols 73-76)."""
    name4 = name[:4] if len(name) >= 4 else (" " + name).ljust(4)
    ic = icode if (icode and icode.strip()) else " "
    return ("ATOM  %5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n"
            % (serial % 100000, name4, " ", resname[:3], (chain or " ")[:1],
               resseq % 10000, ic, xyz[0], xyz[1], xyz[2], 1.0, 0.0,
               segid[:4], element[:2]))


def read_config(path):
    """Read the assembler INI. Returns ``(assembly_dict, [mol_dict, ...])``."""
    cp = configparser.ConfigParser()
    cp.optionxform = str                  # preserve segID case (L2, 23S, PtR)
    if not cp.read(path):
        raise FileNotFoundError(path)
    asm = dict(cp["assembly"])
    trna_segids = {s.strip() for s in asm.get("trna_segids", "PtR,AtR").split(",")}
    mols = []
    for sec in cp.sections():
        if not sec.startswith("mol:"):
            continue
        segid = sec.split(":", 1)[1]
        d = dict(cp[sec])
        d["segid"] = segid
        d["is_trna"] = segid in trna_segids
        mols.append(d)
    return asm, mols


def gen_subunit(config_path, verbose=True):
    """Run stage 2 from an INI config; write the assembled all-atom PDB."""
    cfgdir = os.path.dirname(os.path.abspath(config_path))
    asm, mols = read_config(config_path)

    def rel(p):
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(cfgdir, p))

    main_cif = rel(asm["main_cif"])
    out_path = rel(asm["out"])
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    main_model = gemmi.read_structure(main_cif)[0]
    donor_cache = {}

    def donor_model(cif):
        if cif not in donor_cache:
            donor_cache[cif] = gemmi.read_structure(rel(cif))[0]
        return donor_cache[cif]

    # LSU reference chains (everything from main_cif that is not a tRNA/graft).
    lsu_chains = [m["chain"] for m in mols
                  if not m["is_trna"] and "cif" not in m]

    records = []          # (segid, out_chain, residues)
    for m in mols:
        segid = m["segid"]
        out_chain = m.get("out_chain", segid[:1])
        resrange = _parse_resrange(m["resrange"]) if "resrange" in m else None

        if "cif" in m:                       # grafted from another structure
            dmodel = donor_model(m["cif"])
            residues = extract_molecule(dmodel, m["chain"], resrange)
            mode = m.get("superpose", "").lower()
            R = t = None
            if mode == "lsu":
                R, t = superpose_lsu(dmodel, lsu_chains, main_model)
            elif mode == "via_trna":
                R, t = superpose_via_trna(dmodel, m["ref_donor"],
                                          main_model, m["ref_target"])
            if R is not None:
                residues = [(sid, ic, canon, kind,
                             [(an, el, R @ xyz + t) for an, el, xyz in atoms])
                            for (sid, ic, canon, kind, atoms) in residues]
        else:                                # from the main deposition
            residues = extract_molecule(main_model, m["chain"], resrange)

        # Renumber tRNAs so the 3'-terminal acceptor is residue 76 -- the
        # convention the whole downstream CSP engine assumes (read_anchor /
        # add_trna_tether / optimal_ptc_targets all look up segid:76:{R,P,BR2}).
        # No-op for standard 76-nt tRNAs; shifts e.g. the human 77-nt P-site
        # tRNA (acceptor at 77) down by one so its purine acceptor lands on 76.
        if m["is_trna"] and residues:
            shift = max(sid for sid, *_ in residues) - 76
            if shift:
                residues = [(sid - shift, ic, canon, kind, atoms)
                            for (sid, ic, canon, kind, atoms) in residues]
                print(f"  [{segid}] renumbered by {-shift:+d} so 3'-acceptor = residue 76")

        records.append((segid, out_chain, residues))

    # ---- write ----
    serial = 0
    stats = OrderedDict()
    with open(out_path, "w") as fh:
        fh.write(f"REMARK  TOPO stage-2 assembly from {os.path.basename(config_path)}\n")
        fh.write("REMARK  segID (cols 73-76) labels each biological molecule.\n")
        for segid, out_chain, residues in records:
            n_at = 0
            for sid, ic, canon, kind, atoms in residues:
                for an, el, xyz in atoms:
                    serial += 1
                    n_at += 1
                    fh.write(_pdb_line(serial, an, canon, out_chain, sid, ic,
                                       xyz, el, segid))
            fh.write("TER\n")
            stats[segid] = (len(residues), n_at)
        fh.write("END\n")

    if verbose:
        print(f"Assembled {len(records)} molecules -> {out_path}")
        print(f"  {'segID':6s} {'residues':>9s} {'atoms':>8s}")
        tot_r = tot_a = 0
        for segid, (nr, na) in stats.items():
            print(f"  {segid:6s} {nr:9d} {na:8d}")
            tot_r += nr; tot_a += na
        print(f"  {'TOTAL':6s} {tot_r:9d} {tot_a:8d}")
    return out_path, stats


def main(argv=None):
    """CLI entry point: ``gen_subunit.py -c configs/<organism>.ini``."""
    p = argparse.ArgumentParser(
        description="Stage 2: assemble ribosome large subunit + tRNAs into one "
                    "all-atom PDB with TOPO segIDs.")
    p.add_argument("-c", "--config", required=True, help="organism INI config")
    args = p.parse_args(argv)
    gen_subunit(args.config)


if __name__ == "__main__":
    main()
