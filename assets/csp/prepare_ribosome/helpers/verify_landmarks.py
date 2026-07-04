#!/usr/bin/env python3
"""Verify candidate PTC / exit landmarks for a target ribosome by homology.

Method: superpose the target's P-site tRNA acceptor arm (residues 1-7 + 66-76,
the rigid, universally conserved amino-acyl stem that sits at the PTC) onto the
*validated E. coli oriented* P-site tRNA. That brings the whole target ribosome
into E. coli's oriented frame, where the PTC is at the origin and the tunnel exit
is at roughly (+100, 0, 0). A correct PTC landmark then lands near the origin and
a correct exit landmark lands at large +x on the tunnel axis.

Usage:
  python verify_landmarks.py <cif> <p_tRNA_chain> <bigRNA_chain> <PTC_resid> \
        <uL24_chain> <exit_resid> [--acceptor 1-7,66-76]
"""
import argparse
import numpy as np
import gemmi

ECOLI_ORI = "out/ecoli/4v9d_50S_PtR_5jte_AtR_model_oriented.pdb"


def load_chain_atoms(model, chain_id):
    """Return {(seqid, atomname): xyz} for a gemmi chain."""
    for c in model:
        if c.name == chain_id:
            d = {}
            for r in c:
                for a in r:
                    d[(r.seqid.num, a.name)] = np.array([a.pos.x, a.pos.y, a.pos.z])
            return d
    raise ValueError(f"chain {chain_id} not found")


def load_pdb_seg(path, seg):
    d = {}
    for line in open(path):
        if line[:4] not in ("ATOM", "HETA"):
            continue
        if line[72:76].strip() != seg:
            continue
        d[(int(line[22:26]), line[12:16].strip())] = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return d


def kabsch(P, Q):
    P = np.asarray(P); Q = np.asarray(Q)
    pc = P.mean(0); qc = Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum() / len(P)))
    return R, t, rmsd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("p_trna_chain")
    ap.add_argument("bigrna_chain")
    ap.add_argument("ptc_resid", type=int)
    ap.add_argument("ul24_chain")
    ap.add_argument("exit_resid", type=int)
    ap.add_argument("--acceptor", default="1-7,66-76")
    a = ap.parse_args()

    # acceptor-arm residue set
    accept = set()
    for part in a.acceptor.split(","):
        lo, hi = part.split("-")
        accept |= set(range(int(lo), int(hi) + 1))

    st = gemmi.read_structure(a.cif)
    m = st[0]
    tgt_ptrna = load_chain_atoms(m, a.p_trna_chain)
    eco_ptr = load_pdb_seg(ECOLI_ORI, "PtR")

    # superpose target P-tRNA acceptor arm onto E. coli oriented PtR
    P, Q = [], []
    for (sid, name), xyz in tgt_ptrna.items():
        if sid in accept and name in ("P", "C1'", "C4'") and (sid, name) in eco_ptr:
            P.append(xyz); Q.append(eco_ptr[(sid, name)])
    if len(P) < 6:
        print(f"WARN only {len(P)} acceptor-arm atoms matched; using all P/C1'/C4'")
        P, Q = [], []
        for (sid, name), xyz in tgt_ptrna.items():
            if name in ("P", "C1'", "C4'") and (sid, name) in eco_ptr:
                P.append(xyz); Q.append(eco_ptr[(sid, name)])
    R, t, rmsd = kabsch(P, Q)
    print(f"P-tRNA acceptor-arm superposition onto E. coli PtR: "
          f"{len(P)} atoms, RMSD {rmsd:.2f} A")

    def to_eco(xyz):
        return R @ xyz + t

    # candidate PTC atom
    bigrna = load_chain_atoms(m, a.bigrna_chain)
    ptc_key = None
    for name in ("N6", "N1", "C2"):
        if (a.ptc_resid, name) in bigrna:
            ptc_key = (a.ptc_resid, name); break
    if ptc_key is None:
        print(f"  PTC resid {a.ptc_resid}: NOT FOUND in chain {a.bigrna_chain}")
    else:
        # is it an adenine? check for N6 + N1 + C2 + N9
        atoms_here = {nm for (sid, nm) in bigrna if sid == a.ptc_resid}
        is_a = {"N6", "N1", "C2", "N3", "C4", "C5", "C6", "N9"} <= atoms_here
        p_eco = to_eco(bigrna[(a.ptc_resid, "N6")]) if (a.ptc_resid, "N6") in bigrna \
            else to_eco(bigrna[ptc_key])
        print(f"  PTC {a.bigrna_chain}:{a.ptc_resid} adenine={is_a}  "
              f"N6-in-Ecoli-frame = {p_eco.round(1)}  |dist to origin| = "
              f"{np.linalg.norm(p_eco):.1f} A")

    # candidate exit atom
    ul24 = load_chain_atoms(m, a.ul24_chain)
    if (a.exit_resid, "N") in ul24:
        e_eco = to_eco(ul24[(a.exit_resid, "N")])
        print(f"  EXIT {a.ul24_chain}:{a.exit_resid}:N in-Ecoli-frame = "
              f"{e_eco.round(1)}  (target +x ~ 100, y~0, z~0)")
    else:
        print(f"  EXIT resid {a.exit_resid}:N NOT FOUND in chain {a.ul24_chain}")

    # reference E. coli exit for comparison
    print("  [ref] E. coli exit L24:51:N is at ~(99.9, 0, 0) in this frame")


if __name__ == "__main__":
    main()
