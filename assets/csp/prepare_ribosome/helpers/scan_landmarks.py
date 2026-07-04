#!/usr/bin/env python3
"""Discover PTC / exit landmarks for a target ribosome by homology to E. coli.

Superpose the target P-site tRNA acceptor arm onto the validated E. coli oriented
PtR, then in that frame:
  * report the acceptor-arm fit RMSD (low => this really is the P-site tRNA),
  * list large-rRNA **adenines** whose N6 is closest to the origin (PTC candidates),
  * list uL24-homolog residues whose backbone N is closest to E. coli's exit point
    (~99.9, 0, 0) AND on the +x tunnel side (exit candidates).

Usage:
  python scan_landmarks.py <cif> <p_tRNA_chain> <bigRNA_chain> <uL24_chain>
"""
import sys
import numpy as np
import gemmi

ECOLI_ORI = "out/ecoli/4v9d_50S_PtR_5jte_AtR_model_oriented.pdb"
EXIT_REF = np.array([99.9, 0.0, 0.0])
ACCEPT = set(range(1, 8)) | set(range(66, 77))


def load_chain(model, cid):
    for c in model:
        if c.name == cid:
            out = {}
            for r in c:
                out.setdefault(r.seqid.num, {"name": r.name, "atoms": {}})
                for a in r:
                    out[r.seqid.num]["atoms"][a.name] = np.array(
                        [a.pos.x, a.pos.y, a.pos.z])
            return out
    raise ValueError(f"chain {cid} not found")


def load_pdb_seg(path, seg):
    d = {}
    for line in open(path):
        if line[:4] not in ("ATOM", "HETA") or line[72:76].strip() != seg:
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
    return R, t, float(np.sqrt((((P @ R.T) + t - Q) ** 2).sum() / len(P)))


def main():
    cif, ptrna, bigrna, ul24 = sys.argv[1:5]
    st = gemmi.read_structure(cif); m = st[0]
    ptr = load_chain(m, ptrna)
    eco = load_pdb_seg(ECOLI_ORI, "PtR")
    P, Q = [], []
    for sid, rec in ptr.items():
        if sid not in ACCEPT:
            continue
        for nm in ("P", "C1'", "C4'"):
            if nm in rec["atoms"] and (sid, nm) in eco:
                P.append(rec["atoms"][nm]); Q.append(eco[(sid, nm)])
    R, t, rmsd = kabsch(P, Q)
    print(f"{cif}: P-tRNA {ptrna} acceptor-arm fit RMSD = {rmsd:.2f} A ({len(P)} atoms)")
    to_eco = lambda x: R @ x + t

    # PTC candidates: adenines near origin
    big = load_chain(m, bigrna)
    ade = []
    for sid, rec in big.items():
        at = rec["atoms"]
        if "N6" in at and {"N1", "C2", "N9"} <= set(at):
            d = np.linalg.norm(to_eco(at["N6"]))
            ade.append((d, sid, rec["name"]))
    ade.sort()
    print(f"  PTC candidates (adenines, big rRNA {bigrna}, nearest origin):")
    for d, sid, nm in ade[:5]:
        print(f"    {bigrna}:{sid} {nm}  |N6 - origin| = {d:.1f} A")

    # exit candidates: uL24-homolog residues near the E. coli exit point, +x side
    u = load_chain(m, ul24)
    cand = []
    for sid, rec in u.items():
        if "N" in rec["atoms"]:
            p = to_eco(rec["atoms"]["N"])
            cand.append((np.linalg.norm(p - EXIT_REF), sid, rec["name"], p))
    cand.sort()
    print(f"  EXIT candidates (uL24-homolog {ul24}, nearest E. coli exit ~{EXIT_REF}):")
    for d, sid, nm, p in cand[:6]:
        print(f"    {ul24}:{sid} {nm} N  pos={p.round(1)}  dist={d:.1f} A")


if __name__ == "__main__":
    main()
