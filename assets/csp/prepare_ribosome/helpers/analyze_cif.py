#!/usr/bin/env python3
"""Print per-chain polymer type, length, and entity description for a cif.

Usage: python analyze_cif.py <cif> [--rna] [--grep TEXT]
Helps build the chain->segID map for a new organism.
"""
import sys, argparse
import gemmi

aa3={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
rna={'A','U','G','C'}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--rna",action="store_true",help="only RNA chains")
    ap.add_argument("--grep",default=None,help="filter entity description (case-insensitive)")
    ap.add_argument("--maxlen",type=int,default=None)
    ap.add_argument("--minlen",type=int,default=None)
    a=ap.parse_args()
    st=gemmi.read_structure(a.cif)
    # entity description per subchain / chain
    ent_desc={}
    for ent in st.entities:
        for sub in ent.subchains:
            ent_desc[sub]=ent.name
    m=st[0]
    rows=[]
    for ch in m:
        naa=sum(1 for r in ch if r.name in aa3)
        nnt=sum(1 for r in ch if r.name in rna or (gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_nucleic_acid()))
        typ = "RNA" if nnt>naa else ("PROT" if naa>0 else "other")
        n=naa+nnt
        if n==0: continue
        # entity desc: look up via first residue subchain
        desc=""
        try:
            sub=ch[0].subchain
            desc=ent_desc.get(sub,"")
        except Exception: pass
        # also try polymer entity description
        rows.append((ch.name,typ,n,desc))
    for name,typ,n,desc in sorted(rows,key=lambda x:-x[2]):
        if a.rna and typ!="RNA": continue
        if a.maxlen and n>a.maxlen: continue
        if a.minlen and n<a.minlen: continue
        if a.grep and a.grep.lower() not in desc.lower(): continue
        print(f"  {name:4s} {typ:5s} n={n:5d}  {desc}")

if __name__=="__main__":
    main()
