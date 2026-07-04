#!/usr/bin/env python3
"""List polymer entities (description + auth chains) from a cif.

Usage: python entities.py <cif> [grep]
"""
import sys
import gemmi

def main():
    cif=sys.argv[1]
    filt=sys.argv[2].lower() if len(sys.argv)>2 else None
    doc=gemmi.cif.read(cif)
    b=doc.sole_block()
    # entity_id -> description
    desc={}
    for eid,d in zip(b.find_loop("_entity.id"), b.find_loop("_entity.pdbx_description")):
        desc[eid]=d
    # entity_poly: entity_id -> strand ids, type
    rows=[]
    ids=b.find_loop("_entity_poly.entity_id")
    strands=b.find_loop("_entity_poly.pdbx_strand_id")
    types=b.find_loop("_entity_poly.type")
    for eid,strand,typ in zip(ids,strands,types):
        d=desc.get(eid,"?")
        rows.append((eid,typ,strand,d))
    for eid,typ,strand,d in rows:
        if filt and filt not in d.lower(): continue
        print(f"  ent {eid:3s} [{typ:22s}] chains={strand:20s} {d}")

if __name__=="__main__":
    main()
