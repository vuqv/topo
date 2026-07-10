#!/usr/bin/env python3
"""Generate a stage-2/3 INI for a eukaryotic large subunit from its cif entities.

Selects the large-subunit rRNAs (big rRNA + 5.8S + 5S) and **all** large-subunit
(60S) proteins, plus the P-site tRNA and the A-site tRNA (native chain or a graft
spec), and writes the [assembly] / [mol:*] / [orient] config. The small subunit
(40S/18S), mRNA, factors, RACK1, nascent chain, ions and water are dropped by
simply not whitelisting them.

Per-organism specs are defined in ORGS below; run `python make_euk_config.py <org>`.
"""
import re
import sys
import numpy as np
import gemmi

# ---------------------------------------------------------------------------
# Per-organism specification (chains/landmarks verified in PROVENANCE.md)
# ---------------------------------------------------------------------------
ORGS = {
    "yeast": dict(
        tag="6q8y_60S", cif="../raw/yeast/6q8y.cif", org="yeast",
        big_rrna=("BQ", "25S"), rrna_5_8=("BS", "5.8S"), rrna_5=("BR", "5S"),
        small_rrna="2",
        p_trna="n",
        a_trna_graft=dict(cif="../raw/human/8g61.cif", chain="At",
                          resrange="1-76", ref_donor="Pt", ref_target="n"),
        ptc="25S:2971:N6", exit="L26:91:N",
    ),
    "ncrassa": dict(
        tag="7r81_60S", cif="../raw/ncrassa/7r81.cif", org="ncrassa",
        big_rrna=("A1", "26S"), rrna_5_8=("C1", "5.8S"), rrna_5=("B1", "5S"),
        small_rrna="A2",
        p_trna="u1",                       # label CC = P-site (verified)
        a_trna_native="t1",                # label BC = A-site (verified)
        ptc="26S:2931:N6", exit="L26:91:N",
    ),
    "human": dict(
        tag="8g61_60S", cif="../raw/human/8g61.cif", org="human",
        big_rrna=("L5", "28S"), rrna_5_8=("L8", "5.8S"), rrna_5=("L7", "5S"),
        small_rrna="S2",
        p_trna="Pt",
        a_trna_native="At",
        ptc="28S:4548:N6", exit="L26:93:N",
    ),
}

# Author chain of the exit-tunnel protein per organism. This is the uL24 family
# (universal nomenclature): E. coli protein L24, eukaryotic RPL26. Its segID is
# forced to L26 below so the exit landmark is exact. The name describes the role
# (exit-tunnel protein), not the E. coli-specific "L24".
EXIT_PROTEIN_CHAIN = {"yeast": "AK", "ncrassa": "a1", "human": "LY"}


# Non-ribosomal polypeptides that may sit near the large subunit but must be
# dropped (they are not part of the rigid ribosome scenery).
NON_RIBOSOMAL = ("nascent", "factor", "exoribonuclease", "xrn",
                 "guanine nucleotide", "beta-like", "chaperone", "trigger")


def _entity_desc_map(block):
    """entity_id -> unquoted pdbx_description."""
    out = {}
    for eid, d in zip(block.find_loop("_entity.id"),
                      block.find_loop("_entity.pdbx_description")):
        out[eid] = gemmi.cif.as_string(d)
    return out


def big_protein_chains(cif, big_rrna_chain, small_rrna_chain):
    """Return [(auth_chain, description)] for **large-subunit** proteins.

    A polypeptide is assigned to the large subunit if its atoms are, on the
    whole, closer to the big rRNA (25S/26S/28S) than to the small rRNA (18S) --
    a deposition-agnostic structural test that separates 60S from 40S proteins
    regardless of naming quirks. Explicit non-ribosomal chains (nascent peptide,
    factors, RACK1, ...) are dropped by name.
    """
    doc = gemmi.cif.read(cif)
    b = doc.sole_block()
    desc = _entity_desc_map(b)
    st = gemmi.read_structure(cif)
    m = st[0]

    def chain_atoms(cid, atom_filter=None):
        for c in m:
            if c.name == cid:
                pts = []
                for r in c:
                    for a in r:
                        if atom_filter is None or a.name == atom_filter:
                            pts.append([a.pos.x, a.pos.y, a.pos.z])
                return np.asarray(pts)
        return np.empty((0, 3))

    # subsample rRNA P atoms for a fast min-distance test
    big = chain_atoms(big_rrna_chain, "P")
    small = chain_atoms(small_rrna_chain, "P")
    if len(big) == 0 or len(small) == 0:
        raise RuntimeError("could not load rRNA reference atoms")

    def min_dist(pts, ref):
        # coarse: nearest ref atom to the chain centroid + a few atoms
        c = pts.mean(0)
        return np.linalg.norm(ref - c, axis=1).min()

    # map auth chain -> entity description
    ch_desc = {}
    for eid, strand, typ in zip(b.find_loop("_entity_poly.entity_id"),
                                b.find_loop("_entity_poly.pdbx_strand_id"),
                                b.find_loop("_entity_poly.type")):
        if "polypeptide" not in typ:
            continue
        for ch in strand.split(","):
            ch_desc[ch.strip()] = desc.get(eid, "")

    out = []
    for cid, d in ch_desc.items():
        if any(k in d.lower() for k in NON_RIBOSOMAL):
            continue
        pts = chain_atoms(cid, "CA")
        if len(pts) == 0:
            continue
        if min_dist(pts, big) < min_dist(pts, small):   # large-subunit side
            out.append((cid, d))
    return out


def segid_for(desc, taken):
    """Derive a unique <=4-char segID from a protein description."""
    m = re.search(r"[Ll](\d+)", desc)
    base = f"L{m.group(1)}" if m else "LP"
    seg = base
    i = 0
    while seg in taken:
        i += 1
        seg = (base + chr(ord('a') + i - 1))[:4]
    taken.add(seg)
    return seg


def main():
    org = sys.argv[1]
    s = ORGS[org]
    cif_read = s["cif"].replace("../", "")  # config path is rel to configs/; read rel to prep_rib/
    lines = []
    A = lines.append
    A(f"# {org} 60S large subunit assembly (auto-generated by make_euk_config.py;")
    A(f"# landmarks + tRNA-site assignments verified -- see out/{s['org']}/PROVENANCE.md).")
    A("")
    A("[assembly]")
    A(f"tag = {s['tag']}")
    A(f"main_cif = {s['cif']}")
    A(f"out = ../out/{s['org']}/{s['tag']}_model.pdb")
    A("trna_segids = PtR,AtR")
    A("")
    # rRNAs
    for chain, seg in (s["big_rrna"], s["rrna_5_8"], s["rrna_5"]):
        A(f"[mol:{seg}]"); A(f"chain = {chain}"); A("")
    # P-site tRNA
    A("[mol:PtR]"); A(f"chain = {s['p_trna']}"); A("")
    # A-site tRNA: native or graft
    if "a_trna_native" in s:
        A("[mol:AtR]"); A(f"chain = {s['a_trna_native']}"); A("resrange = 1-76"); A("")
    else:
        g = s["a_trna_graft"]
        A("[mol:AtR]")
        A(f"cif = {g['cif']}")
        A(f"chain = {g['chain']}")
        A(f"resrange = {g['resrange']}")
        A("superpose = via_trna")
        A(f"ref_donor = {g['ref_donor']}")
        A(f"ref_target = {g['ref_target']}")
        A("")
    # 60S proteins
    taken = {s["big_rrna"][1], s["rrna_5_8"][1], s["rrna_5"][1], "PtR", "AtR"}
    exit_chain = EXIT_PROTEIN_CHAIN[org]
    prots_key = big_protein_chains(cif_read, s["big_rrna"][0], s["small_rrna"])
    # force the exit-tunnel protein chain -> L26 first so the landmark segID is exact
    ordered = sorted(prots_key, key=lambda cd: (cd[0] != exit_chain,))
    for chain, desc in ordered:
        if chain == exit_chain:
            seg = "L26"; taken.add(seg)
        else:
            seg = segid_for(desc, taken)
        A(f"[mol:{seg}]"); A(f"chain = {chain}"); A("")
    # orient
    A("[orient]")
    A(f"oriented_out = ../out/{s['org']}/{s['tag']}_model_oriented.pdb")
    A(f"ptc    = {s['ptc']}")
    A(f"exit   = {s['exit']}")
    A("a_site = AtR:last:ribose")   # 'last' = 3'-terminal acceptor (76, or 77 for human Pt)
    A("p_site = PtR:last:ribose")
    A("")
    out = f"configs/{org}.ini"
    open(out, "w").write("\n".join(lines))
    print(f"wrote {out}  ({len(prots_key)} large-subunit proteins + 3 rRNAs + PtR + AtR)")


if __name__ == "__main__":
    main()
