#!/usr/bin/env python
"""Carve one all-atom protein chain out of a ribosome CIF for native-contact use.

The continuous-synthesis ``ribo_free_mask`` feature makes a portion of one
ribosomal protein flexible (a topo-style Go loop). topo builds that loop's
native contacts with :func:`topo.utils.nonbonded.build_nonbonded_interaction`,
which needs an **all-atom** structure of the protein (STRIDE hydrogen bonds +
heavy-atom side-chain contacts) -- the coarse-grained truncated ribosome PDB
(Cα only) cannot supply that. This helper extracts the requested chain's
amino-acid residues, all atoms, from a deposited CIF into a small PDB that
``ribo_free_pdb`` then points at.

Residue numbers (author ``seqid``) are preserved so the carved PDB aligns with
the truncated ribosome's residue ids by residue number (how
``append_flexible_l24_loop`` maps native contacts onto the appended beads).

The exit-tunnel protein is the uL24 family (universal nomenclature): E. coli
protein L24, eukaryotic RPL26. Known author chains in the bundled raw CIFs:

    E. coli  4v9d.cif : DU     (segID L24)
    yeast    6q8y.cif : AK     (segID L26)
    N.crassa 7r81.cif : a1     (segID L26)
    human    8g61.cif : LY     (segID L26)

Example
-------
    python carve_flexible_protein.py ../raw/ecoli/4v9d.cif DU \\
        ../out/ecoli/L24_atomistic.pdb
"""
import argparse
import sys

import gemmi


def carve_chain(cif_path: str, auth_chain: str, out_pdb: str,
                out_chain: str = "A") -> int:
    """Write the amino-acid residues (all atoms) of ``auth_chain`` to ``out_pdb``.

    Parameters
    ----------
    cif_path : str
        Deposited mmCIF structure.
    auth_chain : str
        Author (``label_asym`` / auth) chain id of the target protein.
    out_pdb : str
        Output all-atom PDB path.
    out_chain : str
        Single-character chain id to stamp on the output (default ``"A"``).
        Only used for the PDB chain column; residue numbers are preserved.

    Returns
    -------
    int
        Number of residues written.

    Raises
    ------
    SystemExit
        If the chain is absent from the model (lists the available chains).
    """
    st = gemmi.read_structure(cif_path)
    model = st[0]

    src = next((c for c in model if c.name == auth_chain), None)
    if src is None:
        avail = ", ".join(sorted({c.name for c in model}))
        raise SystemExit(
            f"chain {auth_chain!r} not found in {cif_path}. Available: {avail}")

    new = gemmi.Structure()
    new.add_model(gemmi.Model("1"))
    dst = gemmi.Chain(out_chain)
    n = 0
    for res in src:
        info = gemmi.find_tabulated_residue(res.name)
        if info is not None and info.is_amino_acid():
            dst.add_residue(res)   # keeps atoms + author seqid
            n += 1
    if n == 0:
        raise SystemExit(
            f"chain {auth_chain!r} in {cif_path} has no amino-acid residues "
            f"(is it an RNA/ligand chain?).")
    new[0].add_chain(dst)
    new.setup_entities()
    new.write_pdb(out_pdb)
    return n


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="Carve one all-atom protein chain from a ribosome CIF.")
    ap.add_argument("cif", help="deposited mmCIF structure")
    ap.add_argument("chain", help="author chain id of the target protein "
                                  "(e.g. DU for E. coli L24 in 4v9d.cif)")
    ap.add_argument("out_pdb", help="output all-atom PDB path")
    ap.add_argument("--out-chain", default="A",
                    help="single-char chain id stamped on output (default A)")
    args = ap.parse_args(argv)

    n = carve_chain(args.cif, args.chain, args.out_pdb, args.out_chain)
    print(f"wrote {n} residues (all atoms) from chain {args.chain} "
          f"-> {args.out_pdb}", file=sys.stderr)


if __name__ == "__main__":
    main()
