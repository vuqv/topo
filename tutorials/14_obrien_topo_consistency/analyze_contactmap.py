#!/usr/bin/env python3
"""All-atom residue contact map for a protein PDB.

Two residues i, j are in CONTACT if the minimum distance between any of their
(heavy) atoms is <= CUTOFF (default 4.5 A). This is the standard all-atom contact
definition (the input PDB here is heavy-atom only -- no hydrogens). The map is
symmetric; the near-diagonal band is local secondary structure, off-diagonal
clusters are tertiary (long-range) contacts.

Context (tutorials/14): P0CX28's nascent chain is hard to extrude from the
ribosome exit tunnel. A high density of long-range contacts (and a high contact
order) means the chain forms compact tertiary structure that resists threading
N-first down the narrow tunnel -- this script quantifies that.

Usage:
    python analyze_contactmap.py [structure.pdb] [cutoff_A] [min_seq_sep]
Defaults: P0CX28_clean.pdb, 4.5, 1  (min_seq_sep only affects the printed stats).
"""
import sys

import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pdb = sys.argv[1] if len(sys.argv) > 1 else "P0CX28_clean.pdb"
cutoff = float(sys.argv[2]) if len(sys.argv) > 2 else 4.5
min_sep = int(sys.argv[3]) if len(sys.argv) > 3 else 1   # |i-j| >= this counts as a contact

u = mda.Universe(pdb)
residues = u.residues
n = len(residues)
resids = [r.resid for r in residues]
resnames = [r.resname for r in residues]
atom_pos = [r.atoms.positions for r in residues]   # heavy atoms per residue

# Minimum inter-residue heavy-atom distance (A). n is small (~100), so the
# O(n^2) double loop over small atom blocks is fast and exact.
mindist = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        d = cdist(atom_pos[i], atom_pos[j]).min()
        mindist[i, j] = mindist[j, i] = d
np.fill_diagonal(mindist, 0.0)

contact = mindist <= cutoff
# Sequence separation |i-j| for every pair.
ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
sep = np.abs(ii - jj)

# --- statistics (upper triangle, |i-j| >= min_sep) ---
ut = (jj > ii) & (sep >= min_sep)
contacts_ut = contact & ut
n_contacts = int(contacts_ut.sum())
sep_of_contacts = sep[contacts_ut]
# Relative contact order: mean sequence separation of contacts / N (Plaxco et al.).
rco = float(sep_of_contacts.mean() / n) if n_contacts else 0.0
n_local = int((contacts_ut & (sep <= 4)).sum())          # |i-j| 1..4 (helix/turn)
n_medium = int((contacts_ut & (sep > 4) & (sep <= 11)).sum())
n_long = int((contacts_ut & (sep >= 12)).sum())          # |i-j| >= 12 (tertiary)

print(f"Structure : {pdb}")
print(f"Residues  : {n}  ({resnames[0]}{resids[0]} .. {resnames[-1]}{resids[-1]})")
print(f"Cutoff    : {cutoff} A (min heavy-atom distance); contacts with |i-j| >= {min_sep}")
print(f"Contacts  : {n_contacts} total")
print(f"            local  |i-j| 1-4   : {n_local}")
print(f"            medium |i-j| 5-11  : {n_medium}")
print(f"            long   |i-j| >= 12 : {n_long}   <- tertiary structure (tunnel-relevant)")
print(f"Contacts/residue        : {n_contacts / n:.2f}")
print(f"Relative contact order  : {rco:.4f}  (mean |i-j| of contacts = {sep_of_contacts.mean():.1f})")
print(f"Max contact separation  : {int(sep_of_contacts.max())} residues")

# --- plot: contact map coloured by sequence separation ---
fig, ax = plt.subplots(figsize=(7.2, 6.4))
# background min-distance (light) for context, then overlay contacts coloured by |i-j|.
ax.imshow(mindist, origin="lower", cmap="Greys", vmin=0, vmax=20,
          extent=[resids[0] - 0.5, resids[-1] + 0.5, resids[0] - 0.5, resids[-1] + 0.5],
          alpha=0.25)
ci, cj = np.where(np.triu(contact, k=min_sep))
sc = ax.scatter(np.array(resids)[cj], np.array(resids)[ci], c=sep[ci, cj],
                cmap="viridis", s=10, marker="s", vmin=0, vmax=n)
# mirror (lower triangle) for the full symmetric map
ax.scatter(np.array(resids)[ci], np.array(resids)[cj], c=sep[ci, cj],
           cmap="viridis", s=10, marker="s", vmin=0, vmax=n)
cb = fig.colorbar(sc, ax=ax, label="sequence separation |i-j| (residues)")
ax.set_xlabel("residue i"); ax.set_ylabel("residue j")
ax.set_title(f"{pdb}  all-atom contact map (min heavy-atom dist <= {cutoff} A)\n"
             f"{n_contacts} contacts | {n_long} long-range (|i-j|>=12) | RCO={rco:.3f}")
ax.set_aspect("equal")
fig.tight_layout()
out_png = pdb.rsplit(".", 1)[0] + f"_contactmap_{cutoff:g}A.png"
fig.savefig(out_png, dpi=150)
print(f"Wrote {out_png}")

# Also save the contact pair list (i j |i-j| mindist) for downstream use.
out_dat = pdb.rsplit(".", 1)[0] + f"_contacts_{cutoff:g}A.dat"
with open(out_dat, "w") as fh:
    fh.write(f"# all-atom contacts (min heavy-atom dist <= {cutoff} A), |i-j| >= {min_sep}\n")
    fh.write("# resid_i resname_i resid_j resname_j seq_sep mindist_A\n")
    for i, j in zip(*np.where(np.triu(contact, k=min_sep))):
        fh.write(f"{resids[i]:4d} {resnames[i]:>3s} {resids[j]:4d} {resnames[j]:>3s} "
                 f"{abs(resids[i]-resids[j]):4d} {mindist[i,j]:6.2f}\n")
print(f"Wrote {out_dat}")
