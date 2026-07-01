#!/usr/bin/env python
"""Validation analysis for Tutorial 15 / P0CX28 (claude_fix: equilibrium-PTC + AllBonds).

P0CX28 has **no O'Brien reference run**, so D6 (quantitative dwell/geometry match) does NOT
apply. This script reports the criteria that do:

* D5  -- per-stage potential energy stays finite (no >1e12 kJ/mol blow-ups).
* D5b -- ejection: the released chain diffuses out along +x, does not penetrate the tunnel
         wall, and its CoM moves away from the PTC (prefers ejection_long/ if present).
* internal consistency -- the final nascent chain threads the tunnel (monotonic +x) and folds
         to a sensible radius of gyration for a 106-residue protein.

Run from tutorials/15_claude_fix/P0CX28/:  python analyze_validation.py
"""
from __future__ import annotations

import glob
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
OUT = HERE / "synth_out"
# O'Brien's authentic truncated CG ribosome (.cor + sibling .psf/.prm), as the run used.
from topo.csp.core import optimal_ptc_targets
from topo.csp.ribosome import load_ribosome, anchor_coord
RIBO = load_ribosome(str(HERE / "ribosome_trunc.pdb"))
RIBO_COORDS_A = RIBO.coords_nm * 10.0
# Tunnel wall plane: same equil-PTC rule the runner uses (min of the two PTC target x).
_at, _pt = optimal_ptc_targets(RIBO)
TUNNEL_WALL_X0_NM = float(min(_at[0], _pt[0]))
BLOWUP_LIMIT = 1.0e12


def rg(coords_A: np.ndarray) -> float:
    """Radius of gyration (nm) of an (N,3) coordinate array given in Angstrom."""
    c = coords_A - coords_A.mean(axis=0)
    return float(np.sqrt((c ** 2).sum(axis=1).mean()) / 10.0)


def read_pdb_coords(path, names=("CA",)):
    """Read coordinates (Angstrom) of atoms whose name is in `names` from a PDB."""
    coords = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        nm = line[12:16].strip()
        if names is None or nm in names:
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def scan_energies():
    print("=" * 70)
    print(f"D5  Per-stage potential energy (max |PotE| per stage; limit {BLOWUP_LIMIT:g})")
    print("=" * 70)
    worst, worst_stage, n = 0.0, None, 0
    for f in sorted(glob.glob(str(OUT / "**" / "traj.log"), recursive=True)):
        m = 0.0
        for line in open(f):
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    m = max(m, abs(float(parts[2])))
                except ValueError:
                    continue
        n += 1
        if m > worst:
            worst, worst_stage = m, Path(f).relative_to(OUT)
    ok = worst < BLOWUP_LIMIT
    print(f"  scanned {n} stage logs; worst = {worst:.3g} kJ/mol at {worst_stage}")
    print(f"  D5 {'PASS' if ok else 'FAIL'}")
    return ok


def analyze_ejection():
    import MDAnalysis as mda
    print("\n" + "=" * 70)
    print("D5b  Ejection: egress along +x, no tunnel-wall penetration / ribosome clash")
    print("=" * 70)
    if (OUT / "ejection_long" / "traj.dcd").is_file():
        psf, dcd = OUT / "ejection_long" / "traj.psf", OUT / "ejection_long" / "traj.dcd"
        print("  (using extended ejection_long/)")
    else:
        psf, dcd = OUT / "ejection" / "traj.psf", OUT / "ejection" / "traj.dcd"
        print("  (using in-run ejection/; run eject_demo.py for the full egress demo)")
    if not dcd.is_file():
        print("  no ejection trajectory found -- SKIP")
        return None
    ribo = RIBO_COORDS_A
    u = mda.Universe(str(psf), str(dcd))
    nas = u.select_atoms("all")
    com_x, min_x, min_d = [], [], []
    for _ in u.trajectory:
        p = nas.positions
        com_x.append(float(p[:, 0].mean())); min_x.append(float(p[:, 0].min()))
        d2 = ((p[:, None, :] - ribo[None, :, :]) ** 2).sum(axis=2)
        min_d.append(float(np.sqrt(d2.min())))
    com_x, min_x, min_d = map(np.array, (com_x, min_x, min_d))
    x0_A = TUNNEL_WALL_X0_NM * 10.0
    print(f"  frames: {len(com_x)}")
    print(f"  nascent CoM-x (A):  {com_x[0]:.2f} -> {com_x[-1]:.2f} "
          f"(net {com_x[-1]-com_x[0]:+.2f})")
    print(f"  min nascent x: {min_x.min():.2f} A  (tunnel wall x0 = {x0_A:.2f} A)")
    print(f"  min nascent-ribosome distance: {min_d.min():.2f} A")
    wall_ok = min_x.min() >= x0_A - 1.0
    clash_ok = min_d.min() >= 3.0
    egress_ok = com_x[-1] > com_x[0]
    print(f"  wall not penetrated: {'PASS' if wall_ok else 'FAIL'}; "
          f"no ribosome clash: {'PASS' if clash_ok else 'FAIL'}; "
          f"net +x egress: {'PASS' if egress_ok else 'FAIL'}")
    return dict(wall_ok=wall_ok, clash_ok=clash_ok, egress_ok=egress_ok)


def internal_consistency():
    print("\n" + "=" * 70)
    print("Internal consistency (no O'Brien reference for P0CX28; D6 N/A)")
    print("=" * 70)
    finals = sorted(glob.glob(str(OUT / "L_*" / "stage_3" / "traj_final.pdb")))
    if not finals:
        print("  no final structures found -- SKIP")
        return None
    last = finals[-1]
    nas = read_pdb_coords(last, names=None)
    L = len(nas)
    corr = float(np.corrcoef(np.arange(L), nas[:, 0])[0, 1])
    ribo = RIBO_COORDS_A
    below = int((nas[:, 0] < 0).sum())
    print(f"  final length L = {L}")
    print(f"  threads tunnel: corr(residue index, x) = {corr:.3f} "
          f"(expect strongly negative: N-term extruded, C-term at PTC)")
    print(f"  no collapse into ribosome interior: {below} beads at x<0")
    print(f"  final nascent R_g = {rg(nas):.3f} nm "
          f"(sanity for a {L}-residue chain)")
    return dict(L=L, corr=corr, rg=rg(nas))


if __name__ == "__main__":
    d5 = scan_energies()
    d5b = analyze_ejection()
    ic = internal_consistency()
    print("\n" + "=" * 70 + "\nSUMMARY\n" + "=" * 70)
    print(f"  D5  finite energies : {'PASS' if d5 else 'FAIL'}")
    if d5b:
        print(f"  D5b ejection        : wall {'OK' if d5b['wall_ok'] else 'FAIL'}, "
              f"clash {'OK' if d5b['clash_ok'] else 'FAIL'}, "
              f"egress {'OK' if d5b['egress_ok'] else 'FAIL'}")
    if ic:
        print(f"  internal consistency: L={ic['L']}, threads corr={ic['corr']:.3f}, "
              f"R_g={ic['rg']:.3f} nm")
