#!/usr/bin/env python
"""Extended free-MD ejection demo (FINAL GOAL #3: clean +x egress, no wall leak / clash).

The in-run ``ejection`` phase of ``topo-csp`` is short (20 k steps): far too brief for a
folded nascent chain to diffuse off the truncated ribosome. This standalone demo reloads
the **final synthesized structure** of a completed run and integrates it free (C-terminus
restraint OFF, tunnel wall ON to forbid backward leak) for many steps, writing to
``<outdir>/ejection_long/`` -- the trajectory ``analyze_validation.py`` prefers for the
D5b egress check.

It rebuilds the exact same length-L system + rigid ribosome + equilibrium-PTC geometry as
the run (read from the same INI), so the only difference from synthesis is: no restraint and
a long step budget. Reports C-terminus x, nascent CoM-x and the min nascent--ribosome
distance over the trajectory.

Usage:
    python eject_demo.py -f csp_val.ini --steps 500000 [--L <len>] [--nstout 2000]

``--L`` defaults to the run's L_max (the final length). Use a smaller completed length
(e.g. the debug L=8) to show a short chain fully traversing the tunnel.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit

from topo.csp.protocol import read_csp_config
from topo.csp.core import (precompute_contacts, read_anchor, run_length,
                           optimal_ptc_targets, TUNNEL_AXIS)
from topo.csp.ribosome import load_ribosome, TRNA_TETHER_BOND_NM


def _final_nascent_nm(out_root: Path, L: int) -> np.ndarray:
    """Load the L=`L` stage-3 final nascent coordinates (nm) from a completed run."""
    pdb = out_root / f"L_{L:03d}" / "stage_3" / "traj_final.pdb"
    if not pdb.is_file():
        raise FileNotFoundError(f"no synthesized final structure at {pdb} -- run the "
                                f"synthesis first (topo-csp -f <ini>).")
    coords = []
    for ln in open(pdb):
        if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
            coords.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return np.asarray(coords) / 10.0  # angstrom -> nm


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--config", required=True, help="the same CSP INI the run used")
    ap.add_argument("--steps", type=int, default=500000, help="free-MD ejection steps")
    ap.add_argument("--L", type=int, default=None, help="length to eject (default = L_max)")
    ap.add_argument("--nstout", type=int, default=2000, help="output stride (steps/frame)")
    ap.add_argument("--device", default=None, help="override device (CPU/GPU)")
    args = ap.parse_args()

    cfg = read_csp_config(args.config, verbose=False)
    ep = cfg.params.elong
    ep.trna_tether = False           # match the CSP position-restraint path
    ep.nstout = args.nstout
    if args.device:
        ep.device = args.device

    out_root = Path(cfg.outdir)
    R_full, eps_full = precompute_contacts(cfg.pdb_file, cfg.domain_def,
                                           cfg.stride_output_file)
    N_full = R_full.shape[0]
    L = args.L or (cfg.L_max or N_full)
    final_nm = _final_nascent_nm(out_root, L)
    if final_nm.shape[0] != L:
        raise ValueError(f"final structure has {final_nm.shape[0]} residues but L={L}.")

    ribo = load_ribosome(cfg.ribosome, model="topo")
    p_anchor = read_anchor(cfg.ribosome, "PtR", 76, "R")
    a_anchor = read_anchor(cfg.ribosome, "AtR", 76, "R")

    # Same target / wall geometry the run used (equil-PTC fix path).
    if ep.equil_peptide_geometry:
        a_target, p_target = optimal_ptc_targets(ribo)
        if ep.tunnel_wall:
            ep.tunnel_wall_x0_nm = float(min(a_target[0], p_target[0]))
    else:
        offset = ep.ptc_offset_nm if ep.ptc_offset_nm is not None else TRNA_TETHER_BOND_NM
        p_target = p_anchor + offset * TUNNEL_AXIS
        a_target = a_anchor + offset * TUNNEL_AXIS
        if ep.tunnel_wall:
            ep.tunnel_wall_x0_nm = float(min(p_anchor[0], a_anchor[0]) + offset)

    print(f"=== Extended ejection demo: L={L}, {args.steps} steps, restraint OFF, "
          f"wall x>={ep.tunnel_wall_x0_nm:.3f} nm -> {out_root/'ejection_long'}/ ===")
    run_length(L, full_pdb=cfg.pdb_file, R_full=R_full, eps_full=eps_full,
               p_anchor=p_target, a_anchor=a_anchor, prev_final=None,
               seed_override=final_nm, out_root=out_root, params=ep, ribo=ribo,
               restrain=False, out_subdir="ejection_long",
               n_steps_override=args.steps, label=f"Extended ejection (L={L})")

    # Quick directional summary from the written trajectory.
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings("ignore")
    psf = out_root / "ejection_long" / "traj.psf"
    dcd = out_root / "ejection_long" / "traj.dcd"
    ribo_xyz = ribo.coords_nm * 10.0  # nm -> A
    u = mda.Universe(str(psf), str(dcd))
    nas = u.select_atoms("all")
    cterm_x, com_x, mind = [], [], []
    for _ in u.trajectory:
        p = nas.positions
        cterm_x.append(float(p[-1, 0])); com_x.append(float(p[:, 0].mean()))
        d2 = ((p[:, None, :] - ribo_xyz[None, :, :]) ** 2).sum(2)
        mind.append(float(np.sqrt(d2.min())))
    print(f"\n  frames: {len(cterm_x)}")
    print(f"  C-terminus x (A): {cterm_x[0]:.1f} -> {cterm_x[-1]:.1f} "
          f"(net {cterm_x[-1]-cterm_x[0]:+.1f})")
    print(f"  nascent CoM-x (A): {com_x[0]:.1f} -> {com_x[-1]:.1f} "
          f"(net {com_x[-1]-com_x[0]:+.1f})")
    print(f"  min nascent-ribosome distance (A): start {mind[0]:.2f}, "
          f"min {min(mind):.2f}, end {mind[-1]:.2f}")


if __name__ == "__main__":
    main()
