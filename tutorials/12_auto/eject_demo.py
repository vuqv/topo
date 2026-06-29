#!/usr/bin/env python
"""Extended ejection demo (D5b) from the validated L=10 final structure.

The csp_val.ini ejection phase is short (50 k steps) -- long enough to show the
chain does not collapse into the ribosome, but far too short for a 10-residue chain
to diffuse the ~3-4 nm out of the exit tunnel. This driver continues that final
structure for many more steps with the C-terminus tether released (and the one-sided
tunnel wall still on, so the chain can only move forward, +x), then reports the
egress: nascent CoM-x and the minimum nascent-ribosome distance per frame.

Usage:  python eject_demo.py [n_steps]   (default 2_000_000)
Writes: synth_out/ejection_long/
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import openmm.app as app
from openmm import unit

warnings.filterwarnings("ignore")

from topo.translation.elongate import (read_anchor, TUNNEL_AXIS, TRNA_TETHER_BOND_NM,
                                       precompute_contacts, run_length, ElongationParams)
from topo.translation.ribosome import load_ribosome

HERE = Path(__file__).resolve().parent
N_STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000

full = str(HERE / "4c5c_model_clean.pdb")
ribo_pdb = str(HERE / "ribosome_trunc.pdb")
L = 10

# Match csp_val.ini MD / ribosome knobs.
ep = ElongationParams()
ep.dt_ps = 0.015
ep.ref_t = 310.0
ep.tau_t = 0.01
ep.nstout = 20000
ep.device = "GPU"
ep.constraints = None
ep.buffer_nm = 0.4
ep.minimize = True
ep.rigid_ribosome = True
ep.trna_tether = False          # CSP uses the position-restraint path
ep.tunnel_wall = True
ep.tunnel_wall_x0_nm = 1.05
ep.tunnel_wall_k = 8368.0
ep.nascent_only_output = True

p_anchor = read_anchor(ribo_pdb, "PtR", 76, "R")
a_anchor = read_anchor(ribo_pdb, "AtR", 76, "R")
offset = TRNA_TETHER_BOND_NM
p_target = p_anchor + offset * TUNNEL_AXIS
ribo = load_ribosome(ribo_pdb, model="topo")
R_full, eps_full = precompute_contacts(full, str(HERE / "domain.yaml"),
                                       str(HERE / "4c5c_model_clean_stride.dat"))

# Seed from the validated L=10 final structure.
seed = np.asarray(app.PDBFile(str(HERE / "synth_out" / "L_010" / "stage_3" /
                                   "traj_final.pdb")).getPositions(asNumpy=True
                                   ).value_in_unit(unit.nanometer))[:L]

print(f"=== Extended ejection: L={L}, {N_STEPS} steps, restraint OFF ===")
run_length(L, full_pdb=full, R_full=R_full, eps_full=eps_full,
           p_anchor=p_target, a_anchor=a_anchor, prev_final=None,
           seed_override=seed, out_root=HERE / "synth_out", params=ep, ribo=ribo,
           restrain=False, out_subdir="ejection_long", n_steps_override=N_STEPS,
           label=f"Extended ejection (L={L}, {N_STEPS} steps)")
print("Done -> synth_out/ejection_long/")
