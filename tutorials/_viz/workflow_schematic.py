# Draws the workflow schematic of tutorials/A0_start_here.
# Usage: python tutorials/_viz/workflow_schematic.py tutorials/A0_start_here/img/workflow_schematic.png
import sys
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.8,
                            "axes.spines.top": False, "axes.spines.right": False})

import numpy as np, textwrap, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
S_BASE, S_BODY, S_TICK = 9.5, 7.5, 6.5
MONO = "DejaVu Sans Mono"
C_HEAD="#2f4b6e"; C_BOX="#f4f5f7"; C_EDGE="#b4bac4"; C_RUN="#ffffff"; C_RUNEDGE="#2f4b6e"
C_CHK="#fdf3d8"; C_CHKEDGE="#d2ad4a"; C_TXT="#1c2530"; C_GREY="#5d6673"; C_ARROW="#2f4b6e"
C_A="#0072B2"; C_B="#E69F00"
W, H = 13.0, 9.88
fig = plt.figure(figsize=(W, H))
cv = fig.add_axes([0,0,1,1]); cv.set_xlim(0,W); cv.set_ylim(0,H); cv.axis('off')
boxes = {}; owned = []
def box(name, x, y, w, h, fc, ec, lw=0.8, r=0.05):
    cv.add_patch(FancyBboxPatch((x+r, y+r), w-2*r, h-2*r, boxstyle=f"round,pad={r}", fc=fc, ec=ec, lw=lw, zorder=2))
    boxes[name] = (x, y, w, h)
def txt(name, x, y, s, size=S_BODY, mono=False, color=C_TXT, ha='left', va='top', weight='normal', ls=1.32, **kw):
    t = cv.text(x, y, s, fontsize=size, family=(MONO if mono else None), color=color, ha=ha, va=va,
                fontweight=weight, linespacing=ls, zorder=5, **kw)
    owned.append((t, name)); return t
def bullets(items, width):
    out=[]
    for it in items:
        lines = textwrap.wrap(it, width)
        out.append("\n".join(["• "+lines[0]] + ["  "+l for l in lines[1:]]))
    return "\n".join(out)
X0, CW, GAP = 0.62, 2.92, 0.20
colx = [X0 + i*(CW+GAP) for i in range(4)]
PADX, PADY = 0.09, 0.08
y_title = H - 0.22
band = (H-1.12, 0.50)
rows = {}; ytop = band[0] - 0.14
for name, h in [("head",0.40),("meaning",0.88),("sketch",1.15),("provide",1.70),("run",1.84),("get",0.70),("check",1.04)]:
    rows[name] = (ytop - h, h); ytop -= h + 0.10
cv.text(X0-0.37, y_title, "Refolding simulation of a globular protein with TOPO: what to provide, run and check at each stage",
        fontsize=S_BASE+1.5, color=C_TXT, ha='left', va='top', fontweight='bold', zorder=5)
by, bh = band
box("band", X0-0.37, by, colx[3]+CW-(X0-0.37), bh, "#e9edf3", C_EDGE)
txt("band", X0-0.27, by+bh/2, "0 · Set up once", size=S_BASE, weight='bold', color=C_HEAD, va='center')
chips = ["conda env: openmm parmed\nmdanalysis mdtraj", "pip install -e .\n(in the topo repository)",
         "scripts/install_deps.sh\n(installs STRIDE)", "topo-mdrun, topo-optimize\nprint their help"]
cx = X0 + 1.72; cw_ = [2.35, 2.25, 2.2, 2.35]; cg = 0.36
for i,(c,wd) in enumerate(zip(chips,cw_)):
    box(f"chip{i}", cx, by+0.06, wd, bh-0.12, C_RUN, C_EDGE, lw=0.6)
    txt(f"chip{i}", cx+wd/2, by+bh/2, c, mono=True, ha='center', va='center', ls=1.2)
    if i < 3:
        cv.add_patch(FancyArrowPatch((cx+wd+0.04, by+bh/2), (cx+wd+cg-0.04, by+bh/2), arrowstyle='-|>', mutation_scale=9, color=C_ARROW, lw=1.0, zorder=4))
    cx += wd + cg
for name, lab in [("meaning","What it\ndoes"),("provide","You provide"),("run","You run"),("get","You get"),("check","Check before\nmoving on")]:
    y,h = rows[name]
    cv.text(X0-0.20, y+h/2, lab, rotation=90, fontsize=S_BODY, color=C_GREY, ha='center', va='center', fontweight='bold', linespacing=1.1, zorder=5)
heads = [("1 · Structure and domains","by hand"), ("2 · Optimize nscale","GPU job"),
         ("3 · Temperature quench","GPU job"), ("4 · Analyze native contacts Q","CPU")]
provide = [
 ("sans", "All-atom structure of the protein,\none chain, saved as <PROTEIN_STRUCTURE>.pdb\n\nDomains of this protein (for example\nfrom TED): residue ranges and class of\neach domain (alpha, beta or alpha-beta)\n\nNote: assign every residue to a domain;\ncurate automatic assignments by hand."),
 ("mono", "<PROTEIN_STRUCTURE>.pdb\ndomain.yaml\noptimize.ini:\n  pdb_file   = <PROTEIN_STRUCTURE>.pdb\n  domain_def = domain.yaml\n  ref_t      = ____\n  md_steps   = ____\n  nstxout    = ____\n  ntraj      = ____"),
 ("mono", "<PROTEIN_STRUCTURE>.pdb\nmd.ini:\n  domain_def =\n      opt_out/domain_optimized.yaml\n  anneal = yes      anneal_ramp = jump\n  t_high       = ____\n  anneal_steps = ____\n  ref_t        = ____\n  md_steps     = ____\n  tau_t        = ____\n  n_copies     = ____"),
 ("mono", "<PROTEIN_STRUCTURE>.pdb\nopt_out/domain_optimized.yaml\ntraj.psf\ntraj_quench.dcd\ntraj.dcd"),
]
EXAMPLE_YAML = "n_residues: 139\nintra_domains:\n  A: {residues: [1-90], class: beta,\n      nscale: 1.0}       # placeholder\n  B: {residues: [91-139], class: alpha,\n      nscale: 1.0}\ninter_domains: {A-B: 1.0}"
run = [
 "write domain.yaml:\n" + EXAMPLE_YAML,
 "topo-optimize -f optimize.ini -o opt_out\n\n# watch progress:\ntail -f opt_out/optimization.log\n\n# a stopped optimization cannot resume;\n# it starts again from round 1",
 "topo-mdrun -f md.ini\n\n\n\n\n# job stopped at the walltime during\n# refolding: set restart = yes in md.ini\n# and submit again; the hold is skipped",
 "python -m topo.utils.multichain \\\n  -f traj.dcd -tp traj.psf -o split\n\n\n\n\npython -m topo.analysis.native_contacts \\\n  -d opt_out/domain_optimized.yaml \\\n  -r <PROTEIN_STRUCTURE>.pdb \\\n  -p traj.psf -f split/traj_k.dcd \\\n  -o Q_k.csv\n# then the same for traj_quench.dcd",
]
get = [
 "<PROTEIN_STRUCTURE>.pdb\ndomain.yaml",
 "opt_out/domain_optimized.yaml\nopt_out/optimization.log\nopt_out/round_N/traj/Q_k.csv",
 "traj_quench.dcd  traj_quench.log\ntraj.dcd         traj.log\ntraj.psf  traj_multi.psf  traj.chk",
 "Q_k.csv for each copy k, columns:\nFrame, Q_protein, Q_<domain>,\nQ_<domain1>-<domain2>",
]
check = [
 ["n_residues equals the number of residues in <PROTEIN_STRUCTURE>.pdb", "every residue belongs to a named domain", "each class is alpha, beta or alpha-beta"],
 ["optimization.log ends with CONVERGED and no line of the last round is tagged fallback", "if it ends with a WARNING, the unstable domains or interfaces carry a fallback or untested nscale; record them", "the log shows your ntraj, ref_t and steps/traj"],
 ["mean temperature in traj_quench.log is close to t_high", "mean temperature in traj.log is close to ref_t", "the run ends with 'Finished in ... seconds'"],
 ["Q is near 0 at the end of the hold; if not, raise t_high or anneal_steps and repeat stage 3", "one Q_k.csv per copy for hold and for refolding"],
]
for i,x in enumerate(colx):
    y,h = rows["head"]; box(f"head{i}", x, y, CW, h, C_HEAD, C_HEAD)
    txt(f"head{i}", x+PADX+0.02, y+h/2, heads[i][0], size=S_BASE, weight='bold', color='white', va='center')
    txt(f"head{i}", x+CW-PADX-0.02, y+h/2, heads[i][1], size=S_TICK, color='#d5dde8', va='center', ha='right')
    y,h = rows["sketch"]; box(f"sketch{i}", x, y, CW, h, "#ffffff", C_EDGE, lw=0.6)
    y,h = rows["provide"]; box(f"provide{i}", x, y, CW, h, C_BOX, C_EDGE)
    kind, s = provide[i]; txt(f"provide{i}", x+PADX, y+h-PADY, s, mono=(kind=="mono"))
    y,h = rows["run"]; box(f"run{i}", x, y, CW, h, C_RUN, C_RUNEDGE, lw=1.3)
    txt(f"run{i}", x+PADX, y+h-PADY, run[i], mono=True)
    y,h = rows["get"]; box(f"get{i}", x, y, CW, h, C_BOX, C_EDGE)
    txt(f"get{i}", x+PADX, y+h-PADY, get[i], mono=True)
    y,h = rows["check"]; box(f"check{i}", x, y, CW, h, C_CHK, C_CHKEDGE)
    txt(f"check{i}", x+PADX, y+h-PADY, bullets(check[i], 50))

meaning = [
 "Prepare the native structure and assign every residue to a domain. The native structure defines the contacts that hold the model protein together.",
 "Find the minimum scale factor (nscale) of the sidechain-sidechain contact energies, within each domain and between domains, that keeps the protein stable at ref_t.",
 "Raise the system to a high temperature to generate unfolded conformations, then quench to the target temperature (ref_t) to see how the protein folds from the unfolded state to the native state.",
 "Compute Q, the fraction of native contacts formed in each frame (0 unfolded, 1 native), to confirm that the hold unfolded the protein and to follow each copy as it refolds.",
]
for i,x in enumerate(colx):
    y,h = rows["meaning"]; box(f"meaning{i}", x, y, CW, h, "#e9edf3", "#e9edf3", lw=0.6)
    txt(f"meaning{i}", x+PADX, y+h-PADY, "\n".join(textwrap.wrap(meaning[i], 46)))
cv.text(colx[0], rows["check"][0]-0.09, "Names in <ANGLE BRACKETS> and blanks ____ are placeholders: replace them with your own file name or value.",
        fontsize=S_TICK, color=C_GREY, ha='left', va='top', zorder=5)

# highlighted note in the stage 3 "You run" box
y,h = rows["run"]; LH = S_BODY*1.32/72
txt("run2", colx[2]+PADX+0.04, y+h-PADY-LH-0.13, "Simulates n_copies non-interacting copies\nof the protein in the same simulation.",
    weight='bold', ls=1.25, bbox=dict(boxstyle="round,pad=0.35", fc="#ffe98a", ec="none"))

# highlighted note in the stage 4 "You run" box: what the multichain split does
txt("run3", colx[3]+PADX+0.04, y+h-PADY-2*LH-0.13, "Splits the combined trajectory (traj.dcd)\ninto one trajectory per copy: traj_k.dcd",
    weight='bold', ls=1.25, bbox=dict(boxstyle="round,pad=0.35", fc="#ffe98a", ec="none"))
for i in range(3):
    xa = colx[i]+CW; xb = colx[i+1]; xm = (xa+xb)/2
    yg = rows["get"][0] + rows["get"][1]*0.5; yp = rows["provide"][0] + rows["provide"][1]*0.5
    cv.plot([xa, xm, xm], [yg, yg, yp], color=C_ARROW, lw=1.5, zorder=4, solid_capstyle='round', solid_joinstyle='round')
    cv.add_patch(FancyArrowPatch((xm-0.01, yp), (xb+0.005, yp), arrowstyle='-|>', mutation_scale=10, color=C_ARROW, lw=1.5, zorder=4, shrinkA=0, shrinkB=0))
yc = rows["check"][0] + rows["check"][1]*0.60
cv.add_patch(FancyArrowPatch((colx[3]-0.005, yc), (colx[2]+CW+0.005, yc), arrowstyle='-|>', mutation_scale=10, color=C_CHKEDGE, lw=1.4, ls=(0,(2.5,1.6)), zorder=4, shrinkA=0, shrinkB=0))
def inset(i, l=0.12, r=0.12, b=0.10, t=0.10):
    y,h = rows["sketch"]; x = colx[i]
    return fig.add_axes([(x+l)/W, (y+b)/H, (CW-l-r)/W, (h-b-t)/H])
ax = inset(0, l=0.18, r=0.18, b=0.30, t=0.30); ax.set_xlim(0,140); ax.set_ylim(0,1)
ax.add_patch(Rectangle((1,0.15),89,0.7, fc=C_A, ec='none')); ax.add_patch(Rectangle((90,0.15),49,0.7, fc=C_B, ec='none'))
ax.text(45.5,0.5,"domain A · beta", color='white', ha='center', va='center', fontsize=S_BODY, fontweight='bold')
ax.text(114.5,0.5,"B · alpha", color=C_TXT, ha='center', va='center', fontsize=S_BODY, fontweight='bold')
ax.set_xticks([1,90,139]); ax.set_yticks([]); ax.tick_params(labelsize=S_TICK, length=2.5, pad=1.5)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel("residue number (example: P0A6E6, tutorial A5)", fontsize=S_TICK, labelpad=1.5, color=C_GREY)
ax.set_title("one bead per residue; contacts are scaled per domain", fontsize=S_TICK, color=C_GREY, pad=3, loc='left')
y,h = rows["sketch"]; x = colx[1]
bw, bh2 = 1.18, 0.40; xa_, xb_ = x+0.16, x+CW-0.16-bw; yb_ = y+h-0.14-bh2
box("loopA", xa_, yb_, bw, bh2, C_BOX, C_EDGE, lw=0.6); txt("loopA", xa_+bw/2, yb_+bh2/2, "ntraj copies\nat ref_t", size=S_TICK, ha='center', va='center', ls=1.15)
box("loopB", xb_, yb_, bw, bh2, C_BOX, C_EDGE, lw=0.6); txt("loopB", xb_+bw/2, yb_+bh2/2, "Q of each domain\nand interface", size=S_TICK, ha='center', va='center', ls=1.15)
cv.add_patch(FancyArrowPatch((xa_+bw+0.03, yb_+bh2/2), (xb_-0.03, yb_+bh2/2), arrowstyle='-|>', mutation_scale=8, color=C_ARROW, lw=1.0, zorder=4))
yl = y+0.33
cv.plot([xb_+bw/2, xb_+bw/2, xa_+bw/2], [yb_-0.02, yl, yl], color=C_ARROW, lw=1.0, zorder=4)
cv.add_patch(FancyArrowPatch((xa_+bw/2, yl-0.01), (xa_+bw/2, yb_-0.01), arrowstyle='-|>', mutation_scale=8, color=C_ARROW, lw=1.0, zorder=4, shrinkA=0, shrinkB=0))
txt("sketch1", x+CW/2, yl+0.05, "unstable: next nscale level", size=S_TICK, ha='center', va='bottom', color=C_TXT)
txt("sketch1", x+CW/2, y+0.07, "all stable: write domain_optimized.yaml", size=S_TICK, ha='center', va='bottom', color=C_GREY)
t1 = 0.44
ax3 = inset(2, l=0.52, r=0.14, b=0.27, t=0.22)
ax3.plot([0,t1,t1,1],[1,1,0,0], color=C_HEAD, lw=1.8, solid_joinstyle='miter')
ax3.axvspan(0,t1, color=C_B, alpha=0.16, lw=0)
ax3.set_xlim(0,1); ax3.set_ylim(-0.42,1.30); ax3.set_yticks([0,1]); ax3.set_yticklabels(["ref_t","t_high"], fontsize=S_TICK)
ax3.set_xticks([t1/2,(1+t1)/2]); ax3.set_xticklabels(["anneal_steps","md_steps"], fontsize=S_TICK); ax3.tick_params(length=0, pad=2)
ax3.text(t1/2, 0.42, "traj_quench.dcd", ha='center', va='center', fontsize=S_TICK, family=MONO, color=C_TXT)
ax3.text((1+t1)/2, 0.42, "traj.dcd", ha='center', va='center', fontsize=S_TICK, family=MONO, color=C_TXT)
ax3.text(-0.19, 1.44, "temperature (sketch)", fontsize=S_TICK, color=C_GREY, ha='left', va='bottom', clip_on=False)
ax4 = inset(3, l=0.36, r=0.14, b=0.27, t=0.22)
tt = np.linspace(0,1,500)
ax4.axvspan(0,t1, color=C_B, alpha=0.16, lw=0)
for t0, k in [(0.55,45),(0.64,40),(0.75,50),(0.87,42)]:
    q = np.where(tt<t1, 0.05+0.95*np.exp(-tt/0.04), 0.05+0.87/(1+np.exp(-k*(tt-t0))))
    ax4.plot(tt, q, color=C_A, lw=1.1, alpha=0.85)
ax4.set_xlim(0,1); ax4.set_ylim(-0.08,1.12); ax4.set_yticks([0,1]); ax4.set_yticklabels(["0","1"], fontsize=S_TICK)
ax4.set_xticks([t1/2,(1+t1)/2]); ax4.set_xticklabels(["hold at t_high","refolding at ref_t"], fontsize=S_TICK); ax4.tick_params(length=0, pad=2)
ax4.text(-0.115, 1.21, "Q, one line per copy (sketch)", fontsize=S_TICK, color=C_GREY, ha='left', va='bottom', clip_on=False)
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else "workflow_schematic.png", dpi=200)
