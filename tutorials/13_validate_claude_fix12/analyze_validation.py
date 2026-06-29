#!/usr/bin/env python
"""Validation analysis for Tutorial 12_auto (O'Brien CSP reproduced in topo-style).

Checks the Goal.md Definition-of-Done criteria against the topo-csp output in
``synth_out/`` and the reference run in ``continuous_synthesis/output/``:

* D5   -- per-stage potential energy stays finite (no >1e12 kJ/mol blow-ups).
* D5b  -- ejection: the released chain diffuses out along +x, never penetrates the
          tunnel wall / overlaps a ribosome bead, and its CoM moves away from the PTC.
* D6   -- quantitative consistency vs the reference (length, dwell-time totals,
          per-codon means, final radius of gyration).

Run from tutorials/12_auto/:  python analyze_validation.py
Outputs a human-readable report to stdout (captured into NOTES.md).
"""
from __future__ import annotations

import glob
import re
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
OUT = HERE / "synth_out"
REF = HERE / "continuous_synthesis" / "output"
RIBO_PDB = HERE / "ribosome_trunc.pdb"
TUNNEL_WALL_X0_NM = 1.05            # csp_val.ini tunnel_wall_x0 (nm)
BLOWUP_LIMIT = 1.0e12              # D5 threshold (kJ/mol)


def rg(coords_A: np.ndarray) -> float:
    """Radius of gyration (nm) of an (N,3) coordinate array given in Angstrom."""
    c = coords_A - coords_A.mean(axis=0)
    return float(np.sqrt((c ** 2).sum(axis=1).mean()) / 10.0)


# --------------------------------------------------------------------------
# D5: scan every stage's potential-energy log
# --------------------------------------------------------------------------
def scan_energies():
    print("=" * 70)
    print("D5  Per-stage potential energy (max |PotE| per stage; limit "
          f"{BLOWUP_LIMIT:g} kJ/mol)")
    print("=" * 70)
    worst = 0.0
    worst_stage = None
    n = 0
    for f in sorted(glob.glob(str(OUT / "**" / "traj.log"), recursive=True)):
        m = 0.0
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        v = abs(float(parts[2]))
                    except ValueError:
                        continue
                    m = max(m, v)
        n += 1
        rel = Path(f).relative_to(OUT)
        if m > worst:
            worst, worst_stage = m, rel
        flag = "  <<< BLOWUP" if m > BLOWUP_LIMIT else ""
        if m > BLOWUP_LIMIT or m > 1e4:
            print(f"  {str(rel):40s} max|PotE| = {m:.3g}{flag}")
    print(f"  scanned {n} stage logs; worst = {worst:.3g} kJ/mol at {worst_stage}")
    ok = worst < BLOWUP_LIMIT
    print(f"  D5 {'PASS' if ok else 'FAIL'} (worst {worst:.3g} "
          f"{'<' if ok else '>='} {BLOWUP_LIMIT:g})")
    return ok


# --------------------------------------------------------------------------
# D5b: ejection trajectory analysis
# --------------------------------------------------------------------------
def read_pdb_coords(path, names=("CA",)):
    """Read coordinates (Angstrom) of atoms whose name is in `names` from a PDB."""
    coords = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            nm = line[12:16].strip()
            if names is None or nm in names:
                coords.append([float(line[30:38]), float(line[38:46]),
                               float(line[46:54])])
    return np.array(coords)


def analyze_ejection():
    import MDAnalysis as mda
    print()
    print("=" * 70)
    print("D5b  Ejection: egress along +x, no tunnel-wall penetration / ribosome clash")
    print("=" * 70)
    # Prefer the extended ejection demo (eject_demo.py) if present -- a 10-bead chain
    # needs far longer than the 50 k in-run steps to traverse the ~4 nm tunnel.
    if (OUT / "ejection_long" / "traj.dcd").is_file():
        psf = OUT / "ejection_long" / "traj.psf"
        dcd = OUT / "ejection_long" / "traj.dcd"
        print("  (using extended ejection_long/; run eject_demo.py to (re)generate)")
    else:
        psf = OUT / "ejection" / "traj.psf"
        dcd = OUT / "ejection" / "traj.dcd"
        print("  (using in-run ejection/; 50 k steps -- too short for full egress, "
              "run eject_demo.py for the egress demo)")
    if not dcd.is_file():
        print("  no ejection trajectory found -- SKIP")
        return None
    # Ribosome beads (Angstrom) -- everything in ribosome_trunc.pdb.
    ribo = read_pdb_coords(RIBO_PDB, names=None)
    u = mda.Universe(str(psf), str(dcd))
    nas = u.select_atoms("all")
    n_frames = len(u.trajectory)
    com_x, min_x, min_d = [], [], []
    for ts in u.trajectory:
        p = nas.positions  # Angstrom
        com_x.append(float(p[:, 0].mean()))
        min_x.append(float(p[:, 0].min()))
        # min distance from any nascent bead to any ribosome bead
        d2 = ((p[:, None, :] - ribo[None, :, :]) ** 2).sum(axis=2)
        min_d.append(float(np.sqrt(d2.min())))
    com_x = np.array(com_x); min_x = np.array(min_x); min_d = np.array(min_d)
    x0_A = TUNNEL_WALL_X0_NM * 10.0
    # ribosome excluded-volume bead radius ~ sigma ~ a few A; a clash would be
    # min_d well under ~4 A. We report the min over the whole trajectory.
    print(f"  frames: {n_frames}")
    print(f"  nascent CoM-x (A):  start {com_x[0]:.2f}  ->  end {com_x[-1]:.2f}  "
          f"(net +x displacement {com_x[-1] - com_x[0]:+.2f} A)")
    print(f"  min nascent x over run: {min_x.min():.2f} A  (tunnel wall x0 = {x0_A:.2f} A)")
    print(f"  min nascent-ribosome distance over run: {min_d.min():.2f} A "
          f"(min at frame {int(min_d.argmin())})")
    # CoM-x trend: fraction of frames moving forward, and overall slope.
    dx = np.diff(com_x)
    frac_fwd = float((dx > 0).mean()) if len(dx) else float("nan")
    slope = float(np.polyfit(np.arange(n_frames), com_x, 1)[0]) if n_frames > 1 else 0.0
    print(f"  CoM-x monotonic-ish +x: {100 * frac_fwd:.0f}% of steps advance, "
          f"linear slope {slope:+.3f} A/frame")
    wall_ok = min_x.min() >= x0_A - 1.0    # allow ~1 A numerical slack at the wall
    clash_ok = min_d.min() >= 3.0          # no hard overlap with a ribosome bead
    egress_ok = com_x[-1] > com_x[0]
    print(f"  wall not penetrated: {'PASS' if wall_ok else 'FAIL'}; "
          f"no ribosome clash: {'PASS' if clash_ok else 'FAIL'}; "
          f"net +x egress: {'PASS' if egress_ok else 'FAIL'}")
    return dict(com_x=com_x, min_x=min_x, min_d=min_d, wall_ok=wall_ok,
                clash_ok=clash_ok, egress_ok=egress_ok)


# --------------------------------------------------------------------------
# D6: quantitative comparison vs the reference
# --------------------------------------------------------------------------
def parse_reference_out():
    """Parse continuous_synthesis/output/output/1.out -> per-length dict."""
    path = REF / "output" / "1.out"
    text = path.read_text()
    blocks = re.split(r"--> Elongation at length ", text)[1:]
    data = {}
    for b in blocks:
        L = int(b.split(maxsplit=1)[0])
        def grab(pat):
            m = re.search(pat, b)
            return float(m.group(1)) if m else None
        mean_pt = grab(r"Mean in vivo peptidyl transfer dwell time:\s*([-\d.eE]+)")
        mean_tl = grab(r"Mean in vivo translocation dwell time:\s*([-\d.eE]+)")
        mean_tr = grab(r"Mean in vivo tRNA binding dwell time:\s*([-\d.eE]+)")
        ns1 = grab(r"Sampled in silico peptidyl transfer dwell time:\s*([-\d.eE]+)")
        ns2 = grab(r"Sampled in silico translocation dwell time:\s*([-\d.eE]+)")
        ns3 = grab(r"Sampled in silico tRNA binding dwell time:\s*([-\d.eE]+)")
        s1 = grab(r"steps for in silico dwell time before peptidyl transfer:\s*(\d+)")
        s2 = grab(r"steps for in silico dwell time before translocation:\s*(\d+)")
        s3 = grab(r"steps for in silico dwell time before next tRNA binding:\s*(\d+)")
        data[L] = dict(mean_pt=mean_pt, mean_tl=mean_tl, mean_tr=mean_tr,
                       ns=(ns1, ns2, ns3), steps=(s1, s2, s3))
    return data


def parse_our_dwell():
    rows = {}
    with open(OUT / "dwell_times.dat") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            L = int(p[0])
            rows[L] = dict(codon=p[1], t_total=float(p[2]),
                           t=(float(p[3]), float(p[4]), float(p[5])),
                           ns=(float(p[6]), float(p[7]), float(p[8])),
                           steps=(int(p[9]), int(p[10]), int(p[11])))
    return rows


def analyze_quantitative():
    print()
    print("=" * 70)
    print("D6  Quantitative consistency vs reference (continuous_synthesis/output/)")
    print("=" * 70)
    ref = parse_reference_out()
    ours = parse_our_dwell()
    Ls = sorted(set(ref) & set(ours))
    # (a) length
    print(f"  (a) synthesized length: ours max L = {max(ours)}, ref max L = {max(ref)}")
    # (b) per-codon in-vivo MEAN times (deterministic; should match exactly since
    #     same mRNA + trans_times + stage means). Compare total codon dwell.
    print("  (b) per-residue in-vivo TOTAL dwell (s)  [mean = t1+t2+t3 means]:")
    print(f"      {'L':>2} {'codon':>5} {'ours t_total':>13} {'ref mean_sum':>13}")
    ours_total = 0.0
    ref_total = 0.0
    for L in Ls:
        ref_mean_sum = (ref[L]['mean_pt'] or 0) + (ref[L]['mean_tl'] or 0) + (ref[L]['mean_tr'] or 0)
        # ref mean_tl can be negative (translocation correction); the in-vivo
        # codon total is mean_pt + mean_tl + mean_tr by construction = intrinsic[L].
        o = ours[L]
        ours_sampled_sum = sum(o['t'])
        print(f"      {L:>2} {o['codon']:>5} {o['t_total']:13.4e} {ref_mean_sum:13.4e}")
        ours_total += ours_sampled_sum
        ref_total += sum(x for x in ref[L]['ns'])   # placeholder, replaced below
    # Totals in in-silico ns (sampled) -- the stochastic quantity to compare.
    ours_ns = sum(sum(ours[L]['ns']) for L in Ls)
    ref_ns = sum(sum(ref[L]['ns']) for L in Ls)
    print(f"  total SAMPLED in-silico dwell (ns):  ours {ours_ns:.3f}  ref {ref_ns:.3f}  "
          f"(ratio {ours_ns / ref_ns:.2f}x)")
    ours_invivo = sum(ours[L]['t_total'] for L in Ls)
    # reference in-vivo total = sum of intrinsic means
    ref_invivo = sum((ref[L]['mean_pt'] or 0) + (ref[L]['mean_tl'] or 0) +
                     (ref[L]['mean_tr'] or 0) for L in Ls)
    print(f"  total in-vivo dwell (s):  ours(sampled) {ours_invivo:.4f}  "
          f"ref(mean) {ref_invivo:.4f}  (ratio {ours_invivo / ref_invivo:.2f}x)")
    # (c) final nascent Rg
    print("  (c) final nascent radius of gyration (L=10):")
    our_final = read_pdb_coords(OUT / "L_010" / "stage_3" / "traj_final.pdb", names=None)
    our_rg = rg(our_final)
    ref_rg = ref_nascent_rg()
    print(f"      ours  Rg = {our_rg:.3f} nm ({len(our_final)} beads)")
    print(f"      ref   Rg = {ref_rg:.3f} nm")
    print(f"      ratio = {our_rg / ref_rg:.2f}x")
    return dict(ours_ns=ours_ns, ref_ns=ref_ns, ours_invivo=ours_invivo,
                ref_invivo=ref_invivo, our_rg=our_rg, ref_rg=ref_rg)


def ref_nascent_rg():
    import MDAnalysis as mda
    psf = REF / "traj" / "1" / "rnc_l10.psf"
    cor = REF / "traj" / "1" / "rnc_l10_stage_3_final.cor"
    u = mda.Universe(str(psf), topology_format="PSF")
    u.load_new(str(cor), format="CRD")
    nas = u.select_atoms("segid A")
    return rg(nas.positions)


if __name__ == "__main__":
    d5 = scan_energies()
    d5b = analyze_ejection()
    d6 = analyze_quantitative()
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  D5  finite energies : {'PASS' if d5 else 'FAIL'}")
    if d5b:
        print(f"  D5b ejection        : wall {'OK' if d5b['wall_ok'] else 'FAIL'}, "
              f"clash {'OK' if d5b['clash_ok'] else 'FAIL'}, "
              f"egress {'OK' if d5b['egress_ok'] else 'FAIL'}")
