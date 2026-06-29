"""Continuous Synthesis Protocol (O'Brien), ported to topo (``topo.csp``).

This is the **full O'Brien continuous-synthesis protocol** -- the per-codon, 3-stage
elongation cycle of ``continuous_synthesis_v6.py`` -- expressed in topo style. It is
the kinetic upgrade of :mod:`topo.translation.elongate`: that driver grows the chain
one residue per step at a *fixed* ``n_steps`` (collapsing the elongation cycle into a
single MD segment); **this** driver times every residue from its codon and splits it
into three sub-stages, exactly as O'Brien do.

What is reused vs. new:

- **Reused wholesale** from :mod:`topo.translation.elongate`: the per-length MD
  machinery -- :func:`run_length` (build-once-subset contacts, coordinate seeding,
  rigid-ribosome scenery, tunnel wall, minimize/run/finalize), :func:`read_anchor`,
  :func:`precompute_contacts`, and :class:`ElongationParams` for all the MD/ribosome
  knobs. Nothing about the force field or the OpenMM plumbing is re-implemented here.
- **New** (this module): the O'Brien *kinetics* (:mod:`topo.csp.kinetics`) and the
  outer loop that calls :func:`run_length` **three times per residue**, switching the
  C-terminus restraint target A->P to reproduce translocation.

The 3-stage mapping onto :func:`run_length` (``L`` = nascent length):

==========  =====================================  ==========================
stage       biology                                restraint target / seed
==========  =====================================  ==========================
1           peptidyl transfer / A-site delivery    A-anchor; new residue placed
2           translocation begins                   A-anchor; continue stage 1
3           translocation to P-site / wait         P-anchor; continue stage 2
==========  =====================================  ==========================

Stage 3's final structure seeds the next residue's stage 1. The cold-start segment
(``L == L0``) is laid down the tunnel from the P-anchor (no A-site delivery yet).
Because CSP needs the restraint target to switch A<->P, it drives the **position
restraint** path (``trna_tether`` is forced off); ``rigid_ribosome`` / ``tunnel_wall``
/ excluded volume + electrostatics still work exactly as in build step v2.

Drive it with an INI control file (see :func:`read_csp_config`)::

    topo-csp -f csp.ini
    python -m topo.csp -f csp.ini

**Units:** OpenMM defaults -- length nm, time ps, energy kJ/mol, temperature K,
force constants kJ/mol/nm^2. In-vivo dwell times in the kinetics are **seconds**.
"""
from __future__ import annotations

import argparse
import configparser
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np
import openmm as mm

from topo.translation.elongate import (ElongationParams, TUNNEL_AXIS,
                                       precompute_contacts, read_anchor,
                                       run_length, TRNA_TETHER_BOND_NM)
from topo.translation.ribosome import load_ribosome
from topo.utils.config import strtobool
from topo.csp import kinetics


# --------------------------------------------------------------------------
# Run parameters
# --------------------------------------------------------------------------
@dataclass
class CSPParams:
    """O'Brien continuous-synthesis run parameters (set once from the CLI).

    Composition over duplication: every MD / ribosome knob lives in the reused
    :class:`ElongationParams` (``elong``); only the **kinetic** fields are added
    here. ``elong.trna_tether`` is forced off by the runner (CSP needs the A<->P
    switchable position restraint), and ``elong.n_steps`` is ignored (each stage's
    step count comes from the kinetics).
    """
    elong: ElongationParams = field(default_factory=ElongationParams)

    # --- O'Brien timing ---
    scale_factor: float = 4331293.0     # in-vivo seconds -> in-silico ns compressor
    time_stage_1: float = 0.00034       # mean peptidyl-transfer dwell (s)
    time_stage_2: float = 0.004201      # mean translocation dwell (s)
    uniform_ta: bool = False            # ignore the mRNA; use uniform_mfpt for every codon
    uniform_mfpt: float = 0.05          # uniform mean codon time (s) when uniform_ta
    ribosome_traffic: bool = False      # apply the external traffic correction if available
    initiation_rate: float = 0.083333   # translation initiation rate (1/s), traffic only
    random_seed: Optional[int] = None   # seed for the FPT sampler (reproducibility)

    # --- test clamps (production: leave both at their defaults / None) ---
    max_steps_per_stage: Optional[int] = None  # cap each stage (tutorial: small)
    min_steps_per_stage: int = 1               # floor each stage

    # --- post-synthesis phases (steps; 0 = skip) ---
    ejection_steps: int = 0             # release the restraint; let the chain leave
    dissociation_steps: int = 0         # free run; protein drifts off the ribosome


# --------------------------------------------------------------------------
# The continuous-synthesis loop
# --------------------------------------------------------------------------
def run_continuous_synthesis(full_pdb: str, ribosome_pdb: str, *,
                             L0: int, L_max: Optional[int] = None,
                             out_root: str = "synth_out",
                             mrna: Optional[str] = None,
                             trans_times: Optional[str] = None,
                             domain_def: Optional[str] = None,
                             stride_output_file: Optional[str] = None,
                             params: Optional[CSPParams] = None) -> None:
    """Run the full O'Brien continuous synthesis ``L = L0 .. L_max``.

    Parameters
    ----------
    full_pdb : str
        Full native PDB of the target protein (the nascent chain at full length).
    ribosome_pdb : str
        Truncated CG ribosome PDB -- source of the P-/A-anchors (and, with
        ``params.elong.rigid_ribosome``, the rigid scenery).
    L0, L_max : int
        Nascent-length schedule; ``L_max`` blank -> the full residue count.
    out_root : str
        Root output directory; each residue writes ``L_<L>/stage_<1,2,3>/``.
    mrna, trans_times : str, optional
        mRNA sequence file and per-codon time table -- the codon-resolved kinetics.
        Not needed when ``params.uniform_ta`` (every codon gets ``uniform_mfpt``).
    domain_def, stride_output_file : str, optional
        Passed to the one-time contact precompute (n_scale / STRIDE).
    params : CSPParams, optional
        Run parameters (defaults to the dataclass defaults).
    """
    if params is None:
        params = CSPParams()
    ep = params.elong
    # CSP needs the A<->P switchable position restraint: force the tRNA tether off.
    ep.trna_tether = False

    out_path = Path(out_root)
    out_path.mkdir(parents=True, exist_ok=True)

    # --- anchors (fixed points from the truncated ribosome) -----------------
    p_anchor = read_anchor(ribosome_pdb, "PtR", resid=76, bead="R")
    a_anchor = read_anchor(ribosome_pdb, "AtR", resid=76, bead="R")
    print(f"P-anchor (PtR 76 R): {p_anchor} nm")
    print(f"A-anchor (AtR 76 R): {a_anchor} nm")

    # Hold/seed targets, offset into the tunnel (+x) from the anchor bead so the
    # C-terminus does not sit on top of a ribosome bead (clash). Same auto rule as
    # the elongation driver (tether bond length / 0.4 nm in v2; 0 in v1).
    offset = ep.ptc_offset_nm
    if offset is None:
        offset = TRNA_TETHER_BOND_NM if ep.rigid_ribosome else 0.0
    p_target = p_anchor + offset * TUNNEL_AXIS
    a_target = a_anchor + offset * TUNNEL_AXIS

    # --- rigid ribosome (loaded once; identical every length) ---------------
    ribo = None
    if ep.rigid_ribosome:
        ribo = load_ribosome(ribosome_pdb, model="topo")
        print(f"Rigid ribosome: {ribo.n} beads from {ribosome_pdb} "
              f"(mass-0 scenery; ribosome<->nascent forces on).")

    # --- build-once-subset contacts on the full native structure ------------
    R_full, eps_full = precompute_contacts(full_pdb, domain_def, stride_output_file)
    N_full = R_full.shape[0]
    if L_max is None:
        L_max = N_full
    if not (1 <= L0 <= L_max <= N_full):
        raise ValueError(f"require 1 <= L0 <= L_max <= N_full; got L0={L0}, "
                         f"L_max={L_max}, N_full={N_full}.")

    # --- kinetics: intrinsic / real per-codon mFPT lists --------------------
    # Need intrinsic[L_max] valid -> at least L_max + 1 codons.
    intrinsic, real, codons = kinetics.build_mfpt_lists(
        L_max + 1, uniform_ta=params.uniform_ta, uniform_mfpt=params.uniform_mfpt,
        mrna_path=mrna, trans_times_path=trans_times,
        ribosome_traffic=params.ribosome_traffic,
        initiation_rate=params.initiation_rate)
    rng = random.Random(params.random_seed)

    print()
    print("=" * 66)
    print("[ O'Brien continuous synthesis -- kinetic schedule ]")
    print("=" * 66)
    print(f"  timing mode: {'uniform' if params.uniform_ta else 'per-codon (mRNA)'}; "
          f"scale_factor={params.scale_factor:g}; dt={ep.dt_ps} ps")
    print(f"  stage means (s): peptidyl-transfer={params.time_stage_1:g}, "
          f"translocation={params.time_stage_2:g}, tRNA-binding=remainder")
    if params.max_steps_per_stage is not None:
        print(f"  TEST CLAMP: <= {params.max_steps_per_stage} steps/stage "
              f"(~{3 * params.max_steps_per_stage} steps/residue). Remove for production.")
    print(f"Synthesizing {full_pdb}: L = {L0} .. {L_max} (N_full = {N_full}).")

    # --- per-residue dwell-time log (mirrors v6 output/<id>.out) -------------
    # One machine-parsable row per residue: the sampled 3-stage dwell times (s),
    # their in-silico mapping (ns + integration steps) and the codon. Written
    # incrementally so a partial run still leaves a usable table (D4).
    dwell_log = out_path / "dwell_times.dat"
    dt_ns = ep.dt_ps * 1e-3
    dwell_fh = open(dwell_log, "w")
    dwell_fh.write(
        "# O'Brien continuous-synthesis per-residue dwell times (topo.csp)\n"
        f"#   scale_factor={params.scale_factor:g}  dt={ep.dt_ps} ps  "
        f"time_stage_1={params.time_stage_1:g} s  time_stage_2={params.time_stage_2:g} s\n"
        f"#   timing={'uniform' if params.uniform_ta else 'per-codon'}  "
        f"ribosome_traffic={'on' if params.ribosome_traffic else 'off'}  "
        f"random_seed={params.random_seed}\n"
        "#   t1/t2/t3 = sampled peptidyl-transfer / translocation / tRNA-binding "
        "dwell (s); steps = clamped integration steps actually run\n"
        "# L  codon  t_invivo_total_s  t1_s  t2_s  t3_s  "
        "ns1  ns2  ns3  steps1  steps2  steps3\n")
    dwell_fh.flush()

    # --- main loop: one residue = three sub-stages --------------------------
    prev_final: Optional[np.ndarray] = None
    for L in range(L0, L_max + 1):
        (s1, s2, s3), (t1, t2, t3) = kinetics.stage_steps(
            L, intrinsic, real, time_stage_1=params.time_stage_1,
            time_stage_2=params.time_stage_2, scale_factor=params.scale_factor,
            dt_ps=ep.dt_ps, rng=rng,
            max_steps_per_stage=params.max_steps_per_stage,
            min_steps_per_stage=params.min_steps_per_stage)

        codon = codons[L - 1] if codons is not None else "uniform"
        print()
        print("#" * 66)
        print(f"# Residue L = {L}  codon {codon}  "
              f"(total in-vivo dwell {intrinsic[L]:.4g} s)")
        print(f"#   stage 1 peptidyl transfer : {t1:.4g} s -> {s1} steps")
        print(f"#   stage 2 translocation     : {t2:.4g} s -> {s2} steps")
        print(f"#   stage 3 tRNA binding/wait : {t3:.4g} s -> {s3} steps")
        print("#" * 66)

        # Record the per-residue dwell row (in-vivo seconds + in-silico ns/steps).
        dwell_fh.write(
            f"{L:4d}  {codon:>5s}  {intrinsic[L]:.6e}  "
            f"{t1:.6e}  {t2:.6e}  {t3:.6e}  "
            f"{t1 * 1e9 / params.scale_factor:.6e}  "
            f"{t2 * 1e9 / params.scale_factor:.6e}  "
            f"{t3 * 1e9 / params.scale_factor:.6e}  "
            f"{s1:8d}  {s2:8d}  {s3:8d}\n")
        dwell_fh.flush()

        ldir = f"L_{L:03d}"
        # Stage 1: deliver the new residue at the A-site (cold-start lays the
        # initial segment from the P-anchor instead) and restrain there.
        stage1_anchor = p_target if prev_final is None else a_target
        f1 = run_length(L, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
                        p_anchor=stage1_anchor, a_anchor=a_anchor,
                        prev_final=prev_final, out_root=out_path, params=ep,
                        ribo=ribo, restrain=True,
                        out_subdir=f"{ldir}/stage_1", n_steps_override=s1,
                        label=f"L={L} stage 1 (peptidyl transfer) {s1} steps")

        # Stage 2: continue from stage 1, still held at the A-site.
        f2 = run_length(L, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
                        p_anchor=stage1_anchor, a_anchor=a_anchor,
                        prev_final=None, seed_override=f1, out_root=out_path,
                        params=ep, ribo=ribo, restrain=True,
                        out_subdir=f"{ldir}/stage_2", n_steps_override=s2,
                        label=f"L={L} stage 2 (translocation) {s2} steps")

        # Stage 3: translocate A->P (restrain the C-terminus to the P-anchor).
        f3 = run_length(L, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
                        p_anchor=p_target, a_anchor=a_anchor,
                        prev_final=None, seed_override=f2, out_root=out_path,
                        params=ep, ribo=ribo, restrain=True,
                        out_subdir=f"{ldir}/stage_3", n_steps_override=s3,
                        label=f"L={L} stage 3 (tRNA binding) {s3} steps")
        prev_final = f3

    dwell_fh.close()
    print()
    print(f"Done. Synthesized {L0} -> {L_max}. Per-residue/-stage outputs under {out_path}/")
    print(f"Per-residue dwell-time table: {dwell_log}")

    # --- post-synthesis: ejection then dissociation (both free runs) --------
    if params.ejection_steps > 0:
        print()
        print(f"=== Ejection (L = {L_max}, {params.ejection_steps} steps, "
              f"restraint OFF) -> {out_path / 'ejection'}/ ===")
        prev_final = run_length(
            L_max, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
            p_anchor=p_target, a_anchor=a_anchor, prev_final=None,
            seed_override=prev_final, out_root=out_path, params=ep, ribo=ribo,
            restrain=False, out_subdir="ejection",
            n_steps_override=params.ejection_steps,
            label=f"Ejection (L = {L_max})")

    if params.dissociation_steps > 0:
        print()
        print(f"=== Dissociation (L = {L_max}, {params.dissociation_steps} steps, "
              f"restraint OFF) -> {out_path / 'dissociation'}/ ===")
        run_length(
            L_max, full_pdb=full_pdb, R_full=R_full, eps_full=eps_full,
            p_anchor=p_target, a_anchor=a_anchor, prev_final=None,
            seed_override=prev_final, out_root=out_path, params=ep, ribo=ribo,
            restrain=False, out_subdir="dissociation",
            n_steps_override=params.dissociation_steps,
            label=f"Dissociation (L = {L_max})")


# --------------------------------------------------------------------------
# INI control file
# --------------------------------------------------------------------------
@dataclass
class CSPConfig:
    """Parsed contents of a CSP control file (``csp.ini``)."""
    pdb_file: str
    ribosome: str
    L0: int
    L_max: Optional[int] = None
    outdir: str = "synth_out"
    mrna: Optional[str] = None
    trans_times: Optional[str] = None
    domain_def: Optional[str] = None
    stride_output_file: Optional[str] = None
    params: CSPParams = field(default_factory=CSPParams)
    config_file: Optional[str] = None


def read_csp_config(config_file: str, verbose: bool = True) -> CSPConfig:
    """Parse a CSP control file (INI ``[OPTIONS]``) into a :class:`CSPConfig`.

    Superset of ``elongate.ini``: the same structure / MD / ribosome keys (see
    :func:`topo.translation.elongate.read_elongate_config`) **plus** the O'Brien
    kinetic keys. Required: ``pdb_file``, ``ribosome``, ``L0``. Non-uniform timing
    additionally requires ``mrna`` and ``trans_times``.

    Kinetic keys
    ------------
    - ``mrna`` -- mRNA sequence file (raw nucleotides, wrapped ok); one codon per
      residue + 1 stop.
    - ``trans_times`` -- per-codon mean-time table (``CODON  seconds``).
    - ``scale_factor`` -- in-vivo seconds -> in-silico ns compressor.
    - ``time_stage_1`` / ``time_stage_2`` -- mean peptidyl-transfer / translocation
      dwell (s); stage 3 = codon total minus these.
    - ``uniform_ta`` -- yes: ignore the mRNA, use ``uniform_mfpt`` for every codon.
    - ``uniform_mfpt`` -- the uniform mean codon time (s) when ``uniform_ta``.
    - ``ribosome_traffic`` -- yes: apply the external ``ribosome_traffic`` binary's
      correction if it is on PATH (else falls back to no traffic).
    - ``initiation_rate`` -- translation initiation rate (1/s), traffic only.
    - ``random_seed`` -- seed for the FPT sampler (reproducible schedules).
    - ``max_steps_per_stage`` -- cap each stage's step count (the tutorial uses a
      small value for a ~2000-steps/residue test; blank = uncapped production).
    - ``min_steps_per_stage`` -- floor each stage's step count (default 1).
    - ``ejection_steps`` / ``dissociation_steps`` -- post-synthesis free runs
      (0 = skip).

    MD / ribosome keys (same meaning as ``elongate.ini``): ``n_steps`` is **not**
    used (step counts come from the kinetics); ``dt``, ``ref_t``, ``tau_t``,
    ``nstout``, ``device``, ``ppn``, ``constraints``, ``restraint_k``, ``buffer``,
    ``minimize``, ``rigid_ribosome``, ``tunnel_wall``, ``tunnel_wall_x0``,
    ``tunnel_wall_k``, ``ptc_offset``, ``nascent_only_output``.
    """
    def log(msg: str) -> None:
        if verbose:
            print(msg)

    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read CSP config file: {config_file!r}")
    if "OPTIONS" not in cp:
        raise ValueError(f"{config_file}: missing required [OPTIONS] section.")
    o = cp["OPTIONS"]

    def opt(key: str) -> Optional[str]:
        v = o.get(key, None)
        if v is None:
            return None
        v = v.strip()
        return v if v != "" else None

    def req(key: str) -> str:
        v = opt(key)
        if v is None:
            raise ValueError(f"{config_file}: required option '{key}' is missing or blank.")
        return v

    def as_int(s: str) -> int:
        return int(str(s).replace("_", ""))

    log(f"Reading CSP parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    ribosome = req("ribosome")
    L0 = int(req("L0"))
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None
    outdir = opt("outdir") or "synth_out"
    mrna = opt("mrna")
    trans_times = opt("trans_times")
    domain_def = opt("domain_def")
    stride_output_file = opt("stride_output_file")

    # --- MD / ribosome knobs (reused ElongationParams) ----------------------
    ep = ElongationParams()
    if opt("dt") is not None:
        ep.dt_ps = float(opt("dt"))
    if opt("ref_t") is not None:
        ep.ref_t = float(opt("ref_t"))
    if opt("tau_t") is not None:
        ep.tau_t = float(opt("tau_t"))
    if opt("nstout") is not None:
        ep.nstout = int(opt("nstout"))
    if opt("device") is not None:
        ep.device = opt("device")
    if opt("ppn") is not None:
        ep.ppn = int(opt("ppn"))
    cons = opt("constraints")
    ep.constraints = None if (cons is None or cons.lower() == "none") else cons
    if opt("restraint_k") is not None:
        ep.restraint_k = float(opt("restraint_k"))
    if opt("buffer") is not None:
        ep.buffer_nm = float(opt("buffer"))
    if opt("minimize") is not None:
        ep.minimize = bool(strtobool(opt("minimize")))
    if opt("rigid_ribosome") is not None:
        ep.rigid_ribosome = bool(strtobool(opt("rigid_ribosome")))
    if opt("tunnel_wall") is not None:
        ep.tunnel_wall = bool(strtobool(opt("tunnel_wall")))
    if opt("tunnel_wall_x0") is not None:
        ep.tunnel_wall_x0_nm = float(opt("tunnel_wall_x0"))
    if opt("tunnel_wall_k") is not None:
        ep.tunnel_wall_k = float(opt("tunnel_wall_k"))
    if opt("ptc_offset") is not None:
        ep.ptc_offset_nm = float(opt("ptc_offset"))
    if opt("nascent_only_output") is not None:
        ep.nascent_only_output = bool(strtobool(opt("nascent_only_output")))

    # --- O'Brien kinetic knobs ---------------------------------------------
    p = CSPParams(elong=ep)
    if opt("scale_factor") is not None:
        p.scale_factor = float(opt("scale_factor"))
    if opt("time_stage_1") is not None:
        p.time_stage_1 = float(opt("time_stage_1"))
    if opt("time_stage_2") is not None:
        p.time_stage_2 = float(opt("time_stage_2"))
    if opt("uniform_ta") is not None:
        p.uniform_ta = bool(strtobool(opt("uniform_ta")))
    if opt("uniform_mfpt") is not None:
        p.uniform_mfpt = float(opt("uniform_mfpt"))
    if opt("ribosome_traffic") is not None:
        p.ribosome_traffic = bool(strtobool(opt("ribosome_traffic")))
    if opt("initiation_rate") is not None:
        p.initiation_rate = float(opt("initiation_rate"))
    if opt("random_seed") is not None:
        p.random_seed = as_int(opt("random_seed"))
    if opt("max_steps_per_stage") is not None:
        p.max_steps_per_stage = as_int(opt("max_steps_per_stage"))
    if opt("min_steps_per_stage") is not None:
        p.min_steps_per_stage = as_int(opt("min_steps_per_stage"))
    if opt("ejection_steps") is not None:
        p.ejection_steps = as_int(opt("ejection_steps"))
    if opt("dissociation_steps") is not None:
        p.dissociation_steps = as_int(opt("dissociation_steps"))

    # Validation mirroring v6.
    if not p.uniform_ta and (mrna is None or trans_times is None):
        raise ValueError(f"{config_file}: non-uniform timing needs both 'mrna' and "
                         f"'trans_times' (or set uniform_ta = yes).")

    log(f"  inputs: pdb_file={pdb_file}, ribosome={ribosome}")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}")
    log(f"  timing: {'uniform (mfpt=%g s)' % p.uniform_mfpt if p.uniform_ta else f'per-codon (mrna={mrna}, trans_times={trans_times})'}")
    log(f"          scale_factor={p.scale_factor:g}, time_stage_1={p.time_stage_1:g} s, "
        f"time_stage_2={p.time_stage_2:g} s")
    log(f"          ribosome_traffic={'on' if p.ribosome_traffic else 'off'}"
        + (f" (initiation_rate={p.initiation_rate:g}/s)" if p.ribosome_traffic else ""))
    if p.max_steps_per_stage is not None:
        log(f"          TEST CLAMP max_steps_per_stage={p.max_steps_per_stage} "
            f"(~{3 * p.max_steps_per_stage} steps/residue)")
    log(f"  ribosome forces: {'on (rigid v2)' if ep.rigid_ribosome else 'off (v1)'}"
        + (f"; tunnel wall: {('x>=%.2f nm' % ep.tunnel_wall_x0_nm) if ep.tunnel_wall else 'off'}"
           if ep.rigid_ribosome else ""))
    log(f"  integrator: dt={ep.dt_ps} ps, ref_t={ep.ref_t} K, tau_t={ep.tau_t} /ps, nstout={ep.nstout}")
    if p.ejection_steps or p.dissociation_steps:
        log(f"  post-synthesis: ejection={p.ejection_steps} steps, "
            f"dissociation={p.dissociation_steps} steps")
    log(f"  hardware/output: device={ep.device}, ppn={ep.ppn}, outdir={outdir}")

    return CSPConfig(pdb_file=pdb_file, ribosome=ribosome, L0=L0, L_max=L_max,
                     outdir=outdir, mrna=mrna, trans_times=trans_times,
                     domain_def=domain_def, stride_output_file=stride_output_file,
                     params=p, config_file=config_file)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def csp(argv: Optional[List[str]] = None) -> None:
    """Console entry point: ``topo-csp -f csp.ini``.

    The O'Brien continuous synthesis protocol (per-codon, 3-stage elongation),
    driven by an INI control file (see :func:`read_csp_config`). ``-o`` / ``--device``
    override the output directory / compute device for sweeps.
    """
    parser = argparse.ArgumentParser(
        prog="topo-csp",
        description="O'Brien Continuous Synthesis Protocol (per-codon, 3-stage "
                    "elongation) in topo style. Times every residue from its codon, "
                    "splits it into peptidyl-transfer / translocation / tRNA-binding "
                    "sub-stages, and grows the nascent chain N->C on the ribosome. "
                    "Controlled by an INI file: topo-csp -f csp.ini",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-input", "-f", dest="config", type=str,
                        help="CSP control file (INI, [OPTIONS] section).")
    parser.add_argument("-o", "--outdir", default=None,
                        help="override the output directory from the config file.")
    parser.add_argument("--device", default=None, choices=["CPU", "GPU"],
                        help="override the compute device from the config file.")

    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(argv)
    if not args.config:
        parser.error("a CSP control file is required: -f csp.ini")

    print(f"OpenMM version: {mm.__version__}")

    cfg = read_csp_config(args.config)
    if args.outdir:
        cfg.outdir = args.outdir
    if args.device:
        cfg.params.elong.device = args.device

    run_continuous_synthesis(
        cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max, out_root=cfg.outdir,
        mrna=cfg.mrna, trans_times=cfg.trans_times, domain_def=cfg.domain_def,
        stride_output_file=cfg.stride_output_file, params=cfg.params)


if __name__ == "__main__":
    csp()
