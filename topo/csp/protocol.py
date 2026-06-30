"""Continuous Synthesis Protocol (O'Brien), ported to topo (``topo.csp``).

This is the **full O'Brien continuous-synthesis protocol** -- the per-codon, 3-stage
elongation cycle of ``continuous_synthesis_v6.py`` -- expressed in topo style. It is
the kinetic upgrade of :mod:`topo.csp.core`: that driver grows the chain
one residue per step at a *fixed* ``n_steps`` (collapsing the elongation cycle into a
single MD segment); **this** driver times every residue from its codon and splits it
into three sub-stages, exactly as O'Brien do.

What is reused vs. new:

- **Reused wholesale** from :mod:`topo.csp.core`: the per-length MD
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
restraint** path (``trna_tether`` is forced off); the supplied ribosome is always
rigid scenery and the tunnel wall + excluded volume + electrostatics are on.

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

from topo.csp.core import (ElongationParams, TUNNEL_AXIS,
                                       precompute_contacts, read_anchor,
                                       run_length, TRNA_TETHER_BOND_NM)
from topo.csp.ribosome import load_ribosome
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

    Attributes
    ----------
    elong : ElongationParams
        All MD / ribosome knobs (timestep, temperature, rigid ribosome, tunnel wall,
        restraint constant, output, ...), reused unchanged from the elongation driver.
    scale_factor : float
        In-vivo-seconds -> in-silico-nanoseconds compression factor (default
        ``4331293``). Larger -> fewer MD steps per residue -> faster run.
    time_stage_1 : float
        Mean peptidyl-transfer (peptide-bond) dwell time, seconds (default ``0.00034``).
    time_stage_2 : float
        Mean translocation dwell time, seconds (default ``0.004201``).
    uniform_ta : bool
        If True, ignore the mRNA and use ``uniform_mfpt`` for every codon.
    uniform_mfpt : float
        The uniform per-codon mean time (seconds) used when ``uniform_ta`` (default
        ``0.05``).
    random_seed : int or None
        Seed for the first-passage-time sampler, for reproducible schedules.
    max_steps_per_stage : int or None
        Optional upper clamp on each stage's MD step count (``None`` = uncapped /
        production); clamps MD steps only, not the dwell times in seconds.
    min_steps_per_stage : int
        Lower clamp on each stage's MD step count (default 1).
    ejection_steps : int
        Post-synthesis ejection-phase length in steps (release the C-terminus
        restraint; ``0`` = skip).
    dissociation_steps : int
        Post-synthesis dissociation-phase length in steps (free run; ``0`` = skip).
    """
    elong: ElongationParams = field(default_factory=ElongationParams)

    # --- O'Brien timing ---
    scale_factor: float = 4331293.0     # in-vivo seconds -> in-silico ns compressor
    time_stage_1: float = 0.00034       # mean peptidyl-transfer dwell (s)
    time_stage_2: float = 0.004201      # mean translocation dwell (s)
    uniform_ta: bool = False            # ignore the mRNA; use uniform_mfpt for every codon
    uniform_mfpt: float = 0.05          # uniform mean codon time (s) when uniform_ta
    # ribosome_traffic / initiation_rate: HIDDEN/deferred (off by default; not exposed
    # in the docs or example csp.ini -- see topo/csp/TODO.md). Still parsed if present.
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
                             L0: int = 1, L_max: Optional[int] = None,
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
        Truncated CG ribosome PDB -- source of the P-/A-anchors and the rigid
        (mass-0) scenery (always loaded; providing it is the signal to use it).
    L0 : int, optional
        First nascent length to synthesize (default ``1`` -- start from a single residue).
    L_max : int or None, optional
        Final nascent length (default ``None`` -> the full residue count of the protein).
    out_root : str
        Root output directory; each residue writes ``L_<L>/stage_<1,2,3>/``.
    mrna : str, optional
        mRNA sequence file (one codon per residue) -- the codon-resolved kinetics.
        Required for per-codon timing; not needed when ``params.uniform_ta``.
    trans_times : str, optional
        Per-codon mean-time table. ``None`` -> the bundled E. coli 310 K table
        (organism-universal; see :func:`topo.csp.kinetics.default_trans_times_path`).
    domain_def : str
        Domain-definition file (``domain.yaml``) defining the protein's **native-contact
        strengths** -- per-domain and per-interface Gō well-depth scaling (the
        structure-based analog of O'Brien's ``nscal``). Required via the INI.
    stride_output_file : str, optional
        Precomputed STRIDE file for the contact build (skips re-running STRIDE).
    params : CSPParams, optional
        Run parameters (defaults to the dataclass defaults).

    Returns
    -------
    None
        Side-effecting: writes per-residue/-stage trajectories under
        ``out_root/L_<L>/stage_<1,2,3>/``, a per-residue ``dwell_times.dat`` table,
        and (if requested) ``ejection/`` and ``dissociation/`` phases.

    Raises
    ------
    ValueError
        If the length schedule is invalid (``1 <= L0 <= L_max <= N_full`` fails), or
        if non-uniform kinetics are requested without ``mrna`` / ``trans_times``
        (propagated from :func:`topo.csp.kinetics.build_mfpt_lists`).
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
    # C-terminus does not sit on top of a ribosome bead (clash). The supplied ribosome
    # is always treated as rigid scenery, so the offset defaults to the tether bond
    # length (override via `ptc_offset`).
    offset = ep.ptc_offset_nm
    if offset is None:
        offset = TRNA_TETHER_BOND_NM
    p_target = p_anchor + offset * TUNNEL_AXIS
    a_target = a_anchor + offset * TUNNEL_AXIS

    # Tunnel wall plane (auto-derived from the structure -- never a stale user knob):
    # place the one-sided wall at the LOWER (deeper-in-tunnel, smaller-x) of the two
    # C-terminus hold planes, i.e. min(P,A anchor x) + ptc_offset (the P-site, where
    # the nascent C-terminus is tethered). This blocks the chain from slipping below
    # the synthesis point into the void left where the 50S was truncated, while the
    # held C-terminus still sits at (just on) the plane. Recomputed for whatever
    # ribosome PDB is supplied, so switching structures can never leave it wrong.
    if ep.tunnel_wall:
        ep.tunnel_wall_x0_nm = float(min(p_anchor[0], a_anchor[0]) + offset)
        print(f"Tunnel wall plane: x >= {ep.tunnel_wall_x0_nm:.4f} nm "
              f"(auto: lower P/A-site C-terminus hold plane = min(P.x,A.x)+ptc_offset).")

    # --- rigid ribosome (always: the supplied PDB is rigid scenery) ---------
    # Loaded once; identical at every length. Providing a ribosome PDB *is* the
    # signal to treat it as rigid -- there is no separate on/off flag.
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
        f"{'ribosome_traffic=on  ' if params.ribosome_traffic else ''}"
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
    """Parsed contents of a CSP control file (``csp.ini``).

    Bundles the run inputs (structures, the ``L0..L_max`` schedule, output directory,
    one-time-precompute options) with the kinetic + MD :class:`CSPParams`. Produced
    by :func:`read_csp_config` and consumed by :func:`csp` / passed straight to
    :func:`run_continuous_synthesis`.

    Attributes
    ----------
    pdb_file : str
        Full native PDB of the target protein (the nascent chain at full length).
    ribosome : str
        Truncated CG ribosome PDB (P-/A-anchors; always rigid scenery).
    L0 : int
        Starting nascent-chain length (default ``1``).
    L_max : int or None
        Final length (``None`` -> the full residue count).
    outdir : str
        Root output directory (default ``"synth_out"``).
    mrna, trans_times : str or None
        mRNA sequence file and per-codon time table (required unless
        ``params.uniform_ta``).
    domain_def : str
        Domain-definition file (contact-strength scaling); required via the INI.
    stride_output_file : str or None
        Precomputed STRIDE file for the one-time contact build (optional).
    params : CSPParams
        The kinetic + MD/ribosome run parameters.
    config_file : str or None
        Path of the INI file this config was parsed from (provenance).
    """
    pdb_file: str
    ribosome: str
    L0: int = 1
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

    The structure / MD / ribosome keys configure the shared
    :class:`topo.csp.core.ElongationParams` machinery; the O'Brien **kinetic keys**
    are added on top. Required: ``pdb_file``, ``ribosome``, ``domain_def`` (the protein's
    domain/contact-strength definition). ``L0`` (default ``1``) and ``L_max`` (default =
    full residue count) are optional. Per-codon timing additionally requires ``mrna``
    (``trans_times`` is optional).

    Kinetic keys
    ------------
    - ``mrna`` -- mRNA sequence file (raw nucleotides, wrapped ok); one codon per
      residue + 1 stop. Required for per-codon timing (it is protein-specific).
    - ``trans_times`` -- per-codon mean-time table (``CODON  seconds``). **Optional**:
      defaults to the bundled E. coli 310 K table (the table is organism-universal,
      so it is not protein-specific). See :func:`topo.csp.kinetics.default_trans_times_path`.
    - ``scale_factor`` -- in-vivo seconds -> in-silico ns compressor.
    - ``time_stage_1`` / ``time_stage_2`` -- mean peptidyl-transfer / translocation
      dwell (s); stage 3 = codon total minus these.
    - ``uniform_ta`` -- yes: ignore the mRNA, use ``uniform_mfpt`` for every codon.
    - ``uniform_mfpt`` -- the uniform mean codon time (s) when ``uniform_ta``.
    - ``random_seed`` -- seed for the FPT sampler (reproducible schedules).
    - ``max_steps_per_stage`` -- cap each stage's step count (the tutorial uses a
      small value for a ~2000-steps/residue test; blank = uncapped production).
    - ``min_steps_per_stage`` -- floor each stage's step count (default 1).
    - ``ejection_steps`` / ``dissociation_steps`` -- post-synthesis free runs
      (0 = skip).

    MD / ribosome keys: ``dt``, ``ref_t``, ``tau_t``, ``nstout``, ``device``, ``ppn``,
    ``constraints``, ``restraint_k``, ``buffer``, ``minimize``, ``tunnel_wall``,
    ``ptc_offset``. (The supplied ribosome is always rigid scenery -- there is no
    ``rigid_ribosome`` key -- and output is always nascent-only, so there is no
    ``nascent_only_output`` key either. The wall plane ``tunnel_wall_x0`` and stiffness
    ``tunnel_wall_k`` are **not** keys -- the plane is auto-derived from the ribosome
    structure at run time and the stiffness is a fixed model constant.)

    Parameters
    ----------
    config_file : str
        Path to the INI control file (must contain an ``[OPTIONS]`` section).
    verbose : bool, optional
        If True (default), echo the parsed configuration to stdout.

    Returns
    -------
    CSPConfig
        The parsed configuration (inputs, schedule and a populated
        :class:`CSPParams`), ready to pass to :func:`run_continuous_synthesis`.

    Raises
    ------
    FileNotFoundError
        If ``config_file`` cannot be read.
    ValueError
        If the ``[OPTIONS]`` section is missing, a required key (``pdb_file``,
        ``ribosome``, ``domain_def``) is absent/blank, or per-codon timing is requested
        without ``mrna``.
    """
    def log(msg: str) -> None:
        """Print a status message, but only when ``verbose`` is set.

        Parameters
        ----------
        msg : str
            The message to print to stdout.
        """
        if verbose:
            print(msg)

    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    if not cp.read(config_file):
        raise FileNotFoundError(f"could not read CSP config file: {config_file!r}")
    if "OPTIONS" not in cp:
        raise ValueError(f"{config_file}: missing required [OPTIONS] section.")
    o = cp["OPTIONS"]

    def opt(key: str) -> Optional[str]:
        """Return an optional ``[OPTIONS]`` value, or ``None`` if absent/blank.

        Parameters
        ----------
        key : str
            The option name to look up in the ``[OPTIONS]`` section.

        Returns
        -------
        str or None
            The stripped value, or ``None`` if the key is missing or its value
            is the empty string.
        """
        v = o.get(key, None)
        if v is None:
            return None
        v = v.strip()
        return v if v != "" else None

    def req(key: str) -> str:
        """Return a required ``[OPTIONS]`` value, raising if it is missing.

        Parameters
        ----------
        key : str
            The option name to look up in the ``[OPTIONS]`` section.

        Returns
        -------
        str
            The non-blank option value.

        Raises
        ------
        ValueError
            If the option is missing or blank.
        """
        v = opt(key)
        if v is None:
            raise ValueError(f"{config_file}: required option '{key}' is missing or blank.")
        return v

    def as_int(s: str) -> int:
        """Parse an integer string, tolerating ``_`` digit-group separators.

        Parameters
        ----------
        s : str
            The integer literal, optionally containing underscores (e.g.
            ``"200_000"``).

        Returns
        -------
        int
            The parsed integer value.
        """
        return int(str(s).replace("_", ""))

    log(f"Reading CSP parameters from {config_file} ...")

    pdb_file = req("pdb_file")
    ribosome = req("ribosome")
    L0 = int(opt("L0")) if opt("L0") is not None else 1   # optional; default = 1
    L_max = opt("L_max")
    L_max = int(L_max) if L_max is not None else None      # optional; default = full length
    outdir = opt("outdir") or "synth_out"
    mrna = opt("mrna")
    trans_times = opt("trans_times")
    domain_def = req("domain_def")   # required: defines the protein's contact strengths
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
    # rigid_ribosome is intentionally NOT read: a ribosome PDB is required, and supplying
    # it *is* the signal to treat it as rigid scenery -- there is no separate flag.
    if opt("tunnel_wall") is not None:
        ep.tunnel_wall = bool(strtobool(opt("tunnel_wall")))
    # Neither tunnel_wall_x0 nor tunnel_wall_k is read from the INI: the wall plane is
    # auto-derived from the ribosome structure (see run_continuous_synthesis), and the
    # stiffness is a fixed model constant (O'Brien's 20 kcal/mol/A^2 = 8368 kJ/mol/nm^2,
    # ElongationParams.tunnel_wall_k). Only the `tunnel_wall` on/off toggle is exposed.
    if opt("ptc_offset") is not None:
        ep.ptc_offset_nm = float(opt("ptc_offset"))
    # nascent_only_output is not a key: with a ribosome present (always, for CSP) the
    # output is always nascent-only -- it is the default CSP behavior, not a toggle.

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

    # Validation: per-codon timing needs the (protein-specific) mRNA. The codon-time
    # table is organism-universal, so `trans_times` is optional -- it defaults to the
    # bundled E. coli 310 K table (topo.csp.kinetics.default_trans_times_path).
    if not p.uniform_ta and mrna is None:
        raise ValueError(f"{config_file}: per-codon timing needs an 'mrna' file "
                         f"(or set uniform_ta = yes). 'trans_times' is optional "
                         f"(defaults to the bundled E. coli 310 K table).")

    log(f"  inputs: pdb_file={pdb_file}, ribosome={ribosome}")
    log(f"  schedule: L0={L0}, L_max={L_max if L_max is not None else 'full'}")
    log(f"  timing: {'uniform (mfpt=%g s)' % p.uniform_mfpt if p.uniform_ta else f'per-codon (mrna={mrna}, trans_times={trans_times or 'bundled E. coli 310 K'})'}")
    log(f"          scale_factor={p.scale_factor:g}, time_stage_1={p.time_stage_1:g} s, "
        f"time_stage_2={p.time_stage_2:g} s")
    if p.ribosome_traffic:   # hidden/deferred feature; only mention it when enabled
        log(f"          ribosome_traffic=on (initiation_rate={p.initiation_rate:g}/s)")
    if p.max_steps_per_stage is not None:
        log(f"          TEST CLAMP max_steps_per_stage={p.max_steps_per_stage} "
            f"(~{3 * p.max_steps_per_stage} steps/residue)")
    log(f"  ribosome: rigid scenery (always, from the supplied PDB)"
        f"; tunnel wall: {'on (plane auto-derived from structure)' if ep.tunnel_wall else 'off'}")
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

    Parameters
    ----------
    argv : list of str or None, optional
        Command-line arguments to parse (default ``None`` -> use ``sys.argv``). A
        bare invocation (no args) prints help and exits.

    Returns
    -------
    None
        Runs the synthesis for its side effects (see
        :func:`run_continuous_synthesis`).
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
