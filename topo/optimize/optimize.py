"""
Strength (n_scale) optimizer for TOPO coarse-grained models.

Automatically chooses the per-domain and per-interface ``strength`` in
``domain.yaml`` — the smallest value on a discrete per-class ladder that keeps
every domain and interface folded across ``ntraj`` independent trajectories.

This is the canonical optimizer for the package. Use it as a CLI::

    topo-optimize -f optimize.ini -o opt_out
    python -m topo.optimize -f optimize.ini -o opt_out

or call :func:`run_optimizer` from your own script.

The optimizer owns the search logic and drives the package tools as sub-steps,
one round at a time::

    round loop:
      1. write round_N/domain.yaml with the current strengths
      2. topo.mdrun                 (one multi-copy run -> ntraj chains)
      3. topo.split_chains          (split into per-copy DCDs, in-process)
      4. score Q per domain / per interface  (topo.analysis.native_contacts)
      5. decide: stable units freeze; unstable units climb the ladder
      until all units are stable, or unstable units reach the median fallback.

``optimize.ini`` is a MINIMAL config: a single ``[OPTIONS]`` section. The
optimizer takes the keys it needs (ntraj, q_threshold, frame_fraction,
max_rounds, min_contacts — see :data:`CONTROL_TYPES`); every other key is a
simulation parameter (pdb_file, domain_def, md_steps, sampling, ref_t, ...)
passed through to each round's md.ini. Anything unset uses the optimizer's
implicit protocol defaults (:data:`IMPLICIT_DEFAULTS` / :data:`OPT_DEFAULTS`).
Each round the driver expands ``optimize.ini`` into a full ``round_N/md.ini`` and
writes ``round_N/domain.yaml`` with the current strengths.

Per round it also writes one Q time series per trajectory next to its DCD:
``round_N/traj/Q_<k>.csv`` (paired with ``traj_<k>.dcd``; columns
``frame, Q_<domain>..., Q_<d1>-<d2>...``) for inspection of the Q values behind
each stability decision.

Limitations: the optimization is not resumable — each invocation starts fresh
(decide how to checkpoint ``level[]`` and completed rounds to add this).
"""

from __future__ import annotations

import argparse
import configparser
import copy
import itertools
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import MDAnalysis as mda

# The Q machinery lives in the package. Import the submodule explicitly
# (topo.analysis also re-exports a `native_contacts` function of the same name).
import topo.analysis.native_contacts as ncmod
from topo.utils.multichain import split_chains


# --------------------------------------------------------------------------- #
# The strength ladder (Table 1). Each class: 5 levels + a median fallback.
# --------------------------------------------------------------------------- #
LADDER = {
    "alpha":      ([1.1954, 1.4704, 1.7453, 2.0322, 2.5044], 1.7453),
    "beta":       ([1.4732, 1.8120, 2.1508, 2.5044, 2.5044], 2.1508),
    "alpha_beta": ([1.1556, 1.4213, 1.6871, 1.9644, 2.5044], 1.6871),
    "interface":  ([1.2747, 1.5679, 1.8611, 2.1670, 2.5044], 1.8611),
}
FALLBACK_INDEX = 5   # level index 0..4 -> the 5 levels; 5 -> median fallback

CLASS_ALIASES = {
    "alpha": "alpha", "a": "alpha",
    "beta": "beta", "b": "beta",
    "alpha-beta": "alpha_beta", "alpha/beta": "alpha_beta",
    "alphabeta": "alpha_beta", "ab": "alpha_beta", "c": "alpha_beta",
}


def normalize_class(raw):
    """Map a user-written structural class to a LADDER key."""
    key = str(raw).strip().lower()
    if key not in CLASS_ALIASES:
        raise ValueError(
            f"Unknown structural class {raw!r}; expected one of "
            f"alpha/beta/alpha-beta (or a/b/c)."
        )
    return CLASS_ALIASES[key]


def strength_for(class_key, level):
    """Strength for a given class and ladder level index (>=5 -> fallback)."""
    levels, fallback = LADDER[class_key]
    return fallback if level >= len(levels) else levels[level]


# --------------------------------------------------------------------------- #
# optimize.ini -> per-round md.ini.
#
# The user writes a MINIMAL optimize.ini with only the essentials. The optimizer
# fills in the model-appropriate "implicit" simulation parameters below and, each
# round, expands everything into a full md.ini for topo.mdrun.
#
# IMPLICIT_DEFAULTS are the optimizer's own protocol defaults — NOT the bare
# SimulationConfig dataclass defaults. In particular dt = 0.015 ps (the CG model
# was parameterized at a 15 fs timestep), whereas the dataclass default is 0.01.
# Anything set in optimize.ini [OPTIONS] overrides these.
#
# We use configparser with the SAME settings as the package reader
# (topo.read_simulation_config, config.py) so parsing semantics match exactly.
# read_simulation_config returns a typed SimulationConfig (units + defaults) with
# no serializer back to a file, so we round-trip the text with configparser here.
# --------------------------------------------------------------------------- #
IMPLICIT_DEFAULTS = {
    "dt": "0.015",        # ps; CG model parameterized at a 15 fs timestep
    "model": "topo",
    "tcoupl": "yes",
    "tau_t": "0.01",
    "pcoupl": "no",
    "pbc": "no",
    "minimize": "no",     # native structure is already the model's energy minimum
    "restart": "no",      # every optimization round starts fresh
    "device": "GPU",      # production-scale task; override in optimize.ini for CPU
    "ppn": "4",
    "ref_t": "310",       # K; stability protocol temperature
}

# max_rounds = 6 covers the normal case exactly: 5 ladder levels + the median
# fallback, with every unstable unit climbing one level per round. A run can in
# theory need more only if raising one unit's strength destabilizes a previously
# stable unit (rare); such a protein simply fails to stabilize in 6 rounds and is
# flagged with a WARNING in the report for the user to inspect/exclude.
# min_contacts: a unit (domain or interface) with fewer than this many native
# contacts is considered too weakly structured to fold; it is NOT optimized but
# pinned at the first ladder level and frozen. Default 0 disables the check.
OPT_DEFAULTS = {"ntraj": 10, "q_threshold": 0.6688,
                "frame_fraction": 0.98, "max_rounds": 6, "min_contacts": 0}

# Keys the optimizer consumes itself (with the type to cast them to). Everything
# else in the config section is a simulation parameter passed through to the
# per-round md.ini, so optimize.ini needs only ONE [OPTIONS] section: the
# optimizer takes these keys, topo.mdrun gets the rest.
CONTROL_TYPES = {"ntraj": int, "q_threshold": float, "frame_fraction": float,
                 "max_rounds": int, "min_contacts": int}


def read_optimize_config(path):
    """Read a minimal optimize.ini.

    Returns
    -------
    pdb : str            absolute path to the reference PDB
    domain : str         absolute path to the initial domain.yaml
    sim_options : dict   simulation parameters for the per-round md.ini
                         (file overrides IMPLICIT_DEFAULTS; controls removed)
    controls : dict      optimizer controls (ntraj, thresholds, max_rounds, ...)
    """
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    cp.optionxform = str
    cp.read(path)
    if "OPTIONS" not in cp:
        raise SystemExit(f"{path}: missing required [OPTIONS] section")

    options = dict(cp["OPTIONS"])
    for required in ("pdb_file", "domain_def"):
        if required not in options:
            raise SystemExit(f"{path}: [OPTIONS] must set '{required}'")

    # One flat section: split off the optimizer controls; the rest are simulation
    # parameters for the per-round md.ini. Controls are popped so they never leak
    # into md.ini.
    controls = {}
    for key, cast in CONTROL_TYPES.items():
        controls[key] = cast(options.pop(key)) if key in options else OPT_DEFAULTS[key]

    base = Path(path).resolve().parent          # resolve paths relative to the ini
    pdb = str((base / options["pdb_file"]).resolve())
    domain = str((base / options["domain_def"]).resolve())
    # A user-supplied STRIDE file must be absolute too: each round's md.ini lives in
    # a different round_N/ dir, so a relative path would not resolve at run time.
    if "stride_output_file" in options:
        options["stride_output_file"] = str(
            (base / options["stride_output_file"]).resolve())

    sim_options = {**IMPLICIT_DEFAULTS, **options}   # file overrides defaults
    return pdb, domain, sim_options, controls


def write_round_ini(path, base_options, overrides):
    """Write a per-round md.ini = base [OPTIONS] with `overrides` applied."""
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    cp.optionxform = str
    cp["OPTIONS"] = {**base_options, **{k: str(v) for k, v in overrides.items()}}
    with open(path, "w") as fh:
        cp.write(fh)


# --------------------------------------------------------------------------- #
# domain.yaml per round
# --------------------------------------------------------------------------- #
def write_domain_yaml(path, raw_cfg, interfaces, strength_dom, strength_int):
    """Write a domain.yaml that keeps residues/class from the input and sets the
    current per-domain and per-interface strengths."""
    cfg = copy.deepcopy(raw_cfg)
    for name, vals in cfg["intra_domains"].items():
        vals["strength"] = round(float(strength_dom[name]), 4)
    cfg["inter_domains"] = {
        f"{a}-{b}": round(float(strength_int[(a, b)]), 4) for (a, b) in interfaces
    }
    with open(path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)


# --------------------------------------------------------------------------- #
# Scorer: native-contact lists are built ONCE and reused every round.
# --------------------------------------------------------------------------- #
class Scorer:
    """Per-domain / per-interface native contacts, and per-trajectory folding."""

    def __init__(self, pdb, domain_yaml, cutoff=4.5, local_separation=3,
                 tolerance=1.2):
        self.tolerance = tolerance

        u_ref = ncmod.load_universe(pdb)
        ca_positions, resindex_to_pos = ncmod.reference_residue_geometry(u_ref)
        self.n_res = ca_positions.shape[0]

        heavy = u_ref.select_atoms("protein and not name H*")
        heavy_positions = heavy.positions.copy()
        heavy_res = np.fromiter(
            (resindex_to_pos[ri] for ri in heavy.resindices), dtype=int,
            count=heavy.n_atoms,
        )

        domains, n_residues = ncmod.load_domains(domain_yaml, self.n_res)
        if n_residues != self.n_res:
            raise ValueError(
                f"domain.yaml n_residues={n_residues} but reference has "
                f"{self.n_res} protein residues."
            )
        self.domain_names = list(domains)
        self.interfaces = list(itertools.combinations(self.domain_names, 2))

        # Build the native-contact list (pairs + reference Cα distances) per unit.
        self.units = {}   # key -> (pairs, native_dist)
        for name, idx in domains.items():
            self.units[("domain", name)] = ncmod.build_native_contacts(
                set(idx.tolist()), None, heavy_positions, heavy_res,
                ca_positions, cutoff, local_separation)
        for a, b in self.interfaces:
            self.units[("interface", (a, b))] = ncmod.build_native_contacts(
                set(domains[a].tolist()), set(domains[b].tolist()),
                heavy_positions, heavy_res, ca_positions, cutoff, local_separation)

    def unit_keys(self):
        return list(self.units)

    def n_contacts(self, key):
        """Number of native contacts detected for a unit (domain or interface)."""
        pairs, _ = self.units[key]
        return int(pairs.shape[0])

    @staticmethod
    def label(key):
        kind, name = key
        return f"Q_{name}" if kind == "domain" else f"Q_{name[0]}-{name[1]}"

    def q_per_frame(self, psf, dcd):
        """Return {unit_key: np.array of Q per frame} for one trajectory."""
        # Silence MDAnalysis' DCDReader v3.0 DeprecationWarning (harmless here, and
        # this runs once per chain per round, so it would otherwise spam the log).
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r".*DCDReader currently makes independent timesteps.*",
                category=DeprecationWarning,
            )
            u = mda.Universe(psf, dcd)
        if u.atoms.n_atoms != self.n_res:
            raise ValueError(
                f"CG model has {u.atoms.n_atoms} beads but reference has "
                f"{self.n_res} residues.")
        n_frames = u.trajectory.n_frames
        series = {k: np.empty(n_frames) for k in self.units}
        for i, _ in enumerate(u.trajectory):
            pos = u.atoms.positions
            for k, (pairs, dnat) in self.units.items():
                series[k][i] = ncmod.fraction_native_contacts(
                    pos, pairs, dnat, self.tolerance)
        return series

    @staticmethod
    def folded_fraction(series, q_threshold):
        """Fraction of frames folded per unit, from a q_per_frame() result.
        NaN (no native contacts) is treated as folded (not applicable)."""
        return {k: float(np.mean(np.isnan(v) | (v > q_threshold)))
                for k, v in series.items()}


# --------------------------------------------------------------------------- #
# Round execution
# --------------------------------------------------------------------------- #
def run_subprocess(cmd, log_path, label, cwd=None):
    """Run `cmd`, streaming combined stdout/stderr to `log_path`; raise on failure."""
    print(f"    $ {' '.join(str(c) for c in cmd)}")
    with open(log_path, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{label} failed (exit {proc.returncode}); see {log_path}")


def run_md(round_dir, md_ini, python_exe):
    """Run one multi-copy MD (topo.mdrun) as a subprocess.

    Runs with cwd=round_dir so that, when no stride_output_file is configured, the
    STRIDE file the model build caches ("{pdb_stem}_stride.dat") lands predictably
    inside round_dir, where the caller can pick it up to reuse in later rounds.
    All paths in the md.ini are absolute, so the working directory is otherwise
    irrelevant. A fresh subprocess per round also isolates each OpenMM/GPU
    context."""
    run_subprocess([python_exe, "-m", "topo.mdrun", "-f", str(md_ini)],
                   round_dir / "mdrun.out", "topo.mdrun", cwd=round_dir)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def run_optimizer(config, outdir="opt_out", device=None, md_steps=None,
                  python_exe=None):
    """Run the strength optimization end to end.

    Parameters
    ----------
    config : str or Path
        Path to the minimal ``optimize.ini``.
    outdir : str or Path, optional
        Optimization root directory (created if missing). Default ``opt_out``.
    device : str, optional
        Override the simulation device (``CPU``/``GPU``) for every round.
    md_steps : int, optional
        Override ``md_steps`` for every round (useful for quick test runs).
    python_exe : str, optional
        Python interpreter used to launch ``topo.mdrun`` subprocesses.
        Defaults to the current interpreter (``sys.executable``).

    Returns
    -------
    final_yaml : Path
        Path to the written ``domain_optimized.yaml``.
    converged : bool
        True if every (non-frozen) unit reached stability within ``max_rounds``.
    """
    python_exe = python_exe or sys.executable
    pdb, domain_path, sim_options, controls = read_optimize_config(config)
    ntraj = controls["ntraj"]
    q_threshold = controls["q_threshold"]
    frame_fraction = controls["frame_fraction"]
    max_rounds = controls["max_rounds"]
    min_contacts = controls["min_contacts"]
    raw_cfg = yaml.safe_load(Path(domain_path).read_text())

    out_root = Path(outdir).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    opt_log = out_root / "optimization.log"
    # Line-buffered + explicit flush so the report is readable live during a long
    # run (each line hits disk immediately; tail -f works).
    log_fh = open(opt_log, "w", buffering=1)

    def log(msg):
        print(msg, flush=True)
        log_fh.write(msg + "\n")
        log_fh.flush()

    try:
        return _optimize_loop(
            log, pdb, domain_path, raw_cfg, sim_options, out_root,
            ntraj, q_threshold, frame_fraction, max_rounds, min_contacts,
            device, md_steps, python_exe)
    finally:
        log_fh.close()


def _optimize_loop(log, pdb, domain_path, raw_cfg, sim_options, out_root,
                   ntraj, q_threshold, frame_fraction, max_rounds, min_contacts,
                   device, md_steps, python_exe):
    """Core round loop. Separated from resource setup so the log file is always
    closed (see :func:`run_optimizer`)."""
    # --- build units (domains + interfaces) and their classes ---
    scorer = Scorer(pdb, domain_path)
    unit_class = {}
    for name in scorer.domain_names:
        unit_class[("domain", name)] = normalize_class(
            raw_cfg["intra_domains"][name]["class"])
    for a, b in scorer.interfaces:
        unit_class[("interface", (a, b))] = "interface"

    level = {k: 0 for k in scorer.unit_keys()}   # per-unit ladder index

    # Units with too few native contacts are too weakly structured to fold:
    # pin them at the first ladder level and exclude them from optimization.
    frozen = {k for k in scorer.unit_keys()
              if scorer.n_contacts(k) < min_contacts}

    log(f"# Strength optimization for {Path(pdb).name}")
    log(f"# domains: {scorer.domain_names}  interfaces: {scorer.interfaces}")
    log("# native contacts detected (heavy-atom <= 4.5 A, |i-j| > 3):")
    for key in scorer.unit_keys():
        kind, name = key
        what = f"domain {name}" if kind == "domain" else f"interface {name[0]}-{name[1]}"
        note = f"  (frozen: < min_contacts={min_contacts})" if key in frozen else ""
        log(f"#   {what:<16}: {scorer.n_contacts(key)}{note}")
    if frozen:
        log(f"# {len(frozen)} unit(s) below min_contacts={min_contacts}: "
            f"pinned at level 1, not optimized.")
    log(f"# ntraj={ntraj}  T={sim_options.get('ref_t')} K  "
        f"Q>{q_threshold} for >={frame_fraction:.0%} of frames")

    converged = False
    unstable_labels = []
    rnd = 0
    for rnd in range(1, max_rounds + 1):
        round_dir = out_root / f"round_{rnd}"
        round_dir.mkdir(exist_ok=True)

        # 1. current strengths from each unit's class ladder
        strength_dom = {name: strength_for(unit_class[("domain", name)],
                                            level[("domain", name)])
                        for name in scorer.domain_names}
        strength_int = {(a, b): strength_for("interface", level[("interface", (a, b))])
                        for (a, b) in scorer.interfaces}

        domain_yaml = round_dir / "domain.yaml"
        write_domain_yaml(domain_yaml, raw_cfg, scorer.interfaces,
                          strength_dom, strength_int)

        def level_tag(key):
            if key in frozen:
                return "frozen"
            lv = level[key]
            return "fallback" if lv >= FALLBACK_INDEX else f"level {lv + 1}"

        log(f"\n## Round {rnd}:")
        for name in scorer.domain_names:
            tag = level_tag(("domain", name))
            log(f"   domain {name:<10} {tag:<9} strength = {strength_dom[name]:.4f}")
        for a, b in scorer.interfaces:
            tag = level_tag(("interface", (a, b)))
            log(f"   iface  {a}-{b:<8} {tag:<9} strength = {strength_int[(a, b)]:.4f}")

        # 2-3. produce ntraj independent trajectories (one multi-copy MD run)
        md_ini = round_dir / "md.ini"
        overrides = {
            "pdb_file": pdb,
            "domain_def": str(domain_yaml),
            "n_copies": ntraj,
            "output_dir": str(round_dir / "traj"),
            "outname": "traj",
        }
        if device:
            overrides["device"] = device
        if md_steps is not None:
            overrides["md_steps"] = md_steps
        write_round_ini(md_ini, sim_options, overrides)
        eff_steps = overrides.get("md_steps", sim_options.get("md_steps"))
        log(f"   running {ntraj}-copy MD ({eff_steps} steps/traj) — "
            f"live progress in round_{rnd}/traj/traj.log ...")
        run_md(round_dir, md_ini, python_exe)

        # STRIDE depends only on the (fixed) structure, so it is identical every
        # round. If the user did not supply one, reuse the file the round-1 model
        # build cached so STRIDE runs once instead of once per round.
        if "stride_output_file" not in sim_options:
            cached_stride = round_dir / f"{Path(pdb).stem}_stride.dat"
            if cached_stride.exists():
                sim_options["stride_output_file"] = str(cached_stride)
                log(f"   caching STRIDE output for later rounds: {cached_stride}")

        # Split the combined multi-copy DCD into per-copy DCDs (in-process,
        # memory-bounded streaming — handles large trajectories).
        psf = round_dir / "traj" / "traj.psf"
        chain_dcds = [round_dir / "traj" / f"traj_{k}.dcd" for k in range(ntraj)]
        split_chains(round_dir / "traj" / "traj.dcd", chain_dcds)
        log(f"   MD done; scoring {ntraj} trajectories ...")

        # 4. score: per-frame Q for every trajectory. Save one CSV per trajectory
        #    next to its DCD (round_N/traj/Q_<k>.csv, paired with traj_<k>.dcd) for
        #    easy inspection, and reduce to the folded fraction per unit for the
        #    decision.
        traj_dir = round_dir / "traj"
        per_traj = []
        for k, dcd in enumerate(chain_dcds):
            series = scorer.q_per_frame(str(psf), str(dcd))
            per_traj.append(Scorer.folded_fraction(series, q_threshold))
            table = pd.DataFrame({scorer.label(key): arr
                                  for key, arr in series.items()})
            table.insert(0, "frame", np.arange(len(table)))
            table.to_csv(traj_dir / f"Q_{k}.csv", index=False)

        # 5. decide stability and climb unstable units
        all_stable = True
        progressed = False
        unstable_labels = []
        for key in scorer.unit_keys():
            lbl = scorer.label(key)
            if key in frozen:
                log(f"   {lbl:<14} frozen (< min_contacts={min_contacts}) "
                    f"-> not optimized")
                continue
            fracs = [pt[key] for pt in per_traj]
            stable = all(f >= frame_fraction for f in fracs)
            log(f"   {lbl:<14} folded frac/traj: "
                f"[{', '.join(f'{f:.2f}' for f in fracs)}] -> "
                f"{'STABLE' if stable else 'unstable'}")
            if not stable:
                all_stable = False
                unstable_labels.append(lbl)
                if level[key] < FALLBACK_INDEX:
                    level[key] += 1
                    progressed = True   # at least one unit can still climb

        if all_stable:
            converged = True
            log(f"\n# Converged at round {rnd}: all units stable.")
            break
        if not progressed:
            log(f"\n# Stopped at round {rnd}: remaining unstable units are at the "
                f"median fallback (accepted as-is).")
            break
    else:
        log(f"\n# Stopped: reached max_rounds={max_rounds}.")

    # --- final calibrated domain.yaml ---
    final_dom = {name: strength_for(unit_class[("domain", name)],
                                    level[("domain", name)])
                 for name in scorer.domain_names}
    final_int = {(a, b): strength_for("interface", level[("interface", (a, b))])
                 for (a, b) in scorer.interfaces}
    final_yaml = out_root / "domain_optimized.yaml"
    write_domain_yaml(final_yaml, raw_cfg, scorer.interfaces, final_dom, final_int)

    log("\n## Final strengths:")
    for name in scorer.domain_names:
        log(f"   domain {name:<10} strength = {final_dom[name]:.4f}")
    for a, b in scorer.interfaces:
        log(f"   iface  {a}-{b:<8} strength = {final_int[(a, b)]:.4f}")
    if frozen:
        log(f"# frozen (< min_contacts={min_contacts}, pinned at level 1, "
            f"not optimized): "
            f"{', '.join(scorer.label(k) for k in scorer.unit_keys() if k in frozen)}")

    if converged:
        log(f"\n# CONVERGED in {rnd} rounds. Calibrated model: {final_yaml}")
    else:
        bar = "#" + "=" * 68
        log(f"\n{bar}")
        log(f"# WARNING: {Path(pdb).name} did NOT stabilize after {rnd} rounds "
            f"(max_rounds={max_rounds}).")
        log(f"# Unstable units at the final model: {', '.join(unstable_labels)}")
        log("# Their strengths are left at the highest level / median fallback.")
        log("# Consider excluding this protein or inspecting it manually.")
        log(bar)
        log(f"# Calibrated model written anyway: {final_yaml}")

    return final_yaml, converged


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="topo-optimize",
        description="Per-domain/interface strength (n_scale) optimizer for "
                    "TOPO coarse-grained models.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-f", "--config", required=True,
                   help="optimize.ini (minimal config; see Tutorial 5).")
    p.add_argument("-o", "--outdir", default="opt_out", help="Optimization root dir.")
    p.add_argument("--device", default=None,
                   help="Override device (CPU/GPU) for every round.")
    p.add_argument("--md-steps", type=int, default=None,
                   help="Override md_steps for every round (e.g. quick test runs).")
    p.add_argument("--python", default=sys.executable,
                   help="Python interpreter used to launch topo.mdrun.")
    # A bare `topo-optimize` (no arguments) prints help, like `-h`.
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        p.print_help()
        sys.exit(0)
    return p.parse_args(argv)


def optimize(argv=None):
    """Console entry point (``topo-optimize`` / ``python -m topo.optimize``)."""
    args = parse_args(argv)
    run_optimizer(args.config, outdir=args.outdir, device=args.device,
                  md_steps=args.md_steps, python_exe=args.python)


if __name__ == "__main__":
    optimize()
