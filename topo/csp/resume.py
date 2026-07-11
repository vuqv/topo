"""Resume support for the O'Brien continuous-synthesis driver (``topo.csp``).

A production ``topo-csp`` run is hours to days of wall time and today survives no
interruption. This module makes the run resumable with three small, human-readable
on-disk artifacts and no heavyweight simulation checkpoint (see ``review/F_RESUME.md``
for the full design):

1. **The schedule** (``dwell_times.dat``). The per-residue 3-stage step counts are
   drawn from the seeded generator **once, before the main loop**, and persisted. The
   file header additionally carries the deterministic-but-expensive PTC restraint
   geometry (``a_target`` / ``p_target`` / tunnel-wall plane). On resume the driver
   *reads* this table instead of re-drawing the RNG or re-running the SLSQP PTC solve,
   so the kinetic schedule and the restraint geometry are pinned identical across the
   interruption. This is the immutable **plan**.

2. **The progress log** (``progress.log``). An append-only ``L_XXX RUNNING`` /
   ``L_XXX DONE`` record of how far the run got. The ``DONE`` line is the commit point;
   a unit left ``RUNNING`` at the crash is dropped and redone. This is the mutable
   **status**.

3. **The seed conformation.** The only per-residue state that crosses a residue
   boundary besides the schedule is the previous residue's final coordinates, already
   written to disk every stage; resume reloads them (:func:`load_final_pdb`).

The resume unit is the **residue** (its final structure ``traj_final.pdb`` is the one
recovery target). Per-stage resume is intentionally out of scope. There is no serialized
RNG state anywhere -- the schedule file *is* the materialized RNG output.

**Layout note.** :func:`residue_final_path` / :func:`phase_final_path` centralize the
one piece of output-layout knowledge resume needs (where a completed unit's final
structure lands). Under the consolidated layout (``F_RESUME.md`` §3.5) that is
``L_<L>/traj_final.pdb`` (the single per-residue final written by stage 3).
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit

# Progress-log schema tag (bumped only on an incompatible format change).
PROGRESS_SCHEMA = 1
PROGRESS_FILENAME = "progress.log"


# --------------------------------------------------------------------------
# Output-layout knowledge (the only place resume encodes where finals land)
# --------------------------------------------------------------------------
def residue_final_path(out_root: Path, L: int) -> Path:
    """Path of the final (stage-3) structure for residue length ``L``.

    Consolidated layout (§3.5): ``<out_root>/L_<L>/traj_final.pdb`` -- the single
    per-residue final (written by stage 3). It is the seed for residue ``L+1`` and the
    resume-reload target.
    """
    return Path(out_root) / f"L_{L:03d}" / "traj_final.pdb"


def cylinder_final_path(out_root: Path, L: int) -> Path:
    """Path of the final structure for residue ``L`` in the **cylinder** runner.

    The cylinder runner (:mod:`topo.csp.cylinder`) uses one MD segment per residue and a
    flat layout: ``<out_root>/L_<L>/traj_final.pdb`` (no ``stage_*`` sub-directories).
    """
    return Path(out_root) / f"L_{L:03d}" / "traj_final.pdb"


def phase_final_path(out_root: Path, name: str) -> Path:
    """Path of a post-synthesis phase's final structure (``ejection`` / ``dissociation``).

    ``<out_root>/<name>/traj_final.pdb`` -- the phase runs as a single ``run_length``
    call writing its final directly under ``<name>/``. Same for both runners.
    """
    return Path(out_root) / name / "traj_final.pdb"


# --------------------------------------------------------------------------
# The schedule table (dwell_times.dat): immutable plan
# --------------------------------------------------------------------------
@dataclass
class SchedRow:
    """One residue's persisted schedule row.

    Attributes
    ----------
    L : int
        1-indexed nascent length.
    codon : str
        The codon that adds residue ``L`` (``"uniform"`` for uniform timing).
    t_total : float
        The codon's total in-vivo dwell time ``intrinsic[L]`` (seconds).
    times : tuple of float
        The three sampled sub-stage dwell times ``(t1, t2, t3)`` (seconds).
    steps : tuple of int
        The three clamped integration step counts ``(s1, s2, s3)``.
    """
    L: int
    codon: str
    t_total: float
    times: Tuple[float, float, float]
    steps: Tuple[int, int, int]


def _fmt_float_exact(x: float) -> str:
    """``repr`` of a Python float -- round-trips to the exact same value."""
    return repr(float(x))


def write_schedule(path: Path, rows: List[SchedRow], params,
                   a_target: np.ndarray, p_target: np.ndarray,
                   wall_x: Optional[float]) -> None:
    """Write the immutable schedule table + PTC-geometry header to ``dwell_times.dat``.

    Called **once** at a fresh start (never on resume -- the file is re-read, not
    rewritten). The data rows are byte-compatible with the legacy in-loop
    ``dwell_times.dat`` writer; the ``#PTC`` header block is new and carries the
    restraint geometry at full float precision so resume can skip the SLSQP solve.

    Parameters
    ----------
    path : pathlib.Path
        Output path (``<out_root>/dwell_times.dat``).
    rows : list of SchedRow
        The full per-residue schedule, ascending in ``L``.
    params : RunParams
        Run parameters (for the metadata header: ``scale_factor``, ``dt_ps``, ...).
    a_target, p_target : numpy.ndarray
        The A-site / P-site PTC restraint target points (nm), 3-vectors.
    wall_x : float or None
        The tunnel-wall plane x (nm), or ``None`` if the wall is off.
    """
    dt_ps = params.dt_ps
    timing = "uniform" if params.uniform_codon_time is not None else "per-codon"
    with open(path, "w") as fh:
        fh.write(
            "# O'Brien continuous-synthesis per-residue dwell times (topo.csp)\n"
            f"#   scale_factor={params.scale_factor:g}  dt={dt_ps} ps  "
            f"time_stage_1={params.time_stage_1:g} s  time_stage_2={params.time_stage_2:g} s\n"
            f"#   timing={timing}  "
            f"{'ribosome_traffic=on  ' if params.ribosome_traffic else ''}"
            f"random_seed={params.random_seed}\n"
            "#   t1/t2/t3 = sampled peptidyl-transfer / translocation / tRNA-binding "
            "dwell (s); steps = clamped integration steps actually run\n")
        # PTC restraint geometry (machine-readable, full precision) -- re-read on
        # resume so the SLSQP PTC solve is skipped and the geometry is pinned
        # identical to the residues already on disk.
        fh.write(f"#PTC schema {PROGRESS_SCHEMA}\n")
        fh.write("#PTC a_target " + " ".join(_fmt_float_exact(v) for v in a_target) + "\n")
        fh.write("#PTC p_target " + " ".join(_fmt_float_exact(v) for v in p_target) + "\n")
        fh.write("#PTC wall_x " + (_fmt_float_exact(wall_x) if wall_x is not None else "none") + "\n")
        fh.write(
            "# L  codon  t_invivo_total_s  t1_s  t2_s  t3_s  "
            "ns1  ns2  ns3  steps1  steps2  steps3\n")
        for r in rows:
            t1, t2, t3 = r.times
            s1, s2, s3 = r.steps
            fh.write(
                f"{r.L:4d}  {r.codon:>5s}  {r.t_total:.6e}  "
                f"{t1:.6e}  {t2:.6e}  {t3:.6e}  "
                f"{t1 * 1e9 / params.scale_factor:.6e}  "
                f"{t2 * 1e9 / params.scale_factor:.6e}  "
                f"{t3 * 1e9 / params.scale_factor:.6e}  "
                f"{s1:8d}  {s2:8d}  {s3:8d}\n")


def read_schedule(path: Path) -> Tuple[List[SchedRow], np.ndarray, np.ndarray,
                                       Optional[float]]:
    """Read a ``dwell_times.dat`` table back into rows + PTC geometry.

    Inverse of :func:`write_schedule`. The ``#PTC`` header lines supply the restraint
    geometry; the data rows supply the per-residue schedule (step counts exact, dwell
    times to the file's precision).

    Returns
    -------
    tuple
        ``(rows, a_target, p_target, wall_x)`` -- rows ascending in ``L``, the two PTC
        target points (nm, 3-vectors), and the tunnel-wall plane x (nm) or ``None``.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If a required ``#PTC`` header line is absent or a data row is malformed.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"schedule table not found: {path}")
    a_target = p_target = None
    wall_x: Optional[float] = None
    rows: List[SchedRow] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#PTC"):
                tok = stripped.split()
                if len(tok) >= 5 and tok[1] == "a_target":
                    a_target = np.array([float(x) for x in tok[2:5]])
                elif len(tok) >= 5 and tok[1] == "p_target":
                    p_target = np.array([float(x) for x in tok[2:5]])
                elif len(tok) >= 3 and tok[1] == "wall_x":
                    wall_x = None if tok[2] == "none" else float(tok[2])
                continue
            if stripped.startswith("#"):
                continue
            tok = stripped.split()
            if len(tok) < 12:
                raise ValueError(f"{path}: malformed schedule row: {line!r}")
            rows.append(SchedRow(
                L=int(tok[0]), codon=tok[1], t_total=float(tok[2]),
                times=(float(tok[3]), float(tok[4]), float(tok[5])),
                steps=(int(tok[9]), int(tok[10]), int(tok[11]))))
    if a_target is None or p_target is None:
        raise ValueError(f"{path}: missing '#PTC a_target'/'#PTC p_target' header "
                         f"lines (not a resume-capable schedule file).")
    rows.sort(key=lambda r: r.L)
    return rows, a_target, p_target, wall_x


def schedule_covers(rows: List[SchedRow], L0: int, L_max: int) -> None:
    """Assert the persisted schedule covers exactly ``L0..L_max`` (contiguous).

    A launch config whose ``L0``/``L_max`` disagree with the persisted plan is caught
    here -- extending a run (a larger ``L_max``) needs draws past the materialized
    schedule and is a fresh run, not a resume.

    Raises
    ------
    SystemExit
        With an actionable message if the coverage does not match.
    """
    have = {r.L for r in rows}
    want = set(range(L0, L_max + 1))
    if have != want:
        missing = sorted(want - have)
        extra = sorted(have - want)
        raise SystemExit(
            f"[resume] persisted schedule covers L={sorted(have)[0]}..{sorted(have)[-1]} "
            f"but this run asks for L={L0}..{L_max}"
            + (f"; missing {missing}" if missing else "")
            + (f"; unexpected {extra}" if extra else "")
            + ". Extending a run is a fresh run (the schedule is fixed at first launch).")


# --------------------------------------------------------------------------
# The cylinder schedule table (single MD segment per residue)
# --------------------------------------------------------------------------
# The cylinder runner (topo.csp.cylinder) has no ribosome/PTC geometry to persist -- its
# tunnel geometry is cheap and deterministic from the params -- and one MD segment per
# residue instead of three. So its schedule table needs no #PTC header and carries a
# single (dwell, steps) pair per row. The RNG draw is still the thing that must be
# materialized once and re-read, exactly as for CSP.
@dataclass
class CylSchedRow:
    """One residue's persisted cylinder schedule row (single MD segment).

    Attributes
    ----------
    L : int
        1-indexed nascent length.
    codon : str
        The codon that adds residue ``L`` (``"uniform"`` for uniform timing).
    dwell_s : float
        The sampled codon dwell time (seconds).
    steps : int
        The clamped integration step count for this residue's single MD segment.
    """
    L: int
    codon: str
    dwell_s: float
    steps: int


def write_cylinder_schedule(path: Path, rows: List[CylSchedRow], params) -> None:
    """Write the immutable cylinder schedule to ``dwell_times.dat`` (written once).

    Byte-compatible with the legacy in-loop cylinder ``dwell_times.dat`` writer; drawn up
    front and re-read on resume so the per-residue RNG draw is not repeated.
    """
    timing = "uniform" if params.uniform_codon_time is not None else "per-codon"
    with open(path, "w") as fh:
        fh.write(
            "# cylinder continuous-synthesis per-residue dwell times (topo.csp.cylinder)\n"
            f"#   scale_factor={params.scale_factor:g}  dt={params.dt_ps} ps  "
            f"timing={timing}  random_seed={params.random_seed}\n"
            "#   t_dwell = sampled codon dwell (s); ns = in-silico ns; steps = integration "
            "steps actually run (single MD segment)\n"
            "# L  codon  t_dwell_s  ns  steps\n")
        for r in rows:
            ns = r.dwell_s * 1e9 / params.scale_factor
            fh.write(f"{r.L:4d}  {r.codon:>5s}  {r.dwell_s:.6e}  {ns:.6e}  {r.steps:8d}\n")


def read_cylinder_schedule(path: Path) -> List[CylSchedRow]:
    """Read a cylinder ``dwell_times.dat`` table back into rows (ascending in ``L``).

    Inverse of :func:`write_cylinder_schedule`: step counts exact, dwell times to the
    file's precision.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If a data row is malformed.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"schedule table not found: {path}")
    rows: List[CylSchedRow] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tok = stripped.split()
            if len(tok) < 5:
                raise ValueError(f"{path}: malformed cylinder schedule row: {line!r}")
            rows.append(CylSchedRow(L=int(tok[0]), codon=tok[1],
                                    dwell_s=float(tok[2]), steps=int(tok[4])))
    rows.sort(key=lambda r: r.L)
    return rows


# --------------------------------------------------------------------------
# The progress log (progress.log): mutable status
# --------------------------------------------------------------------------
def progress_path(out_root: Path) -> Path:
    """Path of the progress log under ``out_root``."""
    return Path(out_root) / PROGRESS_FILENAME


def progress_exists(out_root: Path) -> bool:
    """True iff a ``progress.log`` is present under ``out_root``."""
    return progress_path(out_root).is_file()


def write_progress_header(out_root: Path) -> None:
    """Create a fresh ``progress.log`` with its schema header line (truncates any prior)."""
    with open(progress_path(out_root), "w") as fh:
        fh.write(f"# csp progress log -- schema {PROGRESS_SCHEMA}\n")


def append_progress(out_root: Path, unit: str, status: str) -> None:
    """Append one ``<unit> <status>`` line and flush (the crash-safe commit point).

    A single short append is effectively atomic on POSIX; appending ``DONE`` only after
    all of a unit's work is on disk makes every crash recoverable.

    Parameters
    ----------
    out_root : pathlib.Path
        Output root holding ``progress.log``.
    unit : str
        Unit name -- ``L_<L>`` for a residue, or ``ejection`` / ``dissociation``.
    status : str
        ``"RUNNING"`` or ``"DONE"``.
    """
    with open(progress_path(out_root), "a") as fh:
        fh.write(f"{unit} {status}\n")
        fh.flush()


@dataclass
class Progress:
    """Parsed ``progress.log`` -- the last status seen per unit."""
    last_status: Dict[str, str]

    def is_done(self, unit: str) -> bool:
        """True iff ``unit``'s last recorded status is ``DONE``."""
        return self.last_status.get(unit) == "DONE"

    @property
    def last_done_residue(self) -> int:
        """Highest residue length ``L`` whose unit ``L_<L>`` is ``DONE`` (0 if none)."""
        best = 0
        for unit, st in self.last_status.items():
            if st == "DONE" and unit.startswith("L_"):
                try:
                    best = max(best, int(unit[2:]))
                except ValueError:
                    continue
        return best

    def running_units(self) -> List[str]:
        """Units whose last recorded status is ``RUNNING`` (the in-flight-at-crash units)."""
        return [u for u, s in self.last_status.items() if s == "RUNNING"]


def read_progress(out_root: Path) -> Progress:
    """Parse ``progress.log`` into a :class:`Progress` (last status wins per unit)."""
    last: Dict[str, str] = {}
    with open(progress_path(out_root)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()
            if len(tok) != 2:
                continue
            last[tok[0]] = tok[1]
    return Progress(last_status=last)


# --------------------------------------------------------------------------
# Resume actions: verify, drop, reload
# --------------------------------------------------------------------------
def verify_completed_units(out_root: Path, prog: Progress, L0: int,
                           final_path_fn=residue_final_path) -> None:
    """Assert every length ``L0..last_done`` has its final structure on disk.

    ``progress.log`` records intent; the disk records reality, and the two can diverge
    (a ``DONE`` residue's directory deleted or lost to a scratch purge). Because resume
    only reloads the *last* ``DONE`` residue as the seed, a hole earlier would go
    unnoticed and leave a permanently broken assembled trajectory. So verify presence
    of every prior length before continuing.

    Parameters
    ----------
    final_path_fn : callable, optional
        ``(out_root, L) -> Path`` for a residue's final structure. Defaults to the CSP
        per-stage layout (:func:`residue_final_path`); the cylinder runner passes
        :func:`cylinder_final_path`.

    Raises
    ------
    SystemExit
        Naming the first missing length, rather than silently continuing.
    """
    for L in range(L0, prog.last_done_residue + 1):
        fp = final_path_fn(out_root, L)
        if not fp.is_file():
            raise SystemExit(
                f"[resume] L_{L:03d} is marked DONE but its final structure {fp} is "
                f"missing -- the output tree is incomplete. Refusing to resume (re-run "
                f"fresh, or restore the missing length).")


def drop_running_units(out_root: Path, prog: Progress) -> List[str]:
    """Remove the on-disk directory of every ``RUNNING`` unit; return their names.

    At most one unit is ``RUNNING`` at a crash (the unit in flight when the run died);
    its partial output is dropped so the unit re-runs cleanly from its persisted
    schedule row.
    """
    dropped: List[str] = []
    for unit in prog.running_units():
        d = Path(out_root) / unit
        if d.is_dir():
            shutil.rmtree(d)
            dropped.append(unit)
    return dropped


def load_final_pdb(path: Path) -> np.ndarray:
    """Reload a written nascent final structure as an ``(N, 3)`` nm coordinate array.

    The seed for the next residue / post-synthesis phase on resume. Reads the
    nascent-only ``traj_final.pdb`` (all atoms in it are the nascent chain).

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist (surfaced as an actionable resume error upstream).
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"seed structure to resume from not found: {path}")
    pos = mm.app.PDBFile(str(path)).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(pos)


def est_walltime(total_steps: int, params) -> str:
    """Nominal cost annotation for the up-front schedule report.

    Returns a parenthetical string giving the **exact** simulated time
    (``total_steps * dt``) and an explicit caveat that wall-time is only nominal:
    dt-halving stability retries re-run diverged stages at more steps, and per-step
    cost grows as the nascent chain lengthens, so actual wall-time exceeds the naive
    ``steps * const``.
    """
    sim_ns = total_steps * params.dt_ps * 1e-3
    return (f" (~{sim_ns:,.1f} ns simulated at dt={params.dt_ps} ps; wall-time is "
            f"nominal -- dt-halving retries and a growing chain push it higher)")
