"""O'Brien continuous-synthesis kinetics: codon timing and the 3-stage schedule.

This is the timing core of the O'Brien *Continuous Synthesis Protocol*
(``continuous_synthesis_v6.py``), ported to topo as pure, side-effect-free helpers
(no OpenMM here -- just the maths). It answers one question for every residue: **how
many integration steps does each of the three elongation sub-stages run for?**

The pieces (all of them straight out of v6):

1. **Per-codon translation times.** The mRNA is split into codons; a ``trans_times``
   table maps each codon to its mean in-vivo translation time (seconds). That gives a
   per-residue **intrinsic mean first-passage time** list (:func:`codon_mfpt_list`).
2. **Ribosome traffic** (optional). An upstream-queue correction can stretch the
   per-residue time; :func:`ribosome_traffic_times` calls O'Brien's external
   ``ribosome_traffic`` binary if it is on ``PATH``, otherwise it returns
   ``real == intrinsic`` (no traffic) and says so.
3. **The 3-stage split.** For residue ``L`` the total dwell time is partitioned into
   peptidyl transfer (stage 1), translocation (stage 2, + traffic correction) and
   tRNA binding/waiting (stage 3 = remainder). Each stage's dwell is drawn from an
   **exponential** distribution about its mean (:func:`sample_fpt`,
   :func:`stage_dwell_times`) -- O'Brien's first-passage-time sampling.
4. **Time -> steps.** A dwell time in seconds is mapped to in-silico nanoseconds via
   ``scale_factor`` and then to integration steps via the time step
   (:func:`seconds_to_steps`): ``steps = t_s * 1e9 / scale_factor / dt_ns``.

Indexing convention (mirrors v6 exactly): the codon/mFPT lists are **0-indexed**,
``mfpt[i]`` = time of the ``i``-th codon (the codon that makes residue ``i+1``... see
:func:`stage_dwell_times` for how the 1-indexed nascent length ``L`` reads
``mfpt[L]`` / ``mfpt[L-1]``). The mRNA carries ``N+1`` codons (one per residue plus a
stop), so ``mfpt[L]`` is always in range for ``L = 1..N``.

**Units:** times are **seconds** (in-vivo) until :func:`seconds_to_steps`; the time
step is **picoseconds**; ``scale_factor`` is dimensionless.
"""
from __future__ import annotations

import os
import random
import shutil
import subprocess
from typing import Dict, List, Optional, Sequence, Tuple

# --- the genetic-code stop codons (RNA). A stop terminates the codon list. -----
STOP_CODONS = ("UAA", "UAG", "UGA")


# --------------------------------------------------------------------------
# Input tables
# --------------------------------------------------------------------------
def read_trans_times(path: str) -> Dict[str, float]:
    """Read a per-codon mean-translation-time table into ``{codon: seconds}``.

    Format (O'Brien ``trans_times`` file, e.g. the Fluitt *E. coli* table): one
    codon per line, ``CODON<whitespace>TIME`` -- the codon is RNA (``U`` not ``T``),
    the time is the mean in-vivo translation time in **seconds**. Blank lines and
    ``#`` comments are ignored. Codons are upper-cased and ``T`` is normalised to
    ``U`` so a DNA-style table still works.
    """
    table: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}: cannot parse trans_times line: {line!r}")
            codon = parts[0].upper().replace("T", "U")
            table[codon] = float(parts[1])
    if not table:
        raise ValueError(f"{path}: no codon/time rows found.")
    return table


def read_mrna(path: str, stop_at_stop: bool = True) -> List[str]:
    """Read an mRNA sequence file and split it into a list of 3-nt codons.

    The file is raw nucleotides (``A/U/G/C``; ``T`` normalised to ``U``), optionally
    wrapped across several lines (whitespace is stripped and concatenated). The
    sequence length must be a multiple of 3. If ``stop_at_stop`` the list is
    truncated at (and including) the first stop codon -- matching v6, which expects
    ``len(codons) == n_residues + 1`` (one codon per residue plus the terminator).
    """
    seq = ""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            seq += line
    seq = seq.upper().replace("T", "U").replace(" ", "")
    if len(seq) % 3 != 0:
        raise ValueError(f"{path}: mRNA length {len(seq)} is not a multiple of 3.")
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    if stop_at_stop:
        for i, c in enumerate(codons):
            if c in STOP_CODONS:
                return codons[:i + 1]
    return codons


def codon_mfpt_list(codons: Sequence[str],
                    trans_times: Dict[str, float]) -> List[float]:
    """Map a codon list to a 0-indexed list of mean first-passage times (seconds).

    ``mfpt[i] = trans_times[codons[i]]``. Raises if a codon is missing from the
    table (so an incomplete ``trans_times`` is caught early rather than silently
    mis-timing a residue).
    """
    out: List[float] = []
    for i, c in enumerate(codons):
        if c not in trans_times:
            raise KeyError(f"codon #{i + 1} {c!r} not found in trans_times table.")
        out.append(float(trans_times[c]))
    return out


def uniform_mfpt_list(n: int, uniform_mfpt: float) -> List[float]:
    """A constant mFPT list of length ``n`` (the ``uniform_ta = 1`` mode of v6).

    Every codon gets the same mean translation time ``uniform_mfpt`` (seconds);
    used when no mRNA / ``trans_times`` are supplied.
    """
    if uniform_mfpt <= 0:
        raise ValueError("uniform_mfpt must be > 0 when uniform_ta is on.")
    return [float(uniform_mfpt)] * int(n)


# --------------------------------------------------------------------------
# Ribosome traffic (optional external correction)
# --------------------------------------------------------------------------
def ribosome_traffic_times(mrna_path: str, trans_times_path: str,
                           initiation_rate: float,
                           binary: str = "ribosome_traffic",
                           verbose: bool = True) -> Optional[List[float]]:
    """Return per-codon *real* mFPTs from O'Brien's ``ribosome_traffic`` binary.

    The binary models upstream-queue (traffic) effects: given the mRNA, the
    intrinsic per-codon times and the initiation rate it prints one traffic-corrected
    mean first-passage time per codon. We capture that into a list. **If the binary
    is not on ``PATH`` this returns ``None``** (the caller then falls back to
    ``real == intrinsic`` -- no traffic), so the port stays runnable without the
    compiled helper. This mirrors v6's ``ribosome_traffic <mrna> <trans_times>
    <initiation_rate>`` call.
    """
    exe = shutil.which(binary)
    if exe is None:
        if verbose:
            print(f"  [ribosome_traffic] binary {binary!r} not found on PATH; "
                  f"falling back to real == intrinsic (no traffic correction).")
        return None
    cmd = [exe, mrna_path, trans_times_path, str(initiation_rate)]
    if verbose:
        print(f"  [ribosome_traffic] running: {' '.join(cmd)}")
    # The helper is a compiled external binary with its own runtime/library needs;
    # if it cannot run (missing shared libs, bad args, non-zero exit) we must not
    # take the whole synthesis down -- degrade to no traffic and warn.
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, OSError) as exc:
        if verbose:
            print(f"  [ribosome_traffic] {binary} failed to run ({exc}); "
                  f"falling back to real == intrinsic (no traffic correction).")
        return None
    times: List[float] = []
    for tok in proc.stdout.split():
        try:
            times.append(float(tok))
        except ValueError:
            continue
    if not times:
        if verbose:
            print(f"  [ribosome_traffic] {binary} produced no numeric output; "
                  f"falling back to real == intrinsic.")
        return None
    return times


# --------------------------------------------------------------------------
# First-passage-time sampling and the 3-stage split
# --------------------------------------------------------------------------
def sample_fpt(mean_s: float, rng: random.Random) -> float:
    """Draw one first-passage time (seconds) from an exponential of mean ``mean_s``.

    This is v6's ``sample_fpt_dist`` (``random.expovariate(1/mean)``). A non-positive
    mean would be ill-defined, so it is floored to a tiny positive value first.
    """
    mean_s = max(float(mean_s), 1e-12)
    return rng.expovariate(1.0 / mean_s)


def stage_dwell_times(L: int, intrinsic: Sequence[float], real: Sequence[float],
                      time_stage_1: float, time_stage_2: float,
                      rng: random.Random) -> Tuple[float, float, float]:
    """Sample the three sub-stage dwell times (seconds) for nascent length ``L``.

    Reproduces ``run_elongation`` lines 69-86 of v6. ``L`` is the 1-indexed nascent
    chain length; ``intrinsic`` / ``real`` are 0-indexed per-codon mFPT lists (see
    module docstring). The three means are:

    - **stage 1** (peptidyl transfer): ``time_stage_1`` -- a fixed mean.
    - **stage 2** (translocation): ``time_stage_2`` plus the ribosome-traffic
      correction ``real[L-1] - intrinsic[L-1]`` **if that is positive** (else just
      ``time_stage_2`` -- v6's guard against a negative correction from sampling
      noise).
    - **stage 3** (tRNA binding / waiting): the remainder
      ``intrinsic[L] - time_stage_1 - time_stage_2`` -- floored to a tiny positive
      value if a fast codon makes it non-positive.

    Each mean is then passed through :func:`sample_fpt` (exponential sampling).
    Returns the three **sampled** dwell times in seconds.
    """
    # stage 1 -- fixed mean.
    t1 = sample_fpt(time_stage_1, rng)

    # stage 2 -- base + traffic correction (guarded non-negative).
    correction = real[L - 1] - intrinsic[L - 1]
    mean2 = time_stage_2 + correction if correction > 0 else time_stage_2
    t2 = sample_fpt(mean2, rng)

    # stage 3 -- remainder of the codon's total dwell time.
    mean3 = intrinsic[L] - time_stage_1 - time_stage_2
    if mean3 <= 0:
        mean3 = 1e-9  # fast codon: t1+t2 already exceed the codon time.
    t3 = sample_fpt(mean3, rng)

    return t1, t2, t3


def seconds_to_steps(t_s: float, scale_factor: float, dt_ps: float) -> int:
    """Map an in-vivo dwell time (s) to a number of integration steps.

    O'Brien's two-step conversion: ``t_sim_ns = t_s * 1e9 / scale_factor`` (the
    ``scale_factor`` compresses real time into the in-silico timescale), then
    ``steps = t_sim_ns / dt_ns`` with ``dt_ns = dt_ps * 1e-3``. Truncated to an int
    (like v6's ``int(...)``).
    """
    t_sim_ns = t_s * 1e9 / scale_factor
    dt_ns = dt_ps * 1e-3
    return int(t_sim_ns / dt_ns)


def stage_steps(L: int, intrinsic: Sequence[float], real: Sequence[float],
                *, time_stage_1: float, time_stage_2: float,
                scale_factor: float, dt_ps: float, rng: random.Random,
                max_steps_per_stage: Optional[int] = None,
                min_steps_per_stage: int = 1) -> Tuple[Tuple[int, int, int],
                                                       Tuple[float, float, float]]:
    """Full per-residue schedule: sampled dwell times **and** clamped step counts.

    Combines :func:`stage_dwell_times` + :func:`seconds_to_steps` and applies the
    test clamps. ``max_steps_per_stage`` caps each stage (the tutorial uses a small
    cap so a residue runs ~2000 steps total instead of the production ~10^5-10^6);
    ``min_steps_per_stage`` floors it so every stage does at least a little MD.
    ``None`` cap = uncapped (production). Returns ``((s1,s2,s3), (t1,t2,t3))`` --
    the integer step counts and the seconds dwell times they came from (the latter
    handy for logging the O'Brien table).
    """
    t1, t2, t3 = stage_dwell_times(L, intrinsic, real, time_stage_1, time_stage_2, rng)
    steps = []
    for t in (t1, t2, t3):
        s = seconds_to_steps(t, scale_factor, dt_ps)
        if max_steps_per_stage is not None:
            s = min(s, int(max_steps_per_stage))
        s = max(s, int(min_steps_per_stage))
        steps.append(s)
    return (steps[0], steps[1], steps[2]), (t1, t2, t3)


# --------------------------------------------------------------------------
# Convenience: build the intrinsic / real lists from config inputs
# --------------------------------------------------------------------------
def build_mfpt_lists(n_codons_needed: int, *,
                     uniform_ta: bool, uniform_mfpt: float,
                     mrna_path: Optional[str], trans_times_path: Optional[str],
                     ribosome_traffic: bool, initiation_rate: float,
                     verbose: bool = True) -> Tuple[List[float], List[float], Optional[List[str]]]:
    """Assemble the ``(intrinsic, real, codons)`` lists for a run.

    - ``uniform_ta``: every codon gets ``uniform_mfpt`` (no mRNA needed); ``real ==
      intrinsic`` and ``codons`` is ``None``.
    - otherwise: read the mRNA + ``trans_times``, build the intrinsic per-codon list;
      if ``ribosome_traffic`` and the external binary is available, replace ``real``
      with its traffic-corrected output, else ``real == intrinsic``.

    ``n_codons_needed`` is the minimum list length required (``L_max + 1`` so that
    ``intrinsic[L_max]`` is valid). Raises if the supplied mRNA is too short.
    """
    if uniform_ta:
        intrinsic = uniform_mfpt_list(n_codons_needed, uniform_mfpt)
        return intrinsic, list(intrinsic), None

    if not mrna_path or not trans_times_path:
        raise ValueError("non-uniform kinetics require both `mrna` and `trans_times`.")
    trans = read_trans_times(trans_times_path)
    codons = read_mrna(mrna_path)
    intrinsic = codon_mfpt_list(codons, trans)
    if len(intrinsic) < n_codons_needed:
        raise ValueError(
            f"mRNA has {len(intrinsic)} codons but the schedule needs at least "
            f"{n_codons_needed} (L_max + 1). Provide a longer mRNA or lower L_max.")
    real = list(intrinsic)
    if ribosome_traffic:
        traffic = ribosome_traffic_times(mrna_path, trans_times_path,
                                         initiation_rate, verbose=verbose)
        if traffic is not None:
            if len(traffic) < n_codons_needed:
                raise ValueError(
                    f"ribosome_traffic returned {len(traffic)} times but "
                    f"{n_codons_needed} are needed.")
            real = traffic
    return intrinsic, real, codons
