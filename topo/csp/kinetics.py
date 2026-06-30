"""O'Brien continuous-synthesis kinetics: codon timing and the 3-stage schedule.

This is the timing core of the O'Brien *Continuous Synthesis Protocol*
(``continuous_synthesis_v6.py``), ported to topo as pure, side-effect-free helpers
(no OpenMM here -- just the maths). It answers one question for every residue: **how
many integration steps does each of the three elongation sub-stages run for?**

The pieces (all of them straight out of v6):

1. **Per-codon translation times.** The mRNA is split into codons; a ``trans_times``
   table maps each codon to its mean in-vivo translation time (seconds). That gives a
   per-residue **intrinsic mean first-passage time** list (:func:`codon_mfpt_list`).
2. **The 3-stage split.** For residue ``L`` the total dwell time is partitioned into
   peptidyl transfer (stage 1), translocation (stage 2) and tRNA binding/waiting
   (stage 3 = remainder). Each stage's dwell is drawn from an **exponential**
   distribution about its mean (:func:`sample_fpt`, :func:`stage_dwell_times`) --
   O'Brien's first-passage-time sampling.
3. **Time -> steps.** A dwell time in seconds is mapped to in-silico nanoseconds via
   ``scale_factor`` and then to integration steps via the time step
   (:func:`seconds_to_steps`): ``steps = t_s * 1e9 / scale_factor / dt_ns``.

(An optional per-codon ribosome-traffic correction exists in the code but is
**off by default and deferred** -- not exposed in the docs or example configs; see
``topo/csp/TODO.md``. With it off, ``real == intrinsic`` and stage 2's mean is exactly
``time_stage_2``.)

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
from importlib import resources
from typing import Dict, List, Optional, Sequence, Tuple

# --- the genetic-code stop codons (RNA). A stop terminates the codon list. -----
STOP_CODONS = ("UAA", "UAG", "UGA")

# --- bundled default codon-time table (organism/temperature universal) ---------
# The per-codon mean translation times are a property of the *organism* (E. coli)
# at a given temperature (310 K), not of the protein being synthesized, so topo.csp
# ships one and uses it whenever an INI/caller gives no explicit table.
DEFAULT_TRANS_TIMES_FILE = "ecoli_trans_times_310K.txt"  # in topo/csp/data/


def default_trans_times_path() -> str:
    """Return the filesystem path of the bundled E. coli (310 K) codon-time table.

    The table (Fluitt *et al.* 2007) is shipped as package data under
    ``topo/csp/data/`` and is the default used when no ``trans_times`` is supplied.

    Returns
    -------
    str
        Absolute path to the bundled ``ecoli_trans_times_310K.txt`` resource.
    """
    return str(resources.files("topo.csp").joinpath("data", DEFAULT_TRANS_TIMES_FILE))


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

    Parameters
    ----------
    path : str
        Path to the codon-time table file (e.g. ``trans_times.txt``).

    Returns
    -------
    dict of {str: float}
        Mapping ``codon -> mean in-vivo translation time (seconds)``. Keys are
        upper-case RNA codons (``T`` normalised to ``U``).

    Raises
    ------
    ValueError
        If a non-comment line cannot be parsed as ``CODON TIME``, or if the file
        contains no codon/time rows at all.
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

    Parameters
    ----------
    path : str
        Path to the mRNA sequence file (raw nucleotides, blank/``#`` lines ignored).
    stop_at_stop : bool, optional
        If True (default), truncate the codon list at and including the first stop
        codon (``UAA``/``UAG``/``UGA``). If False, return every codon in the file.

    Returns
    -------
    list of str
        Codons in 5'->3' order, each a 3-character upper-case RNA string.

    Raises
    ------
    ValueError
        If the concatenated sequence length is not a multiple of 3.
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

    ``mfpt[i] = trans_times[codons[i]]`` -- i.e. the per-codon **intrinsic** mean
    translation time (no ribosome-traffic correction). Raises if a codon is missing
    from the table, so an incomplete ``trans_times`` is caught early rather than
    silently mis-timing a residue.

    Parameters
    ----------
    codons : sequence of str
        Codons (upper-case RNA), e.g. the output of :func:`read_mrna`.
    trans_times : dict of {str: float}
        Codon -> mean translation time (seconds), e.g. from :func:`read_trans_times`.

    Returns
    -------
    list of float
        The intrinsic mean first-passage time (seconds) for each codon, in the same
        order as ``codons``.

    Raises
    ------
    KeyError
        If any codon in ``codons`` is absent from ``trans_times``.
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

    Parameters
    ----------
    n : int
        Length of the list to build (number of codons needed).
    uniform_mfpt : float
        The single mean translation time (seconds) assigned to every codon.

    Returns
    -------
    list of float
        A list of ``n`` identical values, all equal to ``uniform_mfpt``.

    Raises
    ------
    ValueError
        If ``uniform_mfpt <= 0``.
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
    is not on ``PATH`` (or fails to run) this returns ``None``** (the caller then
    falls back to ``real == intrinsic`` -- no traffic), so the port stays runnable
    without the compiled helper. This mirrors v6's ``ribosome_traffic <mrna>
    <trans_times> <initiation_rate>`` call.

    Parameters
    ----------
    mrna_path : str
        Path to the mRNA sequence file (passed through to the binary).
    trans_times_path : str
        Path to the codon-time table (passed through to the binary).
    initiation_rate : float
        Translation-initiation rate (1/s) -- sets how densely ribosomes load.
    binary : str, optional
        Name (or path) of the external executable to call (default
        ``"ribosome_traffic"``); resolved on ``PATH`` via :func:`shutil.which`.
    verbose : bool, optional
        If True (default), print a message when the binary is missing/fails or
        produces no numeric output.

    Returns
    -------
    list of float or None
        One traffic-corrected mean first-passage time (seconds) per codon, or
        ``None`` if the binary is unavailable, errors out, or yields no numbers
        (signalling the caller to use ``real == intrinsic``).
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

    This is v6's ``sample_fpt_dist`` (``random.expovariate(1/mean)``). A single
    rate-limiting molecular event has an exponentially distributed waiting time, so
    each sub-stage dwell is sampled this way. A non-positive mean would be
    ill-defined, so it is floored to ``1e-12`` s first.

    Parameters
    ----------
    mean_s : float
        Mean of the exponential distribution (seconds); floored to ``1e-12`` if
        non-positive.
    rng : random.Random
        Random generator to draw from (seed it for reproducible schedules).

    Returns
    -------
    float
        A single exponentially distributed dwell time (seconds).
    """
    mean_s = max(float(mean_s), 1e-12)
    return rng.expovariate(1.0 / mean_s)


def stage_dwell_times(L: int, intrinsic: Sequence[float], real: Sequence[float],
                      time_stage_1: float, time_stage_2: float,
                      rng: random.Random) -> Tuple[float, float, float]:
    """Sample the three sub-stage dwell times (seconds) for nascent length ``L``.

    Reproduces the stage-time logic of ``continuous_synthesis_v6.py`` (lines 69-86).
    ``L`` is the 1-indexed nascent
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

    Parameters
    ----------
    L : int
        1-indexed nascent-chain length being synthesized.
    intrinsic : sequence of float
        0-indexed per-codon intrinsic mFPTs (seconds); ``intrinsic[L]`` (the *next*
        codon) sets the stage-3 total. Must have length >= ``L + 1``.
    real : sequence of float
        0-indexed per-codon mFPTs *with* the ribosome-traffic correction (equals
        ``intrinsic`` when traffic is off). Used only via ``real[L-1]``.
    time_stage_1 : float
        Fixed mean peptidyl-transfer dwell (seconds).
    time_stage_2 : float
        Fixed mean translocation dwell (seconds), before any traffic correction.
    rng : random.Random
        Random generator for the exponential sampling.

    Returns
    -------
    tuple of (float, float, float)
        The three **sampled** dwell times ``(t1, t2, t3)`` in seconds, for the
        peptidyl-transfer, translocation and tRNA-binding sub-stages respectively.
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
    (like v6's ``int(...)``). A larger ``scale_factor`` therefore yields fewer steps
    (a faster run) for the same physical dwell time.

    Parameters
    ----------
    t_s : float
        In-vivo dwell time to convert (seconds).
    scale_factor : float
        In-vivo-seconds -> in-silico-nanoseconds compression factor (dimensionless).
    dt_ps : float
        Integration timestep (picoseconds).

    Returns
    -------
    int
        Number of integration steps (truncated toward zero); may be 0 for very short
        dwell times before any min/max clamp is applied.
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
    ``None`` cap = uncapped (production). The clamp limits **MD steps only**, never
    the sampled dwell times in seconds (which are returned unclamped for logging).

    Parameters
    ----------
    L : int
        1-indexed nascent-chain length being synthesized.
    intrinsic, real : sequence of float
        0-indexed per-codon mFPT lists (intrinsic and traffic-corrected); see
        :func:`stage_dwell_times`.
    time_stage_1, time_stage_2 : float
        Fixed mean peptidyl-transfer / translocation dwell times (seconds).
    scale_factor : float
        In-vivo-seconds -> in-silico-nanoseconds compression factor.
    dt_ps : float
        Integration timestep (picoseconds).
    rng : random.Random
        Random generator for the exponential dwell-time sampling.
    max_steps_per_stage : int or None, optional
        Upper clamp on each stage's step count (``None`` = uncapped / production).
    min_steps_per_stage : int, optional
        Lower clamp on each stage's step count (default 1).

    Returns
    -------
    tuple
        ``((s1, s2, s3), (t1, t2, t3))`` -- the clamped integer step counts and the
        (unclamped) sampled dwell times in seconds they were derived from.
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
      with its traffic-corrected output, else ``real == intrinsic``. When
      ``trans_times_path`` is ``None`` the **bundled E. coli (310 K) table**
      (:func:`default_trans_times_path`) is used -- the codon-time table is
      organism-universal, so only the (protein-specific) ``mrna`` is mandatory.

    ``n_codons_needed`` is the minimum list length required (``L_max + 1`` so that
    ``intrinsic[L_max]`` is valid).

    Parameters
    ----------
    n_codons_needed : int
        Minimum required list length (use ``L_max + 1``).
    uniform_ta : bool
        If True, ignore the mRNA and give every codon ``uniform_mfpt`` seconds.
    uniform_mfpt : float
        The per-codon mean time (seconds) used when ``uniform_ta`` is True.
    mrna_path : str or None
        Path to the mRNA file (required unless ``uniform_ta``).
    trans_times_path : str or None
        Path to the codon-time table. If ``None`` (and not ``uniform_ta``), the
        bundled E. coli 310 K table (:func:`default_trans_times_path`) is used.
    ribosome_traffic : bool
        If True, attempt the external ``ribosome_traffic`` correction for ``real``
        (falls back to ``real == intrinsic`` if the binary is unavailable).
    initiation_rate : float
        Translation-initiation rate (1/s), passed to the traffic binary.
    verbose : bool, optional
        Forwarded to :func:`ribosome_traffic_times` for its diagnostic messages.

    Returns
    -------
    tuple
        ``(intrinsic, real, codons)`` -- the intrinsic and (traffic-corrected) real
        per-codon mFPT lists (seconds), and the codon list (or ``None`` in the
        ``uniform_ta`` mode).

    Raises
    ------
    ValueError
        If non-uniform timing is requested without ``mrna_path``, or if the mRNA /
        traffic output has fewer than ``n_codons_needed`` entries.
    """
    if uniform_ta:
        intrinsic = uniform_mfpt_list(n_codons_needed, uniform_mfpt)
        return intrinsic, list(intrinsic), None

    if not mrna_path:
        raise ValueError("non-uniform kinetics require an `mrna` file.")
    if not trans_times_path:
        trans_times_path = default_trans_times_path()
        if verbose:
            print(f"  [kinetics] no trans_times given -- using bundled E. coli "
                  f"310 K table ({DEFAULT_TRANS_TIMES_FILE}).")
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
