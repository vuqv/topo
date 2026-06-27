"""Temperature protocol for the runner: equilibrium and annealing/quenching.

A *protocol* is just a list of ``(temperature, n_steps)`` stages whose step
counts sum to ``cfg.md_steps``. Equilibrium is the degenerate one-stage schedule
``[(ref_t, md_steps)]``; annealing holds the system at ``t_high`` and then either
jumps (delta quench) or linearly ramps down to ``ref_t``. In every case
``ref_t`` is the low / refold temperature -- there is no separate ``t_low`` key.

:func:`run_protocol` drives an existing ``Simulation`` through the schedule,
setting the integrator temperature per stage. It also resumes a restart
mid-schedule by skipping already-completed stages.
"""
from openmm import unit


def temperature_schedule(cfg):
    """Build the ``[(temperature, n_steps), ...]`` schedule from the config.

    Equilibrium (``anneal = no``)::

        [(ref_t, md_steps)]

    Annealing, delta quench (``anneal_ramp = jump``)::

        [(t_high, anneal_steps), (ref_t, md_steps - anneal_steps)]

    Annealing, linear ramp (``anneal_ramp = linear``): the hold at ``t_high`` is
    followed by ``anneal_ramp_increments`` discrete cooling stages spanning
    ``anneal_ramp_steps`` (T decreasing ``t_high -> ref_t``), then a constant
    ``ref_t`` hold for whatever steps remain.

    Step counts always sum to ``cfg.md_steps``.
    """
    if not cfg.anneal:
        return [(cfg.ref_t, cfg.md_steps)]

    pre = cfg.anneal_steps
    ramp = cfg.anneal_ramp_steps if cfg.anneal_ramp == 'linear' else 0
    if pre + ramp > cfg.md_steps:
        raise SystemExit(
            f"annealing schedule needs {pre + ramp} steps (anneal_steps={pre} + "
            f"anneal_ramp_steps={ramp}) but md_steps={cfg.md_steps}; reduce the "
            f"hold/ramp or increase md_steps.")

    stages = []
    if pre > 0:
        stages.append((cfg.t_high, pre))

    if ramp > 0:
        n_inc = max(1, cfg.anneal_ramp_increments)
        t_high = cfg.t_high.value_in_unit(unit.kelvin)
        t_low = cfg.ref_t.value_in_unit(unit.kelvin)
        base = ramp // n_inc
        for i in range(1, n_inc + 1):
            # Temperature at the END of increment i (so the final increment lands
            # exactly on ref_t). Distribute the remainder onto the last increment.
            frac = i / n_inc
            temperature = (t_high + (t_low - t_high) * frac) * unit.kelvin
            n = base if i < n_inc else ramp - base * (n_inc - 1)
            if n > 0:
                stages.append((temperature, n))

    remaining = cfg.md_steps - pre - ramp
    if remaining > 0:
        stages.append((cfg.ref_t, remaining))
    return stages


def describe_schedule(schedule):
    """One-line human summary of a schedule, e.g. ``600 K x 100000 -> 300 K x 900000``."""
    parts = []
    for temperature, n in schedule:
        try:
            t = f"{temperature.value_in_unit(unit.kelvin):g} K"
        except AttributeError:
            t = f"{temperature} K"
        parts.append(f"{t} x {n}")
    return " -> ".join(parts)


def run_protocol(simulation, schedule, done_steps=0):
    """Step ``simulation`` through ``schedule``, one stage at a time.

    ``done_steps`` lets a restart resume mid-schedule: stages already completed
    are skipped, and a stage that was partially done runs only its remaining
    steps. For equilibrium (single stage) this reduces exactly to the original
    ``simulation.step(md_steps - done_steps)``.
    """
    consumed = 0
    for temperature, n in schedule:
        if n <= 0:
            continue
        stage_end = consumed + n
        if stage_end <= done_steps:
            consumed = stage_end
            continue
        run = stage_end - max(consumed, done_steps)
        simulation.integrator.setTemperature(temperature)
        simulation.step(run)
        consumed = stage_end
