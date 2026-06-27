"""Temperature protocol for the runner: equilibrium and annealing/quenching.

A *protocol* is a list of ``(temperature, n_steps)`` stages. The runner builds up
to two of them and runs each into its own output files:

* :func:`quench_schedule`     -- the **quench phase** (only when ``anneal = yes``):
  hold at ``t_high`` and, for a linear ramp, cool down to ``ref_t``. Written to
  ``<outname>_quench.dcd`` / ``.log``. Sums to ``cfg.quench_steps()``.
* :func:`production_schedule` -- the **production phase**: ``md_steps`` at
  ``ref_t``. Written to the usual ``<outname>.dcd`` / ``.log``. This is also the
  whole schedule for a plain equilibrium run.

``anneal_steps`` is therefore *separate* from ``md_steps`` (the grand total is
``quench_steps() + md_steps``), and ``ref_t`` is always the low / refold
temperature -- there is no separate ``t_low`` key.

:func:`run_protocol` drives an existing ``Simulation`` through one schedule,
setting the integrator temperature per stage, and can resume a restart partway
through a stage.
"""
from openmm import unit


def quench_schedule(cfg):
    """Stages for the quench phase, summing to ``cfg.quench_steps()``.

    Returns ``[]`` when annealing is off. Otherwise:

    * ``anneal_ramp = jump``   -> ``[(t_high, anneal_steps)]``. The actual quench
      to ``ref_t`` happens at the phase boundary, when production starts, so the
      delta T-jump lands exactly between the ``_quench`` and production files.
    * ``anneal_ramp = linear`` -> the hold at ``t_high`` followed by
      ``anneal_ramp_increments`` discrete cooling stages spanning
      ``anneal_ramp_steps`` (T decreasing ``t_high -> ref_t``).
    """
    if not cfg.anneal:
        return []

    stages = [(cfg.t_high, cfg.anneal_steps)]

    if cfg.anneal_ramp == 'linear' and cfg.anneal_ramp_steps > 0:
        n_inc = max(1, cfg.anneal_ramp_increments)
        t_high = cfg.t_high.value_in_unit(unit.kelvin)
        t_low = cfg.ref_t.value_in_unit(unit.kelvin)
        base = cfg.anneal_ramp_steps // n_inc
        for i in range(1, n_inc + 1):
            # Temperature at the END of increment i (so the final increment lands
            # exactly on ref_t). Put the division remainder on the last increment.
            frac = i / n_inc
            temperature = (t_high + (t_low - t_high) * frac) * unit.kelvin
            n = base if i < n_inc else cfg.anneal_ramp_steps - base * (n_inc - 1)
            if n > 0:
                stages.append((temperature, n))
    return stages


def production_schedule(cfg):
    """The production phase: a single constant-``ref_t`` stage of ``md_steps``."""
    return [(cfg.ref_t, cfg.md_steps)]


def describe_schedule(schedule):
    """One-line human summary, e.g. ``600 K x 100000 -> 300 K x 900000``."""
    parts = []
    for temperature, n in schedule:
        try:
            t = f"{temperature.value_in_unit(unit.kelvin):g} K"
        except AttributeError:
            t = f"{temperature} K"
        parts.append(f"{t} x {n}")
    return " -> ".join(parts) if parts else "(none)"


def run_protocol(simulation, schedule, done_steps=0):
    """Step ``simulation`` through ``schedule``, one stage at a time.

    ``done_steps`` is measured relative to the start of this schedule. It lets a
    restart resume mid-schedule: stages already completed are skipped, and a
    stage that was partially done runs only its remaining steps. For a single
    stage this reduces to ``simulation.step(n - done_steps)``.
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
