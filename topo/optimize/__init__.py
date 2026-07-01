"""TOPO nscale optimizer.

Exposes :func:`optimize`, the console entry point, available from the shell as
``topo-optimize -f optimize.ini`` or ``python -m topo.optimize -f optimize.ini``,
and :func:`run_optimizer`, the importable core. The reusable building blocks
(:class:`Scorer`, :func:`read_optimize_config`, the nscale ``LADDER``, ...) are
re-exported for programmatic use.
"""
from .optimize import (
    optimize,
    run_optimizer,
    read_optimize_config,
    Scorer,
    normalize_class,
    nscale_for,
    LADDER,
    IMPLICIT_DEFAULTS,
    OPT_DEFAULTS,
    CONTROL_TYPES,
)

__all__ = [
    "optimize",
    "run_optimizer",
    "read_optimize_config",
    "Scorer",
    "normalize_class",
    "nscale_for",
    "LADDER",
    "IMPLICIT_DEFAULTS",
    "OPT_DEFAULTS",
    "CONTROL_TYPES",
]
