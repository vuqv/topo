"""TOPO simulation runner.

Exposes :func:`mdrun`, the canonical runner, available from the shell as
``topo-mdrun -f md.ini`` or ``python -m topo.mdrun -f md.ini``.
"""
from .mdrun import mdrun

__all__ = ["mdrun"]
