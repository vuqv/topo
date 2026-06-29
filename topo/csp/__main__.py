"""Entry point for ``python -m topo.csp``.

Forwards to the O'Brien continuous-synthesis runner.
"""
from .csp import csp

if __name__ == "__main__":
    csp()
