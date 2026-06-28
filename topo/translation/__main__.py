"""Entry point for ``python -m topo.translation.elongate`` style invocation.

``python -m topo.translation`` forwards to the elongation runner.
"""
from .elongate import elongate

if __name__ == "__main__":
    elongate()
