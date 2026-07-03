"""``topo.csp`` -- the O'Brien Continuous Synthesis Protocol, ported to topo.

The per-codon, three-stage co-translational synthesis protocol of
``continuous_synthesis_v6.py`` (Yang Jiang, Dan Nissley, Ed O'Brien), expressed in
topo style. It is the kinetic upgrade of :mod:`topo.csp.core`: it times
every residue from its codon and splits each into peptidyl-transfer / translocation /
tRNA-binding sub-stages, reusing the translation module's per-length MD machinery.

CLI::

    topo-csp -f csp.ini
    python -m topo.csp -f csp.ini

See :mod:`topo.csp.protocol` (runner + INI) and :mod:`topo.csp.kinetics` (timing core).
"""
from topo.csp.core import RunParams
from topo.csp.protocol import (CSPConfig, csp, read_csp_config,
                          run_continuous_synthesis)
from topo.csp import kinetics

__all__ = [
    "CSPConfig",
    "RunParams",
    "csp",
    "read_csp_config",
    "run_continuous_synthesis",
    "kinetics",
]
