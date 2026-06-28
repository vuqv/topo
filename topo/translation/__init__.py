"""topo.translation -- co-translational (protein synthesis) modeling on TOPO.

See DESIGN.md and FILES.md in this folder. Provides the ribosome coarse-graining
tool (:mod:`topo.translation.cg_ribosome`) and the nascent-chain elongation
runner (:mod:`topo.translation.elongate`), available from the shell as
``topo-elongate`` or ``python -m topo.translation``.
"""
# Re-export the high-level runner and its parameters. The CLI entry point
# (the ``elongate`` function) is intentionally NOT re-exported here, so the
# ``elongate`` name keeps referring to the submodule
# (:mod:`topo.translation.elongate`) rather than being shadowed by the function.
from .elongate import (run_elongation, ElongationParams,
                       read_elongate_config, ElongateConfig)

__all__ = ["run_elongation", "ElongationParams",
           "read_elongate_config", "ElongateConfig"]
