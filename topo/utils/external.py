"""Locate external third-party executables (STRIDE, PULCHRA) shipped alongside
or installed for ``topo``.

Resolution order for each tool (first hit wins):

1. Explicit override via environment variable (e.g. ``TOPO_STRIDE``) -- a full
   path to the executable. Lets users point at their own build.
2. The program on ``PATH`` (``shutil.which``) -- the conventional install.
3. A binary vendored inside the package at ``topo/bin/<name>`` -- only present
   if the wheel was built with bundled binaries (off by default; see
   ``scripts/install_deps.sh`` for building/installing the tools instead).

If none resolve, a ``RuntimeError`` with an actionable message is raised.
"""

import os
import shutil
import stat
from pathlib import Path

# Directory where optional vendored binaries would live if bundled.
_BIN_DIR = Path(__file__).resolve().parent.parent / "bin"


def find_executable(name, env_var=None):
    """Return an absolute path to external executable ``name``.

    Parameters
    ----------
    name : str
        Program name, e.g. ``"stride"`` or ``"pulchra"``.
    env_var : str, optional
        Environment variable holding an explicit override path. Defaults to
        ``TOPO_<NAME>`` (e.g. ``TOPO_STRIDE``).

    Raises
    ------
    RuntimeError
        If the executable cannot be located by any method.
    """
    env_var = env_var or "TOPO_{}".format(name.upper())

    # 1. Explicit environment override.
    override = os.environ.get(env_var)
    if override:
        p = Path(override).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
        raise RuntimeError(
            "{env}={val} is set but is not an executable file.".format(
                env=env_var, val=override
            )
        )

    # 2. On PATH.
    on_path = shutil.which(name)
    if on_path:
        return on_path

    # 3. Vendored inside the package (topo/bin/<name>).
    bundled = _BIN_DIR / name
    if bundled.is_file():
        # Ensure it is executable (wheels may strip the bit).
        mode = bundled.stat().st_mode
        if not mode & stat.S_IXUSR:
            bundled.chmod(mode | stat.S_IXUSR)
        return str(bundled)

    raise RuntimeError(
        "'{name}' executable not found. Install it and ensure it is on PATH, "
        "set {env}=/path/to/{name}, or build it via scripts/install_deps.sh.".format(
            name=name, env=env_var
        )
    )
