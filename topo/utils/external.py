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
import subprocess
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
        "set {env}=/path/to/{name}, or build it via "
        "'scripts/install_deps.sh {name}'.".format(
            name=name, env=env_var
        )
    )


def run_stride(pdb_file, out_dir=None, timeout=60, require="LOC"):
    """Run ``stride -h`` on ``pdb_file`` and cache the output to
    ``{stem}_stride.dat``; return that path.

    Shared STRIDE runner used by the model build (hydrogen-bond input) and by
    mirror detection (secondary-structure segments). Locate the executable, run
    it, and validate by output content rather than the exit code -- some STRIDE
    builds return non-zero even on success.

    Parameters
    ----------
    pdb_file : str
        All-atom PDB to run STRIDE on.
    out_dir : str or Path, optional
        Directory to write ``{stem}_stride.dat`` into. Defaults to the PDB's
        own directory.
    timeout : int, optional
        Subprocess timeout in seconds (default 60).
    require : str, optional
        Record marker the output must contain to count as valid. ``"LOC"``
        (secondary-structure records, the default) for SS consumers such as
        mirror detection; ``"DNR"`` (donor hydrogen-bond records) for the model
        build, which parses H-bonds via ``parse_hydrogen_bonds``.

    Raises
    ------
    RuntimeError
        If STRIDE cannot be located, or produces no ``require`` records.
    """
    try:
        stride_exe = find_executable("stride")
    except RuntimeError as exc:
        raise RuntimeError(
            "No STRIDE output supplied and STRIDE could not be located. "
            "Provide a precomputed STRIDE output file or make STRIDE available. "
            + str(exc)
        ) from None
    stem = os.path.splitext(os.path.basename(pdb_file))[0]
    out_dir = Path(out_dir) if out_dir is not None else Path(pdb_file).resolve().parent
    stride_path = out_dir / f"{stem}_stride.dat"
    print(f"Running STRIDE (stride -h {pdb_file} -> {stride_path}).")
    result = subprocess.run(
        [stride_exe, "-h", str(pdb_file)],
        capture_output=True, text=True, timeout=timeout,
    )
    if require not in result.stdout:
        raise RuntimeError(
            f"STRIDE produced no {require} records for {pdb_file} "
            f"(exit code {result.returncode}). stderr: {result.stderr}"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    stride_path.write_text(result.stdout)
    return str(stride_path)
