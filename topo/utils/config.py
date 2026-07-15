"""
Read a simulation control file (``md.ini``) into a structured configuration.

The control file is an INI file with a single ``[OPTIONS]`` section. Parsing it
used to be a ~100-line block copy-pasted into every ``run_simulation.py``. This
module centralizes that logic so a run script can simply do::

    import topo

    cfg = topo.read_simulation_config("md.ini")
    model = topo.models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())
    integrator = openmm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)

All parsed values are returned on a :class:`SimulationConfig` dataclass, with
OpenMM units already applied where appropriate (``dt``, ``ref_t``, ``tau_t``,
``ref_p``).
"""
import configparser
from dataclasses import dataclass, field
from pathlib import Path
from json import loads
from typing import Any, List, Optional

import openmm as mm
from openmm import unit


def strtobool(value):
    """Convert a truthy/falsy string to 0/1.

    Local replacement for ``distutils.util.strtobool`` (distutils was removed in
    Python 3.12). Accepts the same vocabulary: true/yes/on/1 -> 1,
    false/no/off/0 -> 0 (case-insensitive); anything else raises ValueError.
    """
    v = str(value).strip().lower()
    if v in ("y", "yes", "t", "true", "on", "1"):
        return 1
    if v in ("n", "no", "f", "false", "off", "0"):
        return 0
    raise ValueError(f"invalid truth value {value!r}")


@dataclass
class SimulationConfig:
    """Parsed contents of a simulation control file (``md.ini``)."""

    # run control
    md_steps: int = 10_000
    # default_factory: a Quantity is treated as a mutable default by dataclasses
    dt: Any = field(default_factory=lambda: 0.015 * unit.picoseconds)
    nstxout: int = 5000
    nstlog: int = 5000
    # Checkpoint-write frequency. Reading an md.ini that omits ``nstchk`` still
    # falls back to that file's ``nstxout`` (see read_simulation_config), so this
    # default applies only when a SimulationConfig is constructed directly.
    nstchk: Optional[int] = 5000
    # Center-of-mass motion removal frequency. When None (the default), no
    # CMMotionRemover is added. COM removal suits a single chain but couples the
    # drift of multiple independent chains, so it is opt-in via the config.
    nstcomm: Optional[int] = None
    # Decimal places for floating-point columns in the .log (energies, time,
    # temperature, ...). Keeps the log readable; set to None for full precision.
    log_precision: Optional[int] = 4
    # Minimum width (characters) of each .log column, for aligned fixed-width
    # output. Set to None to disable fixed-width formatting.
    log_width: Optional[int] = 14
    model: str = 'topo'

    # bond treatment: 'AllBonds' (rigid, default) or None (flexible harmonic bonds)
    constraints: Any = 'AllBonds'
    # integrator constraint tolerance (relative). Tighter than OpenMM's 1e-5
    # default is unnecessary; only meaningful when constraints = AllBonds.
    constraint_tolerance: float = 1e-5

    # temperature coupling
    tcoupl: bool = True
    ref_t: Any = 300.0
    tau_t: Any = None

    # temperature protocol: constant-T equilibrium (default) vs. annealing.
    # When anneal is on the run has TWO phases: a quench phase (hold at t_high,
    # then quench/cool to ref_t) followed by a production phase (md_steps at
    # ref_t). The quench writes <outname>_quench.dcd/.log; production keeps the
    # usual <outname>.dcd/.log. anneal_steps is SEPARATE from md_steps: the grand
    # total is anneal_steps (+ anneal_ramp_steps for linear) + md_steps. ref_t is
    # the low (refold) temperature, reused as-is (no separate t_low key).
    # See topo.mdrun.protocol and SimulationConfig.quench_steps / total_steps.
    anneal: bool = False
    t_high: Any = None              # high (unfolding) temperature; kelvin when tcoupl on
    anneal_steps: int = 0           # quench-phase steps held at t_high (separate from md_steps)
    anneal_ramp: str = 'jump'       # 'jump' (delta quench) or 'linear' (cooling ramp)
    anneal_ramp_steps: int = 0      # linear only: quench-phase steps to ramp t_high -> ref_t
    anneal_ramp_increments: int = 20  # linear only: number of discrete T steps in the ramp

    # pressure coupling
    pcoupl: bool = False
    ref_p: Any = 1.0
    frequency_p: int = 25

    # periodic boundary conditions
    pbc: bool = False
    box_dimension: Optional[List[float]] = None

    # input
    pdb_file: Optional[str] = None
    # Optional PDB of starting coordinates for a fresh run. If unset, the
    # coordinates of the structure used to build the system (pdb_file) are used.
    init_position: Optional[str] = None
    domain_def: Optional[str] = None
    stride_output_file: Optional[str] = None

    # output: all generated files go to <output_dir>/<outname><suffix>, so a run
    # is one self-contained folder (default traj/traj.dcd, traj.log, traj.psf, ...).
    output_dir: str = 'traj'
    outname: str = 'traj'
    checkpoint: Optional[str] = None     # explicit override; defaults to <output_dir>/<outname>.chk

    # hardware
    device: str = 'CPU'
    ppn: int = 1

    # multi-copy (non-interacting chains in one simulation; see topo.multichain)
    n_copies: int = 1
    copy_shift: float = 5.0   # nm, x-translation between replicated copies

    # restart / minimize
    restart: bool = False
    minimize: bool = True

    # the path this config was read from (for bookkeeping)
    config_file: Optional[str] = None

    # console verbosity: when True, suppress the per-run informational banners
    # ("Running simulation on ...", "[tracking] writing run metadata ...",
    # "Wrote last conformation ...", "--- Finished in ... ---"). Set by the CSP
    # runner, which drives one config per stage; the plain mdrun leaves it False.
    quiet: bool = False

    def build_kwargs(self) -> dict:
        """
        Keyword arguments for
        :func:`topo.models.buildCoarseGrainModel`.

        Always passes ``minimize``, ``model`` and ``box_dimension``; passes
        ``domain_def`` / ``stride_output_file`` only when they are set, so the
        builder's own defaults apply otherwise. On restart the build-time energy
        check is skipped (``check_forces=False``) because the loaded checkpoint
        state, not the input structure, is what gets simulated.
        """
        kwargs = dict(minimize=self.minimize, model=self.model,
                      box_dimension=self.box_dimension, constraints=self.constraints,
                      check_forces=not self.restart)
        if self.domain_def is not None:
            kwargs['domain_def'] = self.domain_def
        if self.stride_output_file is not None:
            kwargs['stride_output_file'] = self.stride_output_file
        return kwargs

    def output_path(self, suffix: str = '') -> str:
        """
        Path for a generated output file: ``<output_dir>/<outname><suffix>``.

        Examples: ``output_path('.dcd')`` -> ``traj/traj.dcd``;
        ``output_path('_multi.psf')`` -> ``traj/traj_multi.psf``.

        Built with :class:`pathlib.Path` but returned as ``str`` so it can be
        passed directly to OpenMM/parmed writers (some of which special-case
        ``str`` and would mishandle a raw ``Path``).
        """
        return str(Path(self.output_dir) / f'{self.outname}{suffix}')

    def checkpoint_path(self) -> str:
        """
        Resolved checkpoint path: the explicit ``checkpoint`` option if given,
        otherwise ``<output_dir>/<outname>.chk``.
        """
        return self.checkpoint if self.checkpoint else self.output_path('.chk')

    def quench_steps(self) -> int:
        """Total steps in the quench phase (``0`` when ``anneal`` is off).

        The quench is a *separate* phase from production: it holds at ``t_high``
        for ``anneal_steps`` and, for a linear ramp, additionally cools over
        ``anneal_ramp_steps``. Production then runs ``md_steps`` at ``ref_t``.
        """
        if not self.anneal:
            return 0
        ramp = self.anneal_ramp_steps if self.anneal_ramp == 'linear' else 0
        return self.anneal_steps + ramp

    def total_steps(self) -> int:
        """Grand total integration steps across the quench + production phases."""
        return self.quench_steps() + self.md_steps

    def prepare_output_dir(self) -> None:
        """
        Ensure ``output_dir`` exists. ``Path.mkdir(parents=True, exist_ok=True)``
        creates any missing parents and is a no-op if the folder already exists
        (no manual "does it exist?" check needed).
        """
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def make_platform(self):
        """
        Build the OpenMM ``(platform, properties)`` pair selected by ``device``.

        ``device = 'GPU'`` -> CUDA (mixed precision, device 0); anything else ->
        CPU using ``ppn`` threads.
        """
        if self.device == 'GPU':
            if not self.quiet:
                print("Running simulation on GPU CUDA")
            return (mm.Platform.getPlatformByName('CUDA'),
                    {'CudaPrecision': 'mixed', 'DeviceIndex': '0'})
        if not self.quiet:
            print(f"Running simulation on CPU using {self.ppn} cores")
        return (mm.Platform.getPlatformByName('CPU'), {'Threads': str(self.ppn)})


def read_simulation_config(config_file: str, verbose: bool = True) -> SimulationConfig:
    """
    Parse a simulation control file into a :class:`SimulationConfig`.

    Parameters
    ----------
    config_file : str
        Path to the control file (e.g. ``md.ini``). Must contain an
        ``[OPTIONS]`` section. Inline comments starting with ``#`` or ``;`` are
        ignored. Underscores in ``md_steps`` (e.g. ``500_000``) are allowed.
    verbose : bool, optional (default: True)
        If True, echo each parsed setting to stdout (the original behaviour of
        the inline parser). Set False for a quiet read.

    Returns
    -------
    SimulationConfig
        All settings with OpenMM units applied where relevant. Options absent
        from the file fall back to the dataclass defaults.

    Notes
    -----
    Behavioural details preserved from the original inline parser:

    - ``ref_t``/``tau_t`` get units only when ``tcoupl`` is on.
    - ``box_dimension`` accepts a scalar ``L`` (cubic ``[L, L, L]``) or an
      explicit ``[x, y, z]`` list, and is ``None`` when ``pbc`` is off.
    - ``pcoupl`` requires ``pbc`` (asserted).
    - ``restart`` forces ``minimize`` off.
    """
    def log(msg: str) -> None:
        """Print ``msg`` only when ``verbose`` is set.

        Parameters
        ----------
        msg : str
            Message echoing a parsed setting.
        """
        if verbose:
            print(msg)

    cfg = SimulationConfig(config_file=config_file)

    log(f"Reading simulation parameters from {config_file} file...")
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(config_file)
    params = config['OPTIONS']

    cfg.md_steps = int(str(params.get('md_steps', cfg.md_steps)).replace('_', ''))
    log(f'Setting number of simulation steps to: {cfg.md_steps}')
    cfg.dt = float(params.get('dt', 0.015)) * unit.picoseconds
    log(f'Setting timestep for integration of equations of motion to: {cfg.dt}')
    cfg.nstxout = int(params.get('nstxout', cfg.nstxout))
    log(f'Setting number of steps to write trajectory frame: {cfg.nstxout}')
    cfg.nstlog = int(params.get('nstlog', cfg.nstlog))
    log(f'Setting number of steps to write logfile: {cfg.nstlog}')
    # Checkpoint frequency: defaults to the trajectory frequency (nstxout) when
    # the config does not set nstchk explicitly.
    cfg.nstchk = int(params.get('nstchk', cfg.nstxout))
    log(f'Setting number of steps to write checkpoint: {cfg.nstchk}')
    nstcomm_val = params.get('nstcomm', None)
    cfg.nstcomm = int(nstcomm_val) if str(nstcomm_val).strip() not in ('None', '') else None
    if cfg.nstcomm is None:
        log('Center-of-mass motion removal is off (nstcomm not set)')
    else:
        log(f'Setting frequency of center of mass motion removal to every {cfg.nstcomm} steps')
    prec_val = params.get('log_precision', None)
    if prec_val is not None:
        # Explicit 'None'/empty disables formatting; otherwise use the integer.
        cfg.log_precision = None if str(prec_val).strip().lower() in ('none', '') else int(prec_val)
    width_val = params.get('log_width', None)
    if width_val is not None:
        cfg.log_width = None if str(width_val).strip().lower() in ('none', '') else int(width_val)
    log(f'Log columns: precision={cfg.log_precision if cfg.log_precision is not None else "full"}, '
        f'width={cfg.log_width if cfg.log_width is not None else "auto"}')
    cfg.model = params.get('model', cfg.model)
    log(f'Setting model: {cfg.model}')

    # Bond treatment: 'AllBonds' (rigid, default) or None/'None'/'' (flexible bonds).
    cfg.constraints = params.get('constraints', cfg.constraints)
    if str(cfg.constraints).strip().lower() in ('none', ''):
        cfg.constraints = None
    log(f'Bond constraints: {"None (flexible bonds)" if cfg.constraints is None else cfg.constraints}')
    cfg.constraint_tolerance = float(params.get('constraint_tolerance', cfg.constraint_tolerance))
    log(f'Constraint tolerance: {cfg.constraint_tolerance}')

    cfg.tcoupl = bool(strtobool(str(params.get('tcoupl', cfg.tcoupl))))
    if cfg.tcoupl:
        cfg.ref_t = float(params.get('ref_t', 300.0)) * unit.kelvin
        cfg.tau_t = float(params.get('tau_t', 0.05)) / unit.picoseconds
        log(f'Turning on temperature coupling with reference temperature: '
            f'{cfg.ref_t} and time constant: {cfg.tau_t}')
    else:
        log("Temperature coupling is off")

    # Temperature protocol. anneal = no (default) -> constant-T equilibrium at
    # ref_t. anneal = yes -> hold at t_high then quench/cool to ref_t.
    cfg.anneal = bool(strtobool(str(params.get('anneal', cfg.anneal))))
    if cfg.anneal:
        assert cfg.tcoupl, "annealing requires temperature coupling (tcoupl = yes)"
        t_high_val = params.get('t_high', None)
        if t_high_val is None or str(t_high_val).strip() == '':
            raise ValueError("anneal = yes requires t_high (the high/unfolding temperature)")
        cfg.t_high = float(t_high_val) * unit.kelvin
        cfg.anneal_steps = int(str(params.get('anneal_steps', cfg.anneal_steps)).replace('_', ''))
        if cfg.anneal_steps <= 0:
            raise ValueError("anneal = yes requires anneal_steps > 0 (the quench-phase hold length)")
        cfg.anneal_ramp = str(params.get('anneal_ramp', cfg.anneal_ramp)).strip().lower()
        if cfg.anneal_ramp not in ('jump', 'linear'):
            raise ValueError(f"anneal_ramp must be 'jump' or 'linear', got {cfg.anneal_ramp!r}")
        if cfg.anneal_ramp == 'linear':
            cfg.anneal_ramp_steps = int(str(params.get('anneal_ramp_steps', cfg.anneal_ramp_steps)).replace('_', ''))
            cfg.anneal_ramp_increments = int(params.get('anneal_ramp_increments', cfg.anneal_ramp_increments))
        log(f'Temperature annealing on -- separate quench phase ({cfg.quench_steps()} steps): '
            f'hold {cfg.t_high} for {cfg.anneal_steps} steps, then {cfg.anneal_ramp} to ref_t = {cfg.ref_t}'
            + (f' over {cfg.anneal_ramp_steps} steps ({cfg.anneal_ramp_increments} increments)'
               if cfg.anneal_ramp == 'linear' else '')
            + f'; then production for {cfg.md_steps} steps at ref_t. '
            f'Quench -> <outname>_quench.dcd/.log, production -> <outname>.dcd/.log. '
            f'Grand total = {cfg.total_steps()} steps.')
    else:
        log("Temperature annealing is off (constant-temperature equilibrium at ref_t)")

    cfg.pbc = bool(strtobool(str(params.get('pbc', cfg.pbc))))
    if cfg.pbc:
        box_val = params.get('box_dimension', '').strip()
        if box_val.startswith('['):
            cfg.box_dimension = loads(box_val)
        else:
            try:
                L = float(box_val)
                cfg.box_dimension = [L, L, L]
            except (ValueError, TypeError):
                cfg.box_dimension = None
        if cfg.box_dimension is not None:
            log(f'Turning on periodic boundary conditions with box dimension: '
                f'{cfg.box_dimension} nm')
        else:
            # Invalid/empty box: keep the pbc flag consistent with the (absent) box,
            # so downstream (enforcePeriodicBox, the pcoupl `assert cfg.pbc`) doesn't
            # act as if PBC were on over a system that has no periodic box vectors.
            cfg.pbc = False
            log('Periodic boundary conditions are off (invalid box_dimension)')
    else:
        cfg.box_dimension = None
        log('Periodic boundary conditions are off')

    cfg.pcoupl = bool(strtobool(str(params.get('pcoupl', cfg.pcoupl))))
    if cfg.pcoupl:
        assert cfg.pbc, "Pressure coupling requires box dimensions and periodic boundary condition is on"
        cfg.ref_p = float(params.get('ref_p', 1.0)) * unit.bar
        cfg.frequency_p = int(params.get('frequency_p', cfg.frequency_p))
        log(f'Pressure is set to reference of {cfg.ref_p} with frequency of coupling {cfg.frequency_p}')
    else:
        log("Pressure coupling is off")

    cfg.pdb_file = params.get('pdb_file', cfg.pdb_file)
    log(f'Input structure: {cfg.pdb_file}')
    cfg.init_position = params.get('init_position', cfg.init_position)
    if cfg.init_position:
        log(f'Initial coordinates from: {cfg.init_position}')
    cfg.domain_def = params.get('domain_def', cfg.domain_def)
    if cfg.domain_def:
        log(f'Domain definition file: {cfg.domain_def}')
    cfg.stride_output_file = params.get('stride_output_file', cfg.stride_output_file)
    if cfg.stride_output_file:
        log(f'STRIDE output file: {cfg.stride_output_file}')

    cfg.output_dir = params.get('output_dir', cfg.output_dir)
    cfg.outname = params.get('outname', cfg.outname)
    cfg.checkpoint = params.get('checkpoint', cfg.checkpoint)
    log(f'Output: {cfg.output_path("")}.* (dir: {cfg.output_dir}/, checkpoint: {cfg.checkpoint_path()})')

    cfg.device = params.get('device', cfg.device)
    log(f'Running simulation on {cfg.device}')
    if cfg.device == "CPU":
        cfg.ppn = int(params.get('ppn', cfg.ppn))
        log(f'Using {cfg.ppn} threads')

    cfg.n_copies = int(params.get('n_copies', cfg.n_copies))
    cfg.copy_shift = float(params.get('copy_shift', cfg.copy_shift))
    if cfg.n_copies > 1:
        log(f'Replicating into {cfg.n_copies} non-interacting copies '
            f'(x-shift {cfg.copy_shift} nm)')

    cfg.restart = bool(strtobool(str(params.get('restart', cfg.restart))))
    log(f'Restart simulation: {cfg.restart}')
    if cfg.restart:
        cfg.minimize = False
    else:
        cfg.minimize = bool(strtobool(str(params.get('minimize', cfg.minimize))))
    log(f'Perform Energy minimization of input structure: {cfg.minimize}')

    return cfg
