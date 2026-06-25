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
from distutils.util import strtobool
from json import loads
from typing import Any, List, Optional

from openmm import unit


@dataclass
class SimulationConfig:
    """Parsed contents of a simulation control file (``md.ini``)."""

    # run control
    md_steps: int = 1000
    # default_factory: a Quantity is treated as a mutable default by dataclasses
    dt: Any = field(default_factory=lambda: 0.01 * unit.picoseconds)
    nstxout: int = 10
    nstlog: int = 10
    nstcomm: int = 100
    model: str = 'topo'

    # bond treatment: 'AllBonds' (rigid, default) or None (flexible harmonic bonds)
    constraints: Any = 'AllBonds'

    # temperature coupling
    tcoupl: bool = True
    ref_t: Any = 300.0
    tau_t: Any = None

    # pressure coupling
    pcoupl: bool = False
    ref_p: Any = 1.0
    frequency_p: int = 25

    # periodic boundary conditions
    pbc: bool = False
    box_dimension: Optional[List[float]] = None

    # input / output
    pdb_file: Optional[str] = None
    protein_code: Optional[str] = None
    domain_def: Optional[str] = None
    stride_output_file: Optional[str] = None
    checkpoint: Optional[str] = None

    # hardware
    device: str = 'CPU'
    ppn: int = 1

    # restart / minimize
    restart: bool = False
    minimize: bool = True

    # the path this config was read from (for bookkeeping)
    config_file: Optional[str] = None

    def build_kwargs(self) -> dict:
        """
        Keyword arguments for
        :func:`topo.models.buildCoarseGrainModel`.

        Always passes ``minimize``, ``model`` and ``box_dimension``; passes
        ``domain_def`` / ``stride_output_file`` only when they are set, so the
        builder's own defaults apply otherwise.
        """
        kwargs = dict(minimize=self.minimize, model=self.model,
                      box_dimension=self.box_dimension, constraints=self.constraints)
        if self.domain_def is not None:
            kwargs['domain_def'] = self.domain_def
        if self.stride_output_file is not None:
            kwargs['stride_output_file'] = self.stride_output_file
        return kwargs


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
        if verbose:
            print(msg)

    cfg = SimulationConfig(config_file=config_file)

    log(f"Reading simulation parameters from {config_file} file...")
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(config_file)
    params = config['OPTIONS']

    cfg.md_steps = int(str(params.get('md_steps', cfg.md_steps)).replace('_', ''))
    log(f'Setting number of simulation steps to: {cfg.md_steps}')
    cfg.dt = float(params.get('dt', 0.01)) * unit.picoseconds
    log(f'Setting timestep for integration of equations of motion to: {cfg.dt}')
    cfg.nstxout = int(params.get('nstxout', cfg.nstxout))
    log(f'Setting number of steps to write checkpoint and coordinate: {cfg.nstxout}')
    cfg.nstlog = int(params.get('nstlog', cfg.nstlog))
    log(f'Setting number of steps to write logfile: {cfg.nstlog}')
    cfg.nstcomm = int(params.get('nstcomm', cfg.nstcomm))
    log(f'Setting frequency of center of mass motion removal to every {cfg.nstcomm} steps')
    cfg.model = params.get('model', cfg.model)
    log(f'Setting model: {cfg.model}')

    # Bond treatment: 'AllBonds' (rigid, default) or None/'None'/'' (flexible bonds).
    cfg.constraints = params.get('constraints', cfg.constraints)
    if str(cfg.constraints).strip().lower() in ('none', ''):
        cfg.constraints = None
    log(f'Bond constraints: {"None (flexible bonds)" if cfg.constraints is None else cfg.constraints}')

    cfg.tcoupl = bool(strtobool(str(params.get('tcoupl', cfg.tcoupl))))
    if cfg.tcoupl:
        cfg.ref_t = float(params.get('ref_t', 300.0)) * unit.kelvin
        cfg.tau_t = float(params.get('tau_t', 0.01)) / unit.picoseconds
        log(f'Turning on temperature coupling with reference temperature: '
            f'{cfg.ref_t} and time constant: {cfg.tau_t}')
    else:
        log("Temperature coupling is off")

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
    cfg.protein_code = params.get('protein_code', cfg.protein_code)
    log(f'Prefix use to write file: {cfg.protein_code}')
    cfg.domain_def = params.get('domain_def', cfg.domain_def)
    if cfg.domain_def:
        log(f'Domain definition file: {cfg.domain_def}')
    cfg.stride_output_file = params.get('stride_output_file', cfg.stride_output_file)
    if cfg.stride_output_file:
        log(f'STRIDE output file: {cfg.stride_output_file}')
    cfg.checkpoint = params.get('checkpoint', cfg.checkpoint)

    cfg.device = params.get('device', cfg.device)
    log(f'Running simulation on {cfg.device}')
    if cfg.device == "CPU":
        cfg.ppn = int(params.get('ppn', cfg.ppn))
        log(f'Using {cfg.ppn} threads')

    cfg.restart = bool(strtobool(str(params.get('restart', cfg.restart))))
    log(f'Restart simulation: {cfg.restart}')
    if cfg.restart:
        cfg.minimize = False
    else:
        cfg.minimize = bool(strtobool(str(params.get('minimize', cfg.minimize))))
    log(f'Perform Energy minimization of input structure: {cfg.minimize}')

    return cfg
