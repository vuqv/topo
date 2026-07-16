#!/usr/bin/env python3
"""
simulated_tempering.py -- enhanced sampling via OpenMM's SimulatedTempering.

NOTE ON METHOD: simulated tempering is NOT parallel tempering / replica exchange.
It runs ONE copy of the system that periodically attempts to jump up/down a
temperature ladder, accepting with a Metropolis rule that uses per-temperature
free-energy *weights*. With one walker there is no exchange between replicas and
no need for a multi-GPU setup -- it is the cheap, single-trajectory cousin of
parallel tempering. Weights are unknown a priori; by default OpenMM learns them
on the fly (Wang-Landau style), which is convenient but means early sampling is
not yet canonical. For the true replica-exchange method see
``examples/parallel_tempering_openmmtools/`` (or the from-scratch
``examples/parallel_tempering_scratch/``).

Built on TOPO's public building blocks: it reuses the ``md.ini`` system and the
``topo.engine`` build/setup/reporter/finalize flow, then wraps the resulting
``Simulation`` in ``openmm.app.SimulatedTempering``.

Outputs (under the md.ini ``output_dir``)
    <outname>.dcd / .log     the single trajectory + standard TOPO log
    <outname>_st.log         SimulatedTempering report: temperature index,
                             temperature, and learned weights over time

Run::

    python examples/simulated_tempering/simulated_tempering.py -f md.ini \
        --tmin 300 --tmax 500 --nreplicas 8 --temp-change-interval 100

``ref_t`` in the md.ini is ignored (the ladder sets the temperatures). The run
length is the md.ini ``md_steps``.
"""
import argparse
import sys
import time

import openmm as mm
from openmm import unit
from openmm.app import SimulatedTempering

from topo import engine
from topo.utils.config import read_simulation_config


def parse_args():
    parser = argparse.ArgumentParser(
        prog="simulated_tempering.py",
        description="Simulated tempering on a TOPO CG model (single walker).")
    parser.add_argument('-input', '-f', required=True,
                        help='simulation config file (md.ini) -- supplies the system')
    parser.add_argument('--tmin', type=float, required=True,
                        help='lowest temperature in the ladder (Kelvin)')
    parser.add_argument('--tmax', type=float, required=True,
                        help='highest temperature in the ladder (Kelvin)')
    parser.add_argument('--nreplicas', type=int, required=True,
                        help='number of temperature rungs (geometric spacing)')
    parser.add_argument('--temp-change-interval', type=int, default=100,
                        help='MD steps between temperature-jump attempts (default 100)')
    parser.add_argument('--report-interval', type=int, default=None,
                        help='steps between _st.log lines (default: nstlog)')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"OpenMM version: {mm.__version__}")

    # --- Base system + single Simulation via the standard TOPO flow ----------
    cfg = read_simulation_config(args.input)
    cfg.prepare_output_dir()

    built = engine.build_system(cfg)
    ctx = engine.setup_simulation(cfg, built, control_file=args.input)
    sim = ctx.simulation
    engine.attach_reporters(cfg, sim, suffix='', append=ctx.restart_active,
                            total_steps=cfg.md_steps)

    if cfg.minimize:
        print("Minimizing starting structure...")
        sim.minimizeEnergy()

    # --- Wrap in SimulatedTempering -----------------------------------------
    # numTemperatures + min/maxTemperature -> a geometric ladder built by OpenMM.
    # weights=None -> weights are learned adaptively during the run.
    report_interval = args.report_interval if args.report_interval is not None else cfg.nstlog
    st_log = cfg.output_path('_st.log')
    st = SimulatedTempering(
        sim,
        numTemperatures=args.nreplicas,
        minTemperature=args.tmin * unit.kelvin,
        maxTemperature=args.tmax * unit.kelvin,
        tempChangeInterval=args.temp_change_interval,
        reportInterval=report_interval,
        reportFile=st_log,
    )
    print(f"Simulated tempering: {args.nreplicas} rungs "
          f"{args.tmin:.1f}-{args.tmax:.1f} K, jump attempt every "
          f"{args.temp_change_interval} steps; report -> {st_log}")

    print("Simulation started")
    start_time = time.time()

    # SimulatedTempering.step drives the integrator AND the temperature moves.
    st.step(cfg.md_steps)

    # Visit histogram per rung is a useful diagnostic: flat-ish occupancy means
    # the learned weights are reasonable. (_histogram counts temperature-change
    # attempts landing on each rung; weights are the learned free energies.)
    print("Occupancy / learned weights per temperature rung:")
    for i, t in enumerate(st.temperatures):
        print(f"  rung {i:2d}  {t.value_in_unit(unit.kelvin):6.1f} K  "
              f"visits={st._histogram[i]:>8d}  weight={st.weights[i]:8.3f}")

    engine.finalize_simulation(cfg, ctx, built.topology, start_time)
    print(f"Tempering report: {st_log}")


if __name__ == '__main__':
    main()
