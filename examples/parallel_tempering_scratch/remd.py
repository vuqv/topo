#!/usr/bin/env python3
"""
remd.py -- serial parallel-tempering (temperature) replica exchange MD.

A standalone, hackable REMD driver built on TOPO's public building blocks
(``topo.engine`` + the ``md.ini`` config), in the same spirit as
``examples/custom_md.py``. It is a *script*, not a package module: REMD is an
occasional experiment here, so everything lives in one readable file you can copy
and edit.

What it does
------------
* Builds the coarse-grained system ONCE and shares it across N OpenMM Contexts
  (one per replica) -- the force field is identical between replicas; only the
  thermostat temperature differs.
* Lays out a GEOMETRIC temperature ladder between ``--tmin`` and ``--tmax``,
  which gives roughly uniform exchange acceptance when the heat capacity is
  roughly constant.
* Runs SERIALLY on a single device: each cycle steps every replica in turn for
  ``--exchange-interval`` steps, then attempts nearest-neighbour swaps.
* Uses the standard Sugita-Okamoto Metropolis criterion. On an accepted swap it
  exchanges the two replicas' COORDINATES and rescales their velocities by
  ``sqrt(T_new / T_old)``. Each Context therefore stays at a FIXED temperature,
  so ``<outname>_rep00.dcd`` is the T=tmin canonical ensemble, ``_rep01`` the
  next rung, and so on -- ready for analysis without demultiplexing.

Outputs (under the md.ini ``output_dir``)
    <outname>_repNN.dcd / .log   one trajectory + log per temperature rung
    <outname>_remd.log           exchange attempts and per-pair acceptance ratios

Run it like the real runner, plus the REMD ladder options::

    python examples/parallel_tempering_scratch/remd.py -f md.ini --tmin 300 --tmax 500 --nreplicas 8 \
        --exchange-interval 1000

The md.ini supplies the system (``pdb_file``, ``domain_def``), timestep, friction
(``tau_t``), output cadence (``nstxout``/``nstlog``), device, etc. ``ref_t`` from
the md.ini is ignored -- the ladder sets each replica's temperature. Total steps
per replica = ``--cycles * --exchange-interval``; if ``--cycles`` is omitted it is
derived from the md.ini ``md_steps``.

Not implemented (keep it simple): restart/checkpointing and multi-GPU/MPI. For a
long production run, raise ``--cycles`` and copy this file to add what you need.
"""
import argparse
import math
import random
import sys
import time

import openmm as mm
from openmm import unit

from topo import engine
from topo.utils.config import read_simulation_config
import topo

# Molar Boltzmann constant, so beta = 1 / (KB_NA * T) is in mol/kJ and pairs with
# potential energies read out in kJ/mol.
KB_NA = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA


def geometric_ladder(t_min, t_max, n):
    """Return ``n`` temperatures (Kelvin, plain floats) spaced geometrically."""
    if n < 2:
        raise SystemExit("need at least 2 replicas for REMD")
    ratio = (t_max / t_min) ** (1.0 / (n - 1))
    return [t_min * ratio ** i for i in range(n)]


def make_replica(cfg, built, temperature, tau_t):
    """Create one OpenMM ``Simulation`` at ``temperature`` (Kelvin float).

    All replicas share ``built.system`` / ``built.topology`` (the force field is
    identical); each gets its own Langevin integrator and Context.
    """
    integrator = mm.LangevinIntegrator(temperature * unit.kelvin, tau_t, cfg.dt)
    integrator.setConstraintTolerance(cfg.constraint_tolerance)
    platform, properties = cfg.make_platform()
    return mm.app.Simulation(built.topology, built.system, integrator,
                             platform, properties)


def attach_replica_reporters(cfg, sim, prefix, total_steps):
    """Attach DCD + aligned .log reporters writing to ``<prefix>.dcd`` / ``.log``.

    Mirrors ``engine.attach_reporters`` but targets a per-replica path prefix
    instead of ``cfg.outname``. No checkpoint reporter (restart is out of scope).
    """
    sim.reporters = []
    sim.reporters.append(
        mm.app.DCDReporter(prefix + '.dcd', cfg.nstxout,
                           enforcePeriodicBox=bool(cfg.pbc)))
    sim.reporters.append(
        topo.topoReporter(prefix + '.log', cfg.nstlog,
                          precision=cfg.log_precision, width=cfg.log_width,
                          step=True, time=True,
                          potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                          temperature=True, remainingTime=True, speed=True,
                          totalSteps=total_steps, separator='  '))


def potential_energy(sim):
    """Potential energy of a replica in kJ/mol (plain float)."""
    state = sim.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def swap_configurations(sim_a, temp_a, sim_b, temp_b, pbc):
    """Exchange coordinates of two replicas, rescaling velocities to each T.

    The config leaving replica a (temperature ``temp_a``) enters replica b
    (``temp_b``), so its velocities are scaled by ``sqrt(temp_b / temp_a)`` to
    match b's thermostat, and vice versa.
    """
    state_a = sim_a.context.getState(getPositions=True, getVelocities=True,
                                     enforcePeriodicBox=pbc)
    state_b = sim_b.context.getState(getPositions=True, getVelocities=True,
                                     enforcePeriodicBox=pbc)
    scale_a2b = (temp_b / temp_a) ** 0.5
    scale_b2a = (temp_a / temp_b) ** 0.5
    # a receives b's config; b receives a's config.
    sim_a.context.setPositions(state_b.getPositions())
    sim_a.context.setVelocities(state_b.getVelocities() * scale_b2a)
    sim_b.context.setPositions(state_a.getPositions())
    sim_b.context.setVelocities(state_a.getVelocities() * scale_a2b)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="remd.py",
        description="Serial temperature replica-exchange MD on TOPO CG models.")
    parser.add_argument('-input', '-f', required=True,
                        help='simulation config file (md.ini) -- supplies the system')
    parser.add_argument('--tmin', type=float, required=True,
                        help='lowest replica temperature (Kelvin)')
    parser.add_argument('--tmax', type=float, required=True,
                        help='highest replica temperature (Kelvin)')
    parser.add_argument('--nreplicas', type=int, required=True,
                        help='number of replicas / temperature rungs')
    parser.add_argument('--exchange-interval', type=int, default=1000,
                        help='MD steps between exchange attempts (default 1000)')
    parser.add_argument('--cycles', type=int, default=None,
                        help='number of exchange cycles (default: md_steps // '
                             'exchange-interval from the md.ini)')
    parser.add_argument('--seed', type=int, default=None,
                        help='RNG seed for the Metropolis accept/reject stream')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    print(f"OpenMM version: {mm.__version__}")

    # --- Base system from the md.ini (ref_t is overridden by the ladder) -----
    cfg = read_simulation_config(args.input)
    cfg.prepare_output_dir()
    # tcoupl may be off in the md.ini (tau_t is then None); REMD always needs a
    # thermostat, so fall back to a sensible Langevin friction.
    tau_t = cfg.tau_t if cfg.tau_t is not None else 0.05 / unit.picoseconds

    built = engine.build_system(cfg)

    # --- Temperature ladder --------------------------------------------------
    temps = geometric_ladder(args.tmin, args.tmax, args.nreplicas)
    betas = [1.0 / (KB_NA * (t * unit.kelvin)).value_in_unit(unit.kilojoule_per_mole)
             for t in temps]
    cycles = args.cycles if args.cycles is not None else max(1, cfg.md_steps // args.exchange_interval)
    total_steps = cycles * args.exchange_interval
    print("Temperature ladder (K): " + ", ".join(f"{t:.1f}" for t in temps))
    print(f"{cycles} cycles x {args.exchange_interval} steps = {total_steps} steps/replica")

    # --- Build replicas (shared system, one Context each) --------------------
    base = cfg.output_path('')
    replicas = []
    for i, t in enumerate(temps):
        sim = make_replica(cfg, built, t, tau_t)
        sim.context.setPositions(built.positions)
        replicas.append(sim)

    # Optional single minimization from the native structure, shared by all
    # replicas (they start from identical coordinates).
    if cfg.minimize:
        print("Minimizing starting structure...")
        replicas[0].minimizeEnergy()
        min_pos = replicas[0].context.getState(getPositions=True).getPositions()
        for sim in replicas:
            sim.context.setPositions(min_pos)

    for i, (sim, t) in enumerate(zip(replicas, temps)):
        sim.context.setVelocitiesToTemperature(t * unit.kelvin)
        attach_replica_reporters(cfg, sim, f"{base}_rep{i:02d}", total_steps)

    # --- Exchange bookkeeping: one neighbour pair (i, i+1) per index i -------
    n_pairs = args.nreplicas - 1
    attempts = [0] * n_pairs
    accepts = [0] * n_pairs

    remd_log = f"{base}_remd.log"
    with open(remd_log, 'w') as fh:
        fh.write("# serial temperature-REMD exchange log\n")
        fh.write("# ladder (K): " + ", ".join(f"{t:.2f}" for t in temps) + "\n")
        fh.write("# columns: cycle  pair(i,i+1)  dU_kJ/mol  accepted\n")

    print("REMD started")
    start_time = time.time()

    for cycle in range(cycles):
        # 1) Propagate every replica independently.
        for sim in replicas:
            sim.step(args.exchange_interval)

        # 2) Attempt neighbour swaps. Alternate the even/odd pair set each cycle
        #    so every neighbour pair is tried on average every other cycle and no
        #    replica appears in two pairs within one cycle.
        offset = cycle % 2
        energies = [potential_energy(sim) for sim in replicas]
        cycle_lines = []
        for i in range(offset, n_pairs, 2):
            j = i + 1
            # Sugita-Okamoto: swapping configs x_i<->x_j is accepted with
            # probability min(1, exp(delta)), delta = (beta_i - beta_j)(U_i - U_j).
            delta = (betas[i] - betas[j]) * (energies[i] - energies[j])
            accepted = delta >= 0.0 or rng.random() < math.exp(delta)
            attempts[i] += 1
            if accepted:
                accepts[i] += 1
                swap_configurations(replicas[i], temps[i],
                                    replicas[j], temps[j], bool(cfg.pbc))
                # Keep our energy cache consistent for any later pair this cycle
                # (disjoint pairs, so this is just tidy bookkeeping).
                energies[i], energies[j] = energies[j], energies[i]
            cycle_lines.append(
                f"{cycle:6d}  ({i:2d},{j:2d})  {energies[j] - energies[i]:12.3f}  "
                f"{int(accepted)}")

        with open(remd_log, 'a') as fh:
            fh.write("\n".join(cycle_lines) + "\n")

        if (cycle + 1) % max(1, cycles // 20) == 0 or cycle == cycles - 1:
            done = (cycle + 1) / cycles
            print(f"  cycle {cycle + 1}/{cycles} ({100 * done:4.1f}%)")

    # --- Per-pair acceptance summary ----------------------------------------
    print("Exchange acceptance ratios (pair -> ratio):")
    summary = []
    for i in range(n_pairs):
        ratio = accepts[i] / attempts[i] if attempts[i] else 0.0
        line = f"  ({i:2d},{i+1:2d})  {temps[i]:6.1f}<->{temps[i+1]:6.1f} K  {ratio:5.3f}  ({accepts[i]}/{attempts[i]})"
        print(line)
        summary.append(line)
    with open(remd_log, 'a') as fh:
        fh.write("\n# per-pair acceptance ratios\n")
        fh.write("\n".join(summary) + "\n")

    # Save final coordinates per replica so a follow-up run could seed from them.
    for i, sim in enumerate(replicas):
        pos = sim.context.getState(getPositions=True,
                                   enforcePeriodicBox=bool(cfg.pbc)).getPositions()
        with open(f"{base}_rep{i:02d}_final.pdb", 'w') as fh:
            mm.app.PDBFile.writeFile(built.topology, pos, fh)

    print(f"--- Finished {cycles} cycles in {time.time() - start_time:.1f} s ---")
    print(f"Exchange log: {remd_log}")


if __name__ == '__main__':
    main()
