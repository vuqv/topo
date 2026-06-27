"""
Replicate a single-chain coarse-grained system into many **non-interacting**
copies inside one OpenMM ``System``.

Coarse-grained protein models are tiny (a few hundred beads), so a single chain
badly underuses a GPU. Packing N independent copies into one simulation fills the
device and yields N trajectories per run for better sampling. The copies must not
interact: bonded terms are duplicated per copy, and every ``CustomNonbondedForce``
is restricted to intra-copy interactions with ``addInteractionGroup`` so copy *i*
never sees copy *j*.

Typical use::

    import topo
    cg = topo.models.buildCoarseGrainModel("protein.pdb", domain_def="domain.yaml")
    system, topology, positions = topo.make_noninteracting_copies(
        cg.system, cg.topology, cg.positions, n_copies=20)
    # ... build a Simulation from system/topology/positions, run, then split the
    # multi-chain DCD into per-chain trajectories for analysis.

After the run, the atoms of copy *k* are the contiguous block
``[k*n_atoms : (k+1)*n_atoms]`` (see :func:`split_indices`), which makes splitting
the trajectory trivial.
"""
import copy as _copy
from pathlib import Path
from typing import List, Tuple

import numpy as np
import openmm as mm
from openmm import unit


def replicate_system_intra_only(template_system: mm.System, n_copies: int) -> mm.System:
    """
    Replicate ``template_system`` into ``n_copies`` independent copies.

    Particles and constraints are duplicated per copy; each supported force is
    rebuilt with per-copy atom-index offsets. ``CustomNonbondedForce`` objects get
    one interaction group per copy so copies never interact with one another.

    Supported forces: ``HarmonicBondForce``, ``CustomAngleForce``,
    ``PeriodicTorsionForce``, ``CustomTorsionForce``, ``CustomNonbondedForce``.
    ``CMMotionRemover`` is skipped (add a fresh one to the full system if wanted);
    any other force type raises ``NotImplementedError``.
    """
    if n_copies < 1:
        raise ValueError(f"n_copies must be >= 1, got {n_copies}")

    n = template_system.getNumParticles()
    full = mm.System()

    # particles
    for _ in range(n_copies):
        for i in range(n):
            full.addParticle(template_system.getParticleMass(i))

    # constraints
    for c in range(n_copies):
        off = c * n
        for i in range(template_system.getNumConstraints()):
            p1, p2, dist = template_system.getConstraintParameters(i)
            full.addConstraint(p1 + off, p2 + off, dist)

    for force in template_system.getForces():
        if isinstance(force, mm.CMMotionRemover):
            # COM removal is global; the caller adds one to the full system.
            continue

        elif isinstance(force, mm.HarmonicBondForce):
            nf = mm.HarmonicBondForce()
            nf.setName(force.getName())
            nf.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
            nf.setForceGroup(force.getForceGroup())
            for c in range(n_copies):
                off = c * n
                for i in range(force.getNumBonds()):
                    p1, p2, length, k = force.getBondParameters(i)
                    nf.addBond(p1 + off, p2 + off, length, k)
            full.addForce(nf)

        elif isinstance(force, mm.CustomAngleForce):
            nf = mm.CustomAngleForce(force.getEnergyFunction())
            nf.setName(force.getName())
            nf.setForceGroup(force.getForceGroup())
            nf.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
            for i in range(force.getNumPerAngleParameters()):
                nf.addPerAngleParameter(force.getPerAngleParameterName(i))
            for i in range(force.getNumGlobalParameters()):
                nf.addGlobalParameter(force.getGlobalParameterName(i),
                                      force.getGlobalParameterDefaultValue(i))
            for c in range(n_copies):
                off = c * n
                for i in range(force.getNumAngles()):
                    p1, p2, p3, params = force.getAngleParameters(i)
                    nf.addAngle(p1 + off, p2 + off, p3 + off, params)
            full.addForce(nf)

        elif isinstance(force, mm.PeriodicTorsionForce):
            nf = mm.PeriodicTorsionForce()
            nf.setName(force.getName())
            nf.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
            nf.setForceGroup(force.getForceGroup())
            for c in range(n_copies):
                off = c * n
                for i in range(force.getNumTorsions()):
                    p1, p2, p3, p4, period, phase, k = force.getTorsionParameters(i)
                    nf.addTorsion(p1 + off, p2 + off, p3 + off, p4 + off, period, phase, k)
            full.addForce(nf)

        elif isinstance(force, mm.CustomTorsionForce):
            nf = mm.CustomTorsionForce(force.getEnergyFunction())
            nf.setName(force.getName())
            nf.setForceGroup(force.getForceGroup())
            nf.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
            for i in range(force.getNumPerTorsionParameters()):
                nf.addPerTorsionParameter(force.getPerTorsionParameterName(i))
            for i in range(force.getNumGlobalParameters()):
                nf.addGlobalParameter(force.getGlobalParameterName(i),
                                      force.getGlobalParameterDefaultValue(i))
            for c in range(n_copies):
                off = c * n
                for i in range(force.getNumTorsions()):
                    p1, p2, p3, p4, params = force.getTorsionParameters(i)
                    nf.addTorsion(p1 + off, p2 + off, p3 + off, p4 + off, params)
            full.addForce(nf)

        elif isinstance(force, mm.CustomNonbondedForce):
            nf = mm.CustomNonbondedForce(force.getEnergyFunction())
            nf.setName(force.getName())
            nf.setForceGroup(force.getForceGroup())
            nf.setNonbondedMethod(force.getNonbondedMethod())
            nf.setCutoffDistance(force.getCutoffDistance())
            nf.setUseSwitchingFunction(force.getUseSwitchingFunction())
            if force.getUseSwitchingFunction():
                nf.setSwitchingDistance(force.getSwitchingDistance())
            nf.setUseLongRangeCorrection(force.getUseLongRangeCorrection())
            for i in range(force.getNumGlobalParameters()):
                nf.addGlobalParameter(force.getGlobalParameterName(i),
                                      force.getGlobalParameterDefaultValue(i))
            for i in range(force.getNumPerParticleParameters()):
                nf.addPerParticleParameter(force.getPerParticleParameterName(i))
            for i in range(force.getNumTabulatedFunctions()):
                nf.addTabulatedFunction(force.getTabulatedFunctionName(i),
                                        _copy.deepcopy(force.getTabulatedFunction(i)))
            # per-particle parameters: repeat the template block for each copy
            for c in range(n_copies):
                for i in range(n):
                    nf.addParticle(force.getParticleParameters(i))
            # exclusions, offset per copy
            for c in range(n_copies):
                off = c * n
                for i in range(force.getNumExclusions()):
                    p1, p2 = force.getExclusionParticles(i)
                    nf.addExclusion(p1 + off, p2 + off)
            # intra-copy only: copy k interacts solely with itself
            for c in range(n_copies):
                group = list(range(c * n, (c + 1) * n))
                nf.addInteractionGroup(group, group)
            full.addForce(nf)

        else:
            raise NotImplementedError(
                f"replicate_system_intra_only does not handle force type "
                f"{type(force).__name__} (name '{force.getName()}'). Add a branch "
                f"for it or strip it from the system before replicating."
            )

    return full


def replicate_topology(template_topology: mm.app.Topology, n_copies: int) -> mm.app.Topology:
    """Replicate an OpenMM ``Topology`` into ``n_copies`` chains."""
    full = mm.app.Topology()
    box = template_topology.getPeriodicBoxVectors()
    if box is not None:
        full.setPeriodicBoxVectors(box)
    for _ in range(n_copies):
        atom_map = {}
        for chain in template_topology.chains():
            new_chain = full.addChain(chain.id)
            for residue in chain.residues():
                new_res = full.addResidue(residue.name, new_chain, id=residue.id,
                                          insertionCode=residue.insertionCode)
                for atom in residue.atoms():
                    atom_map[atom] = full.addAtom(atom.name, atom.element, new_res,
                                                  id=atom.id, formalCharge=atom.formalCharge)
        for bond in template_topology.bonds():
            full.addBond(atom_map[bond.atom1], atom_map[bond.atom2],
                         type=bond.type, order=bond.order)
    return full


def replicate_positions(template_positions, n_copies: int,
                        shift=5.0 * unit.nanometer):
    """
    Replicate positions ``n_copies`` times, translating each copy along x by
    ``shift`` so the chains start in separate regions of space (they do not
    interact, but separation keeps the start clean and aids visualization).
    """
    pos_unit = template_positions.unit
    pos_nm = np.array(template_positions.value_in_unit(unit.nanometer), copy=True)
    shift_nm = shift.value_in_unit(unit.nanometer)
    blocks = []
    for c in range(n_copies):
        b = np.array(pos_nm, copy=True)
        b[:, 0] += c * shift_nm
        blocks.append(b)
    full_nm = np.vstack(blocks)
    return (unit.Quantity(full_nm, unit.nanometer)).value_in_unit(pos_unit) * pos_unit


def split_indices(n_atoms_single: int, n_copies: int) -> List[Tuple[int, int]]:
    """Return [(start, stop), ...] atom-index ranges, one per copy."""
    return [(k * n_atoms_single, (k + 1) * n_atoms_single) for k in range(n_copies)]


def make_noninteracting_copies(system: mm.System, topology: mm.app.Topology,
                               positions, n_copies: int,
                               shift=5.0 * unit.nanometer):
    """
    Convenience wrapper: replicate ``system``, ``topology`` and ``positions``
    into ``n_copies`` non-interacting copies in one call.

    Returns
    -------
    (full_system, full_topology, full_positions)
    """
    full_system = replicate_system_intra_only(system, n_copies)
    full_topology = replicate_topology(topology, n_copies)
    full_positions = replicate_positions(positions, n_copies, shift=shift)
    return full_system, full_topology, full_positions


# --------------------------------------------------------------------------- #
# The inverse operation: split a combined multi-copy trajectory back into one
# DCD per copy. Streaming + chunked so arbitrarily large combined DCDs (that do
# not fit in memory) are handled in a single pass.
# --------------------------------------------------------------------------- #
def split_chains(combined_dcd, output_paths, *, chunk_size: int = 1000,
                         center: bool = True):
    """Split a combined multi-copy DCD into one DCD per copy.

    The combined trajectory's atoms are assumed to be ``n_copies`` contiguous
    equal-size blocks (the layout produced by :func:`make_noninteracting_copies`
    / :func:`split_indices`), so copy *k* is atoms ``[k*m : (k+1)*m]``.

    **Memory-bounded streaming.** The combined DCD is read in chunks of
    ``chunk_size`` frames and each chunk is written straight to the per-copy
    files, so peak memory is ``chunk_size * n_atoms * 3`` floats regardless of how
    long the trajectory is — a combined DCD too large to fit in memory is handled
    fine. A single pass over the input writes all copies.

    Coordinates stay in Å throughout (mdtraj's low-level DCD read **and** write
    are both in Å), so no unit conversion is needed or applied.

    Parameters
    ----------
    combined_dcd : str or Path
        The combined multi-copy DCD.
    output_paths : sequence of str or Path
        One output DCD path per copy; ``n_copies = len(output_paths)``. Parent
        directories are created as needed.
    chunk_size : int, optional
        Frames read/written per chunk (default 1000).
    center : bool, optional
        If True (default), translate each copy to its own centre of geometry per
        frame, so every per-copy trajectory is centred at the origin and the
        inter-copy ``copy_shift`` offset is removed. Set False to preserve the raw
        coordinates exactly.

    Returns
    -------
    (n_copies, atoms_per_copy, n_frames) : tuple[int, int, int]
    """
    from mdtraj.formats import DCDTrajectoryFile

    output_paths = [str(p) for p in output_paths]
    n_copies = len(output_paths)
    if n_copies < 1:
        raise ValueError("output_paths must list at least one copy")
    for p in output_paths:
        Path(p).parent.mkdir(parents=True, exist_ok=True)

    writers = [DCDTrajectoryFile(p, "w") for p in output_paths]
    atoms_per_copy = None
    n_frames = 0
    try:
        with DCDTrajectoryFile(str(combined_dcd), "r") as reader:
            while True:
                xyz, _cell_len, _cell_ang = reader.read(n_frames=chunk_size)
                if xyz.shape[0] == 0:
                    break
                if atoms_per_copy is None:
                    total_atoms = xyz.shape[1]
                    if total_atoms % n_copies != 0:
                        raise ValueError(
                            f"combined DCD has {total_atoms} atoms, not divisible "
                            f"by n_copies={n_copies}")
                    atoms_per_copy = total_atoms // n_copies
                for k in range(n_copies):
                    sub = xyz[:, k * atoms_per_copy:(k + 1) * atoms_per_copy, :]
                    if center:
                        sub = sub - sub.mean(axis=1, keepdims=True)
                    writers[k].write(sub)
                n_frames += xyz.shape[0]
    finally:
        for w in writers:
            w.close()
    return n_copies, atoms_per_copy, n_frames


def _split_cli():
    """CLI: ``python -m topo.utils.multichain`` — split a combined multi-copy DCD."""
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m topo.utils.multichain",
        description="Split a combined multi-copy DCD into one DCD per copy "
                    "(streaming; handles DCDs too large to fit in memory).")
    p.add_argument("-f", "--combined-dcd", required=True,
                   help="Combined multi-copy DCD.")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("-n", "--n-copies", type=int, help="Number of copies.")
    group.add_argument("-tp", "--template-psf",
                       help="Single-chain PSF; n_copies is inferred from the "
                            "combined DCD's atom count.")
    p.add_argument("-o", "--output-dir", default="split_trajs",
                   help="Output directory (default: split_trajs).")
    p.add_argument("--outname", default="traj",
                   help="Output basename -> <outname>_<k>.dcd (default: traj).")
    p.add_argument("--chunk-size", type=int, default=1000,
                   help="Frames per chunk (default: 1000).")
    p.add_argument("--no-center", action="store_true",
                   help="Preserve raw coordinates (default: centre each copy per frame).")
    args = p.parse_args()

    if args.n_copies is not None:
        n_copies = args.n_copies
    else:
        from openmm.app import CharmmPsfFile
        from mdtraj.formats import DCDTrajectoryFile
        atoms_per_copy = CharmmPsfFile(args.template_psf).topology.getNumAtoms()
        with DCDTrajectoryFile(args.combined_dcd, "r") as reader:
            xyz, _, _ = reader.read(n_frames=1)
        total = xyz.shape[1]
        if total % atoms_per_copy != 0:
            raise SystemExit(f"DCD atoms ({total}) not divisible by template "
                             f"atoms ({atoms_per_copy}).")
        n_copies = total // atoms_per_copy

    out_dir = Path(args.output_dir)
    out_paths = [out_dir / f"{args.outname}_{k}.dcd" for k in range(n_copies)]
    nc, apc, nf = split_chains(args.combined_dcd, out_paths,
                                           chunk_size=args.chunk_size,
                                           center=not args.no_center)
    print(f"Split {nf} frames into {nc} copies of {apc} atoms each -> {out_dir}/")


if __name__ == "__main__":
    _split_cli()
