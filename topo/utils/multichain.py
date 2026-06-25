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
