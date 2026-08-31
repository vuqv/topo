"""Tests for the two-force IDR architecture (openspec change ``idr-ashbaugh-hatch``).

Where ``test_idr.py`` covers the *matrices* (radii and well depths), these cover the
*forces*: which OpenMM force evaluates which pair, and what its potential looks like.

Three things are being defended, and all three failure modes are **silent** — the
simulation runs to completion and produces plausible-looking numbers:

1. **The potential's shape.** The Ashbaugh–Hatch split has to put its minimum at
   ``R_ij`` with depth exactly ``-eps_ij``, be C¹ there, vanish beyond ``R_ij`` at
   ``eps_ij = 0``, and keep its core fixed as the well depth varies. The expression is
   extracted from the *shipped source* rather than restated here, so a test cannot pass
   against a formula the model does not actually use.
2. **The partition.** OpenMM takes the **union** of a force's interaction groups, so any
   code that widens a force it did not create silently re-admits pairs that force no
   longer owns — and they are then evaluated by *both* contact forces and summed. No
   exception is raised.
3. **The no-disorder path.** A build with no ``disordered:`` section must take the
   single-force route unchanged.

Fixture: ``tutorials/A1_single_domain_quickstart/P0CX28_clean.pdb`` (106 residues,
resids 1..106, so residue number == 0-based index + 1). STRIDE must be on PATH
(vendored under ``assets/``); the build is deterministic.
"""
import inspect
import re
import textwrap
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("openmm")
import openmm as mm
from openmm import unit

REPO_ROOT = Path(__file__).resolve().parents[1]
PDB = REPO_ROOT / "tutorials" / "A1_single_domain_quickstart" / "P0CX28_clean.pdb"

pytestmark = pytest.mark.skipif(not PDB.exists(), reason=f"fixture PDB missing: {PDB}")

N_RES = 106
MASK_END = 20                       # disorder residues 1..20 (0-based 0..19)
KT_310 = 8.314462618e-3 * 310.0     # kJ/mol


def _write_yaml(path, body):
    path.write_text(textwrap.dedent(body))
    return str(path)


def _idr_yaml(tmp_path, name="idr.yaml", end=MASK_END, **kv):
    """A domain_def whose IDR block is residues 1..end, plus any extra keys."""
    extra = "".join(f"          {k}: {v}\n" for k, v in kv.items())
    return _write_yaml(tmp_path / name, f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{end}]
{extra}""")


def _build(pdb, domain_def=None):
    """Build a coarse-grained model exactly as production does."""
    import topo
    return topo.models.buildCoarseGrainModel(
        str(pdb), domain_def=domain_def, check_forces=False)


def _groups(force):
    """The force's interaction groups as a list of (set, set) index pairs."""
    return [tuple(map(set, force.getInteractionGroupParameters(g)))
            for g in range(force.getNumInteractionGroups())]


def _pairs_covered(force, n):
    """Every unordered pair {i, j} the force's interaction groups admit.

    An unrestricted force (no groups) admits every pair — which is what makes a
    forgotten group so dangerous, and why this returns the full set for one.
    """
    groups = _groups(force)
    if not groups:
        return {(i, j) for i in range(n) for j in range(i + 1, n)}
    out = set()
    for s1, s2 in groups:
        for i in s1:
            for j in s2:
                if i != j:
                    out.add((min(i, j), max(i, j)))
    return out


def _exclusions(force):
    return {tuple(sorted(force.getExclusionParticles(k)))
            for k in range(force.getNumExclusions())}


# --------------------------------------------------------------------------- #
# 6.1 Potential shape — evaluated from the SHIPPED expression
# --------------------------------------------------------------------------- #
def _shipped_idr_expression():
    """The energy expression as it is actually written in the shipped source.

    Parsed out of ``addIDRNonBondedForce`` rather than restated, so this test cannot
    silently drift from (or vacuously agree with) the model that is really built.
    """
    from topo.core.system import system as TopoSystem
    src = inspect.getsource(TopoSystem.addIDRNonBondedForce)
    parts = re.findall(r"energy_function \+?= '([^']*)'", src)
    assert parts, "could not extract the IDR energy expression from the shipped source"
    return "".join(parts)


def _ah_context(expr, R, eps, eps_ev):
    """A 2-particle OpenMM context evaluating ``expr`` for one pair."""
    f = mm.CustomNonbondedForce(expr)
    f.addGlobalParameter("epsEV", eps_ev)
    f.addGlobalParameter("rmin_scale", 2.0 ** (1.0 / 6.0))
    f.addTabulatedFunction("eps_table", mm.Discrete2DFunction(2, 2, [eps] * 4))
    f.addTabulatedFunction("R_table", mm.Discrete2DFunction(2, 2, [R] * 4))
    f.addPerParticleParameter("id")
    f.addParticle((0,))
    f.addParticle((1,))
    f.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    s = mm.System()
    s.addParticle(1.0)
    s.addParticle(1.0)
    s.addForce(f)
    return mm.Context(s, mm.VerletIntegrator(1e-6),
                      mm.Platform.getPlatform("Reference"))


def _u_pair(expr, r_nm, R, eps, eps_ev):
    """U(r) for one pair, in kJ/mol."""
    ctx = _ah_context(expr, R, eps, eps_ev)
    out = []
    for r in np.atleast_1d(r_nm):
        ctx.setPositions([[0, 0, 0], [float(r), 0, 0]] * unit.nanometer)
        out.append(ctx.getState(getEnergy=True).getPotentialEnergy()
                   .value_in_unit(unit.kilojoule_per_mole))
    return np.array(out)


def _dudr(expr, r, R, eps, eps_ev):
    """dU/dr from OpenMM's analytic gradient (not a finite difference)."""
    ctx = _ah_context(expr, R, eps, eps_ev)
    ctx.setPositions([[0, 0, 0], [float(r), 0, 0]] * unit.nanometer)
    fx = ctx.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer)[1][0]
    return -fx


def test_idr_potential_minimum_and_depth():
    """Minimum at ``r = R_ij`` with depth exactly ``-eps_ij``, for any ``eps_EV``.

    This is what makes the split a *decoupling*: the well depth is the parameter
    ``idr_scale * eps_BT`` and nothing else, whatever the core strength.
    """
    expr = _shipped_idr_expression()
    R = 0.60
    for eps_ev in (0.2, 0.8368, 3.0):
        for eps in (0.0, 0.25, 1.0, 2.0):
            rs = np.linspace(0.85 * R, 1.25 * R, 4001)
            u = _u_pair(expr, rs, R, eps, eps_ev)
            assert abs(rs[np.argmin(u)] - R) < 2e-4
            assert _u_pair(expr, [R], R, eps, eps_ev)[0] == pytest.approx(-eps, abs=1e-9)


def test_idr_potential_c1_at_rmin():
    """C¹ at the join: both one-sided analytic derivatives go to zero.

    A discontinuous force at ``r = R`` would be an energy leak the integrator turns
    into heat, so the join is not cosmetic.
    """
    expr = _shipped_idr_expression()
    R, eps_ev = 0.60, 0.8368
    for eps in (0.0, 0.5, 2.0):
        rows = [(d, _dudr(expr, R - d, R, eps, eps_ev),
                 _dudr(expr, R + d, R, eps, eps_ev))
                for d in (1e-3, 1e-4, 1e-5)]
        # Both sides are O(|r - R|) about a smooth minimum, so a 10x smaller offset
        # must shrink each by ~10x -- i.e. both one-sided limits are 0, and equal.
        for (_d0, l0, r0), (_d1, l1, r1) in zip(rows, rows[1:]):
            for a, b in ((l0, l1), (r0, r1)):
                if abs(a) > 1e-12:
                    assert 5 < abs(a) / max(abs(b), 1e-30) < 20


def test_idr_potential_zero_attraction_is_clean_wca():
    """``idr_scale = 0`` gives ``U ≡ 0`` beyond ``R`` and a repulsive core below it.

    Under the old coupled 12-10-6 "no attraction" also meant "no excluded volume":
    the bead went thin because it had no energy. Here the two are independent.
    """
    expr = _shipped_idr_expression()
    R, eps_ev = 0.60, 0.8368
    beyond = np.linspace(R, 3 * R, 200)
    assert np.allclose(_u_pair(expr, beyond, R, 0.0, eps_ev), 0.0, atol=1e-12)
    inside = np.linspace(0.7 * R, 0.999 * R, 60)
    u = _u_pair(expr, inside, R, 0.0, eps_ev)
    assert np.all(u > 0)
    assert np.all(np.diff(u) < 0)          # monotonically repulsive


def test_idr_potential_reduces_to_plain_12_6():
    """``eps_ij = eps_EV`` recovers the plain 12-6 with ``sigma = 2^(-1/6) R``."""
    expr = _shipped_idr_expression()
    R, eps_ev = 0.60, 0.8368
    rs = np.linspace(0.8 * R, 2.5 * R, 120)
    sigma = R * 2 ** (-1 / 6)
    lj = 4 * eps_ev * ((sigma / rs) ** 12 - (sigma / rs) ** 6)
    assert np.allclose(_u_pair(expr, rs, R, eps_ev, eps_ev), lj, atol=1e-9)


def test_idr_core_is_set_by_eps_ev_alone():
    """The core (``U = kT``) moves < 5 % as the well depth runs 0 -> 2 kJ/mol.

    The whole point of the change. Under the coupled 12-10-6 the same sweep moved the
    core ~56 % (~3.7x in excluded volume), which is why adding attraction made the
    chain *expand*.
    """
    expr = _shipped_idr_expression()
    R, eps_ev = 0.60, 0.8368
    cores = []
    for eps in (0.0, 0.5, 1.0, 2.0):
        rs = np.linspace(0.5 * R, R, 20001)
        u = _u_pair(expr, rs, R, eps, eps_ev)
        cores.append(rs[np.argmin(np.abs(u - KT_310))] / R)
    assert (max(cores) - min(cores)) / max(cores) < 0.05
    # ...and the residual runs the benign way: more attraction *softens* the core.
    assert cores == sorted(cores, reverse=True)


# --------------------------------------------------------------------------- #
# 6.2 The force partition
# --------------------------------------------------------------------------- #
def test_no_disorder_builds_one_unrestricted_force(tmp_path, monkeypatch):
    """No ``disordered:`` section -> the pre-IDR path: one force, no interaction group."""
    monkeypatch.chdir(tmp_path)
    m = _build(PDB)
    assert m.idr_non_bonded_force is None
    assert m.custom_non_bonded_force.getNumInteractionGroups() == 0
    assert "IDR Non-Bonded Energy" not in m.forceGroups


def test_forces_partition_the_pairs(tmp_path, monkeypatch):
    """The two forces' domains are disjoint, and together they cover every pair.

    A pair in both would get the Ashbaugh–Hatch well *and* the coupled 12-10-6 summed
    on top of it — reinstating exactly the excluded-volume coupling this change exists
    to remove — with no error raised.
    """
    monkeypatch.chdir(tmp_path)
    m = _build(PDB, _idr_yaml(tmp_path))

    cf, idf = m.custom_non_bonded_force, m.idr_non_bonded_force
    assert idf is not None
    n = m.n_atoms
    idr = set(range(MASK_END))

    go_pairs = _pairs_covered(cf, n)
    ah_pairs = _pairs_covered(idf, n)

    assert not (go_pairs & ah_pairs), "a pair is evaluated by BOTH contact forces"
    all_pairs = {(i, j) for i in range(n) for j in range(i + 1, n)}
    assert go_pairs | ah_pairs == all_pairs, "some pair is evaluated by neither force"

    # The 12-10-6 must not see a disordered bead at all.
    for s1, s2 in _groups(cf):
        assert not (s1 & idr) and not (s2 & idr)
    # The AH force sees exactly the IDR-involving pairs, and nothing else.
    assert all(i in idr or j in idr for i, j in ah_pairs)


def test_forces_share_particles_and_exclusions(tmp_path, monkeypatch):
    """Both forces hold one parameter per particle and an identical exclusion list.

    OpenMM's CPU platform requires every ``CustomNonbondedForce`` in a System to agree
    on exclusions; a mismatch is a hard error at context creation, so this also guards
    the build from failing only at run time.
    """
    monkeypatch.chdir(tmp_path)
    m = _build(PDB, _idr_yaml(tmp_path))
    cf, idf = m.custom_non_bonded_force, m.idr_non_bonded_force
    assert cf.getNumParticles() == idf.getNumParticles() == m.n_atoms
    assert _exclusions(cf) == _exclusions(idf)


def test_eps_ev_reaches_the_force(tmp_path, monkeypatch):
    """``eps_ev_kj`` is honoured, and defaults to the calibrated value when omitted."""
    from topo.utils.nonbonded import DEFAULT_EPS_EV_KJ

    monkeypatch.chdir(tmp_path)

    def eps_ev_of(model):
        f = model.idr_non_bonded_force
        names = [f.getGlobalParameterName(i) for i in range(f.getNumGlobalParameters())]
        return f.getGlobalParameterDefaultValue(names.index("epsEV"))

    assert eps_ev_of(_build(PDB, _idr_yaml(tmp_path, "d.yaml"))) == DEFAULT_EPS_EV_KJ
    assert eps_ev_of(_build(PDB, _idr_yaml(tmp_path, "e.yaml", eps_ev_kj=2.5))) == 2.5


def test_folded_folded_energy_is_unchanged_by_declaring_an_idr(tmp_path, monkeypatch):
    """Folded–folded well positions and depths survive a ``disordered:`` section intact.

    Note this is about the *parameters*, not the total energy: declaring an IDR also
    deletes every native contact that crosses the boundary (a disordered residue's
    crystal contacts are artifacts), so the folded remainder is a less stabilized fold.
    """
    from topo.utils.nonbonded import build_nonbonded_interaction

    monkeypatch.chdir(tmp_path)
    rf, ef, r2f = build_nonbonded_interaction(str(PDB), return_rmin_2=True)
    rm, em, r2m, idr = build_nonbonded_interaction(
        str(PDB), _idr_yaml(tmp_path), return_rmin_2=True, return_idr=True)

    mask = np.zeros(N_RES, dtype=bool)
    mask[list(idr["idx"])] = True
    ff = ~(mask[:, None] | mask[None, :])
    assert np.array_equal(rm[ff], rf[ff])          # exact, not approximate
    assert np.array_equal(em[ff], ef[ff])
    assert np.array_equal(r2m[~mask], r2f[~mask])


def test_fully_disordered_chain_builds(tmp_path, monkeypatch):
    """Every residue disordered: the 12-10-6 gets an empty group and contributes zero."""
    monkeypatch.chdir(tmp_path)
    m = _build(PDB, _idr_yaml(tmp_path, "all.yaml", end=N_RES))
    cf, idf = m.custom_non_bonded_force, m.idr_non_bonded_force
    assert _pairs_covered(cf, m.n_atoms) == set()
    assert len(_groups(idf)) == 1                  # {idr}x{idr} only; no folded set
    ctx = mm.Context(m.system, mm.VerletIntegrator(0.002 * unit.picoseconds),
                     mm.Platform.getPlatform("Reference"))
    ctx.setPositions(m.positions)
    fg = list(m.forceGroups).index("Custom Non-Bonded Energy")
    e = ctx.getState(getEnergy=True, groups={fg}).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    assert e == 0.0


# --------------------------------------------------------------------------- #
# 6.4 Multichain replication preserves the partition
# --------------------------------------------------------------------------- #
def _energy(system, positions):
    ctx = mm.Context(system, mm.VerletIntegrator(0.002 * unit.picoseconds),
                     mm.Platform.getPlatform("Reference"))
    ctx.setPositions(positions)
    return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)


@pytest.mark.parametrize("domain", ["idr", "none"])
def test_replicated_energy_is_exactly_n_times_one_copy(tmp_path, monkeypatch, domain):
    """N non-interacting copies have exactly N x the single-copy energy.

    The cheapest check that no pair is double-counted: any force whose interaction
    groups were widened during replication re-admits pairs it does not own, and the
    identity breaks. Run for both the IDR build and the no-disorder build, so the
    fallback path is covered too.
    """
    from topo.utils.multichain import make_noninteracting_copies

    monkeypatch.chdir(tmp_path)
    dm = _idr_yaml(tmp_path) if domain == "idr" else None
    m = _build(PDB, dm)

    n_copies = 3
    full_sys, _top, full_pos = make_noninteracting_copies(
        m.system, m.topology, m.positions, n_copies, shift=50.0 * unit.nanometer)

    e1 = _energy(m.system, m.positions)
    en = _energy(full_sys, full_pos)
    assert en == pytest.approx(n_copies * e1, rel=1e-9)


def test_replication_preserves_the_partition(tmp_path, monkeypatch):
    """Each force's template groups reappear per copy, offset — still disjoint, no leak."""
    from topo.utils.multichain import replicate_system_intra_only

    monkeypatch.chdir(tmp_path)
    m = _build(PDB, _idr_yaml(tmp_path))
    n, n_copies = m.n_atoms, 3
    full = replicate_system_intra_only(m.system, n_copies)

    nbs = [f for f in full.getForces() if isinstance(f, mm.CustomNonbondedForce)]
    for f in nbs:
        assert f.getNumInteractionGroups() > 0, "replication left a force unrestricted"
        # no interaction group may span two copies
        assert all(i // n == j // n for i, j in _pairs_covered(f, n * n_copies))

    # The two contact forces stay disjoint after replication. (The Yukawa force is
    # also a CustomNonbondedForce and legitimately overlaps both; identify the contact
    # pair by their tabulated functions.)
    contact = [_pairs_covered(f, n * n_copies) for f in nbs
               if f.getNumTabulatedFunctions() > 0]
    assert len(contact) == 2
    assert not (contact[0] & contact[1])


def test_replication_falls_back_for_an_unclaimed_force(tmp_path, monkeypatch):
    """A force with no groups still gets the blanket ``{copy} x {copy}`` per copy."""
    from topo.utils.multichain import replicate_system_intra_only

    monkeypatch.chdir(tmp_path)
    m = _build(PDB)                                   # no disorder -> unrestricted
    assert m.custom_non_bonded_force.getNumInteractionGroups() == 0
    n, n_copies = m.n_atoms, 2
    full = replicate_system_intra_only(m.system, n_copies)
    for f in full.getForces():
        if isinstance(f, mm.CustomNonbondedForce) and f.getNumTabulatedFunctions() > 0:
            assert _groups(f) == [(set(range(n)), set(range(n))),
                                  (set(range(n, 2 * n)), set(range(n, 2 * n)))]


# --------------------------------------------------------------------------- #
# 6.5 Pinned aggregate regression for a fully-IDR build
# --------------------------------------------------------------------------- #
def test_fully_idr_build_pinned(tmp_path, monkeypatch):
    """Fingerprints a fully-disordered build at the calibrated defaults.

    The counterpart of ``test_idr_baseline.test_nonbonded_interaction_baseline`` for
    the IDR path: deterministic given a fixed STRIDE + MDAnalysis, so a change here is
    a real change to the disordered physics and worth reviewing.
    """
    from topo.utils.nonbonded import (
        build_nonbonded_interaction, IDR_RMIN_2_NM, DEFAULT_IDR_SCALE,
        DEFAULT_EPS_EV_KJ, NON_NATIVE_KJ, get_residue_mapping,
    )
    import MDAnalysis as mda
    import warnings

    monkeypatch.chdir(tmp_path)
    rmin, eps, r2, idr = build_nonbonded_interaction(
        str(PDB), _idr_yaml(tmp_path, "all.yaml", end=N_RES),
        return_rmin_2=True, return_idr=True)

    assert rmin.shape == (N_RES, N_RES) and eps.shape == (N_RES, N_RES)
    assert np.array_equal(idr["idx"], np.arange(N_RES))
    assert idr["eps_ev_kj"] == DEFAULT_EPS_EV_KJ
    assert DEFAULT_IDR_SCALE == 0.10

    # Every radius is the per-AA table value, so the sum is a pure sequence property,
    # and every well sits at the plain sum rule.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        u = mda.Universe(str(PDB))
    _, index_to_resname, _ = get_residue_mapping(u)
    expected_r2 = np.array([IDR_RMIN_2_NM[index_to_resname[i]] for i in range(N_RES)])
    assert np.allclose(r2, expected_r2)
    assert np.allclose(rmin, expected_r2[:, None] + expected_r2[None, :])

    # No native contact survives; every well is the IDR-IDR depth at the calibrated scale.
    assert eps.max() < 1.0
    assert eps.min() >= NON_NATIVE_KJ

    # Aggregate fingerprints (nm, kJ/mol).
    assert r2.sum() == pytest.approx(34.942243, rel=1e-6)
    assert rmin.sum() == pytest.approx(7407.7555, rel=1e-6)
    assert eps.sum() == pytest.approx(2724.9555, rel=1e-6)


# --------------------------------------------------------------------------- #
# 6.3 CSP — the partition survives append_ribosome
# --------------------------------------------------------------------------- #
def _toy_ribosome(m_beads=6, x0=5.0):
    """A minimal rigid ribosome: ``m_beads`` Cα scenery beads on a line, well clear
    of the nascent chain. Enough to exercise ``append_ribosome``'s bookkeeping
    without loading a real 10^4-bead structure."""
    from topo.csp.ribosome import Ribosome
    coords = np.array([[x0 + 0.5 * i, 0.0, 0.0] for i in range(m_beads)])
    return Ribosome(coords_nm=coords,
                    Rmin_2_nm=[0.3] * m_beads,
                    charges=[0.0] * m_beads,
                    names=["CA"] * m_beads,
                    resnames=["ALA"] * m_beads,
                    resids=list(range(1, m_beads + 1)),
                    segids=["RB"] * m_beads)


@pytest.mark.parametrize("regime,idr_end", [("fully-disordered", 12), ("mixed", 6)])
def test_append_ribosome_preserves_the_partition(tmp_path, monkeypatch, regime, idr_end):
    """A ribosome append must not widen the contact force it did not create.

    ``append_ribosome`` used to assert ``{nascent} x {nascent}`` on the contact force
    unconditionally. With disorder that force already owns only ``{folded}x{folded}``,
    and OpenMM takes the *union* — so the blanket group would silently re-admit every
    IDR pair, which the Ashbaugh–Hatch force also evaluates. Both regimes are covered
    because synthesis passes through each: while the emerged chain is still inside the
    disordered prefix the folded set is empty, and later it is not.
    """
    from topo.csp.core import build_length_model
    from topo.csp.ribosome import append_ribosome
    from topo.utils.nonbonded import build_nonbonded_interaction

    monkeypatch.chdir(tmp_path)
    L = 12                                     # emerged nascent length
    dm = _idr_yaml(tmp_path, f"{regime}.yaml", end=idr_end)
    R_full, eps_full, _r2, idr_full = build_nonbonded_interaction(
        str(PDB), dm, return_rmin_2=True, return_idr=True)

    idx = np.asarray(idr_full["idx"], dtype=int)
    idx = idx[idx < L]
    idr_L = {"idx": idx, "eps_ev_kj": idr_full["eps_ev_kj"]}

    sub_pdb = str(tmp_path / f"native_1_{L}.pdb")
    from topo.csp.core import write_subset_structure
    write_subset_structure(str(PDB), L, sub_pdb)
    model = build_length_model(sub_pdb, np.ascontiguousarray(R_full[:L, :L]),
                               np.ascontiguousarray(eps_full[:L, :L]), idr=idr_L)

    ribo = _toy_ribosome()
    nascent_idx, ribo_idx = append_ribosome(model, ribo,
                                            nascent_rmin_2=np.full(L, 0.3))
    n_total = L + ribo.n

    cf, idf = model.custom_non_bonded_force, model.idr_non_bonded_force
    assert idf is not None

    # Every CustomNonbondedForce holds one parameter entry per system particle...
    for f in model.system.getForces():
        if isinstance(f, mm.CustomNonbondedForce):
            assert f.getNumParticles() == n_total
    # ...and the contact forces still agree on exclusions.
    assert _exclusions(cf) == _exclusions(idf)

    go_pairs = _pairs_covered(cf, n_total)
    ah_pairs = _pairs_covered(idf, n_total)
    assert not (go_pairs & ah_pairs), "append_ribosome re-admitted pairs to both forces"

    # Neither contact force may reach a ribosome bead, and the 12-10-6 must still see
    # no disordered bead.
    disordered = set(int(i) for i in idx)
    for pairs in (go_pairs, ah_pairs):
        assert all(i < L and j < L for i, j in pairs)
    assert all(i not in disordered and j not in disordered for i, j in go_pairs)

    if regime == "fully-disordered":
        assert go_pairs == set(), "no folded bead has emerged; the 12-10-6 owns nothing"
    else:
        assert go_pairs, "the mixed regime must still have folded-folded pairs"
        assert any((i in disordered) ^ (j in disordered) for i, j in ah_pairs), \
            "the mixed regime must carry {idr}x{folded} cross pairs"

    # The build is usable: finite energy, no NaN.
    positions = np.vstack([np.array(model.positions.value_in_unit(unit.nanometer))[:L],
                           ribo.coords_nm])
    e = _energy(model.system, positions * unit.nanometer)
    assert np.isfinite(e)


def test_append_ribosome_still_restricts_an_unclaimed_force(tmp_path, monkeypatch):
    """With no disorder the contact force is unclaimed, so the append must restrict it.

    The other half of the ownership rule: "only if it has none" must not become
    "never", or a fully-folded nascent chain would see the ribosome beads through its
    L x L contact table (whose entries for those indices are meaningless).
    """
    from topo.csp.core import build_length_model, write_subset_structure
    from topo.csp.ribosome import append_ribosome
    from topo.utils.nonbonded import build_nonbonded_interaction

    monkeypatch.chdir(tmp_path)
    L = 12
    R_full, eps_full, _r2 = build_nonbonded_interaction(str(PDB), return_rmin_2=True)
    sub_pdb = str(tmp_path / f"native_1_{L}.pdb")
    write_subset_structure(str(PDB), L, sub_pdb)
    model = build_length_model(sub_pdb, np.ascontiguousarray(R_full[:L, :L]),
                               np.ascontiguousarray(eps_full[:L, :L]))
    assert model.idr_non_bonded_force is None
    assert model.custom_non_bonded_force.getNumInteractionGroups() == 0

    ribo = _toy_ribosome()
    append_ribosome(model, ribo, nascent_rmin_2=np.full(L, 0.3))
    assert _groups(model.custom_non_bonded_force) == [
        (set(range(L)), set(range(L)))]
