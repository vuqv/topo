"""Verification bundle for the IDR Ashbaugh-Hatch force refactor (lambda / HPS form).

Context
-------
The shipped IDR force is written as

    r <= R:  epsEV*LJ + (epsEV - eps_ij)        r > R:  eps_ij*LJ

which is *exactly* the HPS / Ashbaugh-Hatch lambda form

    r <= R:  eps*LJ + (1 - lambda_ij)*eps       r > R:  lambda_ij*eps*LJ

with ``eps = epsEV`` and ``lambda_ij = eps_ij / epsEV``. Rewriting the force in lambda
form is therefore an **energy-identical refactor** for the current (BT-derived) model;
an HPS lambda table is then just a second way to fill ``lambda_ij``.

System
------
Tutorial A.7: P28566 (365 residues, AlphaFold), one folded domain (nscale 2.5044) plus
two IDRs [1-18, 218-281] that overlap the domain ranges at their boundaries. Three builds:

* ``folded`` -- the tutorial ``domain.yaml`` with its ``disordered:`` block removed;
* ``mixed``  -- the tutorial ``domain.yaml`` as shipped (read, never modified);
* ``full``   -- every residue disordered.

What this file locks
--------------------
Against a snapshot of the pre-refactor build (``tests/data/idr_ah_reference.npz``):

* the non-bonded matrices (rmin, eps, Rmin/2) and the IDR index set;
* the energy of **every force** (bonds, angles, torsions, Yukawa, 12-10-6, IDR) and the
  per-bead forces of the two contact forces,

at three configurations: native, compressed 0.85x about the centroid (probes the
repulsive core), and Gaussian-perturbed (sigma 0.05 nm, fixed seed).

Plus an **independent numpy oracle** that evaluates the IDR force in the lambda form from
the matrices alone and must agree with OpenMM -- the statement of the equivalence above,
and the test to extend for an HPS mode.

Expected status after the refactor, in the default (BT) mode: **every test here stays
green, unchanged.** If the build API changes (e.g. returns a lambda table instead of
eps), adapt ``_oracle_idr_energy``'s *inputs*, not its physics.

Pre-existing tests that pin the current *expression text* and may need deliberate edits:
``test_idr_force.py::test_idr_potential_*`` (they parse ``energy_function`` strings and
set the ``epsEV`` / ``rmin_scale`` globals by name).

Regenerate the snapshot (only from a commit whose physics you trust)::

    python tests/test_idr_ah_verification.py --regen
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest
import yaml

pytest.importorskip("openmm")
import openmm as mm
from openmm import unit

REPO_ROOT = Path(__file__).resolve().parents[1]
A7 = REPO_ROOT / "tutorials" / "A7_idr_mixed_protein"
PDB = A7 / "P28566.pdb"
DOMAIN_YAML = A7 / "domain.yaml"
REF = Path(__file__).resolve().parent / "data" / "idr_ah_reference.npz"

BUILDS = ("folded", "mixed", "full")
CONFIGS = ("native", "compressed", "perturbed")
CONTACT_FORCES = ("Custom Non-Bonded Energy", "IDR Non-Bonded Energy")

# Co-translational case: the A.7 chain at two nascent lengths with a toy ribosome.
# The IDR is N-terminal (1-18), so L=12 is the fully-disordered regime and L=40 the
# mixed one (IDR 1-18 + folded 19-40).
CSP_LENGTHS = (12, 40)
CSP_CONFIGS = ("native", "perturbed")

# Reference platform is double precision and deterministic.
E_RTOL, E_ATOL = 1e-9, 1e-8          # kJ/mol
F_RTOL, F_ATOL = 1e-9, 1e-6          # kJ/mol/nm

pytestmark = pytest.mark.skipif(not (PDB.exists() and DOMAIN_YAML.exists()),
                                reason=f"A7 fixture missing under {A7}")


# --------------------------------------------------------------------------- #
# Builds and evaluation
# --------------------------------------------------------------------------- #
def _domain_def(workdir, build):
    """Path to the domain_def for ``build``; derived files go to ``workdir``."""
    if build == "mixed":
        return str(DOMAIN_YAML)
    cfg = yaml.safe_load(DOMAIN_YAML.read_text())
    if build == "folded":
        cfg.pop("disordered")
    elif build == "full":
        cfg = {"n_residues": cfg["n_residues"],
               "disordered": {"residues": [f"1-{cfg['n_residues']}"]}}
    else:
        raise ValueError(build)
    p = Path(workdir) / f"{build}.yaml"
    with open(p, "w") as fh:
        yaml.safe_dump(cfg, fh)
    return str(p)


def _build_all(workdir):
    """Build every model + matrices with CWD = ``workdir`` (STRIDE caches into CWD)."""
    import topo
    from topo.utils.nonbonded import build_nonbonded_interaction

    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        out = {}
        for b in BUILDS:
            dd = _domain_def(workdir, b)
            model = topo.models.buildCoarseGrainModel(str(PDB), domain_def=dd,
                                                      check_forces=False)
            rmin, eps, r2, idr = build_nonbonded_interaction(
                str(PDB), dd, return_rmin_2=True, return_idr=True)
            out[b] = dict(model=model, rmin=rmin, eps=eps, r2=r2, idr=idr)
        out["csp"] = _build_csp(workdir, out["mixed"])
        return out
    finally:
        os.chdir(cwd)


def _toy_ribosome(nascent_pos):
    """Six rigid, charged scenery beads 0.7 nm off every other nascent bead, so the
    NC<->ribosome excluded-volume and Yukawa terms are both non-zero."""
    from topo.csp.ribosome import Ribosome
    m = 6
    coords = nascent_pos[:2 * m:2] + np.array([0.0, 0.0, 0.7])
    return Ribosome(coords_nm=coords,
                    Rmin_2_nm=[0.3] * m,
                    charges=[1.0, -1.0] * (m // 2),
                    names=["CA"] * m,
                    resnames=["ALA"] * m,
                    resids=list(range(1, m + 1)),
                    segids=["RB"] * m)


def _build_csp(workdir, mixed):
    """Nascent-chain + toy-ribosome systems at ``CSP_LENGTHS``, built the way
    ``topo.csp`` does: full-chain matrices sliced to L, IDR index sliced to L."""
    from topo.csp.core import build_length_model, write_subset_structure
    from topo.csp.ribosome import append_ribosome

    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        out = {}
        for L in CSP_LENGTHS:
            idx = np.asarray(mixed["idr"]["idx"], dtype=int)
            idr_L = dict(mixed["idr"], idx=idx[idx < L])
            sub_pdb = str(Path(workdir) / f"native_1_{L}.pdb")
            write_subset_structure(str(PDB), L, sub_pdb)
            model = build_length_model(
                sub_pdb, np.ascontiguousarray(mixed["rmin"][:L, :L]),
                np.ascontiguousarray(mixed["eps"][:L, :L]), idr=idr_L)
            nascent = np.asarray(model.positions.value_in_unit(unit.nanometer))[:L]
            ribo = _toy_ribosome(nascent)
            append_ribosome(model, ribo, nascent_rmin_2=mixed["r2"][:L])
            out[L] = dict(model=model, nascent=nascent, ribo=ribo.coords_nm)
        return out
    finally:
        os.chdir(cwd)


def _csp_positions(d, config):
    """Ribosome beads stay fixed; only the nascent chain is perturbed."""
    x = d["nascent"]
    if config == "perturbed":
        x = x + np.random.default_rng(20261117).normal(0.0, 0.05, size=x.shape)
    return np.vstack([x, d["ribo"]])


def _evaluate_system(system, pos):
    """Energy (kJ/mol) of every force in ``system``, keyed by index and class, so a
    force that is not in ``model.forceGroups`` (NC<->ribosome terms) is pinned too."""
    for i, f in enumerate(system.getForces()):
        f.setForceGroup(i)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatform("Reference"))
    ctx.setPositions(pos * unit.nanometer)
    return {f"{i:02d}:{type(f).__name__}":
            ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy()
               .value_in_unit(unit.kilojoule_per_mole)
            for i, f in enumerate(system.getForces())}


def _positions(model, config):
    x = np.asarray(model.positions.value_in_unit(unit.nanometer), dtype=float)
    if config == "native":
        return x
    if config == "compressed":
        c = x.mean(axis=0)
        return c + 0.85 * (x - c)
    if config == "perturbed":
        return x + np.random.default_rng(20261115).normal(0.0, 0.05, size=x.shape)
    raise ValueError(config)


def _evaluate(model, pos):
    """Per-force energy (kJ/mol) and per-force per-bead forces (kJ/mol/nm)."""
    ctx = mm.Context(model.system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatform("Reference"))
    ctx.setPositions(pos * unit.nanometer)
    energies, forces = {}, {}
    for name, force in model.forceGroups.items():
        st = ctx.getState(getEnergy=True, getForces=True, groups={force.getForceGroup()})
        energies[name] = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces[name] = st.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer)
    return energies, forces


def _key(*parts):
    return "|".join(parts)


def _snapshot(builds):
    """Flatten everything this bundle pins into an npz-able dict."""
    snap = {}
    for L, d in builds["csp"].items():
        for c in CSP_CONFIGS:
            e = _evaluate_system(d["model"].system, _csp_positions(d, c))
            for name in e:
                snap[_key(f"csp{L}", c, "E", name)] = np.array(e[name])
    for b, d in builds.items():
        if b == "csp":
            continue
        snap[_key(b, "rmin")] = d["rmin"]
        snap[_key(b, "eps")] = d["eps"]
        snap[_key(b, "r2")] = d["r2"]
        snap[_key(b, "idr_idx")] = (np.asarray(d["idr"]["idx"]) if d["idr"]
                                    else np.array([], dtype=int))
        for c in CONFIGS:
            e, f = _evaluate(d["model"], _positions(d["model"], c))
            for name in e:
                snap[_key(b, c, "E", name)] = np.array(e[name])
                if name in CONTACT_FORCES:
                    snap[_key(b, c, "F", name)] = f[name]
    return snap


# --------------------------------------------------------------------------- #
# Independent oracle: the IDR force in HPS lambda form, from the matrices alone
# --------------------------------------------------------------------------- #
def _switch(r, rs=1.8, rc=2.0):
    """OpenMM's CustomNonbondedForce switching polynomial (1 below rs, 0 at rc)."""
    x = np.clip((r - rs) / (rc - rs), 0.0, 1.0)
    return 1.0 - 10 * x**3 + 15 * x**4 - 6 * x**5


def _oracle_idr_energy(pos, rmin, eps_ij, eps_ev, idr_idx, excluded):
    """Ashbaugh-Hatch / HPS lambda form, summed over IDR-involving pairs.

    U = LJ + (1 - lambda)*eps   (r <= R)      U = lambda*eps*LJ   (r > R)
    with eps = eps_ev, lambda = eps_ij / eps_ev, sigma = 2^(-1/6) R, 2.0 nm cutoff
    and 1.8 nm switch. Exclusions are taken from the force (topology, not physics).
    """
    n = len(pos)
    idr = np.zeros(n, dtype=bool)
    idr[list(idr_idx)] = True
    i, j = np.triu_indices(n, k=1)
    keep = (idr[i] | idr[j]) & np.array([(a, b) not in excluded for a, b in zip(i, j)])
    i, j = i[keep], j[keep]

    r = np.linalg.norm(pos[i] - pos[j], axis=1)
    R = rmin[i, j]
    lam = eps_ij[i, j] / eps_ev
    sr6 = ((R / 2.0 ** (1.0 / 6.0)) / r) ** 6
    lj = 4.0 * (sr6**2 - sr6)
    u = np.where(r <= R, eps_ev * lj + (1.0 - lam) * eps_ev, lam * eps_ev * lj)
    u = np.where(r < 2.0, u * _switch(r), 0.0)
    return float(u.sum())


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def builds(tmp_path_factory):
    return _build_all(tmp_path_factory.mktemp("idr_ah"))


@pytest.fixture(scope="module")
def ref():
    if not REF.exists():
        pytest.fail(f"reference snapshot missing: {REF} (run this file with --regen)")
    with np.load(REF) as z:
        return {k: z[k] for k in z.files if not k.startswith("meta|")}


@pytest.fixture(scope="module")
def current(builds):
    return _snapshot(builds)


def test_snapshot_has_the_same_keys(ref, current):
    """No force or build silently appears or disappears."""
    assert set(ref) == set(current)


@pytest.mark.parametrize("build", BUILDS)
def test_matrices_unchanged(ref, current, build):
    for q in ("rmin", "eps", "r2", "idr_idx"):
        np.testing.assert_array_equal(current[_key(build, q)], ref[_key(build, q)],
                                      err_msg=f"{build}:{q}")


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize("build", BUILDS)
def test_every_force_energy_unchanged(ref, current, build, config):
    keys = [k for k in ref if k.startswith(_key(build, config, "E", ""))]
    assert keys
    for k in keys:
        assert float(current[k]) == pytest.approx(float(ref[k]), rel=E_RTOL, abs=E_ATOL), k


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize("build", BUILDS)
def test_contact_forces_unchanged(ref, current, build, config):
    keys = [k for k in ref if k.startswith(_key(build, config, "F", ""))]
    assert keys
    for k in keys:
        np.testing.assert_allclose(current[k], ref[k], rtol=F_RTOL, atol=F_ATOL, err_msg=k)


@pytest.mark.parametrize("config", CSP_CONFIGS)
@pytest.mark.parametrize("length", CSP_LENGTHS)
def test_csp_every_force_energy_unchanged(ref, current, length, config):
    """Co-translational build: every force, NC<->ribosome terms included."""
    keys = [k for k in ref if k.startswith(_key(f"csp{length}", config, "E", ""))]
    assert keys
    for k in keys:
        assert float(current[k]) == pytest.approx(float(ref[k]), rel=E_RTOL, abs=E_ATOL), k


def test_csp_bundle_exercises_the_ribosome_terms(current):
    """Guard on the bundle itself: the IDR force and at least two further non-bonded
    forces (NC<->ribosome excluded volume, Yukawa) must be non-zero at each length."""
    for L in CSP_LENGTHS:
        e = {k.split("|")[-1]: float(v) for k, v in current.items()
             if k.startswith(_key(f"csp{L}", "perturbed", "E", ""))}
        nonzero = [k for k, v in e.items() if "Nonbonded" in k and abs(v) > 1e-6]
        assert len(nonzero) >= 3, (L, e)


def test_bundle_exercises_both_branches(builds):
    """Guard on the bundle itself: across the configs, IDR pairs must sample both
    r <= R (core) and r > R (tail), or one branch of the potential goes unchecked."""
    d = builds["mixed"]
    n = len(d["r2"])
    idr = np.zeros(n, dtype=bool)
    idr[list(d["idr"]["idx"])] = True
    i, j = np.triu_indices(n, k=3)
    sel = idr[i] | idr[j]
    i, j = i[sel], j[sel]
    for c in CONFIGS:
        pos = _positions(d["model"], c)
        r = np.linalg.norm(pos[i] - pos[j], axis=1)
        inside = int(np.sum(r <= d["rmin"][i, j]))
        if c != "native":
            assert inside >= 10, f"{c}: only {inside} IDR pairs inside R"
        assert np.sum((r > d["rmin"][i, j]) & (r < 2.0)) >= 100, c


@pytest.mark.parametrize("config", CONFIGS)
@pytest.mark.parametrize("build", ["mixed", "full"])
def test_idr_force_matches_lambda_form_oracle(builds, build, config):
    """OpenMM's IDR force == HPS lambda-form Ashbaugh-Hatch evaluated in numpy."""
    d = builds[build]
    m = d["model"]
    f = m.idr_non_bonded_force
    excluded = {tuple(sorted(f.getExclusionParticles(k))) for k in range(f.getNumExclusions())}
    pos = _positions(m, config)
    e_omm, _ = _evaluate(m, pos)
    e_ref = _oracle_idr_energy(pos, d["rmin"], d["eps"], d["idr"]["epsilon_kj"],
                               d["idr"]["idx"], excluded)
    assert e_omm["IDR Non-Bonded Energy"] == pytest.approx(e_ref, rel=1e-9, abs=1e-8)


@pytest.mark.parametrize("config", CONFIGS)
def test_idr_force_matches_lambda_form_oracle_under_hps_urry(tmp_path_factory, config):
    """The same oracle for the second model, built on the A.7 system. Not part of the
    snapshot: the oracle is the reference. Covers the negative GLU-GLU lambda."""
    import topo
    from topo.utils.nonbonded import build_nonbonded_interaction

    workdir = tmp_path_factory.mktemp("idr_ah_urry")
    cfg = yaml.safe_load(DOMAIN_YAML.read_text())
    cfg["disordered"]["model"] = "hps_urry"
    dd = str(workdir / "hps_urry.yaml")
    with open(dd, "w") as fh:
        yaml.safe_dump(cfg, fh)
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        m = topo.models.buildCoarseGrainModel(str(PDB), domain_def=dd, check_forces=False)
        rmin, eps, idr = build_nonbonded_interaction(str(PDB), dd, return_idr=True)
    finally:
        os.chdir(cwd)
    assert idr["model"] == "hps_urry"
    f = m.idr_non_bonded_force
    excluded = {tuple(sorted(f.getExclusionParticles(k))) for k in range(f.getNumExclusions())}
    pos = _positions(m, config)
    e_omm, _ = _evaluate(m, pos)
    e_ref = _oracle_idr_energy(pos, rmin, eps, idr["epsilon_kj"], idr["idx"], excluded)
    assert e_omm["IDR Non-Bonded Energy"] == pytest.approx(e_ref, rel=1e-9, abs=1e-8)


@pytest.mark.parametrize("length", CSP_LENGTHS)
def test_idr_force_ignores_ribosome_beads(builds, length):
    """Moving the ribosome with the nascent chain fixed leaves the IDR force energy
    unchanged, while the NC<->ribosome terms do change."""
    d = builds["csp"][length]
    system = d["model"].system
    idr_i = next(i for i, f in enumerate(system.getForces())
                 if f.this == d["model"].idr_non_bonded_force.this)
    pos = _csp_positions(d, "native")
    moved = pos.copy()
    moved[length:] += np.array([0.0, 0.0, 0.15])
    e0, e1 = _evaluate_system(system, pos), _evaluate_system(system, moved)
    key = next(k for k in e0 if k.startswith(f"{idr_i:02d}:"))
    assert e1[key] == e0[key]
    assert any(abs(e1[k] - e0[k]) > 1e-6 for k in e0 if k != key)


# --------------------------------------------------------------------------- #
# Snapshot regeneration
# --------------------------------------------------------------------------- #
def _regen():
    with tempfile.TemporaryDirectory() as td:
        snap = _snapshot(_build_all(td))
    # "-dirty" marks a snapshot taken from uncommitted changes on top of that commit.
    commit = subprocess.run(["git", "-C", str(REPO_ROOT), "describe", "--always",
                             "--dirty", "--abbrev=40"],
                            capture_output=True, text=True).stdout.strip()
    snap["meta|commit"] = np.array(commit)
    REF.parent.mkdir(exist_ok=True)
    np.savez_compressed(REF, **snap)
    print(f"wrote {REF} ({len(snap)} arrays, {commit})")


if __name__ == "__main__":
    if "--regen" in sys.argv:
        _regen()
    else:
        sys.exit(pytest.main([__file__, "-q"]))
