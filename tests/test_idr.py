"""Masked-variant tests for the IDR / disorder feature (spec ``review/DISORDER_IDR_SPEC.md`` §3).

Companion to ``tests/test_idr_baseline.py`` (the no-op / byte-identical floor). These
exercise the *changed* behavior once a ``disordered:`` section is present:

- energy path  — ``build_nonbonded_interaction`` -> ``apply_disorder`` (§2.3)
- analysis / Q — ``build_native_contacts`` + ``load_domains`` disorder mask (§2.7)
- edge case    — fully-IDP protein (§2.8)

Fixture: ``tutorials/A1_single_domain_quickstart/P0CX28_clean.pdb`` (106 residues, single
chain, resids 1..106 so residue number == 0-based index + 1). STRIDE must be on PATH
(vendored under ``assets/``); the build is deterministic.
"""
from pathlib import Path
import textwrap

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PDB = REPO_ROOT / "tutorials" / "A1_single_domain_quickstart" / "P0CX28_clean.pdb"

pytestmark = pytest.mark.skipif(not PDB.exists(), reason=f"fixture PDB missing: {PDB}")

N_RES = 106
MASK_END = 20                         # disorder residues 1..20 (0-based 0..19)


def _write_yaml(path, body):
    path.write_text(textwrap.dedent(body))
    return str(path)


def _resnames():
    from topo.utils.nonbonded import get_residue_mapping
    import MDAnalysis as mda
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        u = mda.Universe(str(PDB))
    _, index_to_resname, _ = get_residue_mapping(u)
    return index_to_resname


def _mask(end=MASK_END, n=N_RES):
    m = np.zeros(n, dtype=bool)
    m[:end] = True
    return m


# --------------------------------------------------------------------------- #
# Range-compression helper for the IDR overlap log line
# --------------------------------------------------------------------------- #
def test_format_residue_ranges():
    """Consecutive residues collapse to a-b so the overlap log stays short."""
    from topo.utils.nonbonded import format_residue_ranges as f
    assert f([]) == "[]"
    assert f([42]) == "[42]"
    assert f([3, 1, 2, 5, 6, 7, 10]) == "[1-3, 5-7, 10]"          # sorted + deduped runs
    assert f([2, 2, 3, 3, 4]) == "[2-4]"                          # duplicates collapse
    # a fully contiguous block (an IDR spanning a domain) is one range, not N numbers
    assert f(list(range(18, 165))) == "[18-164]"


# --------------------------------------------------------------------------- #
# Energy path — apply_disorder (§2.3)
# --------------------------------------------------------------------------- #
def test_contact_removal_self_avoiding(tmp_path, monkeypatch):
    """§3 test 2 + 6: idr_scale=0 -> every IDR-involving pair is the excluded-volume
    floor; folded-folded pairs and folded residues' K-B radius are untouched; IDR
    residues take the per-AA radius."""
    from topo.utils.nonbonded import build_nonbonded_interaction, NON_NATIVE_KJ, IDR_RMIN_2_NM

    monkeypatch.chdir(tmp_path)
    rmin_f, eps_f, r2_f = build_nonbonded_interaction(str(PDB), return_rmin_2=True)

    dm = _write_yaml(tmp_path / "idr0.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 0.0
        """)
    rmin_m, eps_m, r2_m = build_nonbonded_interaction(str(PDB), dm, return_rmin_2=True)

    m = _mask()
    involves = m[:, None] | m[None, :]
    ff = ~involves

    # Every masked pair -> exactly the excluded-volume floor (idr_scale=0, floor guard).
    assert np.allclose(eps_m[involves], NON_NATIVE_KJ)
    # Folded-folded pairs are byte-identical to the folded-only run.
    assert np.allclose(eps_m[ff], eps_f[ff])
    assert np.allclose(rmin_m[ff], rmin_f[ff])

    # IDR residues' Rmin/2 -> per-AA; folded residues' Rmin/2 unchanged.
    resname = _resnames()
    for i in range(MASK_END):
        assert r2_m[i] == pytest.approx(IDR_RMIN_2_NM[resname[i]])
    assert np.allclose(r2_m[MASK_END:], r2_f[MASK_END:])
    # No native contact survives that touches a masked residue.
    assert not np.any(eps_m[involves] > 0.01)


def test_idr_folded_excluded_only(tmp_path, monkeypatch):
    """§3 test 3: with idr_scale>0, IDR-folded pairs (exactly one masked) are still
    excluded-volume only (no attraction leaks across the boundary), at well position
    r2[IDR](per-AA) + r2[folded](K-B)."""
    from topo.utils.nonbonded import build_nonbonded_interaction, NON_NATIVE_KJ

    monkeypatch.chdir(tmp_path)
    dm = _write_yaml(tmp_path / "idr.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 0.03
        """)
    rmin_m, eps_m, r2_m = build_nonbonded_interaction(str(PDB), dm, return_rmin_2=True)

    m = _mask()
    dd = m[:, None] & m[None, :]
    cross = (m[:, None] | m[None, :]) & ~dd     # exactly one masked

    assert np.allclose(eps_m[cross], NON_NATIVE_KJ)
    expected_rmin = r2_m[:, None] + r2_m[None, :]
    assert np.allclose(rmin_m[cross], expected_rmin[cross])


def test_eps_construction(tmp_path, monkeypatch):
    """§3 test 5: IDR-IDR eps == max(NON_NATIVE, idr_scale * ss_interaction_energy),
    with the raw BT identity ss == 4.184*|raw-0.6| and no nscale factor."""
    from topo.utils.nonbonded import (
        build_nonbonded_interaction, get_ss_interaction_energy, load_bt_potential,
        NON_NATIVE_KJ, KCAL_TO_KJ, BT_SHIFT_KCAL,
    )
    import MDAnalysis as mda
    import warnings

    monkeypatch.chdir(tmp_path)
    dm = _write_yaml(tmp_path / "idr.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 0.03
        """)
    _, eps_m, _ = build_nonbonded_interaction(str(PDB), dm, return_rmin_2=True)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        u = mda.Universe(str(PDB))
    ss = get_ss_interaction_energy(u)

    m = _mask()
    dd = m[:, None] & m[None, :]
    expected = np.maximum(NON_NATIVE_KJ, 0.03 * ss)
    assert np.allclose(eps_m[dd], expected[dd])

    # Raw BT identity: the transferable per-pair energy is 4.184*|raw-0.6| from the CSV.
    import pandas as pd
    import topo.utils.nonbonded as nb
    csv = Path(nb.__file__).parent.parent / "parameters" / "data" / "bt_potential.csv"
    raw = pd.read_csv(csv, index_col=0)
    df = load_bt_potential()
    assert np.allclose(df.values, KCAL_TO_KJ * np.abs(raw.values - BT_SHIFT_KCAL))


def test_eps_gen_additive(tmp_path, monkeypatch):
    """Additive generic cohesion: IDR-IDR eps == max(NON_NATIVE, eps_gen_kj + a*ss),
    and eps_gen_kj: 0 is identical to omitting it (idr_scale-only build)."""
    from topo.utils.nonbonded import (
        build_nonbonded_interaction, get_ss_interaction_energy, NON_NATIVE_KJ,
    )
    import MDAnalysis as mda
    import warnings

    monkeypatch.chdir(tmp_path)
    EPS_GEN = 2.5104   # 0.6 kcal/mol in kJ/mol
    dm_gen = _write_yaml(tmp_path / "gen.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 1.0
          eps_gen_kj: {EPS_GEN}
        """)
    dm_zero = _write_yaml(tmp_path / "zero.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 1.0
          eps_gen_kj: 0.0
        """)
    dm_absent = _write_yaml(tmp_path / "absent.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 1.0
        """)
    _, eps_gen, _ = build_nonbonded_interaction(str(PDB), dm_gen, return_rmin_2=True)
    _, eps_zero, _ = build_nonbonded_interaction(str(PDB), dm_zero, return_rmin_2=True)
    _, eps_absent, _ = build_nonbonded_interaction(str(PDB), dm_absent, return_rmin_2=True)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        u = mda.Universe(str(PDB))
    ss = get_ss_interaction_energy(u)

    m = _mask()
    dd = m[:, None] & m[None, :]
    expected = np.maximum(NON_NATIVE_KJ, EPS_GEN + 1.0 * ss)
    assert np.allclose(eps_gen[dd], expected[dd])
    # eps_gen_kj: 0 == omitting it == the original idr_scale-only build.
    assert np.allclose(eps_zero, eps_absent)
    # a positive eps_gen deepens every IDR-IDR well relative to eps_gen_kj = 0.
    assert np.all(eps_gen[dd] >= eps_zero[dd])
    assert np.any(eps_gen[dd] > eps_zero[dd])


def test_default_idr_scale(tmp_path, monkeypatch):
    """§3 test 9: omitting idr_scale is identical to an explicit idr_scale: 0.03."""
    from topo.utils.nonbonded import build_nonbonded_interaction

    monkeypatch.chdir(tmp_path)
    dm_default = _write_yaml(tmp_path / "default.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
        """)
    dm_explicit = _write_yaml(tmp_path / "explicit.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 0.03
        """)
    r1, e1, x1 = build_nonbonded_interaction(str(PDB), dm_default, return_rmin_2=True)
    r2, e2, x2 = build_nonbonded_interaction(str(PDB), dm_explicit, return_rmin_2=True)
    assert np.allclose(r1, r2) and np.allclose(e1, e2) and np.allclose(x1, x2)


def test_overlap_precedence(tmp_path, monkeypatch):
    """§3 test 8: a disordered range overlapping a domain (nscale != 1) is governed
    entirely by the disorder rules -- identical to the same residues disordered but
    in no domain. The domain nscale has no effect on IDR-involving pairs."""
    from topo.utils.nonbonded import build_nonbonded_interaction

    monkeypatch.chdir(tmp_path)
    # Run 1: domain A covers everything at nscale 2.0; 1..20 also disordered (overlap).
    dm_overlap = _write_yaml(tmp_path / "overlap.yaml", f"""
        n_residues: {N_RES}
        intra_domains:
          A: {{ residues: [1-{N_RES}], nscale: 2.0 }}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 0.03
        """)
    # Run 2: domain A only over the folded residues; 1..20 disordered, in no domain.
    dm_disjoint = _write_yaml(tmp_path / "disjoint.yaml", f"""
        n_residues: {N_RES}
        intra_domains:
          A: {{ residues: [{MASK_END + 1}-{N_RES}], nscale: 2.0 }}
        disordered:
          residues: [1-{MASK_END}]
          idr_scale: 0.03
        """)
    r1, e1, _ = build_nonbonded_interaction(str(PDB), dm_overlap, return_rmin_2=True)
    r2, e2, _ = build_nonbonded_interaction(str(PDB), dm_disjoint, return_rmin_2=True)

    m = _mask()
    involves = m[:, None] | m[None, :]
    # IDR-involving pairs are identical regardless of the overlapped domain nscale.
    assert np.allclose(e1[involves], e2[involves])
    assert np.allclose(r1[involves], r2[involves])
    # Folded-folded pairs keep the domain's scaling (both runs put 21..106 in A@2.0).
    assert np.allclose(e1[~involves], e2[~involves])


# --------------------------------------------------------------------------- #
# Analysis path — build_native_contacts / load_domains (§2.7)
# --------------------------------------------------------------------------- #
def _heavy_geometry():
    from topo.analysis.native_contacts import load_universe, reference_residue_geometry
    u = load_universe(str(PDB))
    ca, resindex_to_pos = reference_residue_geometry(u)
    heavy = u.select_atoms("protein and not name H*")
    heavy_pos = heavy.positions.copy()
    heavy_res = np.fromiter(
        (resindex_to_pos[ri] for ri in heavy.resindices), dtype=int, count=heavy.n_atoms
    )
    return ca, heavy_pos, heavy_res


def test_q_excludes_idr_contacts():
    """§3 test 10: build_native_contacts drops every pair touching a disordered
    residue; load_domains subtracts disorder from each domain (effective = dom - dis)."""
    from topo.analysis.native_contacts import build_native_contacts, load_domains

    ca, heavy_pos, heavy_res = _heavy_geometry()
    disorder = set(range(MASK_END))     # 0-based 0..19

    pairs_all, _ = build_native_contacts(None, None, heavy_pos, heavy_res, ca,
                                         cutoff=4.5, local_separation=3)
    pairs_msk, _ = build_native_contacts(None, None, heavy_pos, heavy_res, ca,
                                         cutoff=4.5, local_separation=3, disorder=disorder)
    # No masked residue appears; strictly fewer contacts (the tail has native contacts).
    assert not np.any(np.isin(pairs_msk, list(disorder)))
    assert pairs_msk.shape[0] < pairs_all.shape[0]

    # load_domains: overlapped residues subtracted from the domain; disorder returned.
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".yaml"); os.close(fd)
    Path(path).write_text(textwrap.dedent(f"""
        n_residues: {N_RES}
        intra_domains:
          A: {{ residues: [1-{N_RES}], nscale: 1.0 }}
        disordered:
          residues: [1-{MASK_END}]
        """))
    domains, n_residues, dis = load_domains(path, N_RES)
    os.unlink(path)
    assert dis == disorder
    assert not np.any(np.isin(domains["A"], list(disorder)))     # effective = A - disorder
    assert domains["A"].min() >= MASK_END


def test_native_contacts_no_disorder_regression():
    """§3 test 10 regression: disorder=None reproduces the pre-IDR pair list exactly."""
    from topo.analysis.native_contacts import build_native_contacts

    ca, heavy_pos, heavy_res = _heavy_geometry()
    p0, d0 = build_native_contacts(None, None, heavy_pos, heavy_res, ca, 4.5, 3)
    p1, d1 = build_native_contacts(None, None, heavy_pos, heavy_res, ca, 4.5, 3,
                                   disorder=None)
    assert np.array_equal(p0, p1) and np.allclose(d0, d1)
    assert p0.shape == (146, 2)     # pinned in test_idr_baseline


# --------------------------------------------------------------------------- #
# Fully-IDP edge case (§2.8)
# --------------------------------------------------------------------------- #
def test_fully_idp(tmp_path, monkeypatch):
    """§3 test 11: every residue disordered."""
    from topo.utils.nonbonded import build_nonbonded_interaction, NON_NATIVE_KJ, get_ss_interaction_energy
    from topo.analysis.native_contacts import (
        build_native_contacts, fraction_native_contacts, load_domains,
    )
    import MDAnalysis as mda
    import warnings

    monkeypatch.chdir(tmp_path)
    dm = _write_yaml(tmp_path / "fully_idp.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{N_RES}]
          idr_scale: 0.03
        """)

    # (a) energy build finishes: all pairs IDR-IDR, all radii finite (no divide-by-zero).
    rmin, eps, r2 = build_nonbonded_interaction(str(PDB), dm, return_rmin_2=True)
    assert np.all(np.isfinite(rmin)) and np.all(np.isfinite(eps))
    assert np.all(r2 > 0.0)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        u = mda.Universe(str(PDB))
    ss = get_ss_interaction_energy(u)
    assert np.allclose(eps, np.maximum(NON_NATIVE_KJ, 0.03 * ss))

    # (b) analysis: no native contacts -> empty pairs, Q = NaN (not a crash).
    ca, heavy_pos, heavy_res = _heavy_geometry()
    pairs, dnat = build_native_contacts(None, None, heavy_pos, heavy_res, ca, 4.5, 3,
                                        disorder=set(range(N_RES)))
    assert pairs.shape[0] == 0
    assert np.isnan(fraction_native_contacts(ca, pairs, dnat, tolerance=1.2))

    # (c) load_domains accepts a domain_def with no intra_domains (empty domains dict).
    domains, n_residues, dis = load_domains(dm, N_RES)
    assert domains == {}
    assert dis == set(range(N_RES))


def test_optimizer_nothing_to_optimize(tmp_path, monkeypatch):
    """§3 test 11(d): the Scorer over a fully-IDP domain_def has zero foldable
    contacts -- the guard condition the optimizer uses to exit with 'nothing to
    optimize' (spec §2.8)."""
    from topo.optimize.optimize import Scorer

    monkeypatch.chdir(tmp_path)
    dm = _write_yaml(tmp_path / "fully_idp.yaml", f"""
        n_residues: {N_RES}
        disordered:
          residues: [1-{N_RES}]
          idr_scale: 0.03
        """)
    scorer = Scorer(str(PDB), dm)
    assert scorer.unit_keys() == []
    assert sum(scorer.n_contacts(k) for k in scorer.unit_keys()) == 0
