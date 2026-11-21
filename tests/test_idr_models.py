"""The shipped IDR parameter sets (``topo/parameters/data/idr/<model>.yaml``).

Each model is a precomputed residue-pair table ``lambda_ij`` plus its energy scale.
These tests pin the files to their sources, so a hand edit or a stale table fails.
"""
import importlib.util
import textwrap
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PDB = REPO_ROOT / "tutorials" / "A1_single_domain_quickstart" / "P0CX28_clean.pdb"
N_RES = 106
MASK_END = 20


def _generator():
    """``scripts/make_idr_lambda_tables.py`` as a module."""
    path = REPO_ROOT / "scripts" / "make_idr_lambda_tables.py"
    spec = importlib.util.spec_from_file_location("make_idr_lambda_tables", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("model", ["bt", "hps_urry", "calvados2"])
def test_shipped_table_matches_its_generator(model):
    """The shipped file is exactly what the generator produces from the sources."""
    from topo.utils.nonbonded import load_idr_model

    lam, eps = load_idr_model(model)
    built = _generator().build_tables()[model]
    assert eps == built["epsilon_kj"] == 0.8368
    for a, row in built["lambda"].items():
        for b, value in row.items():
            assert lam.loc[a, b] == value, (a, b)


def test_bt_table_is_the_scaled_bt_potential():
    """``lambda_ij = 0.10 * eps_BT(i, j) / 0.8368``, rounded to 4 decimals."""
    from topo.utils.nonbonded import load_idr_model, load_bt_potential

    lam, eps = load_idr_model("bt")
    bt = load_bt_potential().loc[lam.index, lam.columns].to_numpy()
    assert np.array_equal(lam.to_numpy(), np.round(0.10 * bt / eps, 4))
    assert lam.to_numpy().min() > 0.0


def test_hps_urry_pairs_match_the_publication():
    """``lambda_ij = (lambda_i + lambda_j) / 2`` from the Urry per-residue values
    (Regy et al. 2021): ALA 0.522942, ASP 0.214119, GLU -0.08, TRP 0.92, TYR 0.817059."""
    from topo.utils.nonbonded import load_idr_model

    lam, _ = load_idr_model("hps_urry")
    assert lam.loc["ALA", "ALA"] == 0.5229
    assert lam.loc["TRP", "TYR"] == lam.loc["TYR", "TRP"] == 0.8685
    assert lam.loc["GLU", "ASP"] == 0.0671
    # The only non-positive pair: GLU-GLU keeps its published negative value.
    assert lam.loc["GLU", "GLU"] == -0.08
    assert int((lam.to_numpy() <= 0).sum()) == 1


def test_calvados2_pairs_match_the_publication():
    """``lambda_ij = (lambda_i + lambda_j) / 2`` from the CALVADOS 2 per-residue values
    (Tesei & Lindorff-Larsen 2022): ALA 0.274330, GLU 0.000694, ARG 0.730762,
    TRP 0.989376, TYR 0.977461."""
    from topo.utils.nonbonded import load_idr_model

    lam, eps = load_idr_model("calvados2")
    assert eps == 0.8368
    assert lam.loc["ALA", "ALA"] == 0.2743
    assert lam.loc["TRP", "TYR"] == lam.loc["TYR", "TRP"] == 0.9834
    assert lam.loc["ARG", "GLU"] == 0.3657
    # Every pair is attractive; GLU-GLU is the weakest.
    assert lam.to_numpy().min() == lam.loc["GLU", "GLU"] == 0.0007


def test_negative_lambda_pair_is_purely_repulsive():
    """GLU-GLU under hps_urry has no well: U >= 0 and never increases with r."""
    from test_idr_force import _shipped_idr_expression, _u_pair
    from topo.utils.nonbonded import load_idr_model

    lam, eps = load_idr_model("hps_urry")
    R = 0.60
    rs = np.linspace(0.7 * R, 3.0 * R, 400)
    u = _u_pair(_shipped_idr_expression(), rs, R, lam.loc["GLU", "GLU"] * eps, eps)
    assert np.all(u >= 0.0)
    assert np.all(np.diff(u) <= 1e-12)


def test_unknown_model_is_refused():
    from topo.utils.nonbonded import load_idr_model

    with pytest.raises(ValueError, match="bt, hps_urry, calvados2"):
        load_idr_model("mpipi")


def test_bead_sizes_do_not_depend_on_the_model(tmp_path, monkeypatch):
    """Radii and well positions are identical under both models; only the well
    depths of IDR-involving pairs differ."""
    from topo.utils.nonbonded import build_nonbonded_interaction

    monkeypatch.chdir(tmp_path)
    out = {}
    for model in ("bt", "hps_urry"):
        dm = tmp_path / f"{model}.yaml"
        dm.write_text(textwrap.dedent(f"""
            n_residues: {N_RES}
            disordered:
              residues: [1-{MASK_END}]
              model: {model}
            """))
        out[model] = build_nonbonded_interaction(str(PDB), str(dm), return_rmin_2=True)

    (rmin_a, eps_a, r2_a), (rmin_b, eps_b, r2_b) = out["bt"], out["hps_urry"]
    assert np.array_equal(r2_a, r2_b)
    assert np.array_equal(rmin_a, rmin_b)

    m = np.zeros(N_RES, dtype=bool)
    m[:MASK_END] = True
    involves = m[:, None] | m[None, :]
    assert np.array_equal(eps_a[~involves], eps_b[~involves])   # folded-folded untouched
    assert not np.allclose(eps_a[involves], eps_b[involves])
