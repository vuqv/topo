"""Deterministic unit tests for topo.analysis.mirror.

Covers the acceptance tests in the mirror-detection spec. No simulation is
needed: synthetic coordinates are built from a random native Cα reference R with
a fixed RNG seed, plus a minimal end-to-end CLI smoke test on a 2-frame DCD
(native frame, then reflected frame).

Run with ``pytest tests/test_mirror.py``.
"""
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from topo.analysis import mirror as m


SEED = 20240607
N = 40          # residues in the synthetic reference
M_PTS = 16      # points used for chirality tests (>= 4)


@pytest.fixture
def R():
    """Random native Cα reference, shape (N, 3)."""
    rng = np.random.default_rng(SEED)
    return rng.standard_normal((N, 3)) * 5.0


@pytest.fixture
def proper_rotation():
    """A random proper rotation matrix (det = +1)."""
    rng = np.random.default_rng(SEED + 1)
    Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1.0
    assert np.isclose(np.linalg.det(Q), 1.0)
    return Q


# --------------------------------------------------------------------------- #
# Metric 1 -- RMSD
# --------------------------------------------------------------------------- #
def test_identity_native_frame(R):
    """Test 1: native frame -> RMSD_native ~ 0 and < RMSD_reflected."""
    rmsd_native, rmsd_reflected = m.mirror_rmsd(R, R)
    assert rmsd_native <= 1e-6
    assert rmsd_native < rmsd_reflected


def test_proper_rotation_invariance(R, proper_rotation):
    """Test 2: a proper rotation of R fits the native reference at RMSD ~ 0."""
    rmsd_native, _ = m.mirror_rmsd(R @ proper_rotation, R)
    assert rmsd_native <= 1e-6


def test_mirror_frame(R):
    """Test 3: reflected frame -> RMSD_reflected ~ 0 and < RMSD_native.

    Confirms the reflection is not absorbed by the (proper-rotation) fit.
    """
    reflected = m.reflect_coords(m.center_coords(R), "x")
    rmsd_native, rmsd_reflected = m.mirror_rmsd(reflected, R)
    assert rmsd_reflected <= 1e-6
    assert rmsd_reflected < rmsd_native


def test_kabsch_is_proper(R, proper_rotation):
    """Test 4: the Kabsch fit returns a proper rotation (det ~ +1) for both
    the native and the reflected target."""
    mobile = m.center_coords(R @ proper_rotation)
    for target in (m.center_coords(R),
                   m.reflect_coords(m.center_coords(R), "x")):
        cov = mobile.T @ target
        V, _S, Wt = np.linalg.svd(cov)
        d = np.sign(np.linalg.det(V @ Wt))
        rot = V @ np.diag([1.0, 1.0, d]) @ Wt
        assert np.linalg.det(rot) > 0
        assert np.isclose(np.linalg.det(rot), 1.0)


# --------------------------------------------------------------------------- #
# Metric 2 -- chirality
# --------------------------------------------------------------------------- #
def test_chirality_sign_flip(R):
    """Test 5: local_chirality(reflect(P)) == -local_chirality(P)."""
    P = R[:M_PTS]
    chi = m.local_chirality(P)
    chi_reflected = m.local_chirality(m.reflect_coords(m.center_coords(P), "x"))
    assert np.allclose(chi, -chi_reflected, atol=1e-6)


def test_chirality_rotation_invariance(R, proper_rotation):
    """Test 6: chirality depends only on shape/handedness, not orientation."""
    P = R[:M_PTS]
    assert np.allclose(m.local_chirality(P),
                       m.local_chirality(P @ proper_rotation), atol=1e-6)


def test_K_native_and_mirror(R):
    """K == 1 for the native points, K == 0 for the reflected points."""
    P = R[:M_PTS]
    chi_ref = m.local_chirality(P)
    K_native = m.chirality_agreement(m.local_chirality(P), chi_ref)
    K_mirror = m.chirality_agreement(
        m.local_chirality(m.reflect_coords(m.center_coords(P), "x")), chi_ref)
    assert K_native[0] == 1.0
    assert K_mirror[0] == 0.0


# --------------------------------------------------------------------------- #
# Combined classifier
# --------------------------------------------------------------------------- #
def test_classifier_sanity():
    """Test 7: a folded native frame is not a mirror; a folded reflected frame
    is. (Q high for both; K and the RMSD ratio decide.)"""
    Q = np.array([0.8, 0.8])
    K = np.array([1.0, 0.0])
    rmsd_native = np.array([10.0, 10.0])
    rmsd_reflected = np.array([9.5, 1.0])   # native: ratio 0.95 (>0.9); mirror: 0.1
    ratio, is_mirror = m.classify_mirror(Q, K, rmsd_native, rmsd_reflected)
    assert not is_mirror[0]
    assert is_mirror[1]
    assert np.isclose(ratio[1], 0.1)


def test_classifier_rmsd_ratio_gate():
    """A near-tie RMSD (ratio >= 0.9) is not flagged even with folded + inverted."""
    _ratio, is_mirror = m.classify_mirror(
        Q=[0.8], K=[0.0], rmsd_native=[10.0], rmsd_reflected=[9.5])
    assert not is_mirror[0]


# --------------------------------------------------------------------------- #
# Guards / edge cases
# --------------------------------------------------------------------------- #
def test_too_few_segments_raises(tmp_path):
    """Test 8a: fewer than 2 SS segments (< 4 endpoints) raises."""
    seg = tmp_path / "seg.txt"
    seg.write_text("1 3 8\n")           # only one segment -> 2 endpoints
    with pytest.raises(ValueError):
        m.segment_endpoints_from_file(str(seg))


def test_degenerate_bond_guarded():
    """Test 6 (edge): coincident endpoints must not produce NaN/inf chirality."""
    P = np.zeros((M_PTS, 3))
    P[:, 0] = np.arange(M_PTS)
    P[5] = P[4]                          # coincident points -> zero-length bond
    chi = m.local_chirality(P)
    assert np.all(np.isfinite(chi))


# --------------------------------------------------------------------------- #
# End-to-end CLI smoke test
# --------------------------------------------------------------------------- #
def _write_ca_pdb(path, coords):
    """Minimal all-atom-ish PDB: one CA per ALA residue (enough for the CLI's
    'protein and name CA' selection and Q's heavy-atom pass)."""
    lines = []
    for i, (x, y, z) in enumerate(coords, start=1):
        lines.append(
            f"ATOM  {i:>5d}  CA  ALA A{i:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def test_cli_end_to_end(tmp_path, R):
    """Build a 2-frame DCD (native, reflected), run the CLI, and assert row 0 is
    native-like and row 1 is mirror-like."""
    mda = pytest.importorskip("MDAnalysis")
    import pandas as pd

    ref = tmp_path / "ref.pdb"
    _write_ca_pdb(ref, R)

    # Segment file spanning the chain so chirality is well defined.
    seg = tmp_path / "seg.txt"
    seg.write_text("1 2 10\n2 15 25\n3 30 38\n")

    # 2-frame DCD: frame 0 = native R, frame 1 = reflected R.
    u = m.load_universe(str(ref))
    ca = u.select_atoms("protein and name CA")
    reflected = m.reflect_coords(m.center_coords(R), "x")
    dcd = tmp_path / "traj.dcd"
    with mda.Writer(str(dcd), n_atoms=ca.n_atoms) as w:
        ca.positions = R.astype(np.float32)
        w.write(ca)
        ca.positions = reflected.astype(np.float32)
        w.write(ca)

    out = tmp_path / "mirror.csv"
    subprocess.run(
        [sys.executable, "-m", "topo.analysis.mirror",
         "-r", str(ref), "-f", str(dcd), "--segments", str(seg),
         "-o", str(out)],
        check=True,
    )

    df = pd.read_csv(out)
    assert list(df.columns) == [
        "Frame", "Q", "K", "RMSD_native", "RMSD_reflected", "RMSD_ratio", "is_mirror"]
    assert len(df) == 2
    # Row 0 native-like: RMSD_native < RMSD_reflected, high K, not a mirror.
    assert df.loc[0, "RMSD_native"] < df.loc[0, "RMSD_reflected"]
    assert df.loc[0, "K"] == 1.0
    assert not bool(df.loc[0, "is_mirror"])
    # Row 1 mirror-like: RMSD_reflected < RMSD_native, K == 0, flagged mirror.
    assert df.loc[1, "RMSD_reflected"] < df.loc[1, "RMSD_native"]
    assert df.loc[1, "K"] == 0.0
    assert bool(df.loc[1, "is_mirror"])
