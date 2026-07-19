"""Green-baseline (characterization) tests for the IDR feature — run BEFORE implementing it.

These lock the *current* output of the two code paths that the IDR feature will change
(spec: ``review/DISORDER_IDR_SPEC.md``):

- **Energy** — ``build_nonbonded_interaction`` (gets ``apply_disorder`` added, §2.3 / §2.6).
- **Analysis / Q** — ``build_native_contacts`` (gains a ``disorder`` mask, §2.7).

They pass on the pre-IDR code and MUST still pass afterwards **whenever no ``disordered:``
section is present** — that is the "no-op / byte-identical" guarantee (spec §3 tests 1, 7, 10).
So: run these first for a green baseline, implement IDR, then keep them green. When the
``disordered:`` reader lands, add the *masked* variants (spec tests 2, 3, 5, 6, 8, 9, 10) that
assert the changed behavior; these baseline tests stay as the regression floor.

Fixture: ``tutorials/A1_single_domain_quickstart/P0CX28_clean.pdb`` (106 residues, all-atom).
STRIDE must be importable/on PATH (it is vendored under ``assets/``); the run is deterministic,
so the pinned values below are exact fingerprints of the current build.

Run with ``pytest tests/test_idr_baseline.py``.
"""
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PDB = REPO_ROOT / "tutorials" / "A1_single_domain_quickstart" / "P0CX28_clean.pdb"

# Tolerance for pinned float aggregates. The pipeline is deterministic given a fixed
# STRIDE + MDAnalysis, so a tight relative tolerance is appropriate; loosen only if a
# genuine library upgrade shifts the reference (that is a real change worth reviewing).
RTOL = 1e-6


pytestmark = pytest.mark.skipif(not PDB.exists(), reason=f"fixture PDB missing: {PDB}")


# --------------------------------------------------------------------------- #
# Energy path — build_nonbonded_interaction (apply_disorder will hook in here)
# --------------------------------------------------------------------------- #
def test_nonbonded_interaction_baseline(tmp_path, monkeypatch):
    """Pins the folded-protein (rmin, eps, rmin_2) build — the IDR no-op target.

    With no ``disordered:`` section, the post-IDR ``build_nonbonded_interaction`` must
    reproduce these exactly (``apply_disorder`` is not called).
    """
    from topo.utils.nonbonded import build_nonbonded_interaction

    # STRIDE caches "{stem}_stride.dat" into CWD; keep it out of the repo tree.
    monkeypatch.chdir(tmp_path)
    rmin, eps, rmin_2 = build_nonbonded_interaction(str(PDB), return_rmin_2=True)

    # Shapes & basic structure.
    assert rmin.shape == (106, 106)
    assert eps.shape == (106, 106)
    assert rmin_2.shape == (106,)
    assert np.allclose(rmin, rmin.T), "rmin matrix must be symmetric"
    assert np.allclose(eps, eps.T), "eps matrix must be symmetric"
    assert np.all(np.isfinite(rmin)) and np.all(np.isfinite(eps))
    assert np.all(rmin_2 > 0.0), "per-residue Rmin/2 must be positive"

    # Native contacts: pairs whose well depth sits well above the non-native floor
    # (NON_NATIVE_KJ ~= 0.00055 kJ/mol). Upper triangle only.
    n_native = int(np.triu(eps > 0.01, k=1).sum())
    assert n_native == 144

    # Aggregate fingerprints (kJ/mol and nm). Exact for the deterministic build.
    assert rmin.sum() == pytest.approx(9227.494434, rel=RTOL)
    assert eps.sum() == pytest.approx(1181.164689, rel=RTOL)
    assert rmin_2.sum() == pytest.approx(43.693692, rel=RTOL)   # nm; feeds CSP NC<->ribosome EV


# --------------------------------------------------------------------------- #
# Analysis path — build_native_contacts (gains the IDR mask, §2.7)
# --------------------------------------------------------------------------- #
def test_native_contacts_baseline():
    """Pins the whole-protein native-contact list (the Q denominator).

    With no ``disordered:`` section, the post-IDR ``build_native_contacts`` (with its new
    ``disorder=None`` default) must reproduce these exactly.
    """
    from topo.analysis.native_contacts import (
        load_universe, reference_residue_geometry, build_native_contacts,
    )

    u = load_universe(str(PDB))
    ca, resindex_to_pos = reference_residue_geometry(u)
    heavy = u.select_atoms("protein and not name H*")
    heavy_pos = heavy.positions.copy()
    heavy_res = np.fromiter(
        (resindex_to_pos[ri] for ri in heavy.resindices), dtype=int, count=heavy.n_atoms
    )
    assert ca.shape[0] == 106

    # Whole protein (group1=group2=None), analysis convention: cutoff 4.5 A, sep > 3.
    pairs, dnat = build_native_contacts(
        None, None, heavy_pos, heavy_res, ca, cutoff=4.5, local_separation=3
    )
    assert pairs.shape == (146, 2)
    assert np.all(pairs[:, 0] < pairs[:, 1]), "pairs must be ordered i < j"
    assert dnat.shape == (146,)
    assert dnat.sum() == pytest.approx(971.42517, rel=RTOL)   # Angstrom
    assert dnat.min() == pytest.approx(4.17135, rel=RTOL)
    assert dnat.max() == pytest.approx(11.99229, rel=RTOL)
