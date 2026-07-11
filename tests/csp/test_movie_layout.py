"""Layout-discovery tests for topo.csp.movie under the consolidated §3.5 CSP layout.

Deterministic and MD-free: they build synthetic ``L_<L>/`` trees with marker files and
check the layout detection + per-stage segment discovery (no trajectory parsing, which
happens later in the stitch core). Guards the movie.py blast-radius of the layout change.

Run with ``pytest tests/csp/test_movie_layout.py``.
"""
from topo.csp import movie


def _make_residue(root, L, stages=(1, 2, 3), final=True, empty=()):
    d = root / f"L_{L:03d}"
    d.mkdir()
    (d / "traj.psf").write_text("psf")             # shared topology
    (d / f"native_1_{L}.pdb").write_text("native")
    for s in stages:
        (d / f"traj_s{s}.dcd").write_text("" if s in empty else "dcd-bytes")
    if final:
        (d / "traj_final.pdb").write_text("ATOM\n")
    return d


def test_is_csp_layout_detects_per_stage_dcds(tmp_path):
    _make_residue(tmp_path, 1)
    assert movie._is_csp_layout(str(tmp_path))


def test_flat_layout_is_not_csp(tmp_path):
    d = tmp_path / "L_001"
    d.mkdir()
    (d / "traj.dcd").write_text("d")   # single-segment cylinder layout, no traj_s*
    (d / "traj.psf").write_text("psf")
    assert not movie._is_csp_layout(str(tmp_path))


def test_find_stage_segments_consolidated(tmp_path):
    _make_residue(tmp_path, 1)
    _make_residue(tmp_path, 2)
    segs = movie.find_stage_segments(str(tmp_path))
    assert [s[0] for s in segs] == ["L=1 s1", "L=1 s2", "L=1 s3",
                                    "L=2 s1", "L=2 s2", "L=2 s3"]
    for label, L, psf, traj in segs:
        assert psf.endswith("traj.psf")          # shared per residue
        assert "traj_s" in traj                  # per-stage trajectory
        assert L in (1, 2)


def test_stage3_empty_dcd_falls_back_to_final(tmp_path):
    # s1 present; s2 absent -> dropped; s3 absent -> falls back to the single traj_final.pdb.
    _make_residue(tmp_path, 1, stages=(1,), final=True)
    segs = movie.find_stage_segments(str(tmp_path))
    labels = [s[0] for s in segs]
    assert labels == ["L=1 s1", "L=1 s3"]
    assert segs[1][3].endswith("traj_final.pdb")


def test_zero_byte_dcds_skipped_quietly(tmp_path):
    # An older run: stages 1/2 ran fewer steps than nstout -> 0-byte (headerless) DCDs.
    # They must be skipped (not fed to the reader, which would raise premature-EOF);
    # stage 3 has frames + a final, so the residue is still represented.
    _make_residue(tmp_path, 7, stages=(1, 2, 3), final=True, empty=(1, 2))
    segs = movie.find_stage_segments(str(tmp_path))
    assert [s[0] for s in segs] == ["L=7 s3"]
    # a 0-byte DCD is never chosen as a segment trajectory
    assert all(s[3] == "" or __import__("os").path.getsize(s[3]) > 0 for s in segs)


def test_residue_without_psf_skipped(tmp_path):
    d = tmp_path / "L_001"
    d.mkdir()
    (d / "traj_s1.dcd").write_text("d")   # no traj.psf -> the whole residue is skipped
    assert movie.find_stage_segments(str(tmp_path)) == []
