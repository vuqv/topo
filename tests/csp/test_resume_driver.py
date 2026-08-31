"""Driver-level integration tests for CSP resume.

Exercises the real :func:`topo.csp.protocol.run_continuous_synthesis` control flow --
fresh start, mid-run crash, resume, presence guard, idempotent completed run -- with
the heavy MD primitives (``run_length``, ``optimal_ptc_targets``, ribosome load,
contact precompute) stubbed out. This validates the resume *mechanism* (schedule
persist/re-read, progress log, RUNNING-drop, continuation point, post-phase skip)
deterministically and in seconds, without OpenMM / STRIDE / a ribosome PDB.

Run with ``pytest tests/csp/test_resume_driver.py``.
"""
import types

import numpy as np
import pytest

from topo.csp import protocol
from topo.csp.core import RunParams


# --------------------------------------------------------------------------
# INI parsing of the resume key
# --------------------------------------------------------------------------
_BASE_INI = ("[OPTIONS]\n"
             "pdb_file = p.pdb\n"
             "ribosome = r.pdb\n"
             "domain_def = d.yaml\n"
             "codon_times = 0.005\n")


def _write_ini(tmp_path, extra=""):
    path = tmp_path / "csp.ini"
    path.write_text(_BASE_INI + extra)
    return str(path)


def test_ini_resume_default_auto(tmp_path):
    cfg = protocol.read_csp_config(_write_ini(tmp_path), verbose=False)
    assert cfg.params.resume == "auto"


@pytest.mark.parametrize("val", ["yes", "no", "auto", "YES"])
def test_ini_resume_valid(tmp_path, val):
    cfg = protocol.read_csp_config(_write_ini(tmp_path, f"resume = {val}\n"), verbose=False)
    assert cfg.params.resume == val.lower()


def test_ini_resume_invalid(tmp_path):
    with pytest.raises(ValueError, match="resume"):
        protocol.read_csp_config(_write_ini(tmp_path, "resume = maybe\n"), verbose=False)

N_FULL = 6   # the stubbed "protein" length


def _write_final_pdb(path, natoms):
    """Write a minimal valid CA PDB with ``natoms`` atoms."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for i in range(natoms):
            fh.write(f"ATOM  {i + 1:>5d}  CA  ALA A{i + 1:>4d}    "
                     f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n")
        fh.write("END\n")


def _install_stubs(monkeypatch, crash_at_L=None):
    """Patch protocol's MD primitives with fast, deterministic stubs.

    Returns the list of (L, subdir) run_length calls actually made, for assertions.
    If ``crash_at_L`` is set, the stub raises when a residue of that length starts its
    stage 1 (simulating a node failure mid-residue, leaving that unit RUNNING).
    """
    calls = []

    def fake_load_ribosome(path, model="topo"):
        return types.SimpleNamespace(n=10)

    def fake_anchor_coord(ribo, segid, resid, bead):
        return np.zeros(3)

    def fake_optimal_ptc_targets(ribo, **kw):
        return np.array([1.0, 0.0, 0.0]), np.array([1.381, 0.0, 0.0])

    def fake_precompute_contacts(full_pdb, domain_def, stride_output_file):
        R = np.zeros((N_FULL, N_FULL))
        eps = np.zeros((N_FULL, N_FULL))
        rmin = np.full(N_FULL, 0.3)
        return R, eps, rmin, None      # 4th value: the IDR handle (None = no disorder)

    def fake_run_length(L, **kw):
        subdir = kw["out_subdir"]
        out_root = kw["out_root"]
        outname = kw.get("outname", "traj")
        # Consolidated layout: all three stages share the L_<L>/ dir, distinguished by
        # outname (traj_s1/s2/s3). Simulate a crash mid-residue during stage 1.
        if crash_at_L is not None and L == crash_at_L and outname == "traj_s1":
            raise RuntimeError(f"simulated crash at L={L}")
        calls.append((L, subdir))
        # Only stage 3 / the post-phases persist traj_final.pdb (the resume target).
        if kw.get("persist_final", True):
            _write_final_pdb(out_root / subdir / "traj_final.pdb", L)
        return np.zeros((L, 3))

    monkeypatch.setattr(protocol, "load_ribosome", fake_load_ribosome)
    monkeypatch.setattr(protocol, "anchor_coord", fake_anchor_coord)
    monkeypatch.setattr(protocol, "optimal_ptc_targets", fake_optimal_ptc_targets)
    monkeypatch.setattr(protocol, "precompute_contacts", fake_precompute_contacts)
    monkeypatch.setattr(protocol, "run_length", fake_run_length)
    return calls


def _params(**over):
    p = RunParams()
    p.uniform_codon_time = 0.005     # uniform timing -> no mRNA / codon table needed
    p.random_seed = 42
    p.max_steps_per_stage = 10
    p.device = "CPU"
    for k, v in over.items():
        setattr(p, k, v)
    return p


def _run(out_dir, params, **over):
    protocol.run_continuous_synthesis(
        "prot.pdb", "ribo.pdb", L0=1, L_max=5, out_root=str(out_dir),
        params=params, **over)


# --------------------------------------------------------------------------
def test_golden_run(tmp_path, monkeypatch):
    """A clean run writes the schedule, all-DONE progress, and every residue final."""
    _install_stubs(monkeypatch)
    _run(tmp_path, _params())

    from topo.csp import resume as r
    assert (tmp_path / "dwell_times.dat").is_file()
    prog = r.read_progress(tmp_path)
    assert prog.last_done_residue == 5
    for L in range(1, 6):
        assert r.residue_final_path(tmp_path, L).is_file()


def test_crash_then_resume_equivalence(tmp_path, monkeypatch):
    """Crash after residue 2 (L_003 left RUNNING), resume, reach L_max.

    The resumed schedule is byte-identical to an uninterrupted golden run (same seed,
    re-read not redrawn); the partial L_003 dir is dropped and rebuilt.
    """
    from topo.csp import resume as r

    # Golden reference in a separate directory (same seed/config).
    golden = tmp_path / "golden"
    _install_stubs(monkeypatch)
    _run(golden, _params())
    golden_schedule = (golden / "dwell_times.dat").read_bytes()

    # Interrupted run: crash when residue 3 starts.
    work = tmp_path / "work"
    monkeypatch.undo()
    _install_stubs(monkeypatch, crash_at_L=3)
    with pytest.raises(RuntimeError, match="simulated crash"):
        _run(work, _params())
    prog = r.read_progress(work)
    assert prog.last_done_residue == 2
    assert prog.running_units() == ["L_003"]
    # The schedule was persisted before the crash and matches golden byte-for-byte.
    assert (work / "dwell_times.dat").read_bytes() == golden_schedule
    (work / "L_003" / "stage_1").mkdir(parents=True, exist_ok=True)
    (work / "L_003" / "stage_1" / "junk.dcd").write_text("partial")

    # Resume (auto): drops L_003, continues to L_max.
    monkeypatch.undo()
    calls = _install_stubs(monkeypatch)
    _run(work, _params())

    # Only residues 3,4,5 were (re)run -- 1,2 were already DONE.
    rerun_L = sorted({L for (L, sub) in calls if sub.startswith("L_")})
    assert rerun_L == [3, 4, 5]
    # The dropped partial file is gone (dir cleanly rebuilt by the rerun).
    assert not (work / "L_003" / "stage_1" / "junk.dcd").exists()
    prog2 = r.read_progress(work)
    assert prog2.last_done_residue == 5
    # Schedule still byte-identical (never rewritten on resume).
    assert (work / "dwell_times.dat").read_bytes() == golden_schedule


def test_presence_guard_aborts(tmp_path, monkeypatch):
    """A DONE length whose final was deleted aborts the resume, naming that length."""
    _install_stubs(monkeypatch)
    _run(tmp_path, _params())

    from topo.csp import resume as r
    r.residue_final_path(tmp_path, 2).unlink()   # simulate a scratch purge of L_002

    monkeypatch.undo()
    _install_stubs(monkeypatch)
    with pytest.raises(SystemExit, match="L_002"):
        _run(tmp_path, _params(resume="yes"))


def test_resume_completed_is_noop(tmp_path, monkeypatch):
    """Resuming a finished run runs no MD and exits cleanly (idempotent)."""
    _install_stubs(monkeypatch)
    _run(tmp_path, _params())

    monkeypatch.undo()
    calls = _install_stubs(monkeypatch)
    _run(tmp_path, _params())   # resume=auto, all DONE
    assert calls == []          # run_length never called


def test_resume_yes_without_run_errors(tmp_path, monkeypatch):
    """resume = yes on a pristine directory is an error (nothing to resume)."""
    _install_stubs(monkeypatch)
    with pytest.raises(SystemExit, match="nothing to resume"):
        _run(tmp_path, _params(resume="yes"))


def test_post_synthesis_resume(tmp_path, monkeypatch):
    """Crash during ejection: resume reruns ejection (not the residues)."""
    from topo.csp import resume as r

    # Custom stub: crash when the ejection phase starts.
    calls = []

    def make_stubs(crash_ejection):
        def fake_run_length(L, **kw):
            sub = kw["out_subdir"]
            if crash_ejection and sub == "ejection":
                raise RuntimeError("crash in ejection")
            calls.append(sub)
            _write_final_pdb(kw["out_root"] / sub / "traj_final.pdb", L)
            return np.zeros((L, 3))
        monkeypatch.setattr(protocol, "load_ribosome", lambda *a, **k: types.SimpleNamespace(n=10))
        monkeypatch.setattr(protocol, "anchor_coord", lambda *a, **k: np.zeros(3))
        monkeypatch.setattr(protocol, "optimal_ptc_targets",
                            lambda ribo, **k: (np.array([1.0, 0, 0]), np.array([1.381, 0, 0])))
        monkeypatch.setattr(protocol, "precompute_contacts",
                            lambda *a: (np.zeros((N_FULL, N_FULL)), np.zeros((N_FULL, N_FULL)),
                                        np.full(N_FULL, 0.3), None))
        monkeypatch.setattr(protocol, "run_length", fake_run_length)

    make_stubs(crash_ejection=True)
    with pytest.raises(RuntimeError, match="ejection"):
        _run(tmp_path, _params(ejection_steps=10))
    prog = r.read_progress(tmp_path)
    assert prog.last_done_residue == 5
    assert prog.running_units() == ["ejection"]

    calls.clear()
    monkeypatch.undo()
    make_stubs(crash_ejection=False)
    _run(tmp_path, _params(ejection_steps=10))
    # No residue reruns; ejection completed.
    assert "ejection" in calls
    assert not any(s.startswith("L_") for s in calls)
    prog2 = r.read_progress(tmp_path)
    assert prog2.is_done("ejection")
