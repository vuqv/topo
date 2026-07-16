"""The ``[OPTIONS]`` header is optional in topo control files.

``[OPTIONS]`` is the only section a topo control file ever has, so it carries no
information -- it exists because configparser demands a header. These tests pin
that both spellings parse identically, and that the leniency does not change how
files that *do* carry the header behave.
"""
import configparser

import pytest

from topo.utils.config import read_ini, read_simulation_config


WITH_HEADER = """\
[OPTIONS]
md_steps = 3_000
dt = 0.015
tcoupl = yes
ref_t = 300
nstxout = 50
"""

WITHOUT_HEADER = """\
md_steps = 3_000
dt = 0.015
tcoupl = yes
ref_t = 300
nstxout = 50
"""


def _write(tmp_path, text, name="md.ini"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_header_optional_parses_same_as_header(tmp_path):
    """A bare key = value file and a headered one give identical settings."""
    with_h = read_simulation_config(_write(tmp_path, WITH_HEADER, "a.ini"),
                                    verbose=False)
    without_h = read_simulation_config(_write(tmp_path, WITHOUT_HEADER, "b.ini"),
                                       verbose=False)
    for field in ("md_steps", "nstxout", "tcoupl"):
        assert getattr(with_h, field) == getattr(without_h, field)
    assert with_h.md_steps == 3000
    assert with_h.dt == without_h.dt
    assert with_h.ref_t == without_h.ref_t


def test_headerless_file_with_leading_comments(tmp_path):
    """Comments before the first key must not defeat the fallback."""
    text = ("# a leading comment\n"
            "; another, in the other comment style\n"
            "\n"
            "md_steps = 1234\n"
            "nstxout = 7   ; inline comment still stripped\n")
    cfg = read_simulation_config(_write(tmp_path, text), verbose=False)
    assert cfg.md_steps == 1234
    assert cfg.nstxout == 7


def test_headerless_preserves_inline_comment_stripping(tmp_path):
    """Both `#` and `;` inline comments are stripped without the header too."""
    text = "md_steps = 500  # hash comment\nmodel = topo  ; semicolon comment\n"
    cp = read_ini(_write(tmp_path, text))
    assert cp["OPTIONS"]["md_steps"] == "500"
    assert cp["OPTIONS"]["model"] == "topo"


def test_preserve_case_option(tmp_path):
    """preserve_case keeps key spelling; the default lower-cases it."""
    text = "MD_Steps = 10\n"
    assert "md_steps" in read_ini(_write(tmp_path, text))["OPTIONS"]
    assert "MD_Steps" in read_ini(_write(tmp_path, text),
                                  preserve_case=True)["OPTIONS"]


def test_missing_file_raises_filenotfound(tmp_path):
    """A missing file is a FileNotFoundError, not a confusing KeyError."""
    with pytest.raises(FileNotFoundError):
        read_ini(str(tmp_path / "does_not_exist.ini"))


@pytest.mark.parametrize("header", ["[OPTIONS]", "[options]", "[Options]",
                                    "[SETTINGS]", "[MD]", "[DEFAULT]", ""])
def test_any_header_name_is_accepted(tmp_path, header):
    """Header names are decoration: every spelling reads the same settings."""
    text = f"{header}\nmd_steps = 4242\nnstxout = 11\n"
    cfg = read_simulation_config(_write(tmp_path, text), verbose=False)
    assert cfg.md_steps == 4242
    assert cfg.nstxout == 11


def test_keys_under_extra_headers_are_not_dropped(tmp_path):
    """Regression: keys under a second header used to be silently ignored.

    A file that set nstxout under any header but [OPTIONS] parsed "fine" and ran
    with the *default* nstxout -- a wrong run with no error. Every key in the
    file must be applied regardless of which header it sits under.
    """
    text = ("[OPTIONS]\nmd_steps = 7777\n"
            "[EXTRA]\nnstxout = 9\n"
            "[MORE]\nnstlog = 3\n")
    cfg = read_simulation_config(_write(tmp_path, text), verbose=False)
    assert cfg.md_steps == 7777
    assert cfg.nstxout == 9, "key under [EXTRA] was dropped"
    assert cfg.nstlog == 3, "key under [MORE] was dropped"


def test_conflicting_duplicate_key_across_headers_raises(tmp_path):
    """Contradictory values must not resolve silently to last-one-wins."""
    text = "[A]\nref_t = 300\n[B]\nref_t = 310\n"
    with pytest.raises(ValueError, match="set more than once"):
        read_simulation_config(_write(tmp_path, text), verbose=False)


def test_duplicate_key_raises_even_within_one_header(tmp_path):
    """A repeated key is an error wherever it is written, not just across headers."""
    text = "md_steps = 50\nmd_steps = 60\n"
    with pytest.raises(ValueError, match="set more than once"):
        read_simulation_config(_write(tmp_path, text), verbose=False)


def test_box_dimension_list_value_survives(tmp_path):
    """Regression: a `[...]` *value* must not be mistaken for a section header.

    An unanchored `\\[.*\\]` strip would blank out box_dimension, and an empty
    box makes read_simulation_config silently switch pbc off -- a physics change
    with no error.
    """
    text = ("[OPTIONS]\n"
            "pbc = yes\n"
            "box_dimension = [10, 10, 10]\n"
            "[EXTRA]\n"
            "nstxout = 9\n")
    cfg = read_simulation_config(_write(tmp_path, text), verbose=False)
    assert cfg.pbc is True, "pbc was silently turned off"
    assert cfg.box_dimension == [10, 10, 10]
    assert cfg.nstxout == 9


def test_indented_bracket_value_is_not_a_header(tmp_path):
    """An indented `[...]` is a value continuation, not a header."""
    text = "pbc = yes\nbox_dimension =\n    [12, 12, 12]\n"
    cfg = read_simulation_config(_write(tmp_path, text), verbose=False)
    assert cfg.box_dimension == [12, 12, 12]


def test_header_with_trailing_comment(tmp_path):
    """A header line with an inline comment is still a header."""
    text = "[OPTIONS]  ; the old mandatory header\nmd_steps = 21\n"
    cfg = read_simulation_config(_write(tmp_path, text), verbose=False)
    assert cfg.md_steps == 21


def test_parse_error_line_number_points_at_the_real_line(tmp_path):
    """Headers are blanked, not deleted, so reported line numbers stay honest."""
    text = ("[OPTIONS]\n"       # line 1
            "md_steps = 10\n"   # line 2
            "[EXTRA]\n"         # line 3
            "garbage-with-no-separator\n")  # line 4  <- the offending line
    with pytest.raises(configparser.ParsingError) as exc:
        read_simulation_config(_write(tmp_path, text), verbose=False)
    assert "[line  4]" in str(exc.value), f"expected line 4, got: {exc.value}"


def test_headerless_parse_error_line_number_is_also_honest(tmp_path):
    """The headerless path reuses a comment line, so numbering holds there too."""
    text = ("# a leading comment\n"          # line 1  <- [OPTIONS] is hosted here
            "md_steps = 10\n"                # line 2
            "garbage-with-no-separator\n")   # line 3  <- the offending line
    with pytest.raises(configparser.ParsingError) as exc:
        read_simulation_config(_write(tmp_path, text), verbose=False)
    assert "[line  3]" in str(exc.value), f"expected line 3, got: {exc.value}"


def test_duplicate_reports_a_line_number(tmp_path):
    text = "md_steps = 1\nmd_steps = 2\n"
    with pytest.raises(ValueError, match=r"line 2"):
        read_simulation_config(_write(tmp_path, text), verbose=False)


def test_empty_file_is_reported(tmp_path):
    """A file with no settings is an error, not a silent all-defaults run."""
    with pytest.raises(ValueError, match="no settings found"):
        read_simulation_config(_write(tmp_path, "# just a comment\n"), verbose=False)


def test_percent_in_value_is_literal(tmp_path):
    """`%` is data (e.g. in a path), not configparser %(interpolation)s."""
    cp = read_ini(_write(tmp_path, "outname = run_50%_done\n"))
    assert cp["OPTIONS"]["outname"] == "run_50%_done"


def test_optimize_ini_outdir_resolves_against_the_ini(tmp_path, monkeypatch):
    """A relative `outdir` follows the .ini, not the working directory.

    pdb_file/domain_def already resolve against the file; outdir must too, or the
    same optimize.ini would write its rounds somewhere else depending on where it
    was launched from.
    """
    from topo.optimize.optimize import read_optimize_config

    proj = tmp_path / "proj"
    proj.mkdir()
    (proj / "p.pdb").touch()
    (proj / "d.yaml").touch()
    (proj / "optimize.ini").write_text(
        "pdb_file = p.pdb\ndomain_def = d.yaml\noutdir = my_run\n")

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)          # launch from a different cwd

    _, _, sim_options, controls = read_optimize_config(str(proj / "optimize.ini"))
    assert controls["outdir"] == str(proj / "my_run")
    assert "outdir" not in sim_options, "outdir leaked into the per-round md.ini"


def test_optimize_ini_outdir_unset_is_none(tmp_path):
    """Unset must stay None, so -o/--outdir can be told apart from a file value."""
    from topo.optimize.optimize import read_optimize_config

    (tmp_path / "p.pdb").touch()
    (tmp_path / "d.yaml").touch()
    _, _, _, controls = read_optimize_config(
        _write(tmp_path, "pdb_file = p.pdb\ndomain_def = d.yaml\n", "optimize.ini"))
    assert controls["outdir"] is None


def test_optimize_ini_header_optional(tmp_path):
    """read_optimize_config accepts a headerless optimize.ini too."""
    from topo.optimize.optimize import read_optimize_config

    (tmp_path / "p.pdb").touch()
    (tmp_path / "d.yaml").touch()
    text = "pdb_file = p.pdb\ndomain_def = d.yaml\nntraj = 3\n"
    pdb, domain, sim_options, controls = read_optimize_config(
        _write(tmp_path, text, "optimize.ini"))
    assert pdb.endswith("p.pdb")
    assert domain.endswith("d.yaml")
    assert controls["ntraj"] == 3
    # implicit defaults still fill in
    assert sim_options["md_steps"] == "10_000"
