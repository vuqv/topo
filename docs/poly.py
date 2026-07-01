#!/usr/bin/env python3
"""sphinx-polyversion driver configuration for the TOPO docs.

Builds every whitelisted branch/tag of ``docs/`` into its own subfolder of a
single site (``<out>/<branch>/``) plus a landing ``index.html`` and a
``versions.json`` manifest, so the GitHub Pages site can host several branches
at once with an in-page version switcher (see ``_templates/versions.html``).

Usage
-----
Build every published version (used by .github/workflows/static.yml)::

    sphinx-polyversion docs/poly.py docs/_build/html

Preview only the local working tree (no git checkout, fast)::

    sphinx-polyversion docs/poly.py docs/_build/html --local

Override a setting from the CLI (e.g. add a branch, or build local branches)::

    sphinx-polyversion docs/poly.py docs/_build/html -o BRANCH_REGEX='^(main|dev)$'
    sphinx-polyversion docs/poly.py docs/_build/html -o REMOTE=          # local branches
"""
from datetime import datetime
from pathlib import Path

from sphinx_polyversion import apply_overrides
from sphinx_polyversion.driver import DefaultDriver
from sphinx_polyversion.environment import Environment
from sphinx_polyversion.git import Git, GitRef, GitRefType, file_predicate
from sphinx_polyversion.sphinx import SphinxBuilder

# -- Configuration (each overridable from the CLI via -o KEY=VALUE) -----------
#: Branches to publish (full-match regex). Add branches here as needed.
BRANCH_REGEX = r"^(main|translation|tut15-claude-fix)$"
#: Tags to publish (none by default; use r".*" to publish every tag).
TAG_REGEX = r"^$"
#: Which git remote the branches live on. Left empty so we enumerate *local*
#: branches: sphinx-polyversion 2.0.0 mis-parses remote refs (its ref regex lets
#: ``\w+`` swallow "remotes", so refs/remotes/* are dropped), so `remote="origin"`
#: finds nothing. The CI workflow creates a local branch for each origin branch
#: before building; locally your branches already exist. Override with -o REMOTE=... .
REMOTE = ""
#: Sphinx source dir, relative to the repo root.
SOURCE_DIR = "docs"
#: Output dir (overridden by the positional OUT argument from the workflow).
OUTPUT_DIR = "docs/_build/html"
#: Extra arguments passed to sphinx-build.
SPHINX_ARGS = ["-b", "html"]

#: Mock metadata for `--local` preview builds (build the working tree only).
MOCK_DATA = {
    "revisions": [
        GitRef("local", "", "", GitRefType.BRANCH, datetime.fromtimestamp(0)),
    ],
    "current": GitRef("local", "", "", GitRefType.BRANCH, datetime.fromtimestamp(0)),
}
MOCK = False
SEQUENTIAL = False

# Apply -o overrides / --local / positional OUT from the command line.
apply_overrides(globals())

root = Git.root(Path(__file__).parent)
src = Path(SOURCE_DIR)

DefaultDriver(
    root,
    OUTPUT_DIR,
    vcs=Git(
        branch_regex=BRANCH_REGEX,
        tag_regex=TAG_REGEX,
        remote=REMOTE or None,
        buffer_size=1 * 10**9,  # 1 GB git-archive buffer for large ref checkouts
        predicate=file_predicate([src]),  # skip refs that lack a docs/ dir
    ),
    builder=SphinxBuilder(src, args=SPHINX_ARGS),
    env=Environment.factory(),  # reuse the current (CI) environment; no per-version venv
    static_dir=src / "_polyversion_root",  # copied to the site root (landing index.html)
    mock=MOCK_DATA,
).run(MOCK, SEQUENTIAL)
