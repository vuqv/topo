# Docs build — consistency check & warning cleanup (2026-07-01)

Branch `tut15-claude-fix`. Toolchain: Sphinx 9.1.0, furo, napoleon, myst-parser (env
`bioenv`). Build command: `cd docs && make clean && sphinx-build -b html . _build/html`
(the repo's `docs/build_docs.sh` additionally runs `sphinx-apidoc` first).

## Are the docs consistent with the current source?

**Yes**, with one stale reference (now fixed). The API pages are `automodule`/`autoclass`
based, so they import the live package and track the current code — the build imports every
module with **no import errors or missing-symbol errors**, including the recently changed CSP
modules. Checks run:

- No references to removed symbols (`load_obrien_ribosome`, `load_ribosome_auto`,
  `calculate_sigma_values`, `_parse_prm_radii`, `ribosome_psf/prm`) anywhere in `docs/`.
- One **stale hand-written reference** found and fixed: `docs/usage/python_api.rst` said
  `cg.distance_matrix`; the attribute was renamed to `cg.rmin_matrix` (CHANGELOG §4,
  `core/models.py`). autodoc can't catch prose like this — fixed manually.

## Warning count: before → after

| Build | Warnings | Notes |
|--|--:|--|
| Committed baseline (`build_docs.sh`, apidoc regen) | **160** | as recorded in `docs/build.log` / CHANGELOG §8-2 |
| Clean `sphinx-build` baseline | **112** | 109 duplicate-object + 2 myst header + 1 orphan toctree |
| **After this cleanup** | **0** | `build succeeded.` (no warning count printed) |

Reduction path across the fixes: **112 → 42 → 20 → 0**.

## What caused the warnings, and the fix

1. **109 × "duplicate object description"** — the docs carried **two parallel API trees**:
   curated pages (`docs/modules/*.rst`, wired into `index.rst`) *and* auto-generated flat
   stubs (`modules.rst → topo.rst → topo.*.rst`, **not** in any toctree). Four modules
   (system, models, parameters, reporter) were documented in **both** → one duplicate per
   member. Fixed by making each module document once:
   - Deleted the fully-redundant flat stubs `docs/topo.parameters.rst`,
     `docs/topo.reporter.rst` (their modules are covered by the curated pages).
   - Trimmed `docs/topo.core.rst` to the `geometry` submodule only (system/models are
     curated).
   - Trimmed `docs/topo.rst` subpackage list (dropped `topo.parameters`, `topo.reporter`)
     and removed the package-level `automodule:: topo`.
   - Verified nothing was orphaned: `dihedral.load_dihedral_params`, the `models` class,
     the `system` class, and the whole `topo.reporter` module are each still documented by
     the curated pages (checked module-level defs — no other public members exist).

2. **1 × "document isn't included in any toctree" (`modules.rst`)** — the flat tree was
   orphaned. Fixed by wiring `modules` into `index.rst` under a new "Full module index"
   toctree, so `csp / mdrun / optimize / analysis / utils / core.geometry / engine` are all
   reachable from the nav (they have no curated page).

3. **20 × duplicate dataclass fields (`CSPConfig`/`CSPParams`)** — these dataclasses have a
   napoleon **Attributes** docstring section, which napoleon rendered as standalone
   `.. attribute::` directives that collided with the fields autodoc already emits via
   `:undoc-members:`. Fixed globally with **`napoleon_use_ivar = True`** in `conf.py` (renders
   Attributes as inline `:ivar:` fields — the exact option CHANGELOG §8-2 flagged). Also
   gave the `topo.{csp,mdrun,optimize}` package-level "Module contents" `automodule` a
   `:no-index:` (they re-export their submodule classes via `__all__`).

4. **2 × "Non-consecutive header level H1→H3" (myst)** — `tutorials/README.md` jumped from
   `#` to `###` for "Part A/B". Promoted both to `##`.

5. **1 × "Definition list ends without a blank line"** — a malformed docstring in
   `topo/csp/cg_ribosome.py::_parse_pdb`; converted the inline dict example to a literal
   block (`::` + blank lines).

## Files changed

Docs infra: `docs/conf.py` (napoleon_use_ivar), `docs/index.rst` (wire modules),
`docs/topo.rst`, `docs/topo.core.rst`, `docs/topo.csp.rst`, `docs/topo.mdrun.rst`,
`docs/topo.optimize.rst`; **deleted** `docs/topo.parameters.rst`, `docs/topo.reporter.rst`.
Prose/source: `docs/usage/python_api.rst` (rmin_matrix), `tutorials/README.md` (headers),
`topo/csp/cg_ribosome.py` (docstring formatting).

## Content preserved (spot-checked in rendered HTML)

`CSPParams`/`CSPConfig` fields still documented (once each); curated System page keeps its
methods; `topo.core.geometry` intact; `csp`/`mdrun`/`optimize` fully covered. All 14 top-level
+ 5 curated HTML pages present under `docs/_build/html/`.

## To rebuild

```bash
cd docs
make clean && sphinx-build -b html . _build/html   # -> build succeeded. (0 warnings)
# open docs/_build/html/index.html
```

To fail the build on any future warning (CI), add `-W`:
`sphinx-build -W -b html . _build/html`.
