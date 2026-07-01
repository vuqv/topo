# Refactor: `strength` → `nscale` (native-contact scaling factor)

**Rationale.** "strength" was misleading: the per-domain / per-interface value is a
*scaling factor* on the sidechain–sidechain native-contact well depths, not an absolute
interaction strength — at `nscale = 1.0` the contacts still interact. Renamed to `nscale`
(the term already used by the optimizer: "the nscale optimizer", `topo-optimize`).

**Decisions (user-approved).** Rename + keep `strength` as a **deprecated YAML alias**
(non-breaking); update docs/tutorials too.

## What changed

**Python identifiers** (`topo/`)
- `topo/utils/nonbonded.py`: `read_yaml_config` returns `intra_nscales`/`inter_nscales`;
  reads the `nscale` key, **falls back to legacy `strength`** with a one-time
  `DeprecationWarning` + stderr notice. `get_scaling_ss_matrix` updated. Module docstring
  clarified (nscale = scaling factor, not absolute energy).
- `topo/optimize/optimize.py`: `strength_for`→`nscale_for`, `strength_dom`→`nscale_dom`,
  `strength_int`→`nscale_int`; `write_domain_yaml` writes the `nscale` key and **pops any
  legacy `strength`** so optimizer output is always clean. Docstrings/logs/CLI help updated.
- `topo/optimize/__init__.py`: export `nscale_for` (was `strength_for`).
- `topo/csp/protocol.py`, `topo/analysis/native_contacts.py`: docstring/comment wording.
- Normalized `n_scale` → `nscale` spelling across code/docs (kept the LaTeX
  `n_\mathrm{scale}` math symbol and O'Brien's own `nscal` term).

**Config (user-facing)**
- YAML key `strength:` → `nscale:` in **18** `domain.yaml` files (repo-wide). Legacy key
  still accepted by the reader.
- INI comments referencing "strength"/"n_scale" updated (6 files).

**Docs / tutorials**
- ~25 `.md`/`.rst` files updated (114 tokens); `docs/usage/domain_definition.rst` gains a
  note that `nscale` is the key and `strength` is a deprecated alias. Sphinx rebuild: **0
  warnings**.
- Top-level `README.md`, `pyproject.toml` comment, `topo/csp/{DESIGN,FILES,README}.md`.

## Verification
- `python -m compileall topo/` OK; all submodules import; `nscale_for` importable.
- `read_yaml_config` + `get_scaling_ss_matrix` build correctly from a migrated repo YAML.
- Alias: new `nscale` key silent; legacy `strength` key works + emits deprecation notice.
- Optimizer round-trip (`write_domain_yaml` → `read_yaml_config`), incl. a legacy-keyed
  seed: output contains only `nscale`, reads back correctly, no warning.
- Global sweep: no `strength`/`n_scale` remains anywhere except the intentional
  legacy-alias handling in `nonbonded.py` and the deprecation notes.

## Follow-up (optional)
- The deprecated `strength` alias can be removed in a future release once external
  `domain.yaml` files are migrated. Search `nonbonded.py` for `used_legacy_key`.
