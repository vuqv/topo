# TOPO — Code & Documentation Review

**Date:** 2026-07-03
**Scope:** full `topo` package source + all docs (`docs/`, `topo/csp/README.md`, tutorials).
**Method:** two parallel review agents (non-CSP modules; all docs-vs-code) + direct review of
`topo/csp/*` and independent verification of the high-impact findings against source.

**Overall:** the codebase is solid — careful math, edge cases handled, no severe correctness
bugs on the main single-chain path. The recent `codon_times` refactor is clean and consistent.
Findings cluster in (a) **documentation drift** — one broken example and pervasive stale
tutorial numbers in `continuous_synthesis.md`, plus a few wrong defaults — and (b) a handful of
**latent edge-case bugs** in `topo/core` / `topo/utils`.

Findings are IDed (`D#` docs, `S#` source) for easy reference.

Verification status: `tau_t` default, `pbc`-on-invalid-box, `checkBondDistances >=`, the broken
API example, and the `codon_times` refactor consistency were **verified directly against source**
during this review. Items marked *(agent-reported)* were found by the review agents reading the
code but not independently re-verified here — quick to confirm before acting.

---

## Documentation findings

### HIGH

**D1 — Broken Python API example.** ✅ **FIXED (2026-07-03).**
`docs/usage/continuous_synthesis.md` example (a) did `run_continuous_synthesis(**kwargs)` on the
`CSPConfig` **dataclass** returned by `read_csp_config` — a `TypeError`, and the field names don't
match the params anyway. Replaced with the correct unpack: `cfg = read_csp_config("csp.ini")` then
`run_continuous_synthesis(cfg.pdb_file, cfg.ribosome, L0=cfg.L0, L_max=cfg.L_max,
out_root=cfg.outdir, mrna=cfg.mrna, trans_times=cfg.trans_times, domain_def=cfg.domain_def,
stride_output_file=cfg.stride_output_file, params=cfg.params)`. Field names + signature verified
against `protocol.py`; docs rebuilt clean.

**D2 — Pervasive stale tutorial references** in `docs/usage/continuous_synthesis.md`.
✅ **FIXED (2026-07-03).** The page pointed at tutorial dirs that no longer exist. All references
retargeted to **Tutorial 8** (`tutorials/08_ribosome_synthesis/`, subdirs `4c5c/` + `P0CX28/`,
INIs `csp_debug.ini` smoke / `csp_val.ini` full-length): worked-example bullet, quick-start block
(`cd tutorials/08_ribosome_synthesis/4c5c`, `topo-csp -f csp_debug.ini`), and the in-text mentions
at the tunnel-wall, worked-example, scale_factor, stability-guard, and `dwell_times.dat` /
"See also" spots. Two references to files that never existed in the current tree
(`eject_demo.py`, a bilingual `THEORY.md`) were reworded out (egress is now described via the
`ejection_steps`/`dissociation_steps` INI knobs). Verified: zero residual `Tutorial 12/13`,
`12_auto`, `13_validate`, `eject_demo`, or `THEORY.md` strings; docs rebuilt clean.

### MEDIUM

**D3 — Stale tutorial-12 links** in `topo/csp/README.md:15-17,122`. ✅ **FIXED (2026-07-03).**
Both `../../tutorials/12_auto/` + `PLAN.md` links (dirs/files that no longer exist) retargeted:
the example link → `../../tutorials/08_ribosome_synthesis/` (with its real INIs), and the
"done/remaining matrix" links → `../../review/TODO.md`. Both targets confirmed to exist.

**D4 — Wrong `ref_t` default (and cross-page inconsistency).** ✅ **FIXED (2026-07-03).**
`continuous_synthesis.md:363` listed `310` while `RunParams.ref_t` was `300` and
`synthesis_control.rst` said `300`. **Resolution:** per user decision the *code* default was
changed **300 → 310** (`core.py:826`, so kinetics and thermostat agree by default, matching the
bundled 310 K codon-time table). Docs reconciled to `310`: `synthesis_control.rst:243` (and its
note reworded), `cylinder_synthesis.md:208`; `continuous_synthesis.md:363` was already `310` and
is now correct. Isolated-protein `SimulationConfig.ref_t = 300` left unchanged (separate runner,
correctly documented as 300 in `simulation_control.rst`). Verified `RunParams`/`CylinderParams`
report `ref_t = 310`; docs rebuilt clean.

**D5 — Wrong `tau_t` default.** ✅ **FIXED (2026-07-03).** `simulation_control.rst:168` said
`0.01`; the code default is **0.05** everywhere (`config.py:308` parser default and
`RunParams.tau_t`), so this was a **doc-only** error — no code change was needed (the requested
"tau_t default 0.05 in code" was already the case). Corrected the doc to `0.05`, now consistent
with the same page's example (L45) and annealing note (L371). Docs rebuilt clean.

**D6 — Imprecise runner claim.** ✅ **FIXED (2026-07-03).** `overview.rst:60-65`: "Both tutorials
use … `topo-csp`." Tutorial 7 uses **`topo-cylinder`**; only Tutorial 8 uses `topo-csp`. Rewritten
to name both runners explicitly (Tutorial 7 → `topo-cylinder`, Tutorial 8 → `topo-csp`); docs
rebuilt clean.

**D7 — Stale "only synthesis runner" claim.** ✅ **FIXED (2026-07-03).** `continuous_synthesis.md:6-8`
said "CSP is the only synthesis runner." Reworded: CSP is now framed as the **explicit
coarse-grained ribosome** runner, dropped the false claim, and added a pointer to its
analytic-tunnel counterpart ({doc}`cylinder_synthesis`, `topo-cylinder`). Docs rebuilt clean.

### LOW

**D8 — Wrong `strtobool` attribution.** ✅ **FIXED (2026-07-03).** `simulation_control.rst:294`
cited `distutils.util.strtobool`; code uses `topo.utils.config.strtobool` (distutils was removed).
Updated the cross-reference to `:func:topo.utils.config.strtobool`; docs rebuilt clean.

**D9 — `nstout`/`device`/`ppn` shown as "—" (no default).** ✅ **FIXED (2026-07-03).**
`continuous_synthesis.md:364-366` showed no defaults. Filled in `device = CPU`, `ppn = 1`, and
`nstout = 5000`. **Note:** the `RunParams.nstout` code default was also changed **50 → 5000**
(`core.py:828`) per user request, and the other two doc references updated to match
(`synthesis_control.rst:251`, `cylinder_synthesis.md:208`). Docs rebuilt clean; both runners
(`RunParams`/`CylinderParams`) confirmed to report `nstout = 5000`.

**D10 — Incomplete Part-B reference list.** ✅ **FIXED (2026-07-03).** `overview.rst:76-80` Part-B
Reference rubric omitted `synthesis_control` and `cylinder_synthesis`. Added both (cylinder page +
`Synthesis control options (csp.ini)`); docs rebuilt clean.

**Docs verified clean:** toctree (every entry resolves, no orphans, all 3 synthesis pages linked);
tutorial include-stubs 01–08 all resolve; the entire `codon_times` surface (synthesis_control,
cylinder_synthesis, csp/README) matches the parsers with no live references to the removed keys
(`trans_times`/`uniform_ta`/`uniform_mfpt`); defaults spot-checked correct (`scale_factor`,
`time_stage_1/2`, `restraint_k`, `buffer`, `ptc_offset`, `nascent_ev_radii`, `trna_tether`
INI-default, cylinder `tunnel_*`, `post_elongation`).

---

## Source-code findings

### MEDIUM

**S1 — `checkBondDistances` rejects a valid nucleic bond.** ✅ **FIXED (2026-07-03).**
`topo/core/system.py:832` used `>=` with default `threshold = 0.5` nm while
`bond_length_nucleic = 0.5`, so an RNA/nucleic bond built at equilibrium *equals* the threshold and
raised on a valid structure. Changed the comparison to strict `>` (a bond exactly at the threshold
is now allowed; only strictly-longer bonds raise) and replaced the stale `//TODO` at 816-817 with a
note explaining the strict comparison. Verified: a 0.5 nm bond no longer triggers the raise, a
0.6 nm bond still does; called at build time (`system.py:792`); module imports clean.

**S2 — `pbc` flag disagrees with box state on malformed input.** ✅ **FIXED (2026-07-03).**
`topo/utils/config.py`: when `pbc = yes` but `box_dimension` wouldn't parse, it set
`box_dimension = None` but left `cfg.pbc = True`, so `enforcePeriodicBox` acted over a system with
no box vectors and the `pcoupl` `assert cfg.pbc` passed spuriously. Added `cfg.pbc = False` in the
invalid-box branch (with an explanatory comment). Verified: invalid box → `pbc=False`/`box=None`;
valid box → `True`/`[30,30,30]`; and `pcoupl=yes` + invalid box now correctly raises
`AssertionError` (the `pcoupl` assert at config.py:368 runs after the pbc block, so it sees the
corrected flag).

**S3 — Q-analysis and the force field use different contact sets.** ✅ **RESOLVED — DOCUMENTED
(2026-07-03).** Decision (per user): this is **intentional**, not a bug — `|i−j| ≥ 3`
(`LOCAL_SEPARATION = 2`) is the **model parameterization** (which native pairs carry energy), while
the Q-analysis default `--local-separation 3` (`|i−j| ≥ 4`) is a deliberate **analysis convention**
for the folding metric. No behavior change. Documented the difference in both places: the
`--local-separation` argparse default in `native_contacts.py` now carries an explanatory comment +
help text (noting `--local-separation 2` matches the model), and the `docs/usage/native_contacts.rst`
note was expanded to state the model-vs-convention framing and the actionable knob.
*(Also fixed an unrelated Sphinx substitution warning introduced by the S4 docstring edit: `|i-j|`
→ `abs(i - j)`, since RST reads `|…|` as a substitution reference. Docs now build with zero
warnings.)*

### LOW

**S4 — Stale class-attribute docstring.** ✅ **FIXED (2026-07-03).** `topo/core/system.py:60-62`
said `bonded_exclusions_index … = 3`; the actual value is `2` (set from `model_parameters.py:46`,
consumed at system.py:675/753 — correct and consistent with `LOCAL_SEPARATION`). Docstring-only
fix (no runtime effect): corrected to `=2` and expanded to note it excludes 1-2/1-3 neighbours,
preserving native contacts at `|i-j| >= 3`. Module imports clean.

**S5 — `checkLargeForces` skips its check when `minimize=False`.** ✅ **FIXED (2026-07-03), option
(b) — enforce.** Previously, with `minimize=False` the method only reported energy and never
compared any force to `threshold`, so an exploding/clashing input passed silently. Restructured so
the max per-particle force is computed unconditionally; `minimize=True` keeps the iterative
minimization, and `minimize=False` now emits a `RuntimeWarning` ("Largest force … above the
threshold … minimize=False so the structure was not relaxed …") when a force exceeds `threshold`.
Chose a **warning, not a raise**, so the diagnostic-only callers (build/restart at `system.py:810`)
aren't broken. Docstring updated to describe both paths. Module imports clean; structure verified.

**S6 — Unclosed file handle.** ✅ **FIXED (2026-07-03).** `topo/core/system.py` `dumpStructure`
passed `file=open(output_file,'w')` to `structure.writeFile`, which doesn't close it (fragile off
CPython refcounting). Wrapped in `with open(output_file, 'w') as ff:` so the handle is deterministically
closed/flushed. Confirmed it's the only such `file=open(...)` pattern in the module; import clean.

**S7 — Dead guard.** ✅ **FIXED (2026-07-03).** `topo/core/system.py:1064` `if self.atoms !=
OrderedDict():` was always True (`self.atoms = []` is a `list`, never equal to an `OrderedDict`).
Changed to `if self.atoms:` (the intended non-empty check). The bonds/angles/torsions guards at
1091/1110/1121 legitimately use `!= OrderedDict()` (those really are `OrderedDict`s) and were left
as-is, so the `OrderedDict` import stays used. In `dumpForceFieldData` (diagnostic writer), so
behavior only changes for the never-hit empty-atoms case. Module imports clean.

**S8 — Cryptic error on empty domain.** ✅ **FIXED (2026-07-03).**
`topo/analysis/native_contacts.py`: an empty residue list for a domain in `domain.yaml`
(`residues: []`) made `idx` an empty array, so `idx.min()`/`idx.max()` raised the cryptic
`zero-size array to reduction operation minimum which has no identity`. Added an `idx.size == 0`
guard *before* the bounds check that raises `"<domain_yaml>: domain 'X' has no residues listed"`.
Verified with a crafted empty-domain YAML — now emits the clear, actionable message.

**S9 — API/file naming lag.** ✅ **FIXED (2026-07-03) — full rename for systematic consistency.**
The file-path layer now speaks `codon_time_table_path` everywhere, matching the `codon_times` INI
key and the `codon_time_table` dict:
- `run_continuous_synthesis(..., codon_time_table_path=…)` and
  `run_cylinder_synthesis(..., codon_time_table_path=…)`
- `CSPConfig.codon_time_table_path` / `CylinderConfig.codon_time_table_path`
- `build_codon_time_lists(..., codon_time_table_path=…)`, `ribosome_traffic_times(..., codon_time_table_path=…)`
- `default_trans_times_path()` → `default_codon_time_table_path()`;
  `DEFAULT_TRANS_TIMES_FILE` → `DEFAULT_CODON_TIME_TABLE_FILE`

Intentionally **kept** (these name the literal on-disk artifact / old INI keys, not API):
the bundled data file `ecoli_trans_times_310K.txt` (constant *value* unchanged), example filenames
like `trans_times.txt`, and the legacy-guard key strings `("trans_times", "uniform_ta",
"uniform_mfpt")`. Docs/README examples updated to the new param name. Verified: both runners'
signatures + config fields use `codon_time_table_path`, the legacy `trans_times` INI key still
errors, `default_codon_time_table_path()` resolves the bundled table, zero non-keep-case
`trans_times` remain anywhere, docs build clean.

**Source verified clean:** `mdrun` (quench/production step accounting + restart resume),
`reporter`, `utils/multichain` (force replication + interaction-group isolation + streaming DCD
split), `runinfo`, `read_parms`, `parameters/dihedral` + `model_parameters` (0.756 scaling applied
once), `geometry`, the non-bonded force math (12-10-6 minimum, Yukawa, Gaussian angles), and
`optimize` ladder-climb/convergence. **CSP module:** the `codon_times` refactor is fully
consistent — no stale `uniform_ta` / `uniform_mfpt` / `codon_mfpt_list` / `build_mfpt_lists`
identifiers; imports OK; `parse_codon_times`, `stage_dwell_times`, `build_codon_time_lists` and the
ribosome-traffic guard are correct; the legacy-key error fires on old configs.

---

## Recommended fix order

**✅ ALL DOCUMENTATION FINDINGS (D1–D10) FIXED — 2026-07-03.**

Remaining (source code):

**✅ ALL FINDINGS RESOLVED — 2026-07-03.** Docs D1–D10 fixed; source S1, S2 bugs fixed; S4–S8
cleanups fixed; S3 documented as intentional; S9 API rename completed for full `codon_times`
naming symmetry. Nothing outstanding.
2. **S9** — optional API rename, only if you want full `codon_times` symmetry.
