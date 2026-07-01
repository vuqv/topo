# CHANGELOG — Standalone ribosome migration + naming/docs cleanup

Branch: `tut15-claude-fix`. This session made the CSP ribosome **standalone** (no O'Brien
`.cor/.psf/.prm` at runtime), standardized non-bonded parameter naming, removed the O'Brien
loaders, revalidated 4c5c, and brought docstrings/docs to 100% coverage. All changes below are
**numerically identical** to before unless noted (they are renames / plumbing / removals).

Not yet committed — everything is in the working tree. Resume points are in **§8 Open items**.

---

## 1. Standalone ribosome geometry (C5′ ribose bead)

**Root cause found:** topo's CG ribose (`R`) bead was the centroid of the 5-atom ribose ring
`C1'–C4' + O4'`; O'Brien uses the centroid of the five ribose **carbons** `C1'–C5'` (no O4′)
— confirmed in his `create_cg_ribosome_model.py:326`. On the **same** all-atom structure, the
5-carbon set reproduces his `.cor` R beads to **median 0.0000 Å** (mean 0.003), vs 0.485 Å for
the old set. `C5'` *alone* was the worst (2.15 Å) — the fix is the atom **set**, not one atom.

- `topo/csp/cg_ribosome.py` — `RIBOSE_RING` now `["C1'","C2'","C3'","C4'","C5'"]`.
- Regenerated `topo/csp/structures/4v9d_50S_PtR_5jte_AtR_model_cg.pdb` and `…_cg_trunc.pdb`
  (via `cg_ribosome` + `truncate_ribosome`, `--r-cyl 30 --x-lo -8 --x-exit 58`).
- Result: the topo CG truncated ribosome (4575 beads) reproduces O'Brien's `.cor` (4577) to
  **median 0.0003 Å** for coords, **exact** for RNA radii/charges. Only exception: **PtR:76
  A76** (the P-site peptidyl-tRNA acceptor) is 3.45 Å off — a genuine source-structure
  difference; **kept as topo's** per decision (PTC-seeding absorbs it). See §8.

## 2. Parameters folded into `model_parameters` (single source of truth)

- `model_parameters["topo"]` key **`radii` → `Rmin_2`** (it was always Rmin/2, never a σ/diameter).
- RNA `P/R/BR` `Rmin_2` set to O'Brien per-type: **0.644766 / 0.523140 / 0.534244 nm**.
- Protein `Rmin_2` set to O'Brien per-AA sidechain values (== `OBRIEN_SC_RMIN_2_NM`).
- These are the **fixed per-type radii for rigid ribosome scenery** (RNA + ribosomal protein),
  read by `load_ribosome`. `OBRIEN_SC_RMIN_2_NM` (`ribosome.py`) is now **derived** from
  `model_parameters` (+ HSD/HSE/HSP His aliases) — no duplicated numbers.

**Key semantic (documented in `model_parameters.py`):** a mobile protein chain (nascent, or a
folded-protein sim) does **not** use these — its Rmin/2 is the **per-residue, structure-derived
Karanicolas–Brooks** value. So `ALA`-in-nascent = K-B (varies per residue) while `ALA`-in-
ribosome = fixed table. This mirrors O'Brien, who separates them by atom type (`A_i` vs `S<aa>`)
while keeping standard residue names. Verified: `load_ribosome` reproduces O'Brien's ribosome
radii (RNA exact; protein exact except L24) and **charges exact** (0 mismatches / 4454 beads).

## 3. Nascent radius sourced from the K-B array (not the fixed table)

`rf_sigma` (per-particle excluded-volume radius) only feeds `dumpForceFieldData` (never called),
so this is a clarity/consistency change with no simulation effect:
- `csp/core.py build_length_model` takes `nascent_rmin_2` and calls `setParticlesRadii(...)`
  with the K-B array (Option A) instead of `setCARadiusPerResidueType()`.
- `core/models.py` (general builder) requests `build_nonbonded_interaction(return_rmin_2=True)`
  and sets `rf_sigma` from the K-B array.
- `setCARadiusPerResidueType` kept but documented as fixed-per-type/scenery-only (no longer
  called by any builder).

## 4. Non-bonded naming standardized

- `calculate_sigma_values` → **`calculate_rmin_2_values`** — now returns **Rmin/2** directly
  (was the collision *diameter* Rmin; the `0.5×` moved inside). Numerically identical.
- returned matrix `distance_matrix` → **`rmin_matrix`** (per-pair well position `R_ij`: native
  distance for contacts, `Rmin/2_i + Rmin/2_j` sum rule for non-native). Also the
  `topo_model.distance_matrix` attribute and `addCustomNonBondedForce` param.
- **`rmin2` → `rmin_2`** everywhere (identifier `_2` reads as subscript, not squared):
  `return_rmin_2`, `nascent_rmin_2`, `OBRIEN_SC_RMIN_2_NM`, `rmin_2_full`, `aa_rmin_2_nm`, …
  (80 replacements across 10 files incl. tutorials).
- Vocabulary: **`Rmin_2`/`rmin_2`** = per-particle collision radius; **`rmin_matrix`** = per-pair
  well position; `rmin_pair = Rmin_2_i + Rmin_2_j`.

## 5. O'Brien loaders removed (fully standalone)

- Removed **`load_obrien_ribosome`**, `_parse_prm_radii`, `_OBRIEN_RNA_NAME` (`.cor/.psf/.prm`).
- Removed **`load_ribosome_auto`** (redundant once `.cor` was gone).
- Removed `ribosome_psf`/`ribosome_prm` plumbing from `protocol.py` (signature, `CSPConfig`,
  `read_csp_config`, pass-throughs).
- Single loader now: **`load_ribosome(pdb, model="topo")`**. All callers updated (protocol,
  standalone/tut15/P0CX28 `analyze_validation.py` + `eject_demo.py`).

## 6. tut15 / P0CX28 repointed to the standalone ribosome

- INIs (`csp.ini`, `csp_val.ini`, `csp_tether.ini`, P0CX28 `csp.ini`/`csp_val.ini`) and
  `analyze_validation.py`: `ribosome_obrien.cor` → **`ribosome_trunc.pdb`**.
- Refreshed the stale local `ribosome_trunc.pdb` copies (were pre-C5′) with the current
  canonical `topo/csp/structures/…_cg_trunc.pdb`.

## 7. Docstrings + docs

- **100% docstring coverage** (186/186 functions/methods, all classes). Added numpy-style
  docstrings to the 16 that were missing.
- Fixed 2 malformed docstrings in `core/system.py` (`getCAlphaOnly` pseudo-code literal block;
  `system` class Parameters/Attributes indentation).
- Rebuilt Sphinx docs (`docs/`, furo + napoleon): **build succeeds**, HTML in `docs/_build/html`.
  Now also covers `csp`/`mdrun`/`optimize` (new flat stubs `docs/topo.{csp,mdrun,optimize}.rst`).
  See §8 for the remaining duplicate-object warnings.

---

## Verification

- **Debug smoke (L=1→8)** on `tutorials/16_csp_standalone/4c5c/` after every major step: `Done 1→8`,
  standalone `ribosome_trunc.pdb` (4575 beads via `load_ribosome`), fix active (|A−P|=0.3810 nm),
  **0** stability/dt-halving/leak/clash/NaN.
- **Full 4c5c validation (L=1→306 + ejection)** — see `tutorials/16_csp_standalone/4c5c/analyze_standalone.py`:

  | Check | Standalone | tut15 baseline (O'Brien `.cor`) |
  |---|---|---|
  | completes, 0 dt-halving | ✓ | ✓ |
  | worst \|PotE\| | 1.05e3 kJ/mol | 1.48e3 |
  | chain threading corr(idx,x) | −0.870 | −0.926 |
  | wall (min x ≥ x0) | 10.4 ≥ 9.1 Å PASS | 8.37 ≥ 8.71 PASS |
  | clash (min NC–ribo) | 3.36–4.82 Å | 2.2–2.9 |
  | D6 dwell (in-vivo total) | 1.01× ref | 1.01× |

  Differences trace to the two accepted deviations (topo A76 P-anchor shifts the wall/PTC;
  L24 per-AA approximation). Energy/clash slightly better; dwell matches the O'Brien reference.
- Non-bonded numerics unchanged across all renames: native contacts 819, K-B Rmin/2 range
  0.251–0.531 nm, `rmin_matrix[i,j] = rmin_2[i]+rmin_2[j]` verified.

---

## 8. Open items (resume here)

1. **Commit** — nothing committed yet; working tree only. Suggested logical commits: (a) C5′ +
   regenerated structures; (b) model_parameters Rmin_2; (c) K-B nascent radius wiring; (d) naming
   renames; (e) remove O'Brien loaders; (f) repoint tut15/P0CX28; (g) docstrings + docs.
2. **Docs duplicate warnings (160)** — all "duplicate object description": the docs carry BOTH
   flat apidoc stubs (`docs/topo.*.rst`) AND curated pages (`docs/modules/*.rst`:
   system/parameters/models/topo.reporter) for the same modules. Pre-existing structural overlap;
   build still succeeds. Fix options: drop the overlapping curated pages, or `:no-index:` one set,
   or `napoleon_use_ivar = True` (kills the class-attribute half). Not yet decided.
3. **PtR:76 A76 P-anchor** — kept as topo's (3.45 Å from O'Brien). PTC-seeding absorbs it; if a
   bit-exact P-site is ever wanted, splice O'Brien's single A76 R coord into the CG PDB.
4. **L24** — approximated per-AA (like other ribosomal proteins) instead of O'Brien's per-residue
   B-types (102 beads). Accepted; mean protein Rmin/2 error 0.013 nm.
5. **Egress D5b** — the in-run 20k-step ejection is too short for a 306-mer; run `eject_demo.py`
   for the extended egress (apples-to-apples with tut15's PASS).
6. **`interaction_details.md`** (`tutorials/15_claude_fix/`) documents the NC↔ribosome model
   (single CustomNonbondedForce, sum rule, 12-10-6, nascent K-B vs ribosome per-AA/per-type).
7. **Performance note** — CSP throughput is CPU-bound (per-stage OpenMM context rebuild + JIT),
   not GPU-bound. ~900 steps/s effective on a contended node; a node with free CPU cores is the
   win, not a faster GPU.

## Key paths

- Standalone run dir: `tutorials/16_csp_standalone/4c5c/` (`csp_debug.ini` L=1→8, `csp_val.ini` L=1→306,
  `ribosome_trunc.pdb`, `analyze_standalone.py`).
- Ribosome pipeline: `topo/csp/cg_ribosome.py`, `topo/csp/truncate_ribosome.py`,
  structures in `topo/csp/structures/`.
- Loader/params: `topo/csp/ribosome.py` (`load_ribosome`, `OBRIEN_SC_RMIN_2_NM`),
  `topo/parameters/model_parameters.py`, `topo/utils/nonbonded.py`.
