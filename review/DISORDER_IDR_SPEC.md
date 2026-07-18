# Disordered / IDR regions in the Cα model — implementation spec

> **Date:** 2026-07-17
> **Reference (read-only):**
> `sandbox/original_code/perl/cg_model/create_cg_protein_model_v34_0.03.pl` (v1.34, O'Brien).
> **topo code touched:** `topo/core/models.py` (`buildCoarseGrainModel`),
> `topo/utils/nonbonded.py` (`build_nonbonded_interaction`), `topo/utils/config.py`,
> `topo/csp/core.py`.
> **Scope:** the α-carbon (Cα) model only. Proposes a minimal implementation for marking
> protein regions as disordered/unstructured (IDRs). Companion:
> [`DIFFERENCES.md`](DIFFERENCES.md) (CSP-side deviations).

---

## TL;DR

The original defines disorder in the Cα model through two cooperating, per-segment
mechanisms:
1. **Non-Go local backbone** (3.81 Å bonds, double-well angle, Karanicolas dihedral).
   **topo already applies all of this globally** — so an IDR in topo needs nothing extra
   here.
2. **`GENERIC` non-bonded potential** — the segment gets **no STRIDE and no native/Go
   contacts**; instead every non-local pair gets a **weak, uniform, non-specific**
   attraction, with excluded volume from `2·rvdw`.

Because topo already supplies the generic backbone and self-avoidance, **the whole IDR
feature reduces to per-residue contact masking inside `build_nonbonded_interaction`**
(optionally plus the weak generic attraction). No changes to bonds, angles, dihedrals, or
the OpenMM force objects are required. This mirrors the existing `ribo_free_mask` pattern.

---

## 1. How the original builds an unstructured/IDR segment (Cα)

The Perl builds a mixed folded+disordered protein by **stitching multiple PDB segments**
N→C (`@pdb`, `@pot`, per-segment `@dihedral_go`, `@bondlength_go`, `@angle_dw`), bonded at
3.81 Å at each junction. A segment is made disordered by:

**(a) Non-Go local backbone** — set for that segment:
- `bondlength_go = 0` → Cα–Cα bonds = 3.81 Å (line 1418),
- `angle_dw = 1` → double-well backbone angle (line 1445),
- `dihedral_go = 0` → Karanicolas transferable dihedral from `karanicolas_dihe_parm.dat`,
  scaled 0.756 (lines 682, 1531–1541).

**(b) `GENERIC` non-bonded** — `pot = generic-bt` (or `-mj`, `-kgs`):
- **No STRIDE, no native contacts** for the segment (line 1121 guards STRIDE on
  `pot !~ GENERIC`; the whole native-contact block is skipped).
- **Excluded volume**: `sigmin[i] = 2·rvdw(resname_i)` (line 1638) instead of the K–B
  nearest-neighbour distance.
- **Weak uniform attraction** in `NBFIX`, for every non-local pair |i−j| ≥ 3
  (lines 1866–1874, `casm = 0` branch):

  ```
  ε_ij  = (0.3 / 10) · ε_BT(res_i, res_j)      # = 0.3·nscal scaled, /10 CHARMM factor
  r_ij  = rvdw(res_i) + rvdw(res_j)            # sum of side-chain vdW radii
  ```

  i.e. a shallow, **non-specific** attraction between all residues in the segment — a
  weakly collapsing self-avoiding chain, *not* a Go fold.
- Residue/atom naming switches to `U<n>` / `UA<n>` (unstructured) vs `G<n>` for structured
  (lines 1229–1235); cross-segment interactions optionally via `interact_scale`.

The `%rvdw` side-chain van-der-Waals radii table (Perl lines 720–745) is the only extra
parameter needed for (b); it is **not currently shipped in topo**.

> **What topo already has (so this feature is small).** topo's global backbone is exactly
> the original's non-Go choice: the **double-well angle** (`angle_dw = 1`; constants
> 106.4/91.7/26.3/130.0/0.1/4.3) and the **Karanicolas sequence-dependent dihedral**
> (`dihedral_go = 0`, 0.756 scale). The K–B excluded-volume term already gives
> self-avoidance. Only the Go-contact removal (1b) is missing.

---

## 2. Proposed implementation for topo

### 2.1 Two levels (pick per study)

- **Level A — self-avoiding IDR (recommended default).** For any residue in the disorder
  mask, **drop all native contacts** (H-bond, BS, SS) in which it participates; those
  pairs fall back to the existing non-native repulsive/excluded-volume term. Result: a
  generic self-avoiding chain with transferable backbone — the cleanest "no fold"
  baseline, and *zero* new parameters.

- **Level B — weakly-collapsing IDR (O'Brien `generic-bt` parity).** In addition to A,
  add the shallow uniform attraction `ε = 0.3·nscale·ε_BT(i,j)` at `r = rvdw_i + rvdw_j`
  for pairs where **both** residues are in the disorder mask (and, optionally, for
  disorder↔folded cross pairs via a separate scale). Requires shipping the `%rvdw` table.

> **Why masking (not a `nscale = 0` domain) is the right lever.** The existing
> `domain.yaml` scaling multiplies **only the SS energy** (`scaling_matrix * ss_*`).
> Setting a domain's `nscale = 0` would zero SS contacts but **leave H-bond and BS
> contacts intact**, so the region would still fold via backbone H-bonds. A true IDR must
> remove H-bond + BS + SS together — hence a dedicated mask applied to all three, not a
> domain scale.

### 2.2 API / plumbing

New optional argument threaded exactly like `domain_def`:

- `topo/utils/nonbonded.py`:
  `build_nonbonded_interaction(pdb_file, domain_def=None, stride_output_file=None,
  disorder_def=None, *, generic_attraction=False, return_rmin_2=False)`.
- `topo/core/models.py::buildCoarseGrainModel(..., disorder_def=None)` → forwards to the
  call above.
- `topo/utils/config.py`: add `disorder_def: Optional[str] = None` to the config
  dataclass + `build_kwargs()` (guarded like `domain_def`, lines 154–155) and read from
  the INI (`disorder_def = ...`, mirroring line 485). Also expose a `generic_attraction`
  boolean toggle (default `no`).
- `topo/csp/core.py::build_contacts_from_pdb` (line ~374): accept and forward
  `disorder_def` so continuous-synthesis nascent chains can mark regions disordered too.

**Disorder-mask file format** — reuse the residue-range syntax already parsed by
`parse_residue_list` (accepts ints, `"i"`, and `"start-end"`). Minimal YAML:

```yaml
# disorder.yaml
n_residues: 283
disordered:
  - 1-24          # N-terminal IDR
  - 150-165       # internal loop
# optional (Level B): scale for disorder<->folded cross attraction (default 0.0 = none)
cross_scale: 0.0
generic_scale: 0.3   # multiplies epsilon_BT for disorder<->disorder generic attraction
```

A bare newline/whitespace list of ranges would also be acceptable; YAML keeps it uniform
with `domain.yaml`.

### 2.3 Algorithm (inside `build_nonbonded_interaction`)

After the existing `eps_ij` / `binary_contact_matrix` are assembled and **before** the
non-native fill loop:

```python
if disorder_def is not None:
    dis_res = parse_disorder(disorder_def)                 # -> set of 1-based resids
    dis_idx = np.array([resid_to_index[k] for k in keys    # map to 0-based matrix idx
                        if k[1] in dis_res])
    m = np.zeros(n_residues, bool); m[dis_idx] = True
    involves_dis = m[:, None] | m[None, :]                 # pair touches a disordered res

    # Level A: remove native contacts touching any disordered residue.
    eps_ij[involves_dis] = 0.0
    binary_contact_matrix[involves_dis] = 0                # -> non-native / excluded-vol

    # Level B (optional): shallow uniform attraction among disordered residues.
    if generic_attraction:
        dd = m[:, None] & m[None, :]                       # both disordered
        # + optional cross term via cross_scale (disorder<->folded)
        # eps set to generic_scale * eps_BT; rmin set to rvdw_i + rvdw_j (needs %rvdw table)
```

The subsequent non-native loop (`binary_contact_matrix[i,j] == 0 → NON_NATIVE_KJ`,
sum-rule Rmin) then automatically handles every de-contacted pair for Level A. Level B
overwrites those specific entries with the weak well afterwards.

> **Rmin/2 subtlety for Level A.** `calculate_rmin_2_values` derives each residue's
> collision radius from the nearest **non-contact** Cα distance. Since masking *increases*
> the number of non-contact pairs, it must run on the **post-mask**
> `binary_contact_matrix` (already the case if masking precedes it). No extra work — just
> ordering.

### 2.4 What does *not* change (important)

- Bonds, angles, dihedrals, Yukawa: untouched — topo's global transferable backbone is
  already the disordered-appropriate choice (§1).
- OpenMM force construction (`addCustomNonBondedForce`): unchanged; it consumes the same
  two matrices.
- Folded-only runs (`disorder_def=None`): byte-for-byte identical to today.

### 2.5 New parameter to ship (Level B only)

Add `topo/parameters/data/rvdw.csv` = the Perl `%rvdw` table (21 rows, Å). Loaded like
`bt_potential.csv` (path-resolved under `parameters/data/`). Not needed for Level A.

---

## 3. Validation plan

1. **Regression (no-op):** a folded PDB with `disorder_def=None` and with an *empty*
   disorder list must produce identical `rmin_matrix` / `energy_matrix`. → `np.allclose`.
2. **Contact removal (Level A):** mask an α-helix's residues; assert every native contact
   (H-bond/BS/SS) touching a masked residue is gone from `energy_matrix`, and those pairs
   now carry `NON_NATIVE_KJ`. Assert non-masked–non-masked contacts are unchanged.
3. **Behavioural:** short MD (a small folded domain + a masked N-terminal tail). Confirm
   energies finite, the tail's radius of gyration expands relative to the fully-Go run,
   and the folded core RMSD stays low. Use the Tutorial-13 short-debug pattern.
4. **Level B parity (optional):** build a fully-`GENERIC-bt` chain in the Perl and the
   same chain fully-masked with `generic_attraction=yes` in topo; compare the per-pair
   ε and Rmin tables (expect match up to the CHARMM `/10` and unit conventions).

---

## 4. Open questions for the user

1. **Level A vs B default.** Is the intended IDR a *self-avoiding* chain (Level A, no new
   params) or a *weakly-collapsing* chain (Level B, needs `%rvdw` + generic attraction)?
   Recommendation: ship **Level A** first (small, parameter-free), add B behind the
   `generic_attraction` toggle if a study needs collapse.
2. **Disorder↔folded interactions.** When a masked tail passes near the folded core,
   should it feel (i) only excluded volume (Level A), or (ii) the weak generic attraction
   (Level B `cross_scale > 0`)? O'Brien's `interact_scale` is the analog.
3. **Input geometry.** Masking a region of a *folded* PDB keeps that region's native
   (folded) starting coordinates. If you want a genuinely extended IDR start, we should
   also generate an extended-chain segment for those residues (analog of the Perl
   `create_unstructured_pdb.pl`). Is an extended-start needed, or is masking-in-place
   enough for the intended use?
4. **Scope.** Single mask file per system (this spec), or per-chain/per-domain masks?
