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
   contacts**; instead every non-local pair gets a **weak, non-specific** attraction
   (`0.03·nscal·ε_BT(i,j)` — the depth still varies by residue-pair type), with excluded
   volume from `2·rvdw`.

Because topo already supplies the generic backbone and self-avoidance, **the whole IDR
feature reduces to per-residue contact masking inside `build_nonbonded_interaction`**
(optionally plus the weak generic attraction). No changes to bonds, angles, dihedrals, or
the OpenMM force objects are required. This mirrors the existing `ribo_free_mask` pattern.

---

## 1. How the original builds an unstructured/IDR segment (Cα)

The Perl builds a mixed folded+disordered protein by **stitching multiple PDB segments**
N→C (`@pdb`, `@pot`, per-segment `@dihedral_go`, `@bondlength_go`, `@angle_dw`), bonded at
3.81 Å at each junction. "Which region is IDR" is chosen **per PDB file** — there is no
residue-range keyword; the granularity is the whole segment (see §1.3).

### 1.1 Structured vs IDR differ on a single axis: `pot`

The two model choices are **orthogonal**:

- **Local backbone geometry** — `dihedral_go`, `angle_dw`, `bondlength_go`.
- **Tertiary structure (native contacts)** — `pot`.

In the O'Brien production configuration that topo is ported from (e.g. the 1GB1 control
file: `casm=0`, `pot=bt`, `dihedral_go=0`, `angle_dw=1`, `bondlength_go=0`), the
**transferable backbone is used for *every* segment, folded or not**. So when you compare
a structured segment to a disordered one *in the real config*, **only `pot` changes**:

| Keyword | Structured segment (1GB1 config) | Disordered (IDR) segment |
|--|--|--|
| `pot` | `bt` | **`generic-bt`** |
| `dihedral_go` | `0` (Karanicolas transferable) | `0` (Karanicolas transferable) |
| `bondlength_go` | `0` (3.81 Å) | `0` (3.81 Å) |
| `angle_dw` | `1` (double-well transferable) | `1` (double-well transferable) |
| `casm` | `0` (Cα-only) | `0` (Cα-only) |

> **Caution — not the code defaults.** The Perl's *built-in* defaults are the opposite
> (`dihedral_go=1` Go, `angle_dw=0` Go-harmonic, `bondlength_go=1` native); O'Brien
> overrides all three globally. The table above reflects how the model is **actually run**,
> not the code defaults.

### 1.2 What `pot = bt → generic-bt` flips (the only real difference)

| Behavior | `pot=bt` (structured) | `pot=generic-bt` (IDR) |
|--|--|--|
| STRIDE / H-bonds | run; native H-bonds | **skipped** (line 1121) |
| Native SS / BS contacts | built from structure | **none** (line 1028) |
| Non-local pair interaction | Go wells at native Cα distance | **weak, non-specific** attraction on *every* pair (below) |
| Excluded volume | K–B nearest-neighbour | `2·rvdw(resname)` (line 1638) |

The weak attraction, written to `NBFIX` for every non-local pair |i−j| ≥ 3
(lines 1866–1874, `casm = 0` branch):

```
ε_ij  = (0.3 / 10) · ε_BT(res_i, res_j)      # = 0.03 · nscal · ε_BT  (nscal is inside eps)
r_ij  = rvdw(res_i) + rvdw(res_j)            # sum of side-chain vdW radii
```

> **Effective factor is 0.03, not 0.3.** The literal prefactor is `0.3/10 = 0.03`. A
> **native** SS contact of the same residue pair is written as `ε_BT` with *no* prefactor
> (line 1762), so a generic pair attracts at **3% of a native contact**, not 30%. The
> author's comment (line 104) says "0.3*nscal" but the code applies the extra `/10`. When
> porting to topo, **match the code: use 0.03.** The `/10` is not a unit conversion (both
> NBFIX entries are kcal/mol) — it is a genuine, generic-term-only 10× weakening.

It is **non-specific in coverage** (it acts on *every* non-local pair, not just native
contacts) with a **uniform 0.03·nscal prefactor**, but its **depth is pair-type dependent**
via `ε_BT(res_i, res_j)` — so the well is chemically heterogeneous, not the same for all
pairs. Physically: a weak, sequence-modulated, non-fold-encoding attraction — a weakly
collapsing self-avoiding chain, *not* a Go fold. The `%rvdw` side-chain vdW-radii table
(Perl lines 720–745) is the only extra parameter it needs; it is **not currently shipped in
topo**.

### 1.3 Region granularity

The Perl keys disorder on the **PDB file**: it iterates `for($np=0;$np<=$#pdb;$np++)` and
tests `uc($pot[$np]) =~ m/GENERIC/`. So an IDR is a *separate input PDB* (Cα coords only,
generated by the companion `create_unstructured_pdb.pl`) tagged `pot=generic-*`. An
internal IDR (a disordered loop between two folded domains) therefore requires **three**
PDBs (folded-N, loop, folded-C). Residue/atom names switch to `U<n>`/`UA<n>` (unstructured)
vs `G<n>` (structured) (lines 1229–1235); cross-segment interactions optionally via
`interact_scale`.

> **What topo already has (so this feature is small).** Because structured and IDR differ
> *only* in `pot` (native contacts), and topo already applies the transferable backbone
> globally — the **double-well angle** (`angle_dw = 1`; 106.4/91.7/26.3/130.0/0.1/4.3) and
> the **Karanicolas dihedral** (`dihedral_go = 0`, 0.756 scale) — with K–B excluded volume
> for self-avoidance, the *only* thing left to localize is **native-contact removal**. And
> because topo reads one structure, it can localize by **residue range** (a mask) rather
> than by splitting PDBs — so an internal IDR is one `disordered:` entry, not three PDBs.

---

## 2. Proposed implementation for topo

### 2.1 Two levels (pick per study)

- **Level A — self-avoiding IDR (recommended default).** For any residue in the disorder
  mask, **drop all native contacts** (H-bond, BS, SS) in which it participates; those
  pairs fall back to the existing non-native repulsive/excluded-volume term. Result: a
  generic self-avoiding chain with transferable backbone — the cleanest "no fold"
  baseline, and *zero* new parameters.

- **Level B — weakly-collapsing IDR (O'Brien `generic-bt` parity).** In addition to A,
  add the shallow, pair-type-dependent attraction `ε = generic_scale·nscale·ε_BT(i,j)` at
  `r = rvdw_i + rvdw_j` for pairs where **both** residues are in the disorder mask (and,
  optionally, for disorder↔folded cross pairs via `cross_scale`). **O'Brien parity =
  `generic_scale = 0.03`** (the effective code factor `0.3/10`; see §1.2 — the "0.3" in
  his comment omits the `/10`). `generic_scale` is defined as the true depth ratio to a
  native contact of the same pair, so 0.03 reads as "3% of a native contact." Requires
  shipping the `%rvdw` table.

> **Why masking (not a `nscale = 0` domain) is the right lever.** The existing
> `domain.yaml` scaling multiplies **only the SS energy** (`scaling_matrix * ss_*`).
> Setting a domain's `nscale = 0` would zero SS contacts but **leave H-bond and BS
> contacts intact**, so the region would still fold via backbone H-bonds. A true IDR must
> remove H-bond + BS + SS together — hence a dedicated mask applied to all three, not a
> domain scale.

### 2.2 API / plumbing — one file, no new knobs

> **Naming note (deferred).** The key/arg stays `domain_def` for now. Once the file also
> carries a `disordered:` section it is really a *protein-architecture* definition (folded
> domains + interfaces + disordered regions), so a clearer name (evoking the protein's
> architecture/layout) is planned — to be decided later. `domain_def` will be kept as a
> deprecated alias when that rename happens, so nothing here breaks.

Disorder is defined as an **optional `disordered:` section inside the existing
`domain_def` YAML** — *not* a separate file or argument. `domain_def` is already threaded
end-to-end (INI key → config dataclass → `buildCoarseGrainModel` →
`build_nonbonded_interaction`), so this adds **no new input file, no new function
argument, and no new INI key**. The only code changes are the YAML reader and the masking
step:

- `topo/utils/nonbonded.py`: `read_yaml_config` also parses the optional `disordered:`
  section (returning the residue set + `generic_scale` + `cross_scale`);
  `build_nonbonded_interaction` applies the mask (§2.3). Signature is unchanged except the
  already-planned `return_rmin_2` — **no** `disorder_def` / `generic_attraction` params.
- `read_yaml_config`: relax so `intra_domains` is **optional**. A file with only a
  `disordered:` section = single domain (scale 1.0 everywhere) + mask; a file with only
  domains = today's behavior; both = domains + IDR.
- `topo/core/models.py`, `topo/utils/config.py`, `topo/csp/core.py`: **unchanged** — they
  already pass `domain_def` through. Continuous-synthesis nascent chains get IDR support
  for free once the reader understands the section.

**Precedence — disorder wins.** The mask is applied *after* domain scaling, so a residue
listed in both a domain and `disordered:` has its contacts zeroed regardless of its
`nscale`. The two sections therefore never conflict.

**Unified file format** — residue ranges reuse the syntax already parsed by
`parse_residue_list` (ints, `"i"`, `"start-end"`). Only `n_residues` is required; all
three sections are optional:

```yaml
# domain_def.yaml  — the SAME file already passed as domain_def
n_residues: 283

# --- domain scaling of native side-chain contacts (optional, unchanged) ---
intra_domains:
  A: { residues: [1-50],   nscale: 1.0 }
  B: { residues: [60-283], nscale: 1.0 }
inter_domains:
  A-B: 0.5

# --- disordered / IDR regions (optional) ---
disordered:
  residues: [1-24, 150-165]   # native contacts removed for these residues
  generic_scale: 0.0          # 0 = Level A (self-avoiding); O'Brien parity = 0.03 (x nscale x eps_BT)
  cross_scale:   0.0          # disorder<->folded attraction (Level B); 0 = excl-vol only
```

Level A vs B is expressed entirely in-file: `generic_scale == 0` → Level A;
`generic_scale > 0` → Level B (with `cross_scale` optionally adding the disorder↔folded
term).

### 2.3 Algorithm (inside `build_nonbonded_interaction`)

After the existing `eps_ij` / `binary_contact_matrix` are assembled (with domain scaling
already applied) and **before** the non-native fill loop — so the mask overrides any
domain `nscale`:

```python
# dis_res, generic_scale, cross_scale come from the `disordered:` section of domain_def
if dis_res:                                                # empty/absent -> skip entirely
    dis_idx = [resid_to_index[k] for k in resid_to_index if k[1] in dis_res]
    m = np.zeros(n_residues, bool); m[dis_idx] = True
    involves_dis = m[:, None] | m[None, :]                 # pair touches a disordered res

    # Level A: remove native contacts touching any disordered residue.
    eps_ij[involves_dis] = 0.0
    binary_contact_matrix[involves_dis] = 0                # -> non-native / excluded-vol

    # Level B (generic_scale > 0): shallow attraction, depth = generic_scale * eps_BT(i,j).
    if generic_scale > 0:
        dd = m[:, None] & m[None, :]                       # both disordered
        # + optional disorder<->folded term when cross_scale > 0
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
- Folded-only runs (no `disordered:` section): byte-for-byte identical to today.

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
4. **Scope.** The `disordered:` section is one flat residue set per system (this spec).
   Do you need per-chain masks (e.g. residue 1-24 of chain A only)? If so, the `residues:`
   entries would need chain-qualified keys; otherwise a bare residue list is enough.
