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
   volume via the 12-10-6 sum rule `Rmin/2(i) + Rmin/2(j)`.

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
| Excluded volume (non-native repulsion) | per-atom `Rmin/2 = (K–B nearest non-contact)·2^(1/6)/2` | per-atom `Rmin/2 = rvdw·2^(1/6)` (from `sigmin = 2·rvdw`, L1638→L1697) |

The nonbonded potential is **12-10-6** (O'Brien's custom CHARMM), and CHARMM combines the
per-atom half-radii by the **sum rule**, so the excluded-volume *pair* distance is
`R_ij = Rmin/2(i) + Rmin/2(j)` in **both** columns (the `2·rvdw` at L1638 is a per-atom
`sigmin`, never the pair distance). eps for this term = 0.000132 kcal/mol.

The weak attraction is a separate `NBFIX` override, written for every non-local pair
|i−j| ≥ 3 (lines 1866–1874, `casm = 0` branch):

```
ε_ij  = 0.03 · nscal · ε_BT(res_i, res_j)    # 3% of a native contact of the same pair
R_ij  = rvdw(res_i) + rvdw(res_j)            # attractive-well position (NOTE: no 2^(1/6))
```

Note the two distances differ: the excluded-volume repulsion uses `Rmin/2(i)+Rmin/2(j) =
2^(1/6)·(rvdw_i+rvdw_j)`, while this attractive NBFIX well sits at `rvdw_i+rvdw_j` (no
2^(1/6)) — slightly inside the pure excluded-volume minimum.

> **Strength = 0.03.** A **native** SS contact of the same residue pair has depth
> `nscal·ε_BT` (line 1762), so a generic pair attracts at **3% of a native contact**. Use
> `0.03` when porting to topo.

It is **non-specific in coverage** (it acts on *every* non-local pair, not just native
contacts) with a **uniform 0.03·nscal prefactor**, but its **depth is pair-type dependent**
via `ε_BT(res_i, res_j)` — so the well is chemically heterogeneous, not the same for all
pairs. Physically: a weak, sequence-modulated, non-fold-encoding attraction — a weakly
collapsing self-avoiding chain, *not* a Go fold. In the original this uses the `%rvdw`
table (Perl lines 720–745) for the well position; the topo port instead reuses the
existing per-AA `Rmin_2` radii (equal to `rvdw·2^(1/6)` within ~1%), so **no new
parameter file is needed** — see §2.1 / §2.5.

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

- **Level A — self-avoiding IDR.** For any residue in the disorder mask, **drop all native
  contacts** (H-bond, BS, SS) in which it participates; those pairs fall back to the
  existing non-native repulsive/excluded-volume term (well depth `NON_NATIVE ≈ 0.000132
  kcal/mol` — a negligible, *non-zero* floor; essentially pure steric repulsion). Result: a
  generic self-avoiding chain with transferable backbone — the cleanest "no fold" baseline,
  and *zero* new parameters.

- **Level B — weakly-collapsing IDR (O'Brien `generic-bt`, energy parity).** In addition
  to A, add the shallow, pair-type-dependent attraction between masked residues (and,
  optionally, disorder↔folded cross pairs via `cross_scale`):

  ```
  ε_ij = generic_scale · nscale · ε_BT(i,j)        # kJ/mol;  generic_scale = 0.03 for O'Brien
  R_ij = Rmin_2[type_i] + Rmin_2[type_j]           # nm;  per-AA, from model_parameters
  ```

  `generic_scale` is the depth ratio to a native contact of the same pair, so **0.03**
  reads as "3% of a native contact." Two **design decisions** (2026-07-18):

  1. **Well position uses the Rmin/2 sum rule** `Rmin_2_i + Rmin_2_j`, *not* O'Brien's
     `rvdw_i + rvdw_j`. This makes the attractive well consistent with the excluded-volume
     term (same Rmin/2 convention) and fixes O'Brien's missing-`2^(1/6)` inconsistency
     (see §1.2); the well sits ~12% further out — a deliberate, accepted deviation.
  2. **Reuse the existing `model_parameters` per-AA `Rmin_2`** (the transferable radii) —
     **no `rvdw.csv` to ship.** These equal `rvdw·2^(1/6)` to ~1% (mean 0.85%, max 1.40%),
     a negligible difference for a sub-kT well.

  So Level B is **energy-parity** with `generic-bt` (`0.03·nscale·ε_BT`) but uses topo's
  consistent Rmin/2 well position rather than O'Brien's `rvdw` sum.

#### Which level for which use case — **Level B is the recommended default for a real IDP**

Decided (2026-07-18): for modelling an actual intrinsically disordered protein/region,
**use Level B.** The reasoning:

- **Flexibility is *not* a reason to prefer A.** Both levels share the identical flexible
  backbone, and Level B's attraction is **sub-kT (~0.03 RT/pair) and non-specific**, so it
  cannot lock in a fold or a persistent contact. B samples the same broad, flexible
  ensemble as A — it only *reweights* it toward more compact configurations. An IDP is
  flexible because it has **no stable fold**, not because it has **no interactions**; B
  captures exactly that (no fold, but realistic transient contacts).
- **Real IDPs are more compact than a self-avoiding walk.** SAXS/smFRET place most IDPs at
  scaling exponent ν ≈ 0.5–0.55 (between theta and good solvent), vs ν ≈ 0.588 for a pure
  self-avoiding chain. **Level A (pure excluded volume) systematically over-expands** a
  typical IDP; Level B's weak `0.03·ε_BT` attraction pulls the ensemble into the observed
  window without folding it. (This is what O'Brien's `generic-bt` was calibrated for.)
- **topo's always-on Yukawa makes B the balanced model.** Charged IDPs already get
  Debye–Hückel charge–charge repulsion (expansion). The physical picture is *repulsion
  (electrostatics) balanced by weak attraction (hydrophobic/transient)* = **B + Yukawa**;
  A + Yukawa keeps only the repulsive half and biases further toward over-expansion.

**When Level A is the better call instead:** a disordered **linker** whose role is reach /
entropic tethering between folded domains (compaction not the observable); a **strongly
charged / highly expanded** IDP that genuinely approaches self-avoiding-walk statistics; or
a deliberately minimal, assumption-free reference ensemble.

**Calibration.** `generic_scale = 0.03` (O'Brien) is a sensible starting point but is
tunable — if SAXS/smFRET Rg or ν is available for the target sequence, calibrate
`generic_scale` to match it.

*(Implementation note: Level A is a strict subset of B — B = A's contact removal + the
attractive well — so building A's masking first and then adding B's well is a natural build
order, even though B is the recommended run mode for IDPs.)*

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
  generic_scale: 0.03         # Level B, RECOMMENDED for a real IDP. 0 -> Level A: NO added
                              #   attraction (masked pairs keep the non-native excluded-volume
                              #   floor, ~0.000132 kcal/mol -- NOT zero energy)
  cross_scale:   0.0          # disorder<->folded attraction (Level B); 0 = excl-vol only
```

Level A vs B is expressed entirely in-file: `generic_scale == 0` → Level A;
`generic_scale > 0` → Level B (with `cross_scale` optionally adding the disorder↔folded
term).

### 2.3 Algorithm (inside `build_nonbonded_interaction`)

**Ordering is critical** (three ordered phases). The masking runs before the non-native
fill; the Level B overwrite runs **after** it. Getting this order wrong either loses the
excluded volume or clobbers the attraction — see the two guard rails below.

**Phase 1 — mask** (before the non-native fill loop; overrides any domain `nscale`):

```python
# dis_res, generic_scale, cross_scale come from the `disordered:` section of domain_def
if dis_res:                                                # empty/absent -> skip entirely
    dis_idx = [resid_to_index[k] for k in resid_to_index if k[1] in dis_res]
    m = np.zeros(n_residues, bool); m[dis_idx] = True
    involves_dis = m[:, None] | m[None, :]                 # pair touches a disordered res
    eps_ij[involves_dis] = 0.0                             # transient; refilled in Phase 2
    binary_contact_matrix[involves_dis] = 0               # -> treated as non-native
```

**Phase 2 — existing non-native fill loop (UNCHANGED):** every `binary_contact_matrix[i,j]
== 0` pair (which now includes all masked pairs) gets `eps_ij = NON_NATIVE_KJ` and
`rmin = Rmin/2_i + Rmin/2_j`. **This is what makes `generic_scale = 0` safe:** the masked
pairs come out of Phase 2 with the ~0.000132 kcal/mol excluded-volume floor, *not* zero.

**Phase 3 — Level B overwrite (only if `generic_scale > 0`, and only AFTER Phase 2):**

```python
if dis_res and generic_scale > 0:
    dd = m[:, None] & m[None, :]                           # both residues disordered
    # for the dd pairs (and disorder<->folded pairs when cross_scale > 0):
    #   eps_ij[pair]  = generic_scale * eps_BT[pair]       # nscale already in scaling_matrix
    #   rmin_matrix[pair] = Rmin_2[i] + Rmin_2[j]          # per-AA model_parameters
    # This REPLACES the Phase-2 floor for those specific pairs (one 12-10-6 per pair,
    # not an added term). Optional safety: eps = max(NON_NATIVE_KJ, generic_scale*eps_BT)
    # so Level B is never shallower than the floor.
```

> **Guard rail 1 — `generic_scale = 0` must NOT zero ε.** Never compute
> `ε = generic_scale·ε_BT` for masked pairs *unconditionally*: at `generic_scale = 0` that
> gives `ε = 0`, i.e. **no excluded volume — the chain would pass through itself.** The
> floor must come from Phase 2; Phase 3 only *overwrites* it when `generic_scale > 0`.
>
> **Guard rail 2 — Phase 3 after Phase 2.** If the Level B overwrite runs *before* the
> non-native loop, Phase 2 (which fills every `binary_contact == 0` pair) clobbers it back
> to the floor. Level B must be applied **after** Phase 2 (or set those pairs'
> `binary_contact = 1` so Phase 2 skips them).

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

### 2.5 No new parameter file

**Decided: reuse the existing `model_parameters[self.model]` per-AA `Rmin_2`** for the
Level B well position (`R_ij = Rmin_2_i + Rmin_2_j`). No `rvdw.csv` is shipped. These
transferable radii equal the Perl `%rvdw · 2^(1/6)` to ~1% (mean 0.85%, max 1.40%);
topo's table is O'Brien's rounded ribosome `S<aa>` set (a few residues share a value,
e.g. GLN/HIS/ILE/LEU/MET), but the difference is negligible for a sub-kT well. Level A
needs no radii at all (it reuses the existing K–B excluded-volume term).

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
   same chain fully-masked with `generic_scale=0.03` in topo; compare the per-pair
   ε and Rmin tables (expect match up to unit conventions).

---

## 4. Open questions for the user

1. ~~**Level A vs B default.**~~ **RESOLVED (2026-07-18): Level B is the default for a real
   IDP** (`generic_scale = 0.03`); Level A is for flexible linkers / strongly-expanded IDPs
   / a minimal reference. Rationale in §2.1 ("Which level for which use case"). Build order
   may still start with A's masking (A ⊂ B), but B is the recommended run mode.
2. **Excluded-volume radius source for masked residues (RECOMMENDED, not yet locked).**
   Today Level A leaves masked residues on the **structure-derived K–B** `Rmin/2` (from the
   arbitrary IDR input coordinates), while Level B's attractive well uses the **per-AA
   `model_parameters` `Rmin_2`**. Recommendation: give masked residues the per-AA `Rmin_2`
   for excluded volume in **both** levels — it is conformation-independent (K–B is
   meaningless for a disordered chain) and makes A and B differ **only in well depth ε**.
   Trade-off: masked residues then use a different radius *source* than the folded chain
   (per-AA vs K–B), which is intentional. Confirm at implementation.
3. **Disorder↔folded interactions.** When a masked tail passes near the folded core,
   should it feel (i) only excluded volume (Level A), or (ii) the weak generic attraction
   (Level B `cross_scale > 0`)? O'Brien's `interact_scale` is the analog.
4. **Input geometry.** Masking a region of a *folded* PDB keeps that region's native
   (folded) starting coordinates. If you want a genuinely extended IDR start, we should
   also generate an extended-chain segment for those residues (analog of the Perl
   `create_unstructured_pdb.pl`). Is an extended-start needed, or is masking-in-place
   enough for the intended use?
5. **Scope.** The `disordered:` section is one flat residue set per system (this spec).
   Do you need per-chain masks (e.g. residue 1-24 of chain A only)? If so, the `residues:`
   entries would need chain-qualified keys; otherwise a bare residue list is enough.
