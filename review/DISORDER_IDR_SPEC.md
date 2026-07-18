# Disordered / IDR regions in the Cα model — implementation spec

> **Date:** 2026-07-17
> **Reference (read-only):**
> `sandbox/original_code/perl/cg_model/create_cg_protein_model_v34_0.03.pl` (v1.34, O'Brien).
> **topo code CHANGED (3 files):**
> • `topo/utils/nonbonded.py` — `read_yaml_config` parses `disordered:` (residue set +
>   `idr_scale`, default 0.03; `intra_domains` now optional); **new `apply_disorder`**;
>   `build_nonbonded_interaction` calls it at the end (energy path, §2.3).
> • `topo/analysis/native_contacts.py` — `load_domains` parses `disordered:` + subtracts it
>   from each domain; `build_native_contacts` gains optional `disorder=None` and drops IDR pairs;
>   the `main` driver threads it (Q path, §2.7).
> • `topo/optimize/optimize.py` — `Scorer` threads the disorder set into its
>   `build_native_contacts` calls (optimizer Q, §2.7).
> **NOT changed (pass-through / read-only):** `topo/core/models.py`, `topo/utils/config.py`,
> `topo/csp/core.py` (already thread `domain_def` → get IDR for free, §2.2/§2.6);
> `topo/parameters/model_parameters.py` (per-AA `Rmin_2` read only); other `build_native_contacts`
> callers (e.g. `topo/analysis/mirror.py`) keep the `disorder=None` default → byte-identical.
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
feature reduces to per-residue contact masking inside `build_nonbonded_interaction`**, plus
a per-pair-class energy fill. No changes to bonds, angles, dihedrals, or the OpenMM force
objects are required. This mirrors the existing `ribo_free_mask` pattern.

**One model, three pair classes, one knob** (§2.1). Every pair falls into exactly one class:
**folded–folded** (unchanged: Go / non-native), **IDR–IDR** (`max(NON_NATIVE,
idr_scale·ε_BT)`), **IDR–folded** (**excluded-volume only** — plain non-native). Well
position is a per-residue effective radius (IDR residues → transferable per-AA `Rmin_2`; folded
residues keep K–B). The sole knob is `idr_scale`: `0` is the pure self-avoiding chain (no
separate "Level A"), default `0.03` reproduces O'Brien's default configuration. (IDR↔folded
attraction is a documented future extension, not an active knob — §2.1.)

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

### 1.4 How the original treats IDR↔folded (cross) interactions — **default is steric-only**

The intra-IDR generic attraction (§1.2, `0.03·ε_BT`) fires **only within a single `generic`
PDB** (the loop runs over one `$np`, pairs `|i−j|≥3`, lines 1862–1874). Interactions
**between two different PDB files** — which is how O'Brien represents an IDR next to a folded
domain — are a **separate, opt-in mechanism**, written to a separate file and governed by the
`interact_scale` keyword (lines 1885–1908):

```perl
if($interact_scale ne "undefined") {         # OFF by default ($interact_scale = "undefined", L164)
  for each cross-PDB pair (i in np, j in np2):
    $ene  = $interact_scale * $eps[i][j];    # interact_scale · ε_BT
    $temp = $rvdw{i} + $rvdw{j};             # well at rvdw sum
}
```

Three properties matter for the port:

1. **Default OFF → excluded volume only.** `$interact_scale` defaults to `"undefined"`, and the
   header is explicit (lines 111–113): *without the keyword these cross parameters "will not be
   written out even if there are two PDB files listed."* So **by default an IDR and a folded
   domain feel only excluded volume** across the boundary — no attraction.
2. **Separate scale from the intra 0.03.** The within-IDR attraction is hardcoded `0.03·ε_BT`;
   the cross attraction is `interact_scale·ε_BT` with a user-chosen scale. O'Brien deliberately
   did *not* reuse 0.03 for cross pairs — they are independent knobs over disjoint pair sets.
3. **Keyed on "different PDB," not on disorder.** `interact_scale` scales *every* inter-PDB pair
   (folded↔folded, folded↔IDR, IDR↔IDR alike); it errors with a single PDB (line 349). O'Brien
   has no residue mask, so he cannot single out disorder↔folded — the granularity is the whole
   PDB. topo, reading one structure with a mask, **could** target disorder↔folded specifically
   if this cross attraction is ever wanted — a strict improvement (see the extension note in §2.1).

> **Divergence (multiple IDRs).** topo's flat `disordered:` set treats *all* disordered–disordered
> pairs with `idr_scale`, including residues in two *different* IDRs. In O'Brien those would be
> inter-PDB → `interact_scale`, not the intra `0.03`. For a single contiguous IDR the two are
> identical; the distinction only appears with multiple separate IDRs. Accepted for now (one flat
> set, §4 Q5); revisit if per-IDR cross scaling is ever needed.

---

## 2. Proposed implementation for topo

### 2.1 One model, three pair classes (one knob)

There is a **single** IDR model, not two levels. Every residue pair falls into exactly one of
three classes by how many of its residues are in the disorder mask, and the model just picks
the energy for each class. The old "Level A / Level B" split collapses into **one continuous
knob** (`idr_scale`): the pure self-avoiding chain is simply `idr_scale = 0`.

| Pair class | Native contacts | Well depth ε_ij | Well position R_ij | O'Brien analog |
|--|--|--|--|--|
| **folded–folded** (neither in mask) | kept (Go) | today's HB+BS+scaled_SS, else `NON_NATIVE` | native Cα dist, else K–B sum | Go segment (unchanged) |
| **IDR–IDR** (both in mask) | removed | `max(NON_NATIVE, idr_scale·ε_BT(i,j))` | per-AA + per-AA | intra-`generic` `0.03·ε_BT` |
| **IDR–folded** (exactly one in mask) | removed | `NON_NATIVE` (**excluded-volume only**) | per-AA + K–B | cross `interact_scale`, **off by default** |

- **`idr_scale`** (default **0.03**, the *only* knob) controls IDR–IDR attraction. `0` →
  pure excluded volume (self-avoiding chain); `0.03` → O'Brien-like weak collapse.
- **IDR–folded is excluded-volume only** — plain non-native. Those pairs need *no* special
  energy handling: the `apply_disorder` transform (§2.3) sets them to `NON_NATIVE` at the
  per-residue radius sum (per-AA on the IDR side, K–B on the folded side). This reproduces
  O'Brien's default (`interact_scale` undefined → steric only, §1.4).

> **Extension point (not implemented) — transient IDR↔folded (fuzzy) binding.** If a study ever
> needs weak IDR↔core attraction, promote the IDR–folded pairs to `max(NON_NATIVE,
> cross_scale·ε_BT)` with a `cross_scale` knob (default 0) — the direct analog of O'Brien's
> `interact_scale`, kept separate from `idr_scale` just as O'Brien kept `interact_scale`
> separate from the intra `0.03` (§1.4). It is deliberately **left out now**: the current use
> case defines IDR–folded as steric-only, and an always-zero knob is unrequested surface. Adding
> it later is a small change to `apply_disorder` (§2.3): add a `cross = involves_dis & ~dd` mask
> and fill `eps_ij[cross] = max(NON_NATIVE, cross_scale·ε_BT)`, not a redesign.
>
> **Only `eps_ij` changes — the well position is already correct.** `apply_disorder` sets
> `rmin_matrix` for *all* `involves_dis` pairs (cross included) to the sum-rule position
> `rmin_2_nm[IDR](per-AA) + rmin_2_nm[folded](K–B)`. Since the model uses one 12-10-6 term per
> pair, whose minimum is at `rmin` independent of depth, turning on the attraction only *deepens*
> that existing well — `rmin_matrix` and `rmin_2_nm` are untouched. (This yields a topo-consistent
> cross attraction, not a bit-exact O'Brien `interact_scale`, whose well sat at the `rvdw` sum
> with a per-AA folded side — matching that would also require changing the folded radius, which
> contradicts the "folded keeps K–B" decision.)

The IDR–IDR attraction uses the **same 12-10-6 well** the model already has — one term per
pair, no new force. The `max(NON_NATIVE, …)` floor guarantees excluded volume even at
`idr_scale = 0` (or for the handful of pairs whose `ε_BT ≈ 0` near the BT reference), so the
chain can never pass through itself.

#### ε: use the raw BT pair energy, `nscale = 1`

```
ε_BT(i,j) = ss_interaction_energy[i,j] = 4.184 · |raw − 0.6|      # kJ/mol, from bt_potential.csv
ε_ij      = max(NON_NATIVE, idr_scale · ε_BT(i,j))            # IDR-IDR pairs only
```

**Design decision — `nscale = 1` (no domain ladder for IDR pairs, 2026-07-18).** In topo,
`nscale` is the per-structural-class LADDER factor (`topo/optimize/optimize.py`, ~1.15–2.5) —
*the smallest multiplier that makes a folded domain marginally stable*, a folding-thermodynamics
calibration. An IDR has no fold and no stability target, so there is no principled ladder value
to inherit; a masked residue is also not in any `intra_domain`, so its `scaling_matrix` entry is
ill-defined. We therefore multiply the **raw** per-pair BT matrix `ss_interaction_energy` (from
`get_ss_interaction_energy`) — **not** `scaled_ss_interaction_energy`, **not** `scaling_matrix·…`,
and **not** a hand-rolled `(raw − 0.6)` re-read of `bt_potential.csv` (that would drop the
`abs()`, flipping the sign for 394/400 pairs, and drop the kcal→kJ ×4.184 factor). Reusing
`ss_interaction_energy` keeps abs, shift, and units identical to every other energy term, from a
single source of truth. `idr_scale` thus becomes the physically-interpretable coupling —
calibrate to SAXS/smFRET Rg/ν rather than reading 0.03 as "3% of a native contact"
(a counterfactual for a region with no native contacts).

This is **structurally** O'Brien's `generic-bt` (the intra-IDR attraction) but **deliberately
decoupled** from topo's folding ladder (`nscale = 1`). It is *not* bit-exact parity with O'Brien's
`0.03·nscal·ε_BT`: his `nscal` was a single global scale, whereas topo's is per-domain and
inapplicable to a fold-free region. If bit-exact reproduction were ever required, his global
`nscal` would be folded into `idr_scale` explicitly.

#### R: override the per-residue radius for IDR residues, keep folded K–B

The excluded-volume radius is a property of the **residue**, so `apply_disorder` overrides the
per-residue radius array **in place** (§2.3) and lets the existing sum rule combine pairs:

```
rmin_2_nm[dis] = Rmin_2_perAA[dis]     # IDR residues -> transferable per-AA (model_parameters)
                                       #   (folded residues keep their calculate_rmin_2_values K-B value)
R_ij           = rmin_2_nm_i + rmin_2_nm_j    # sum rule, as today
```

This yields **IDR–IDR** = per-AA + per-AA, **IDR–folded** = per-AA + K–B, **folded–folded** =
K–B + K–B (byte-identical to today). Two design decisions ride along:

1. **Override only IDR residues** (not both sides of an IDR-involving pair). A folded residue
   keeps its own collision radius regardless of partner — physically correct, and it makes the
   folded core's excluded volume *exactly* unchanged. This is why the override is on the
   **per-residue `rmin_2_nm` array**, not a per-pair choice (and why the same array stays
   consistent across the NC↔ribosome channel in CSP, §2.6).
2. **IDR residues use the transferable per-AA `Rmin_2`** (from `model_parameters`, = O'Brien's
   `S<aa>` set ≈ `rvdw·2^(1/6)` to ~1%) rather than their **structure-derived K–B** value — the
   K–B radius is meaningless for arbitrary/disordered input coordinates. **No `rvdw.csv` to
   ship.** This also fixes O'Brien's missing-`2^(1/6)` inconsistency between his attractive well
   (`rvdw` sum) and excluded volume (`rvdw·2^(1/6)` sum): topo uses one Rmin/2 convention for
   both, so the well sits ~12% further out — a deliberate, accepted deviation.

#### Recommended defaults for a real IDP

- **`idr_scale = 0.03`** (weakly-collapsing IDR). Rationale:
  - **`idr_scale = 0` (self-avoiding) systematically over-expands a typical IDP.** SAXS/smFRET
    place most IDPs at scaling exponent ν ≈ 0.5–0.55 (between theta and good solvent), vs ν ≈ 0.588
    for a pure self-avoiding chain. The sub-kT (~0.03 RT/pair), non-specific `0.03·ε_BT` attraction
    pulls the ensemble into the observed window **without** locking in a fold — it *reweights* the
    same broad, flexible ensemble toward compaction. An IDP is flexible because it has **no stable
    fold**, not because it has **no interactions**.
  - **topo's always-on Yukawa makes this the balanced model.** Charged IDPs already get
    Debye–Hückel repulsion (expansion); the physical picture is *repulsion balanced by weak
    attraction* = `idr_scale > 0` + Yukawa. Setting `idr_scale = 0` keeps only the
    repulsive half and biases further toward over-expansion.
  - `idr_scale` is **tunable** — if SAXS/smFRET Rg or ν is available, calibrate to match.
- **`idr_scale = 0` (self-avoiding) is the better call** for: a disordered **linker** whose
  role is reach / entropic tethering (compaction not the observable); a **strongly charged /
  highly expanded** IDP that genuinely approaches self-avoiding-walk statistics; or a deliberately
  minimal, assumption-free reference ensemble.
- **IDR↔folded is steric-only** (the faithful O'Brien default, §1.4) — no knob. See the §2.1
  extension note if transient/fuzzy IDR–domain binding is ever needed.

> **Why masking (not a `nscale = 0` domain) is the right lever.** The existing `domain.yaml`
> scaling multiplies **only the SS energy** (`scaling_matrix * ss_*`). Setting a domain's
> `nscale = 0` would zero SS contacts but **leave H-bond and BS contacts intact**, so the region
> would still fold via backbone H-bonds. A true IDR must remove H-bond + BS + SS together — hence
> a dedicated mask applied to all three, not a domain scale.

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
argument, and no new INI key**. The only code changes are the YAML reader and the
`apply_disorder` transform:

- `topo/utils/nonbonded.py`: `read_yaml_config` also parses the optional `disordered:`
  section (returning the residue set + `idr_scale`). **`idr_scale` is optional and defaults to
  `0.03` when the key is omitted** (only `residues:` is required); an explicit `0` selects the
  self-avoiding chain. `build_nonbonded_interaction` runs its **folded build unchanged**, then
  calls a new **`apply_disorder`** transform on the three outputs (§2.3). Signature is unchanged
  except the already-planned `return_rmin_2` — **no** `disorder_def` / `generic_attraction` params.
- `read_yaml_config`: relax so `intra_domains` is **optional**. A file with only a
  `disordered:` section = single domain (scale 1.0 everywhere) + mask; a file with only
  domains = today's behavior; both = domains + IDR.
- `topo/core/models.py`, `topo/utils/config.py`, `topo/csp/core.py`: **unchanged** — they
  already pass `domain_def` through. Continuous-synthesis nascent chains get IDR support
  for free once the reader understands the section.

**Precedence — disorder wins (overlap is allowed and well-defined).** Domain and `disordered:`
ranges *may* overlap. `apply_disorder` runs **after** the whole folded build — including domain
scaling (`scaling_matrix · ss_*`) — and **unconditionally overwrites** `eps_ij`, `rmin_matrix`,
and `rmin_2_nm` for every pair touching a disordered residue (§2.3). So if a residue is listed in
both a domain and `disordered:`, its domain `nscale` is computed and then discarded; the pair is
governed entirely by the disorder rules. This holds for **inter-domain** pairs too: **if either
residue in a pair is disordered, the disorder rules govern it — domain membership of either side
is irrelevant.** The two sections therefore never conflict. This makes overlap a *feature*: define
a domain broadly (e.g. `A: 1-100`) and carve a disordered loop out of it (`disordered: 40-50`)
without splitting the domain definition.

> **Implementation note — log detected overlap.** Because disorder wins *silently*, a residue
> accidentally left in both sections is disordered with no warning — a typo could quietly
> disorder part of a domain. `read_yaml_config` should emit an **info** line (not an error —
> overlap is legal) listing the overlapping residues, e.g.
> `"IDR overlap: residues 40-50 are in both domain A and disordered:; treating as disordered."`
> Cheap insurance; keep it a one-line message, not a per-residue dump.

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
  residues: [1-24, 150-165]   # native contacts removed for these residues (only required key)
  idr_scale: 0.03             # OPTIONAL, defaults to 0.03 if omitted. IDR-IDR attraction (the only
                              #   knob). 0 = self-avoiding chain: masked pairs keep only the
                              #   excluded-volume floor (~0.000132 kcal/mol -- NOT zero energy).
                              # IDR-folded pairs are always excluded-volume only (no knob, §2.1).
```

The IDR model is a single continuous knob: `idr_scale = 0` → self-avoiding chain;
`idr_scale > 0` → weakly-collapsing IDR. IDR↔folded pairs are excluded-volume only. See the
three-class table in §2.1.

### 2.3 Algorithm — build folded, then `apply_disorder` transform

**`build_nonbonded_interaction` builds the normal fully-folded parameters exactly as today** —
the pair matrices `rmin_matrix` (nm), `eps_ij` (kJ/mol), and the per-residue `rmin_2_nm` (nm,
from `calculate_rmin_2_values`). Disorder is then a **single transform on those three outputs**,
called at the very end (after the nm conversion, so everything is in nm). This keeps the folded
build free of IDR logic and makes the transform independently testable. The model consumes only
the *combined* pair matrices (plus the per-residue `rmin_2_nm` for CSP, §2.6), never the
individual HB/BS/SS terms — so a pure output transform is sufficient.

```python
# Called at the end of build_nonbonded_interaction, only if the `disordered:` section is present.
# Auxiliary inputs are already in local scope from the folded build:
#   index_to_resname       -- from get_residue_mapping (per-AA Rmin_2 + resname lookup)
#   ss_interaction_energy  -- the RAW per-pair BT matrix from get_ss_interaction_energy (line 934),
#                             already = 4.184*|raw - 0.6| kJ/mol (abs + shift + units applied once)
def apply_disorder(rmin_matrix, eps_ij, rmin_2_nm,      # the three folded outputs (nm, kJ/mol, nm)
                   dis_idx, idr_scale,
                   index_to_resname, ss_interaction_energy):
    m = np.zeros(len(rmin_2_nm), bool); m[dis_idx] = True
    involves_dis = m[:, None] | m[None, :]              # pair touches a disordered residue
    dd           = m[:, None] & m[None, :]              # both residues disordered

    # 1. IDR residues' radius -> transferable per-AA (K-B is meaningless for disordered coords).
    #    Folded residues are left on their exact folded-run K-B value. This per-residue array
    #    also feeds the NC<->ribosome excluded volume in CSP (§2.6), so overriding it here keeps
    #    the intra-chain and NC<->ribosome excluded volume consistent for IDR beads.
    rmin_2_nm[dis_idx] = [Rmin_2_perAA[index_to_resname[i]] for i in dis_idx]   # nm

    # 2. every IDR-involving pair: drop the native contact -> reposition the well at the sum rule
    #    of the (now-updated) per-residue radii. IDR-folded comes out per-AA + K-B automatically.
    rmin_matrix[involves_dis] = (rmin_2_nm[:, None] + rmin_2_nm[None, :])[involves_dis]
    eps_ij[involves_dis]      = NON_NATIVE_KJ           # excluded-volume floor (covers IDR-folded)

    # 3. IDR-IDR only: the weak generic attraction, never below the floor.
    eps_ij[dd] = np.maximum(NON_NATIVE_KJ, idr_scale * ss_interaction_energy[dd])
    return rmin_matrix, eps_ij, rmin_2_nm
```

Three points make this correct and clean:

- **No `Rmin_2_eff` intermediate.** The override is applied **in place to `rmin_2_nm`** — the
  per-residue array the model and CSP actually consume — and `rmin_matrix` for IDR pairs is
  rebuilt from it. IDR–IDR = per-AA+per-AA, IDR–folded = per-AA+K–B, folded–folded = K–B+K–B.
- **Folded core exactly unchanged.** `calculate_rmin_2_values` runs on the **unmasked** structure
  (as today), so every folded residue keeps its exact folded-run K–B radius — the mask cannot
  perturb it, and arbitrary IDR input coordinates cannot leak into a folded residue's radius.
  Only `involves_dis` entries of the matrices are overwritten; folded–folded entries are byte-
  identical to a folded-only run.
- **The `max(NON_NATIVE_KJ, …)` floor** makes `idr_scale = 0` (and any pair whose `ε_BT ≈ 0`
  near the BT reference) safe: no pair drops below the excluded-volume floor, so the chain can
  never pass through itself. One 12-10-6 well per pair — no added term.

> **Why the raw `ss_interaction_energy`, not `scaling_matrix`.** `nscale` lives in
> `scaling_matrix` and is the per-domain LADDER stability factor (§2.1); the IDR well must not
> inherit it, so the transform reads the *unscaled* BT matrix. Reusing `ss_interaction_energy`
> (rather than re-reading `bt_potential.csv` and recomputing `(raw − 0.6)`) also preserves the
> `abs()` — without it `(raw − 0.6)` goes negative for 394/400 pairs and flips the
> hydrophobicity ordering — and the kcal→kJ ×4.184 factor.

### 2.4 What does *not* change (important)

- Bonds, angles, dihedrals, Yukawa: untouched — topo's global transferable backbone is
  already the disordered-appropriate choice (§1).
- OpenMM force construction (`addCustomNonBondedForce`): unchanged; it consumes the same
  two matrices.
- The folded build inside `build_nonbonded_interaction`: unchanged — `apply_disorder` runs
  *after* it and only rewrites IDR-involving entries.
- Folded-only runs (no `disordered:` section): byte-for-byte identical to today (`apply_disorder`
  is not called). Even *with* a `disordered:` section, folded–folded pairs and every folded
  residue's K–B radius are untouched.

### 2.5 No new parameter file

**Decided: reuse the existing `model_parameters[self.model]` per-AA `Rmin_2`** as the IDR
residues' override radius (`rmin_2_nm[dis_idx]`; the sum rule `R_ij = rmin_2_nm_i + rmin_2_nm_j`
then combines pairs). No `rvdw.csv` is shipped. These transferable radii equal the Perl
`%rvdw · 2^(1/6)` to ~1% (mean 0.85%, max 1.40%); topo's table is O'Brien's rounded ribosome
`S<aa>` set (a few residues share a value, e.g. GLN/HIS/ILE/LEU/MET), but the difference is
negligible for a sub-kT well. Folded residues need no new radii (they keep the existing K–B
excluded-volume term).

### 2.6 CSP / ribosome interaction (nascent-chain IDR)

Continuous synthesis gets IDR support **for free** and stays consistent, because the whole
feature rides the per-residue `rmin_2_nm` array that CSP already threads. The data path:

```
build_nonbonded_interaction(full_pdb, domain_def, return_rmin_2=True)   # apply_disorder runs here
   -> rmin_2_full  (per-residue Rmin/2, nm; IDR residues already overridden to per-AA)
precompute_contacts  returns rmin_2_full                                 (csp/core.py:391)
   -> per length L:  nascent_rmin_2 = rmin_2_full[:L]                    (csp/core.py:442)
   -> build_length_model: setParticlesRadii(nascent_rmin_2)             (csp/core.py:522)
   -> the {nascent}x{ribosome} CustomNonbondedForce combines each nascent particle's Rmin/2
      with each ribosome bead's Rmin/2 by the sum rule.
```

Because `apply_disorder` overrides `rmin_2_nm` (§2.3, operation 1) — **not** just the pair matrix —
the same per-residue radius reaches **both** excluded-volume channels, so they cannot disagree:

| Bead | Rmin/2 for excluded volume | Source |
|--|--|--|
| Ribosome bead | **per-AA transferable** | `Ribosome.Rmin_2_nm` (fixed scenery, unchanged) |
| Nascent **folded** residue | **K–B structure-derived** | `calculate_rmin_2_values` (unchanged) |
| Nascent **IDR** residue | **per-AA transferable** | `apply_disorder` override |

So an IDR nascent bead meets the ribosome with **per-AA on both sides** — the correct convention
(its K–B value is meaningless for disordered coordinates, and it matches the ribosome scenery).

> **Requirement.** `apply_disorder` **must** override the per-residue `rmin_2_nm`, not only the
> `rmin_matrix` pair block. Overriding only the pair matrix would give an IDR bead per-AA against
> other nascent beads but K–B against the ribosome — an inconsistency. Overriding the array (and
> rebuilding the pair block from it, §2.3) keeps both channels from one source.

**Safety / independence:**

- **No `disordered:` section → CSP byte-identical.** `precompute_contacts` passes `domain_def`
  straight through (`csp/core.py:391`); `apply_disorder` is not called when the section is absent.
- **The L24 free-loop path is independent.** `append_flexible_l24_loop` calls
  `build_nonbonded_interaction(atomistic_pdb, return_rmin_2=True)` with **no** `domain_def`
  (`ribosome.py:703`), so `apply_disorder` never fires there; that path keeps its own K–B radius
  override for the freed *ribosomal-protein* loop (`ribosome.py:756-759`). No cross-contamination.
  (An IDR *inside* a freed ribosomal loop would be a separate, future concern.)

### 2.7 Native-contact analysis (Q) & the nscale optimizer

topo has **two independent native-contact definitions**, and the IDR mask must reach **both**:

| Path | Function | Consistent today? |
|--|--|--|
| **Energy** (the model) | `build_nonbonded_interaction` → `apply_disorder` (§2.3) | ✅ IDR contacts removed |
| **Analysis (Q)** | `build_native_contacts` / `load_domains`, `topo/analysis/native_contacts.py` | ❌ built from the reference, mask-blind |

The analysis path derives contacts purely from the all-atom reference (heavy-atom cutoff) and
`load_domains` reads only `intra_domains` — it never sees `disordered:`. Left unfixed, Q counts
IDR-involving native contacts that the energy function removed and that therefore **never form**,
sitting permanently in the denominator and **deflating** `Q_protein` and any overlapping
`Q_domain`. The nscale optimizer uses this Q as its `folded_fraction` stability metric, so with a
`disordered:` section present it would never reach target and would drive `nscale` to the ceiling.

**Rule (decided 2026-07-18): a pair touching *any* IDR residue is not a native contact.** The
analysis and optimizer must apply the **same mask as the energy path** — exclude every native
contact involving a disordered residue from `build_native_contacts` and from the optimizer's Q.
`load_domains` / the Q driver read the `disordered:` section from the *same* `domain_def` (it is
already the required `-d/--domain` input) and drop those pairs — the identical one-liner as the
energy path ("drop pairs where either residue is in `dis_res`").

**Two call sites reuse these functions**, so both get the fix once `native_contacts.py` is updated:
the standalone Q driver (`native_contacts.main`) and the optimizer's `Scorer`
(`topo/optimize/optimize.py`), which builds Q via `load_domains` + `build_native_contacts`.
Concretely: `build_native_contacts` gains an optional `disorder=None` set (default = today's
behavior, so other callers like `mirror.py` are byte-identical); `load_domains` returns the
disorder set; and both `main` and `Scorer` pass it through. The optimizer's **energy** side needs
no change — its `topo-mdrun` subprocess already builds through `build_nonbonded_interaction` and so
gets `apply_disorder` for free.

**Effective domain membership = domain residues − disordered residues.** Overlap is allowed for
convenience (§2.2), but **in reality an overlapped residue is disordered only**. So `load_domains`
also **subtracts** disordered residues from each domain's residue set, so a residue listed in both
domain A and `disordered:` no longer contributes to `Q_A` (nor to any interface Q). This mirrors
the energy path exactly, where "disorder wins" already makes such a residue disorder-only (§2.2).

> **Worked example.** `intra_domains: {A: 1-100}` with `disordered: {residues: 40-50}` gives
> **effective A = {1-39, 51-100}** — a discontinuous set but still **one domain with a hole**, not
> two pieces. The optimizer calibrates `nscale_A` over all native contacts among {1-39, 51-100},
> **including the cross-loop 1-39 ↔ 51-100 contacts** (both residues folded → retained), which are
> exactly what holds the two halves together across the excised loop. Only pairs touching 40-50 are
> dropped from Q / rescaled in energy (40-50↔40-50 → `idr_scale·ε_BT`; 40-50↔folded → excluded
> volume). The loop stays tethered — the 39-40 and 50-51 backbone bonds remain — so A folds as one
> unit joined by a flexible (non-Go) loop.

**The optimizer must SEE the IDR — optimize in the production (IDR-present) condition.** Masking a
region removes its native contacts, *including* any it made with the domain (interface / contiguous
contacts) and the folded scaffold it provided, so **the domain is genuinely less stable when the
IDR is present.** `nscale` must therefore be calibrated with the IDR active, not on the full fold.
There are two *opposite* ways to get this wrong:

| Optimizer setup | Energy | Q metric | `nscale` error | Domain at production |
|--|--|--|--|--|
| Folded-first (optimize with **no** `disordered:` section) | full Go — IDR region folds | full contacts | over-counts the IDR scaffold → **too low** | **under-stable** (may unfold) |
| Q-blind (IDR masked in energy, but Q counts IDR contacts) | `apply_disorder` | counts IDR pairs | deflated Q → **too high** | **over-stabilized** |
| **Correct** | `apply_disorder` | **excludes** IDR pairs (this §) | measures foldable core under production energy | **marginally stable** |

So the required setup is: run the optimizer with the **`disordered:` section active** (energy via
`apply_disorder`) **and** the **masked Q** above. The optimizer then observes the true, IDR-reduced
stability of the foldable core and raises `nscale` to compensate for the contacts the IDR removed —
neither under- nor over-shooting. The ladder `nscale` tunes only the **folded** domains; IDR
residues use `idr_scale` (`nscale = 1`, decoupled) and are outside its scope. *(Folded-first is a
valid shortcut only in the special case where the disordered region shares no native contacts with
the optimized domain — but including the IDR is always the robust choice, so don't rely on it.)*

**Net effect (two guarantees):**
1. **Disorder does not block convergence.** With the masked Q, never-forming IDR contacts are out
   of the denominator, so `Q` can reach target when the foldable core folds — no artificial cap.
   *(If the core genuinely cannot fold without the IDR's contacts even at the ladder ceiling
   `nscale = 2.5044`, the optimizer lands on the ceiling/fallback — a true physical result,
   correctly reported, not a masking artifact.)*
2. **Disorder is never scaled by `nscale`.** `apply_disorder` overwrites every pair touching a
   disordered residue *after* domain scaling, so its `scaling_matrix` (`nscale`) factor is
   discarded; the residue depends only on `idr_scale`. Varying `nscale` during the search has zero
   effect on disordered residues — the search acts solely on the folded core.

---

## 3. Validation plan

1. **Regression (no-op):** a folded PDB with `disorder_def=None` and with an *empty*
   disorder list must produce identical `rmin_matrix` / `energy_matrix`. → `np.allclose`.
   Also assert **folded–folded pairs are untouched** even *with* a `disordered:` section
   present (only IDR-involving pairs change).
2. **Contact removal + self-avoiding (`idr_scale=0`):** mask an α-helix's residues; assert
   every native contact (H-bond/BS/SS) touching a masked residue is gone from `energy_matrix`
   and those pairs now carry `NON_NATIVE_KJ`. Assert non-masked–non-masked contacts are
   unchanged, and that masked residues' overridden `rmin_2_nm` equals the per-AA
   `model_parameters` value (not their K–B value), while folded residues' `rmin_2_nm` is
   byte-identical to the folded-only run.
3. **IDR–folded is excluded-only:** with `idr_scale>0`, assert every IDR–folded pair (exactly
   one residue masked) carries exactly `NON_NATIVE_KJ` (no attraction leaks across the boundary)
   at well position `rmin_2_nm[IDR] (per-AA) + rmin_2_nm[folded] (K–B)`.
4. **Behavioural:** short MD (a small folded domain + a masked N-terminal tail). Confirm
   energies finite, the tail's radius of gyration expands relative to the fully-Go run,
   and the folded core RMSD stays low. Use the Tutorial-13 short-debug pattern.
5. **ε construction (unit test):** for a chain with a masked region and `idr_scale=0.03`,
   assert `energy_matrix[i,j] == max(NON_NATIVE_KJ, idr_scale·ss_interaction_energy[i,j])`
   for IDR–IDR pairs, and the raw BT identity `ss_interaction_energy == 4.184·|raw−0.6|` from
   `bt_potential.csv` — with **no** `scaling_matrix`/`nscale` factor (the `nscale=1` decoupling,
   §2.1). A Perl `GENERIC-bt` comparison is *not* bit-exact: O'Brien's `0.03·nscal·ε_BT` carries
   his global `nscal`, deliberately dropped here — the per-pair *shape* (which pairs are deepest)
   should match, but absolute depths differ by his `nscal`.
6. **Floor guard (`idr_scale=0`):** assert every masked pair still carries exactly
   `NON_NATIVE_KJ` (excluded volume preserved — the chain cannot pass through itself).
7. **CSP consistency (§2.6):** with a nascent-chain IDR, assert the per-residue `rmin_2_full`
   returned by `precompute_contacts` has per-AA values at IDR indices and K–B elsewhere, and
   that the value fed to `setParticlesRadii` (NC↔ribosome channel) matches the radius used in
   the NC↔NC pair matrix for the same residue. Regression: with **no** `disordered:` section,
   a CSP contact build is byte-identical to today.
8. **Overlap precedence (§2.2):** define a domain with `nscale ≠ 1` and a `disordered:` range
   that overlaps it. Assert every pair touching an overlapped residue is governed by the disorder
   rules (`NON_NATIVE` or `idr_scale·ε_BT`), i.e. **identical** to a run where the same
   residues are disordered but *not* in any domain — the domain `nscale` has no effect on them.
   Non-overlapped domain residues keep their scaled contacts.
9. **Default `idr_scale` (§2.2):** a `disordered:` section with `residues:` but **no** `idr_scale`
   key must produce `energy_matrix` / `rmin_matrix` identical to the same file with an explicit
   `idr_scale: 0.03`. → `np.allclose`.
10. **Q analysis excludes IDR contacts (§2.7):** with a `disordered:` section in `domain_def`,
    assert `build_native_contacts` yields **no** pair touching a disordered residue (in
    `Q_protein`, any `Q_domain`, and interfaces), and that a residue listed in both a domain and
    `disordered:` is **absent** from that domain's Q set (effective membership = domain − disorder).
    Regression: with no `disordered:` section, the native-contact lists are byte-identical to today.

---

## 4. Open questions for the user

1. ~~**Level A vs B default.**~~ **RESOLVED (2026-07-18): one model, one knob.** `idr_scale
   = 0.03` (weakly-collapsing) is the default for a real IDP; `idr_scale = 0` (self-avoiding)
   is for flexible linkers / strongly-expanded IDPs / a minimal reference. There are no longer
   two "levels" — see the three-class table in §2.1.
2. ~~**Excluded-volume radius source for masked residues.**~~ **RESOLVED (2026-07-18): per-AA
   for IDR residues, K–B for folded.** `apply_disorder` overrides the per-residue `rmin_2_nm`
   array in place (§2.3): IDR residues take the conformation-independent per-AA `model_parameters`
   `Rmin_2` (K–B is meaningless for disordered coords); folded residues keep their K–B value
   regardless of partner. Intentional mixed radius *source*; the folded core is exactly unchanged,
   and the same array keeps NC↔NC and NC↔ribosome excluded volume consistent (§2.6).
3. ~~**Disorder↔folded interactions.**~~ **RESOLVED (2026-07-18): excluded-volume only, no knob.**
   IDR–folded pairs are plain non-native (steric), reproducing O'Brien's `interact_scale`-undefined
   default (§1.4). No `cross_scale` is implemented — the current use case defines IDR–folded as
   steric-only, and an always-zero knob is unrequested surface. Transient/fuzzy IDR↔core attraction
   is a documented one-line **extension point** (§2.1) — the `interact_scale` analog — to add only
   if a study needs it.
4. ~~**Input geometry.**~~ **RESOLVED (2026-07-18): mask a single folded structure in place — no
   extended-chain segment needed.** The equilibrium IDR ensemble is set by the *potential* (no
   native contacts + flexible transferable backbone + `idr_scale` attraction), **not** by the
   starting coordinates. With its native contacts removed, the region is no longer held in the
   folded conformation and relaxes toward the disordered ensemble, forgetting its folded start; the
   initial relaxation is discarded in the analysis step, as is routine for any MD run. So **no
   `create_unstructured_pdb.pl` analog is required.** (Moot for CSP anyway: the nascent chain's
   simulated coordinates come from the seeding scheme (`cold_start_positions`), not from folded PDB
   coordinates — see `write_subset_structure`.)
5. ~~**Scope.**~~ **RESOLVED (2026-07-18): single-chain only — a bare residue list, no chain
   qualifier.** The `disordered:` section is one flat residue set per system (`residues:` reuses
   `parse_residue_list`). Chain-qualified keys (e.g. "1-24 of chain A only") are a future extension
   if multi-chain masking is ever needed; not implemented now.
