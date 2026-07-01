# NC ↔ Ribosome interaction — how O'Brien's model works, and how topo reproduces it

Reference notes for the ribosome-nascent-chain (RNC) non-bonded interactions, verified
against O'Brien's actual code and parameter files. Edit freely; annotate what you want changed.

**Sources (O'Brien, read directly):**
- Simulation driver: `…/Continuous_synthesis_protocol/continuous_synthesis_v6.py`
- Force-field/xml generator: `…/CG_protein_parameterization/parse_cg_prm.py`
- CG ribosome builder: `…/CG_ribosome_parameterization/create_cg_ribosome_model.py`
- Example setup: `tutorials/11_reproduce_csp/setup/`
  (`rnc.prm`, `rnc.xml`, `4c5c_model_clean_nscal1_fnn1_go_bt.prm`,
  `50S_tRNA_cg_truncated.psf`, `combine_ribo_L24_Yang.prm`)

---

## 1. Setup: what actually matters

The RNC is one OpenMM system = nascent chain + truncated 50S/tRNA ribosome. The ribosome
beads are **rigid (mass = 0)**, so all intra-ribosome interactions (RNA–RNA, RNA–protein,
protein–protein) are constant and irrelevant. The only forces that matter for the dynamics:

1. **nascent ↔ nascent** (the protein folding on itself), and
2. **nascent ↔ ribosome** (RNA and ribosomal protein).

The whole force field is built by `forcefield.createSystem(...)` from `rnc.xml`
(`continuous_synthesis_v6.py:1548`, `create_elongation_system:639`). `rnc.xml` is generated
from `rnc.prm` by `parse_cg_prm.py`.

---

## 2. One force for all non-bonded pairs

All non-bonded pairs — nascent↔nascent, nascent↔rRNA, nascent↔ribosomal-protein — go through
a **single** `CustomNonbondedForce` (from `rnc.xml`):

```
U = ke·q1·q2/ep/r·exp(−r/ld)   +   kv·( a/r¹² + b/r¹⁰ + c/r⁶ )
    └── Debye–Hückel electrostatics ──┘   └──── 12-10-6 van der Waals ────┘
a = acoef(index1, index2),  b = bcoef(index1, index2),  c = ccoef(index1, index2)
```

The `a/b/c` are **tabulated per type-pair** (Discrete2D functions indexed by atom-type index).
They encode the 12-10-6 well:

```
U_vdw(r) = ε·[ 13(R/r)¹² − 18(R/r)¹⁰ + 4(R/r)⁶ ]      # minimum −ε at r = R
```

## 3. The combination rule (`parse_cg_prm.py:232–246`)

For **every** type pair (the default), the coefficients are built from each type's own
`ε` and `Rmin/2`:

```python
epsilon = sqrt(epsilon1 * epsilon2)     # all types share ε = 0.000132 kcal/mol
R_min   = R_min1 + R_min2               # SUM rule (each Ri is that type's r_min/2)
a = 13*epsilon*R_min**12
b = -18*epsilon*R_min**10
c =  4*epsilon*R_min**6
```

**Only NBFIX pairs override this** (`parse_cg_prm.py:247–260`) — and NBFIX = the **native
(Go) contacts, which exist only inside the nascent chain.** Verified: there are **no NBFIX
entries between nascent types and any ribosome type**, so every nascent↔ribosome pair uses
the default sum rule (pure excluded volume).

**Takeaway:** a nascent bead always uses **its own `Rmin/2`** against any partner; the pair
distance is `R = Rmin/2_nascent + Rmin/2_partner`, `ε = 0.000132 kcal/mol`, 12-10-6 form.

### 3a. Same amino-acid name, two different `Rmin/2` (the source of the confusion)

topo uses **standard residue names** (`ALA`, `MET`, …) — not O'Brien's unique per-residue type
names (`A1`, `A2`, …). The per-residue-ness of the nascent chain is **not** carried by the name;
it is carried by a separate structure-derived array (`nascent_rmin2`). Consequently the mapping
**name → `Rmin/2` is context-dependent**:

- an **`ALA` in the nascent chain** → **K–B, per-residue** value (varies residue to residue;
  e.g. one ALA 0.29 nm, another 0.31 nm), derived from the nascent reference PDB;
- an **`ALA` in a ribosomal protein** → **O'Brien table** value (`OBRIEN_SC_RMIN2_NM["ALA"]`
  = 0.2862 nm, fixed, structure-independent).

Same name, two `Rmin/2`, because one is structure-derived and the other is a fixed table.
**This is exactly why a single `Rmin/2` cannot be stored under `ALA` in `model_parameters`** —
the name alone does not determine it (nascent → K–B from structure; ribosomal → O'Brien table).

Important: a given nascent residue `i` still has **one** `Rmin/2` (its K–B value); it uses that
same value against **every** partner (nascent or ribosome). It is not that particle `i` has two
values — it is the amino-acid *name* that maps to two different values depending on which
molecule (nascent vs ribosome) the bead belongs to.

So, restating the two partner cases:

- **intra-nascent** `i–j`: `Rmin/2(i)` and `Rmin/2(j)` **both** from the nascent native structure
  (K–B);
- **nascent `i` ↔ ribosomal-protein `j`**: `Rmin/2(i)` from the nascent native structure (K–B),
  `Rmin/2(j)` from the O'Brien table (structure-independent per-AA).

### 3b. Intra-nascent has two regimes (native vs non-native)

The "both from the native structure" statement above is precise for the **non-native** part.
Intra-nascent pairs split into:

- **native (Go) contacts** → **NBFIX**: the well minimum is the **actual native pair distance**
  measured in the reference PDB (not a sum of two `Rmin/2`);
- **non-native pairs** → **sum rule** `Rmin/2(i) + Rmin/2(j)`, both K–B (excluded volume).

Either way the values are structure-derived. Nascent↔ribosome pairs are **always** non-native
(no NBFIX), so they only ever use the sum rule.

---

## 4. Nascent-protein Rmin/2 = per-residue, structure-derived (K–B) = "Option A"

Nascent beads have **per-residue** atom types `A1..An` (`4c5c…go_bt.prm` NONBONDED), each with
its own `Rmin/2`:

```
A1  →  4.565666 Å = 0.45657 nm
A2  →  4.204688 Å = 0.42047 nm
A3  →  3.710752 Å = 0.37108 nm   ...
```

These are **structure-derived (Karanicolas–Brooks)**, not a per-AA lookup — verified by mapping
4c5c's 306 residues to their `A_i`: identical amino acids get **different** radii (22 ALA → 22
distinct values, 39 LEU → 39 values, etc.).

**In topo this is "Option A":** `build_nonbonded_interaction(return_rmin2=True)` →
`precompute_contacts` → `nascent_rmin2` → `append_ribosome`. It is computed from the nascent
structure and is **independent** of `model_parameters`.

---

## 5. Ribosomal-RNA Rmin/2 = per-type

Three RNA bead types (`rnc.prm` NONBONDED), one `Rmin/2` each:

| type | Rmin/2 (Å) | Rmin/2 (nm) | note |
|------|-----------|-------------|------|
| `P`  | 6.447660  | 0.644766    | phosphate (q = −1) |
| `R`  | 5.231399  | 0.523140    | ribose = **centroid of C1′,C2′,C3′,C4′,C5′** (5 carbons, no O4′) |
| `BR` | 5.342436  | 0.534244    | base bead (PU1/PU2/PY all map to `BR`) |

The ribose bead position is O'Brien's `create_cg_ribosome_model.py:326` (average of the five
ribose carbons). topo's `cg_ribosome.py` now uses the same set (was `C1′..C4′,O4′` before).

**In topo:** `model_parameters["topo"]` P/R/BR `Rmin_2` (set to the values above); read by
`load_ribosome` for the ribosome RNA beads.

---

## 6. Ribosomal-protein Rmin/2 = per-AA (generic), except L24

Each ribosomal-protein bead is typed by its **amino acid**:

- **All proteins except L24** (L2, L3, L4, L17, L22, L23, L29, L32, L34, …): type **`S<aa>`
  (per-AA)**. e.g. L23 beads = `SM, SI, SR, SE, SL, SK, SV, SA, SP, SH, …`. Every ALA → `SA`
  = 0.2862 nm, every LYS → `SK` = 0.3536 nm (20 S-types, 14 distinct values). This per-AA
  table = topo's `OBRIEN_SC_RMIN2_NM`.
- **L24 only**: **per-residue `B`-types** (`B1..B102`, 81 distinct values / 102 beads),
  structure-derived (file `combine_ribo_L24_Yang.prm`).

**In topo:** `load_ribosome` routes protein beads (name `CA`) to `OBRIEN_SC_RMIN2_NM` — exact
for every ribosomal protein **except L24**, which topo currently **approximates per-AA**
(decision: treat L24 like any other ribosomal protein; not per-residue).

---

## 7. Concrete examples (their files)

All pairs: `ε = 0.000132 kcal/mol`, `R = Rmin/2_nascent + Rmin/2_partner`, 12-10-6.

| nascent | partner | Rmin/2 nascent (nm) | Rmin/2 partner (nm) | R (nm) |
|---------|---------|---------------------|---------------------|--------|
| `A1` | rRNA phosphate `P`   | 0.45657 | 0.64477 | 1.10134 |
| `A1` | rRNA ribose `R`      | 0.45657 | 0.52314 | 0.97971 |
| `A1` | ribosomal Ala `SA`   | 0.45657 | 0.28623 | 0.74279 |
| `A1` | L24 bead `B1`        | 0.45657 | 0.42422 | 0.88079 | (topo approximates L24 per-AA) |

---

## 8. Summary — nascent-side vs ribosome-side parameter source

| bead population | Rmin/2 source | per-residue? | topo source |
|-----------------|---------------|--------------|-------------|
| nascent protein | K–B (structure) `A_i` | yes | Option A (`nascent_rmin2`, K–B) |
| ribosomal RNA (P/R/BR) | per-type | no | `model_parameters` P/R/BR |
| ribosomal protein (most) | per-AA `S<aa>` | no | `OBRIEN_SC_RMIN2_NM` |
| ribosomal protein L24 | per-residue `B_i` | yes | approximated per-AA (`OBRIEN_SC_RMIN2_NM`) |

ε is 0.000132 kcal/mol for every non-bonded pair; combination is the **sum rule**
`R = Rmin/2_i + Rmin/2_j`; native (Go) contacts are NBFIX and exist **only within the nascent
chain** (never nascent↔ribosome).

---

## 9. Open item (not yet done)

`model_parameters["topo"]` still carries a **static per-AA protein `Rmin_2`** that is **not**
used by any force (it only fed `rf_sigma` → `dumpForceFieldData`, which is never called). The
real nascent radius is the per-residue K–B value (§4). Proposed: delete the static protein
`Rmin_2`, and have the dump (if ever used) report the per-residue K–B `Rmin_2`. **Awaiting your
direction.**
