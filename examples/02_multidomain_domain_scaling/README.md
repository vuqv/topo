# Tutorial 2 — Multidomain proteins & per-domain contact scaling

**Goal:** simulate a multidomain protein and control the contact-energy strength
**within** each domain and **across** domain interfaces using a `domain.yaml`
file. We use **adenylate kinase** (`1AKE` chain A, 214 residues), a textbook
multidomain enzyme.

**Time:** ~4 seconds on a CPU.

**Prerequisite:** do [Tutorial 1](https://vuqv.github.io/topo/tutorials/01_single_domain.html) first.

---

## Files in this folder

| File | Role |
|------|------|
| `1AKE_A.pdb` | Input structure (adenylate kinase, chain A). |
| `domain.yaml` | **Domain definition** — the new ingredient in this tutorial. |
| `md.ini` | Simulation config; note the added `domain_def = domain.yaml` line. |
| `run_simulation.py` | The runner script. |

## Why domains matter

In a structure-based model, every native contact gets an attractive well. By
default all wells have the same relative strength. But a multidomain protein
often needs **different stabilities** for each domain and for the interface
between them (e.g. to reproduce experimental melting temperatures, or to study
domain-wise unfolding). `domain.yaml` lets you set a **scale factor** on the
sidechain–sidechain contact energy:

- **within** each domain (`intra_domains[...].strength`), and
- **between** domains (`inter_domains`).

## The key idea: discontiguous domains

Adenylate kinase's CORE domain is **not a single stretch of sequence**. It is
made of two segments — residues **1–117** and **166–214** — that fold together
around the inserted NMP-binding domain (**118–165**). TOPO handles this by
letting one domain list **multiple residue ranges**.

Open `domain.yaml`:
```yaml
n_residues: 214
intra_domains:
  A:
    residues: [1-117, 166-214]   # ONE domain, TWO sequence segments
    strength: 1.1556
  B:
    residues: [118-165]
    strength: 1.6871
inter_domains:
  A-B: 1.8611
```

Reading it out:
- **Domain A** = residues 1–117 *and* 166–214; all A–A native contacts
  (including contacts *between* the two segments) scaled by **1.1556**.
- **Domain B** = residues 118–165; B–B contacts scaled by **1.6871**.
- **A–B interface** contacts scaled by **1.8611**.

> **Note:** an unspecified interface defaults to scale **1.0** (native contacts
> kept, unscaled) — the scaling only *modulates* contacts that already exist in
> the structure. To intentionally **decouple** two domains, set their pair to
> **`0.0`** explicitly. Full rules:
> [Domain definition file](https://vuqv.github.io/topo/usage/domain_definition.html).

## Step-by-step

### 1. Note the one change in `md.ini`
Compared with Tutorial 1, the only addition is:
```ini
domain_def = domain.yaml
```
That single line switches on domain-aware contact scaling.

### 2. Run it
```bash
python run_simulation.py -f md.ini
```
In the build log you'll see `Domain definition file: domain.yaml` and
`Built non-bonded interaction matrices: (214, 214), (214, 214)`, confirming the
214×214 contact matrices were built with your scaling applied.

### 3. Outputs
Same set as Tutorial 1, prefixed with `protein_code = 1AKE_A`
(`1AKE_A.log`, `1AKE_A.dcd`, `1AKE_A.chk`, `1AKE_A.psf`,
`1AKE_A_init.pdb`, `1AKE_A_final.pdb`, `1AKE_A_stride.dat`).

## Build your own `domain.yaml`

1. Decide how many domains and which residues belong to each (a discontiguous
   domain just gets several ranges in its `residues` list).
2. Set `n_residues` to the true residue count (must cover all residues; any
   residue you forget is auto-assigned to a fallback domain `X` with strength 1.0).
3. Give every domain a numeric `strength` (a blank value errors).
4. Add an `inter_domains` entry only for interfaces whose strength you want to
   change from the default of 1.0 (set `0.0` to remove an interface; omitted
   pairs stay at 1.0).

The reference doc has six worked scenarios (single domain, contiguous,
discontiguous, decoupled, 3+ domains, partial assignment):
[Domain definition file](https://vuqv.github.io/topo/usage/domain_definition.html).

## Try next

- Make domain B much weaker (e.g. `strength: 0.5`), run, and watch B unfold
  first at elevated temperature while A stays folded.
- Set `inter_domains: A-B: 0.0` to decouple the domains and observe them moving
  independently.
- Continue to **Tutorial 3** to learn restarting and the output files in depth.
