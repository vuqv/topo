# Tutorial A.7 — A protein with intrinsically disordered regions (IDRs)

**Goal:** simulate a protein that is **part folded, part disordered**, and see the
disordered regions behave as flexible, non-folding chains while the folded core
keeps its structure. This is the only new ingredient beyond Tutorial A.2: a
`disordered:` section in `domain.yaml`.

**Time:** ~2 minutes on a CPU for the 100 000-step demo (365 residues).

**Prerequisite:** do [Tutorial A.1](A1_single_domain.md) and
[Tutorial A.2](A2_multidomain.md) first — this tutorial reuses the same runner and
the same `domain.yaml` file, adding one section.

---

## The system

`P28566` (5-HT<sub>1E</sub> receptor, 365 residues) modelled from its **AlphaFold**
structure. We treat it as **one folded domain plus two intrinsically disordered
regions (IDRs)**:

- **folded core** — most of the chain (domain `A`), held together by native contacts;
- **N-terminal tail** (residues 1–18) — disordered;
- **long intracellular loop** (residues 218–281) — disordered.

An IDR is not a fold, so its equilibrium shape is set by the *potential*, not by the
AlphaFold coordinates: TOPO removes the native contacts of the disordered residues
and replaces them with a weak, non-specific attraction, so those regions relax into
a flexible ensemble while the core stays folded.

## Files in this folder

| File | Role |
|------|------|
| `P28566.pdb` | Input structure (AlphaFold model; TOPO keeps only the CA atoms). |
| `domain.yaml` | Folded domain `A` **plus** a `disordered:` section marking the IDRs. |
| `md.ini` | Simulation configuration (steps, temperature, I/O, hardware). |
| `P28566_stride.dat` | Precomputed STRIDE H-bonds (so the run does not call STRIDE). |

## The one new thing — the `disordered:` section

Open `domain.yaml`. Compared with a normal folded `domain.yaml`, the only addition
is the `disordered:` block:

```yaml
n_residues: 365
intra_domains:
  A:
    residues: [18-164, 177-219, 279-363]
    nscale: 2.5044
disordered:
  residues: [1-18, 218-281]     # N-terminal tail + long intracellular loop
  eps_gen_kj: 2.25              # generic cohesion, kJ/mol (default)
  idr_scale: 1.0                # sequence-dependent BT scale (default)
```

- **`residues`** — the residues to treat as disordered. Their native contacts
  (H-bond, backbone–sidechain, and sidechain–sidechain) are removed.
- **`eps_gen_kj`** — the sequence-independent generic cohesion in kJ/mol (optional,
  default `2.25`). This is the knob to turn first when tuning compaction.
- **`idr_scale`** — the scale on the sequence-dependent BT energy, which varies by
  residue-pair type (optional, default `1.0`).
  Both defaults are calibrated against SAXS radii of gyration for 24 disordered
  proteins. Set **both** to `0` for a pure self-avoiding chain — zeroing one alone
  leaves the other channel on. Disordered↔folded pairs feel only excluded volume
  (steric) regardless.

> **You choose the disordered residues.** TOPO applies exactly the set you list —
> it does not detect disorder. Sources to inform the choice include
> [MobiDB](https://mobidb.org/) and the AlphaFold **pLDDT < 70** proxy (read from
> the B-factor column). See the full write-up:
> [Disordered / IDR regions](../usage/disordered_regions.rst).

> **Overlap is allowed — disorder wins.** Here domain `A` and the IDR ranges touch
> at a few residues (18, 218–219, 279–281). The run prints an info line
> (`IDR overlap: residues [...] are in both a domain and 'disordered:'; treating as
> disordered.`) and those residues are treated as disordered. This is a feature:
> you can define a domain broadly and carve a disordered loop out of it.

## Run it

```bash
topo-mdrun -f md.ini            # == python -m topo.mdrun -f md.ini
```

Early in the build you will see the overlap info line and the contact count. With
the IDRs masked, `P28566` has **823 native contacts** (the disordered residues
contribute none). Outputs land in `traj/` exactly as in Tutorial A.1.

> The disordered regions start off the model's energy minimum (their native
> contacts are gone), so they carry some force at step 0 and relax over the first
> steps. This is expected — discard the initial relaxation in analysis, as for any
> MD run. (`minimize = no`; see the note in `md.ini`.)

## What to look for

The folded core keeps its shape while the tail and loop swing freely. Aligning each
frame on the folded core and measuring the per-residue fluctuation (RMSF) over the
100-frame demo makes this quantitative:

| Region | mean RMSF |
|--------|-----------|
| Folded core (19–164) | ~0.7 Å |
| N-terminal IDR (1–18) | ~17 Å |
| Intracellular-loop IDR (218–281) | ~26 Å |

i.e. the disordered regions are **~20–35× more mobile** than the core — flexible,
but still attached and (at the default attraction) compact rather than fully
extended. Watch it in VMD:

```bash
vmd traj/traj.psf traj/traj.dcd
```

## Try next

- **Self-avoiding vs collapsing.** Set **both** `eps_gen_kj: 0` and `idr_scale: 0`
  in `domain.yaml` and rerun; the IDRs expand further. Compare the radius of gyration
  of a disordered region against the default run.
- **A fully disordered protein.** Drop `intra_domains` and list every residue under
  `disordered:` — a whole-IDP run (see the *fully-IDP* note in the
  [Disordered / IDR regions](../usage/disordered_regions.rst) page).
- **Calibrate the core.** The `nscale = 2.5044` here is a placeholder; use the
  optimizer ([Tutorial A.5](A5_opt_nscal.md)) to calibrate it **with the
  `disordered:` section active**, so the folded core is tuned under the real
  (IDR-reduced) energy.
- **Longer sampling.** The demo is short; a converged IDR ensemble needs far more
  steps (the validation run used `md_steps = 1_000_000`) and is best on a GPU
  (`device = GPU`).

> Full theory and the exact energy/well-position rules for each pair class
> (folded–folded, IDR–IDR, IDR–folded) are in
> [Disordered / IDR regions](../usage/disordered_regions.rst).
