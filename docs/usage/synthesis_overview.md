# Protein synthesis: overview

TOPO can grow a protein **vectorially, N-terminus first**, out of a ribosome while
the nascent chain folds *as it emerges* — reproducing co-translational folding
under codon-resolved kinetics. This page is the map: the biology being modelled,
the two ribosome models TOPO offers, and which to use. The detailed references are
linked at the end.

## The biology in one minute

A ribosome synthesizes a protein one residue at a time, N→C, threading the growing
("nascent") chain out through the ~80 Å **exit tunnel**, where it begins to fold
before the full sequence exists. How long the ribosome dwells on each codon sets how
much time each segment has to sample conformations before the next residue arrives;
slow (rare) codons act as translational pauses that can change folding outcomes.

Each residue is added through one **elongation cycle** — (1) codon-dependent
aminoacyl-tRNA **decoding** at the A site (usually rate-limiting), (2) **peptidyl
transfer** at the peptidyl-transferase centre (PTC), (3) **translocation** (A→P→E).
TOPO times each residue from its mRNA codon and, in the explicit model, splits the
cycle into these three sub-stages.

## Two ribosome models

Both grow the same structure-based (Gō-like) Cα nascent chain — see
{doc}`model_theory` for the base force field — and both use the same codon kinetics.
They differ only in **how the exit tunnel is represented**:

| | Explicit CG ribosome | Analytic cylinder |
|---|---|---|
| **Runner** | `topo-csp` ({doc}`continuous_synthesis`) | `topo-cylinder` ({doc}`cylinder_synthesis`) |
| **Tunnel** | a real truncated, coarse-grained large subunit (~4,600 rigid beads) | an analytic cylindrical bore through an infinite wall |
| **Interactions** | ribosome↔nascent **excluded volume + electrostatics**; tRNA tether; A→P translocation | **geometric confinement only** |
| **Needs a structure?** | yes — a `ribosome_trunc.pdb` ({doc}`ribosome_preparation`) | no |
| **Per residue** | three MD sub-stages (transfer / translocation / wait) | one MD segment |
| **Cost / robustness** | heavier; realistic tunnel wall | fast; never jams |

## Which to use

- **Cylinder** — for fast exploration of how *tunnel geometry + codon kinetics* shape
  folding, or when you have no ribosome structure. Simplest starting point.
- **Explicit** — when the **tunnel-wall charge, tunnel shape** (constriction site /
  vestibule), or **translocation-coupled forces** matter to your question.

```{warning}
**The two models are comparable only in the *mean* per-residue dwell time, not in
confinement chemistry.** The cylinder omits the ribosome's electrostatics and
surface excluded volume and has a uniform straight bore. Do not compare folding
observables (folding order, Q-vs-length, radius of gyration) across the two models
without accounting for those missing terms. See {doc}`cylinder_synthesis` §"Physical
scope".
```

## Where to go next

* {doc}`ribosome_preparation` — get or build the truncated CG ribosome the explicit
  model needs (four ready-made organisms ship with TOPO).
* {doc}`codon_dwell_times` — the per-codon dwell-time tables that set the timing
  (E. coli, yeast, human, N. crassa ship with TOPO), and the `fastest`/`slowest` mRNA.
* {doc}`continuous_synthesis` — the explicit-ribosome runner: the RNC force field,
  the tRNA tether/PTC anchors, the 3-stage mechanics, the codon→MD-step kinetics.
* {doc}`cylinder_synthesis` — the analytic-tunnel runner.
* {doc}`synthesis_control` — the `csp.ini` / `cylinder.ini` parameter reference.
* {doc}`../tutorials/07_translation_cylinder` and
  {doc}`../tutorials/08_ribosome_synthesis` — hands-on, ready-to-run examples.
