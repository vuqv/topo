## Why

The disordered-region model does not work, and the reason is structural rather than
a matter of tuning.

Every IDR pair is evaluated by the same Karanicolas–Brooks 12-10-6 the folded Gō
model uses, with a single `eps` multiplying *both* the repulsive `r^-12` wall and
the attractive well. Turning on any attraction therefore inflates the bead: the
effective hard core (where `U = kT`) moves from **0.58 R to ~0.91 R**, roughly
**3.7× in excluded volume**, in exchange for a well shallower than `kT`. Excluded
volume wins, and the chain **expands** when attraction is added.

A 222-run benchmark over the 18 genuinely disordered proteins of a 24-protein SAXS
set measured the consequence. Fitting `Rg = R0 N^ν`:

- `ν` was pinned at **0.68–0.71** across `eps_gen` 0 → 5.02 kJ/mol — the knob moves
  only the prefactor `R0` (0.146 → 0.098), never the exponent.
- `ν` *rose* monotonically with `idr_scale` (0.605 → 0.723), i.e. every increment of
  sequence-specific attraction moved the scaling **away** from the experimental
  0.551.
- Both optima therefore sat at the corner `idr_scale = 0`, `eps_gen = 0` — not
  because attraction is unphysical, but because the functional form cannot express
  it.

Two further defects compound this:

1. **The 12-10-6 carries a desolvation barrier** (`+0.143·eps` at `1.45 R`). It
   models solvent expulsion when a *native* contact forms and has no counterpart
   between disordered residues. Because it sits at larger `r` than the well it
   dominates the second virial coefficient, pushing the θ point from `0.35 kT` out
   to `1.94 kT` and fighting the attraction across the whole usable range.
2. **`eps_gen_kj` is degenerate with `idr_scale`.** Both feed one well depth
   `eps_gen + idr_scale·eps_BT`, so any value of one trades for the other. Its only
   non-degenerate effect is to raise the *mean* depth without its *spread*, diluting
   the sequence contrast — the opposite of what the model needs.

## What Changes

**Change the functional form for IDR-involving pairs — two independent fixes, both
required** (derived in full in `design.md`):

1. *Decouple excluded volume from attraction.* An Ashbaugh–Hatch split holds the
   repulsive branch fixed and scales only the attractive branch, so `eps_ev_kj` sets
   the core alone and `idr_scale · eps_BT` the depth alone. The core then moves 3.2 %
   across the usable range instead of 56 %.
2. *Drop the desolvation barrier — use a 12-6, not a 12-10-6.* Decoupling alone is
   **not sufficient**: the barrier sits beyond the well, carries more `4πr²` weight
   in `B₂`, and makes the chain swell as attraction is added. It keeps the θ point at
   1.94 kT where a 12-6 puts it at 0.35 kT. The barrier represents solvent expulsion
   on forming a *native* contact — real physics for a Gō contact, meaningless between
   two disordered residues.

Together these make `idr_scale` a true solvent-quality dial: raising it *lowers* `ν`,
as attraction physically should. The prediction is directional and was tested — same
benchmark, same 18 proteins, `ν` rose with `idr_scale` under the old form and falls
under the new one.

**Partition the pairs across two forces** with disjoint interaction groups — 12-10-6
on `{folded}×{folded}`, Ashbaugh–Hatch on everything touching an IDR bead — each
force setting its own groups, and **no downstream code permitted to widen them**.

**Remove `eps_gen_kj`** and recalibrate the defaults.

**Fix the two downstream sites that violate group ownership** — `append_ribosome`
and the non-interacting-copy replication — both of which currently assert a
chain-wide interaction group over a force they did not create.

## Impact

Measured on the 18 IDPs, `idr_scale = 0.10`, `eps_ev_kj = 0.8368`:

| | ν | R0 | RMS | r | slope |
|---|---|---|---|---|---|
| this change | **0.566** | **0.223** | **12.0 %** | **0.89** | **0.81** |
| before | 0.591 | 0.202 | 12.9 % | 0.89 | 0.86 |
| HPS-Urry | 0.490 | 0.301 | 19.7 % | 0.70 | 0.68 |
| experiment | 0.551 | 0.244 | — | — | — |

`ν` becomes tunable and lands on target — no cell in the prior 222-run grid could do
this. A 3-seed confirmation at `idr_scale = 0.12` gives `ν = 0.559 ± 0.009`.

**Non-goals.** This does not close the sequence-discrimination gap: RMS 12.0 % is
still above the 9.5 % residual of a pure `Rg = 0.244 N^0.551` fit, so the model does
not yet beat chain length alone. Three chains at `N = 185` span 36 % experimentally
and 7 % in the model, in the wrong rank order. That is a separate problem, believed
to lie in the bonded terms or charge patterning rather than the contact channel.
