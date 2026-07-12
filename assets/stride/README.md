# STRIDE (vendored)

This directory ships the exact **STRIDE** binary that TOPO's contact/H-bond
model was built and validated against. STRIDE performs secondary-structure and
backbone hydrogen-bond assignment; TOPO runs it (`stride -h`) when it builds the
native-contact list and when it locates helix/strand segments for mirror-image
detection.

## When to use this binary

This is the **fallback** copy — a guaranteed-to-work binary for when you can't
build STRIDE from source. The preferred ways to obtain STRIDE, in order, are:

1. **Build from the upstream source** —
   `https://webclu.bio.wzw.tum.de/stride/stride.tar.gz` (often unreachable).
2. **Compile from the GitHub mirror** —
   [`MDAnalysis/stride`](https://github.com/MDAnalysis/stride) (`git clone`, `make`).
3. **Use this vendored binary** — if neither source build works.

`scripts/install_deps.sh` walks exactly this order automatically.

## Provenance

| | |
|---|---|
| Program | STRIDE — knowledge-based secondary-structure assignment |
| Algorithm | D. Frishman & P. Argos, *Proteins* **23**, 566–579 (1995) |
| Binary built | 2026-07 (mtime `Nov 10 2022`) from the upstream STRIDE source |
| Platform | Linux **x86-64** ELF, dynamically linked (glibc + libm only) |

This is the unmodified academic STRIDE, not a fork or patched build.

**Platform note:** this is a Linux x86-64 binary. On any other OS/arch you must
build STRIDE yourself from source and point TOPO at it (see below).

## Making TOPO use it

TOPO resolves STRIDE at run time (see `topo/utils/external.py`) in this order:
`$TOPO_STRIDE` → the program on `$PATH` → a copy vendored in `topo/bin/`.

Pick whichever you prefer:

**Add this directory to your `PATH`:**

```bash
export PATH="/path/to/topo/assets/stride:$PATH"
```

**…or point TOPO straight at the binary:**

```bash
export TOPO_STRIDE="/path/to/topo/assets/stride/stride"
```

Confirm resolution:

```bash
python -c "from topo.utils.external import find_executable; print(find_executable('stride'))"
```

## License

STRIDE is distributed for academic use with restrictive redistribution terms.
It is included here for internal/academic reproducibility only. Cite Frishman &
Argos (1995) in any work that uses STRIDE assignments.
