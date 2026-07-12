# PULCHRA (vendored)

This directory ships a prebuilt **PULCHRA** binary — the *PowerfUL CHain
Restoration Algorithm* (Rotkiewicz & Skolnick), which backmaps a coarse-grained
(Cα) structure to an all-atom model. TOPO uses it only when you reconstruct
atomistic coordinates from a coarse-grained trajectory; it is **optional**.

## When to use this binary

This is the **fallback** copy — a guaranteed-to-work binary for when you can't
build PULCHRA from source. The preferred ways to obtain PULCHRA, in order, are:

1. **Build from the upstream source** — `pulchra_306.tgz` from
   `http://www.pirx.com/downloads/pulchra_306.tgz` (often unreachable).
2. **Compile from the GitHub mirror** —
   [`euplotes/pulchra`](https://github.com/euplotes/pulchra) (`git clone`, then
   `cc -O3 pulchra.c pulchra_data.c -lm -o pulchra`).
3. **Use this vendored binary** — if neither source build works.

`scripts/install_deps.sh` walks exactly this order automatically.

## Provenance

| | |
|---|---|
| Program | PULCHRA 3.06 — all-atom reconstruction from reduced protein models |
| Reference | P. Rotkiewicz & J. Skolnick, *J. Comput. Chem.* **29**, 1460–1465 (2008) |
| Copyright | © 2000–2009 Piotr Rotkiewicz |
| Platform | **Linux x86-64** ELF, dynamically linked (glibc + libm only) |

**Platform note:** this is a **Linux x86-64** binary and will not run on any
other OS/architecture. On macOS, Windows, or non-x86-64 Linux you must build
PULCHRA from source (see above) and point TOPO at it.

## License

PULCHRA is distributed under the **MIT License** (see `LICENSE` in this
directory), so redistribution is unrestricted — this binary may be shipped in
the repository freely. Cite Rotkiewicz & Skolnick (2008) in work that uses it.

## Making TOPO use it

TOPO resolves PULCHRA at run time (see `topo/utils/external.py`):
`$TOPO_PULCHRA` → the program on `$PATH` → a copy vendored in `topo/bin/`.

```bash
export PATH="/path/to/topo/assets/pulchra:$PATH"
# ...or point TOPO straight at the binary:
export TOPO_PULCHRA="/path/to/topo/assets/pulchra/pulchra"
```
