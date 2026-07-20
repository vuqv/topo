#!/usr/bin/env bash
#
# install_deps.sh -- build/install the external binaries topo relies on:
#   * STRIDE  (secondary-structure assignment; required for contact/H-bond builds)
#   * PULCHRA (coarse-grained -> all-atom backmapping; optional)
#
# topo does NOT bundle these binaries in the wheel (compiled C, platform-specific,
# and STRIDE's redistribution terms are restrictive). For each tool this script
# tries, in order: (1) build from the upstream source, (2) build from a GitHub
# source mirror, (3) fall back to the binary vendored under assets/ (Linux x86-64).
#
# topo locates the results at runtime in this order (see topo/utils/external.py):
#   1. $TOPO_STRIDE / $TOPO_PULCHRA  (explicit path)
#   2. the program on $PATH
#   3. topo/bin/<name>               (only if you install there; see PREFIX below)
#
# Usage:
#   scripts/install_deps.sh                 # STRIDE only, into $HOME/.local/bin
#   scripts/install_deps.sh pulchra         # PULCHRA only (opt-in)
#   scripts/install_deps.sh stride pulchra  # both
#   PREFIX=/path/to/topo/bin scripts/install_deps.sh   # vendor into the package
#
# NOTE: verify the source URLs below against the upstream projects before use --
# they are the canonical distribution points at time of writing but may change.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-$HOME/.local/bin}"

# Upstream source archives. Override via env if these move.
STRIDE_URL="${STRIDE_URL:-https://webclu.bio.wzw.tum.de/stride/stride.tar.gz}"
# Maintained STRIDE source mirror (secondary, used if the upstream site is down).
STRIDE_GIT="${STRIDE_GIT:-https://github.com/MDAnalysis/stride.git}"
PULCHRA_URL="${PULCHRA_URL:-http://www.pirx.com/downloads/pulchra_306.tgz}"
# Maintained PULCHRA source mirror (secondary, used if the upstream site is down).
PULCHRA_GIT="${PULCHRA_GIT:-https://github.com/euplotes/pulchra.git}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$PREFIX"

# Build STRIDE from a source tarball URL ($1). Returns non-zero on any failure
# (e.g. the upstream site being unreachable) so callers can fall through.
build_stride_from_tarball() {
    echo ">> Building STRIDE from source ($1) ..."
    ( cd "$WORK" \
        && curl -fSL "$1" -o stride.tar.gz \
        && mkdir -p stride_tar && tar -xzf stride.tar.gz -C stride_tar \
        && src="$(dirname "$(find stride_tar -name Makefile -print -quit)")" \
        && [ -n "$src" ] && make -C "$src" \
        && cp "$src/stride" "$PREFIX/stride" )
}

# Build STRIDE from the MDAnalysis/stride GitHub mirror ($1). Returns non-zero
# on failure.
build_stride_from_git() {
    command -v git >/dev/null 2>&1 || { echo ">> git not available; skipping mirror."; return 1; }
    echo ">> Building STRIDE from the MDAnalysis/stride mirror ($1) ..."
    ( cd "$WORK" \
        && git clone --depth 1 "$1" stride_git \
        && src="$(dirname "$(find stride_git -name Makefile -print -quit)")" \
        && [ -n "$src" ] && make -C "$src" \
        && cp "$src/stride" "$PREFIX/stride" )
}

install_stride() {
    # STRIDE resolution preference:
    #   1. Build from the upstream STRIDE source ($STRIDE_URL).
    #   2. Build from the MDAnalysis/stride GitHub mirror ($STRIDE_GIT).
    #   3. Fall back to the validated binary vendored at assets/stride/.
    if build_stride_from_tarball "$STRIDE_URL"; then
        echo ">> STRIDE installed to $PREFIX/stride (upstream source)"; return
    fi
    echo ">> Upstream STRIDE source unavailable; trying the GitHub mirror ..."
    if build_stride_from_git "$STRIDE_GIT"; then
        echo ">> STRIDE installed to $PREFIX/stride (GitHub mirror)"; return
    fi
    echo ">> Source builds failed; falling back to the vendored STRIDE binary ..."
    vendored="$REPO_ROOT/assets/stride/stride"
    if [[ -x "$vendored" ]]; then
        cp "$vendored" "$PREFIX/stride" && chmod +x "$PREFIX/stride"
        echo ">> STRIDE installed to $PREFIX/stride (vendored fallback)"; return
    fi
    echo ">> ERROR: could not build STRIDE from source and no vendored binary found." >&2
    return 1
}

# Compile PULCHRA from a source tree rooted at $1 (the dir holding pulchra.c and
# its *.h data files). Builds inside that dir so the headers resolve. The glob
# picks up pulchra.c + pulchra_data.c per PULCHRA's documented build command.
build_pulchra_in() {
    local src
    src="$(dirname "$(find "$1" -name 'pulchra.c' -print -quit)")"
    [ -n "$src" ] || return 1
    ( cd "$src" && cc -O3 pulchra*.c -lm -o "$PREFIX/pulchra" )
}

install_pulchra() {
    # PULCHRA resolution preference:
    #   1. Build from the upstream source ($PULCHRA_URL).
    #   2. Build from the euplotes/pulchra GitHub mirror ($PULCHRA_GIT).
    #   3. Fall back to the (MIT-licensed) binary vendored at assets/pulchra/.
    echo ">> Building PULCHRA from source ($PULCHRA_URL) ..."
    if ( cd "$WORK" && curl -fSL "$PULCHRA_URL" -o pulchra.tgz \
            && mkdir -p pulchra_tar && tar -xzf pulchra.tgz -C pulchra_tar ) \
            && build_pulchra_in "$WORK/pulchra_tar"; then
        echo ">> PULCHRA installed to $PREFIX/pulchra (upstream source)"; return
    fi
    echo ">> Upstream PULCHRA source unavailable; trying the GitHub mirror ..."
    if command -v git >/dev/null 2>&1 \
            && ( cd "$WORK" && git clone --depth 1 "$PULCHRA_GIT" pulchra_git ) \
            && build_pulchra_in "$WORK/pulchra_git"; then
        echo ">> PULCHRA installed to $PREFIX/pulchra (GitHub mirror)"; return
    fi
    echo ">> Source builds failed; falling back to the vendored PULCHRA binary ..."
    vendored="$REPO_ROOT/assets/pulchra/pulchra"
    if [[ -x "$vendored" ]]; then
        cp "$vendored" "$PREFIX/pulchra" && chmod +x "$PREFIX/pulchra"
        echo ">> PULCHRA installed to $PREFIX/pulchra (vendored fallback)"; return
    fi
    echo ">> ERROR: could not build PULCHRA from source and no vendored binary found." >&2
    return 1
}

# Default target is STRIDE alone: it is the only binary topo needs to build a
# contact map. PULCHRA is opt-in -- only backmapping runs use it.
targets=("$@")
[ ${#targets[@]} -eq 0 ] && targets=(stride)
for t in "${targets[@]}"; do
    case "$t" in
        stride)  install_stride ;;
        pulchra) install_pulchra ;;
        *) echo "Unknown target: $t (expected 'stride' or 'pulchra')" >&2; exit 2 ;;
    esac
done

echo
echo "Done. Ensure $PREFIX is on your PATH, or point topo at the binaries:"
echo "    export PATH=\"$PREFIX:\$PATH\""
echo "    # or:"
echo "    export TOPO_STRIDE=\"$PREFIX/stride\"        # or TOPO_PULCHRA for pulchra"
