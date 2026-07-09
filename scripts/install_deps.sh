#!/usr/bin/env bash
#
# install_deps.sh -- build/install the external binaries topo relies on:
#   * STRIDE  (secondary-structure assignment; required for contact/H-bond builds)
#   * PULCHRA (coarse-grained -> all-atom backmapping; optional)
#
# topo does NOT bundle these binaries (compiled C, platform-specific, and STRIDE's
# redistribution terms are restrictive). This script fetches and builds them from
# source, or installs STRIDE from bioconda if conda is available.
#
# topo locates the results at runtime in this order (see topo/utils/external.py):
#   1. $TOPO_STRIDE / $TOPO_PULCHRA  (explicit path)
#   2. the program on $PATH
#   3. topo/bin/<name>               (only if you install there; see PREFIX below)
#
# Usage:
#   scripts/install_deps.sh                 # install into $HOME/.local/bin
#   PREFIX=/path/to/topo/bin scripts/install_deps.sh   # vendor into the package
#   scripts/install_deps.sh stride          # just STRIDE
#   scripts/install_deps.sh pulchra         # just PULCHRA
#
# NOTE: verify the source URLs below against the upstream projects before use --
# they are the canonical distribution points at time of writing but may change.

set -euo pipefail

PREFIX="${PREFIX:-$HOME/.local/bin}"

# Upstream source archives. Override via env if these move.
STRIDE_URL="${STRIDE_URL:-https://webclu.bio.wzw.tum.de/stride/stride.tar.gz}"
PULCHRA_URL="${PULCHRA_URL:-https://www.pirx.com/downloads/pulchra/pulchra_306.tgz}"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$PREFIX"

install_stride() {
    # Prefer bioconda if conda is available -- no redistribution/build concerns.
    if command -v conda >/dev/null 2>&1; then
        echo ">> Installing STRIDE from bioconda ..."
        conda install -y -c bioconda stride
        return
    fi
    echo ">> Building STRIDE from source ($STRIDE_URL) ..."
    cd "$WORK"
    curl -fSL "$STRIDE_URL" -o stride.tar.gz
    mkdir -p stride && tar -xzf stride.tar.gz -C stride
    # The tarball may or may not have a top-level dir; find the Makefile.
    src="$(dirname "$(find stride -name Makefile -print -quit)")"
    make -C "$src"
    cp "$src/stride" "$PREFIX/stride"
    echo ">> STRIDE installed to $PREFIX/stride"
}

install_pulchra() {
    echo ">> Building PULCHRA from source ($PULCHRA_URL) ..."
    cd "$WORK"
    curl -fSL "$PULCHRA_URL" -o pulchra.tgz
    mkdir -p pulchra && tar -xzf pulchra.tgz -C pulchra
    src="$(dirname "$(find pulchra -name 'pulchra*.c' -print -quit)")"
    # PULCHRA ships a single C file; compile with the standard flags from its README.
    cc -O3 -o "$PREFIX/pulchra" "$src"/pulchra*.c -lm
    echo ">> PULCHRA installed to $PREFIX/pulchra"
}

targets=("${@:-stride pulchra}")
for t in ${targets[@]}; do
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
echo "    export TOPO_STRIDE=\"$PREFIX/stride\""
echo "    export TOPO_PULCHRA=\"$PREFIX/pulchra\""
