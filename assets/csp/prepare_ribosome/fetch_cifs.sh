#!/usr/bin/env bash
# Download the raw mmCIF depositions needed to (re)generate the reference
# ribosomes. These files are large (18-35 MB each) and are NOT shipped with the
# package -- they land in ./raw/ (git-ignored). Only the processed CG structures
# under ./structures/ are shipped.
#
# Usage:  bash fetch_cifs.sh            # all organisms
#         bash fetch_cifs.sh ecoli      # one organism
set -euo pipefail
cd "$(dirname "$0")"

# organism -> "PDBID:localname [PDBID:localname ...]"
declare -A CIFS=(
  [ecoli]="4V9D:ecoli/4v9d 5JTE:ecoli/5jte"   # 50S (4V9D) + A-site tRNA donor (5JTE)
  [yeast]="6Q8Y:yeast/6q8y"                    # 60S + P/E tRNA (A-site grafted from 8G61)
  [ncrassa]="7R81:ncrassa/7r81"                # 60S + native A/P tRNA (drop nascent chain)
  [human]="8G61:human/8g61"                    # 60S + native A/P tRNA
)
# yeast graft also needs the human cif as the A-site tRNA donor
[[ "${1:-}" == "yeast" ]] && CIFS[yeast]="6Q8Y:yeast/6q8y 8G61:human/8g61"

organisms=("${@:-ecoli yeast ncrassa human}")
for org in ${organisms[@]}; do
  for spec in ${CIFS[$org]}; do
    pdb="${spec%%:*}"; dest="raw/${spec##*:}.cif"
    mkdir -p "$(dirname "$dest")"
    if [[ -f "$dest" ]]; then echo "have $dest"; continue; fi
    echo "downloading $pdb -> $dest"
    curl -sfL "https://files.rcsb.org/download/${pdb}.cif.gz" | gunzip > "$dest"
  done
done
echo "done. raw cifs are in ./raw/ (git-ignored)."
