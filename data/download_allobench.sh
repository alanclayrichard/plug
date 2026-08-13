#!/bin/sh
# download allobench (proteins with known allosteric + active sites) into data/allobench.
# we only need AlloBench.csv — its `sequence` column holds each test protein.
set -e
d="$(dirname "$0")/allobench"
mkdir -p "$d"
curl -fsSL https://raw.githubusercontent.com/djmaity/allobench/HEAD/AlloBench.csv -o "$d/AlloBench.csv"
# its site residues are numbered like the pdb file, so we need sifts to renumber them
sh "$(dirname "$0")/download_sifts.sh"
