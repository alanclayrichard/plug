#!/bin/sh
# download passerrank allosteric proteins: pull the asd table, then fetch each unique
# pdb_uniprot accession's sequence from uniprot into one fasta.
set -e
d="$(dirname "$0")/passerrank"
mkdir -p "$d"
curl -fsSL https://raw.githubusercontent.com/smu-tao-group/PASSerRank/HEAD/data/source_data/ASD_Release_201909_AS.txt -o "$d/ASD_Release_201909_AS.txt"
awk -F'\t' 'NR>1 && $4 ~ /^[A-Z0-9]+$/ {print $4}' "$d/ASD_Release_201909_AS.txt" | sort -u > "$d/accessions.txt"
: > "$d/passerrank.fasta"
while read -r acc; do
  curl -fsSL "https://rest.uniprot.org/uniprotkb/$acc.fasta" >> "$d/passerrank.fasta" || true
done < "$d/accessions.txt"
# its site residues are numbered like the pdb file, so we need sifts to renumber them
sh "$(dirname "$0")/download_sifts.sh"
