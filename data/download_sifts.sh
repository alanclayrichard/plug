#!/bin/sh
# download the sifts table that lines pdb residue numbers up with uniprot ones, into
# data/sifts. the allosteric benchmarks need it to renumber their site labels.
set -e
d="$(dirname "$0")/sifts"
mkdir -p "$d"
curl -fsSL https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/uniprot_segments_observed.csv.gz -o "$d/uniprot_segments_observed.csv.gz"
