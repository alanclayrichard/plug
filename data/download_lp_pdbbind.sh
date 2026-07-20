#!/bin/sh
# download LP-PDBBind (protein-ligand binding affinity, leak-proof split) into data/lp-pdbbind.
# we only need LP_PDBBind.csv — its `seq` column holds each protein, `value` the affinity.
set -e
d="$(dirname "$0")/lp-pdbbind"
mkdir -p "$d"
curl -fsSL https://raw.githubusercontent.com/THGLab/LP-PDBBind/master/dataset/LP_PDBBind.csv -o "$d/LP_PDBBind.csv"
