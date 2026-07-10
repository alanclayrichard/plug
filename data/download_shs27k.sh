#!/bin/sh
# download the SHS27k STRING PPI subset (GNN-PPI / PIPR) into data/shs27k.
# SHS27k = 1,690 human proteins (<40% mutual identity) sampled from STRING v10.5, with the
# multi-type interaction network between them. two files: the action network (one row per
# ordered pair x interaction mode, with a STRING confidence score) + the sequence dictionary.
set -e
d="$(dirname "$0")/shs27k"
mkdir -p "$d"
base="https://raw.githubusercontent.com/lvguofeng/GNN_PPI/main/data"
curl -fsSL "$base/protein.actions.SHS27k.STRING.txt" -o "$d/protein.actions.SHS27k.STRING.txt"
curl -fsSL "$base/protein.SHS27k.sequences.dictionary.tsv" -o "$d/protein.SHS27k.sequences.dictionary.tsv"
