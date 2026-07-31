#!/bin/sh
# download the megascale ddG benchmark (tsuboyama et al. 2023, as curated + split by thermompnn)
# into data/megascale. we pull the three csvs directly instead of cloning the whole repo.
# train/val/test are thermompnn's published split, held out by WT-protein cluster.
set -e
d="$(dirname "$0")/megascale"
mkdir -p "$d"
base="https://raw.githubusercontent.com/Kuhlman-Lab/ThermoMPNN/main/data_all"
curl -fsSL "$base/training/mega_train.csv" -o "$d/train.csv"
curl -fsSL "$base/training/mega_val.csv"   -o "$d/val.csv"
curl -fsSL "$base/testing/mega_test.csv"   -o "$d/test.csv"
