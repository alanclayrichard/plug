#!/bin/sh
# download the fireprot ddG benchmark (fireprotdb, as curated + split by thermompnn) into
# data/fireprot. we pull the three csvs directly instead of cloning the whole repo.
# train/val/test are thermompnn's published split, held out by sequence cluster.
set -e
d="$(dirname "$0")/fireprot"
mkdir -p "$d"
base="https://raw.githubusercontent.com/Kuhlman-Lab/ThermoMPNN/main/data_all"
curl -fsSL "$base/training/fireprot_train.csv" -o "$d/train.csv"
curl -fsSL "$base/training/fireprot_val.csv"   -o "$d/val.csv"
curl -fsSL "$base/testing/fireprot_test.csv"   -o "$d/test.csv"
