#!/bin/sh
# download the cafa5 function test set (temporal holdout, as packaged by bioreason-pro)
# into data/cafa5
set -e
here="$(dirname "$0")"
d="$here/cafa5"
hf="$here/../.venv/bin/hf"; [ -x "$hf" ] || hf=hf   # prefer repo venv, else PATH
"$hf" download wanglab/bioreason-pro-test-data --repo-type dataset --local-dir "$d/test"
