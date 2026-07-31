# LP-PDBBind protein and ligand binding affinity. the published split is the leak proof
# one from the paper, in the new_split column. note it judges leakage on the protein and
# the ligand together, which is far looser on the protein than the rule used here.
import csv
import sys

from . import Benchmark
from .. import config as c

# the sequences in this file are long enough to trip the default csv field limit
csv.field_size_limit(sys.maxsize)


class LpPdbbind(Benchmark):
    name = "lp-pdbbind"
    # every item is about one pdb entry
    key = "pdb"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        with open(c.data / "lp-pdbbind" / "LP_PDBBind.csv", newline="") as f:
            for r in csv.DictReader(f):
                s = (r.get("seq") or "").strip().upper()
                v = (r.get("value") or "").strip()
                if not s or not v: continue
                yield {"sequence": s,
                       # the pdb code sits in the unnamed first column
                       "pdb": r.get("", ""),
                       "value": float(v),
                       "smiles": (r.get("smiles") or "").strip(),
                       # a few hundred rows are in none of the three splits
                       "split": r.get("new_split", "")}
