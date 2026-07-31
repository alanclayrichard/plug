# megascale ddG (tsuboyama et al. 2023, curated by thermompnn). one item per single
# substitution, and the ddg sign follows thermompnn. comes with train/val/test, split
# by wild type protein.
import csv

from . import Benchmark
from .. import config as c


class Megascale(Benchmark):
    name = "megascale"
    # every item is about one wild type protein
    key = "pdb"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        for split in ("train", "val", "test"):
            path = c.data / "megascale" / f"{split}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    mut, ddg = r.get("mut_type", "").strip(), r.get("ddG_ML", "").strip()
                    if not mut or mut == "wt" or ddg in ("", "-"): continue
                    # keep single substitutions only
                    if any(x in mut for x in ("ins", "del", ":")): continue
                    yield {"sequence": r["wt_seq"].strip().upper(),
                           "mutation": mut,
                           "ddg": -float(ddg),
                           "pdb": r.get("WT_name", "").replace(".pdb", ""),
                           "split": split}

    # the mutant sequences are part of the benchmark too, so hold those out as well
    def sequences(self):
        yield from super().sequences()
        for split in ("train", "val", "test"):
            path = c.data / "megascale" / f"{split}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    s = (r.get("aa_seq") or "").strip().upper()
                    if s: yield s
