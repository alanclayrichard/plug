# fireprot ddG (fireprotdb, curated by thermompnn). one item per substitution, numbered
# against the uniprot sequence. comes with train/val/test, split by sequence.
import csv

from . import Benchmark
from .. import config as c


class Fireprot(Benchmark):
    name = "fireprot"
    # every item is about one wild type protein
    key = "pdb"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        for split in ("train", "val", "test"):
            path = c.data / "fireprot" / f"{split}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    if not (r.get("ddG") or "").strip(): continue
                    yield {"sequence": r["sequence"].strip().upper(),
                           "mutation": f"{r['wild_type']}{r['position']}{r['mutation']}",
                           "ddg": float(r["ddG"]),
                           "pdb": r.get("pdb_id_corrected", ""),
                           "split": split}

    # the benchmark also carries the pdb sequence of each protein, so hold that out too
    def sequences(self):
        yield from super().sequences()
        for split in ("train", "val", "test"):
            path = c.data / "fireprot" / f"{split}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    s = (r.get("pdb_sequence") or "").strip().upper()
                    if s: yield s
