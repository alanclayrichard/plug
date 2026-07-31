# allobench: proteins with a known allosteric site and a known active site, one item per
# entry. ships no split of its own.
import ast
import csv

from . import Benchmark
from .. import config as c


# a residue column arrives either as a list or as a string that looks like one
def terms(v):
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []


class Allobench(Benchmark):
    name = "allobench"
    # every item is about one target protein
    key = "target_id"
    published = ()

    @classmethod
    def rows(cls):
        with open(c.data / "allobench" / "AlloBench.csv", newline="") as f:
            for r in csv.DictReader(f):
                s = (r.get("sequence") or "").strip().upper()
                if not s: continue
                yield {"sequence": s,
                       "target_id": r.get("target_id", ""),
                       "gene": r.get("target_gene", ""),
                       "organism": r.get("organism", ""),
                       "uniprot": r.get("pdb_uniprot", ""),
                       # the residues making up each site
                       "allosteric_site": terms(r.get("allosteric_site_residue", "")),
                       "active_site": terms(r.get("active_site_residue", ""))}
