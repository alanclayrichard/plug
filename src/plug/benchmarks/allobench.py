# allobench: proteins with a known allosteric site and a known active site, one item per
# entry. ships no split of its own.
import ast
import csv
import re

from . import Benchmark
from .. import config as c
from ..sifts import THREE, sites


# a residue column arrives either as a list or as a string that looks like one
def terms(v):
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []


# "B-THR-7" -> ("B", "T", 7), ready to be renumbered against the uniprot sequence
def labels(v):
    out = []
    for lab in terms(v):
        m = re.match(r"^([A-Za-z0-9]+)-([A-Z]{3})-(-?\d+)", str(lab).strip())
        if m and m.group(2) in THREE: out.append((m.group(1), THREE[m.group(2)], int(m.group(3))))
    return out


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
                pdb, acc = r.get("allosteric_pdb", ""), r.get("pdb_uniprot", "")
                yield {"sequence": s,
                       "target_id": r.get("target_id", ""),
                       "gene": r.get("target_gene", ""),
                       "organism": r.get("organism", ""),
                       "uniprot": acc,
                       "pdb": pdb,
                       # the residues making up each site, counted along the sequence above
                       "allosteric_site": sites(labels(r.get("allosteric_site_residue", "")),
                                                pdb, acc, s),
                       "active_site": [p for p in terms(r.get("active_site_residue", ""))
                                       if isinstance(p, int) and 1 <= p <= len(s)]}
