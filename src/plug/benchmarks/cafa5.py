# cafa5 function prediction, as packaged by bioreason-pro. the labels are gene ontology
# terms; the proteins were held back by date rather than by similarity. ships no split
# of its own.
import ast

import pandas as pd

from . import Benchmark
from .. import config as c


# a go column arrives either as a list or as a string that looks like one
def terms(v):
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []


class Cafa5(Benchmark):
    name = "cafa5"
    # every item is about one protein
    key = "protein_id"
    published = ()

    @classmethod
    def rows(cls):
        path = c.data / "cafa5" / "test" / "data" / "test-00000-of-00001.parquet"
        for r in pd.read_parquet(path).itertuples():
            yield {"sequence": str(r.sequence).strip().upper(),
                   "protein_id": r.protein_id,
                   "organism": r.organism,
                   # the three go branches: biological process, molecular function,
                   # cellular component
                   "go_bp": terms(r.go_bp),
                   "go_mf": terms(r.go_mf),
                   "go_cc": terms(r.go_cc)}
