# cafa5 gene ontology function prediction, one item per protein
'''
The CAFA5 function test set as packaged by bioreason-pro. The labels are gene
ontology terms across the three branches — biological process, molecular
function and cellular component — and the proteins were held back by date
rather than by similarity, which is exactly the kind of split a homology split
is meant to replace. It ships no split of its own. There are thousands of go
terms but only the ones with enough positives are ever scored, so only those
become strata: the rest would be three constraints apiece over a handful of
proteins, which costs solver time and buys nothing.
'''

import ast

import pandas as pd

from . import Benchmark
from .. import config as c

def terms(v):
    '''
    a go column arrives either as a list or as a string that looks like one.

    inputs:
        v: str | list -- the raw cell
    outputs:
        terms: list[str] -- the go terms it holds
    '''
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []

class Cafa5(Benchmark):
    '''go term prediction on a temporal holdout.'''

    name = "cafa5"
    # every item is about one protein
    key = "protein_id"
    published = ()

    label = "go"
    # a term with fewer positives than this is not stratified on
    min_positives = 20

    @classmethod
    def strata(cls, rows):
        '''
        one stratum per go term a protein carries, keeping only the terms with
        enough positives to be worth constraining.

        inputs:
            rows: list[dict] -- the benchmark's items
        outputs:
            strata: list[tuple[str, ...]] -- the terms each protein belongs to
        '''
        from collections import Counter
        n = Counter(t for it in rows for a in ("bp", "mf", "cc") for t in it[f"go_{a}"])
        keep = {t for t, c in n.items() if c >= cls.min_positives}
        return [tuple(sorted({t for a in ("bp", "mf", "cc") for t in it[f"go_{a}"]} & keep))
                for it in rows]

    @classmethod
    def rows(cls):
        '''
        one item per protein, with its go terms split by branch.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, protein_id, organism and the
                  three go branches
        '''
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
