# lp-pdbbind protein-ligand binding affinity, one item per pdb entry
'''
The published split is the leak-proof one from the paper, in the new_split
column; note it judges leakage on the protein and the ligand together, which is
far looser on the protein than the rule used here, and a few hundred rows sit
in none of its three splits. The label is about a protein and a ligand
together, so this benchmark declares the ligand view: the homology split has to
keep the chemistry apart as well, or a test ligand comes back as a train
molecule and the model has seen the answer.
'''

import csv
import sys

from . import Benchmark
from .. import config as c

# the sequences in this file are long enough to trip the default csv field limit
csv.field_size_limit(sys.maxsize)

class LpPdbbind(Benchmark):
    '''binding affinity for a protein and a ligand, with the paper's split.'''

    name = "lp-pdbbind"
    # every item is about one pdb entry
    key = "pdb"
    # what this benchmark is scored on, for the stratified split
    label = "value"
    published = ("train", "val", "test")
    # the ligand leaks too, not just the protein
    views = ("ligand",)

    @classmethod
    def rows(cls):
        '''
        one item per pdb entry with a sequence and a measured affinity.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, pdb, value, smiles and split
        '''
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

    @classmethod
    def ligands(cls):
        '''
        the molecules each pdb entry was measured against, so the split can
        keep the chemistry apart as well as the proteins.

        inputs:
            none
        outputs:
            ligands: dict[str, list[str]] -- pdb -> its smiles
        '''
        out = {}
        for it in cls.rows():
            if it["smiles"]: out.setdefault(it[cls.key], []).append(it["smiles"])
        return out
