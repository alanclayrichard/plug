# megascale ddG, one item per single substitution
'''
Tsuboyama et al. 2023, curated by ThermoMPNN, which is also where the sign
convention comes from. It ships train/val/test already split by wild type
protein, so a protein and its mutants never span two splits. Only single
substitutions are kept: insertions, deletions and multiples are dropped,
because the item a scorer is asked about is one wild type sequence plus one
mutation. The mutant sequences are part of the benchmark too even though no
item carries them, so sequences() adds them on top, which is what keeps the
unlabeled dataset clear of them as well.
'''

import csv

from . import Benchmark
from .. import config as c

class Megascale(Benchmark):
    '''ddG for single substitutions, split by wild type protein.'''

    name = "megascale"
    # every item is about one wild type protein
    key = "pdb"
    # what this benchmark is scored on, for the stratified split
    label = "ddg"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        '''
        one item per single substitution, carrying the wild type sequence, the
        mutation and its ddg.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, mutation, ddg, pdb and split
        '''
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

    def sequences(self):
        '''
        the mutant sequences are part of the benchmark too, so hold those out
        as well. megascale spells each one out, so nothing has to be applied.

        inputs:
            none
        outputs:
            seqs: iterator[str] -- wild type and mutant sequences
        '''
        yield from super().sequences()
        for split in ("train", "val", "test"):
            path = c.data / "megascale" / f"{split}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    s = (r.get("aa_seq") or "").strip().upper()
                    if s: yield s
