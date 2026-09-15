# proteingym deep mutational scanning, one item per measured variant
'''
ProteinGym's substitution benchmark, one item per variant and grouped by assay.
Every assay targets one protein, so the assay id is what names the protein here
rather than the sequence — two assays on the same protein still measure
different things. It ships no split of its own. The measured variants are
sequences the benchmark carries too and ProteinGym spells every one of them
out, so sequences() yields them straight from the per-assay files with nothing
to apply. Pass assay to keep one kind of measurement or dms_id to keep a single
assay.
'''

import csv
import sys

from . import Benchmark
from .. import config as c

# target sequences are long enough to trip the default csv field limit
csv.field_size_limit(sys.maxsize)

class Dms(Benchmark):
    '''deep mutational scanning scores, grouped by assay.'''

    name = "dms"
    # every item is about the protein its assay targets
    key = "dms_id"
    # what this benchmark is scored on, for the stratified split
    label = "assay"
    published = ()

    def __init__(self, assay=None, dms_id=None, split="test", how=None):
        '''
        load the benchmark, optionally narrowed to one kind of measurement or
        one assay.

        inputs:
            assay: str | None -- Activity, Binding, Expression,
                   OrganismalFitness or Stability
            dms_id: str | None -- a single assay
            split: str -- "train", "val", "test" or "all"
            how: str | None -- "published", "homology" or "whole"
        outputs:
            none
        '''
        super().__init__(split, how)
        if assay: self.items = [it for it in self.items if it["assay"] == assay]
        if dms_id: self.items = [it for it in self.items if it["dms_id"] == dms_id]

    @classmethod
    def rows(cls):
        '''
        one item per measured variant, with the assay's target sequence and
        what that assay measured.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, mutation, score, assay, dms_id
                  and uniprot
        '''
        folder = c.data / "proteingym"
        # one row per assay: its target sequence, what it measured, and the protein
        with open(folder / "DMS_substitutions.csv", newline="") as f:
            meta = {r["DMS_id"]: (r["target_seq"].strip().upper(),
                                  r["coarse_selection_type"], r["UniProt_ID"])
                    for r in csv.DictReader(f)}
        for dms_id, (seq, assay, uniprot) in meta.items():
            path = folder / "DMS_ProteinGym_substitutions" / f"{dms_id}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    if not (r.get("DMS_score") or "").strip(): continue
                    yield {"sequence": seq, "mutation": r["mutant"],
                           "score": float(r["DMS_score"]), "assay": assay,
                           "dms_id": dms_id, "uniprot": uniprot}

    def sequences(self):
        '''
        the measured variants are sequences this benchmark carries too, so hold
        those out as well. proteingym ships every one already spelled out, so
        nothing has to be applied.

        inputs:
            none
        outputs:
            seqs: iterator[str] -- target and mutated sequences
        '''
        yield from super().sequences()
        folder = c.data / "proteingym"
        with open(folder / "DMS_substitutions.csv", newline="") as f:
            ids = [r["DMS_id"] for r in csv.DictReader(f)]
        for dms_id in ids:
            path = folder / "DMS_ProteinGym_substitutions" / f"{dms_id}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    s = (r.get("mutated_sequence") or "").strip().upper()
                    if s: yield s
