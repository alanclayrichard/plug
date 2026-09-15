# fireprot ddG, one item per substitution on the uniprot sequence
'''
FireProtDB, curated by ThermoMPNN. It ships train/val/test split by sequence
and its positions are numbered against the uniprot sequence rather than the pdb
one. The benchmark also carries the pdb sequence of each protein and the mutant
of every measured substitution, so sequences() holds those out too. FireProt
spells no mutant out, so each one is applied here — and only where the residue
it names really is at that position, which is the same check the scoring makes
before it trusts a row.
'''

import csv

from . import Benchmark
from .. import config as c

class Fireprot(Benchmark):
    '''ddG for substitutions, numbered against uniprot and split by sequence.'''

    name = "fireprot"
    # every item is about one wild type protein
    key = "pdb"
    # what this benchmark is scored on, for the stratified split
    label = "ddg"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        '''
        one item per substitution with a measured ddg.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, mutation, ddg, pdb and split
        '''
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

    def sequences(self):
        '''
        the pdb sequence of each protein and the mutant of every measured
        substitution are in the benchmark too. the mutants aren't spelled out,
        so each is applied here, and only where the residue it names really is
        at that position.

        inputs:
            none
        outputs:
            seqs: iterator[str] -- wild type, pdb and mutant sequences
        '''
        yield from super().sequences()
        for split in ("train", "val", "test"):
            path = c.data / "fireprot" / f"{split}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    s = (r.get("pdb_sequence") or "").strip().upper()
                    if s: yield s
                    seq = (r.get("sequence") or "").strip().upper()
                    wt, new = (r.get("wild_type") or "").strip(), (r.get("mutation") or "").strip()
                    try: i = int(float(r.get("position") or ""))
                    except ValueError: continue
                    if seq and len(new) == 1 and 1 <= i <= len(seq) and seq[i - 1] == wt:
                        yield seq[:i - 1] + new + seq[i:]
