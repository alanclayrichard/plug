# SHS27k: the string human interaction subset used by GNN-PPI (lv et al. 2021), built
# with the recipe from PIPR (chen et al. 2019). 1,690 proteins that are less than 40%
# alike, and the pairs between them. ships no split of its own.
import csv

from . import Benchmark
from .. import config as c

# the kinds of interaction string records for a pair
MODES = ("reaction", "binding", "catalysis", "activation", "inhibition", "ptmod", "expression")


class Shs27k(Benchmark):
    name = "shs27k"
    # an item is about two proteins, so both are named and both have a sequence
    key = ("id_a", "id_b")
    seq = ("sequence_a", "sequence_b")
    published = ()

    # pass mode to keep only the pairs that interact that way
    def __init__(self, mode=None, split="test", how=None):
        super().__init__(split, how)
        if mode: self.items = [it for it in self.items if mode in it["types"]]

    @classmethod
    def rows(cls):
        folder = c.data / "shs27k"
        seqs = {}
        with open(folder / "protein.SHS27k.sequences.dictionary.tsv") as f:
            for line in f:
                acc, _, s = line.partition("\t")
                if s.strip(): seqs[acc] = s.strip().upper()
        # string lists one row per interaction kind, so gather the kinds per pair first
        pairs = {}
        with open(folder / "protein.actions.SHS27k.STRING.txt", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                a, b = r["item_id_a"], r["item_id_b"]
                if a not in seqs or b not in seqs: continue
                e = pairs.setdefault(tuple(sorted((a, b))), {"types": set(), "score": 0})
                if r["mode"]: e["types"].add(r["mode"])
                # keep the best confidence string gives the pair
                e["score"] = max(e["score"], int(r["score"]) if r["score"].strip() else 0)
        for (a, b), e in pairs.items():
            yield {"sequence_a": seqs[a], "sequence_b": seqs[b], "id_a": a, "id_b": b,
                   "types": sorted(e["types"]), "score": e["score"]}
