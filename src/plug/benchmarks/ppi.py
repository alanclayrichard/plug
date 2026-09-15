# human protein-protein interaction gold standard, one item per pair
'''
Bernett et al. 2023. Each item is a pair of proteins labelled 1 if they
interact and 0 if they don't, balanced 50/50 on every split. The paper divides
its proteins into three blocks and keeps only the pairs inside a block, so
neither a protein nor a close relative spans two splits — which is the same
property a homology split here builds, arrived at differently. An item is about
two proteins, so both are named and both carry a sequence, and a pair survives
a split only when both of its ends agree; a pair is skipped entirely unless
both proteins have a sequence in the fasta.
'''

from . import Benchmark
from .. import config as c

# the paper's blocks, in the order it trains on them
BLOCKS = {1: "train", 0: "val", 2: "test"}

class Ppi(Benchmark):
    '''interacting and non-interacting human protein pairs.'''

    name = "ppi"
    # an item is about two proteins, so both are named and both have a sequence
    key = ("id_a", "id_b")
    # what this benchmark is scored on, for the stratified split
    label = "label"
    seq = ("sequence_a", "sequence_b")
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        '''
        one item per pair, positive and negative, over the paper's three
        blocks.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- both sequences, both ids, label and split
        '''
        folder = c.data / "ppi"
        seqs = dict(c.iter_fasta(folder / "human_swissprot_oneliner.fasta"))
        for block, split in BLOCKS.items():
            for label, sign in ((1, "pos"), (0, "neg")):
                path = folder / f"Intra{block}_{sign}_rr.txt"
                if not path.exists(): continue
                for line in open(path):
                    p = line.split()
                    # skip a pair unless both of its proteins have a sequence
                    if len(p) >= 2 and p[0] in seqs and p[1] in seqs:
                        yield {"sequence_a": seqs[p[0]], "sequence_b": seqs[p[1]],
                               "id_a": p[0], "id_b": p[1],
                               "label": label, "split": split}
