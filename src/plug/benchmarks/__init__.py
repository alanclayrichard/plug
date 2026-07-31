# one benchmark = a list of items. each benchmark file only says how to read its own
# rows; this class does the rest: which split to use, and which one to hand back.
from torch.utils.data import Dataset
import csv
from collections import Counter

from .. import config as c
from .. import homology

class Benchmark(Dataset):
    # what this benchmark is called, also the name of its saved split file
    name = ""
    # the field that names the protein an item is about
    # (two fields, as a tuple, when an item is about a pair of proteins)
    key = "sequence"
    # the field holding that protein's sequence (also two, for pairs)
    seq = "sequence"
    # the splits this benchmark came with, if it came with any
    published = ()
    # worked out the first time a homology split is asked for, then kept here
    _split = None

    # each benchmark writes this one and nothing else: yield one dict per item.
    # if the benchmark came with splits, give every item the split it belongs to.
    @classmethod
    def rows(cls):
        raise NotImplementedError

    # every sequence in this benchmark, for the searches that keep the training set clean.
    # a benchmark holding sequences its items don't carry adds them on top of these.
    def sequences(self):
        for it in self.items:
            for f in self.seqfields(): yield it[f]

    # split is the split you want: "train", "val", "test", or "all".
    # how is where that split comes from:
    #   published — the split the benchmark came with (if any)
    #   homology  — a split worked out here, so no two splits share similarity
    #   whole     — no split at all, the whole benchmark is one test set (for self/unsurpervised)
    def __init__(self, split="test", how=None):
        # if not specified, use the benchmark's own split when it has one or the whole thing
        how = how or ("published" if self.published else "whole")
        # read every item first — the split has to be decided over all of them
        items = list(self.rows())

        if how == "published":
            # if published requested but there is no published throw error
            if not self.published:
                raise ValueError(f"{self.name} ships no published split — "
                                 "use how='homology' to create a split or how='whole'"
                                 " to use data as one big test set")
            # if published requested but published doesnt have split throw error
            if split not in self.published + ("all",):
                raise ValueError(f"{self.name} publishes {'/'.join(self.published)}"
                                 f", {split!r} is not allowed")

        # label every item from the homology split, working it out if it isn't cached yet
        elif how == "homology":
            items = self.stamp(items, self.split_map())

        # no splits, so everything counts as test
        elif how == "whole":
            for it in items: it["split"] = "test"

        else:
            raise ValueError(f"how: 'published', 'homology' or 'whole', not {how}")

        # keep only the splits that were asked for
        self.items = items if split == "all" else [it for it in items if it["split"] == split]

    # the number of items in this benchmark
    def __len__(self): return len(self.items)

    # get an item by index
    def __getitem__(self, i): return self.items[i]

    # key and seq are one field for most benchmarks and two for pairs. these always
    # give back a tuple, so nothing further down has to check which it is.
    @classmethod
    def keyfields(cls): return (cls.key,) if isinstance(cls.key, str) else tuple(cls.key)

    @classmethod
    def seqfields(cls): return (cls.seq,) if isinstance(cls.seq, str) else tuple(cls.seq)

    # location of cached splots. the settings are part of the file name, so a split made
    # at 20% identity is never reused after you change that setting to 30%.
    @classmethod
    def split_file(cls):
        return c.splits / f"{cls.name}__{c.split_tag()}.csv"

    # which split each protein is in, as {protein: "train"/"val"/"test"}. worked out
    # once, cached, and read back from then on — grouping the sequences takes
    # minutes on the big benchmarks, so it should only ever happen once per setting.
    @classmethod
    def split_map(cls):
        if cls._split is None:
            path = cls.split_file()
            # read it if cached, build it if not
            cls._split = _read(path) if path.exists() else _build(cls, path)
        return cls._split

    # label each item with the split its protein is in. for a pair, both proteins have
    # to be in the same split — if one is in train and the other in test, the pair
    # belongs to neither, so it is left out.
    @classmethod
    def stamp(cls, items, split_map):
        out = []
        for it in items:
            names = {split_map[it[k]] for k in cls.keyfields()}
            if len(names) == 1:
                it["split"] = names.pop()
                out.append(it)
        return out

# work out fresh homology splits: collect every protein and its sequence, count how many items
# each protein has, group the similar ones together, then deal out those groups out 80/10/10
# so that a protein and everything like it always end up in the same split.
def _build(cls, path):
    seqs, weight = {}, Counter()
    for it in cls.rows():
        for k, s in zip(cls.keyfields(), cls.seqfields()):
            seqs.setdefault(it[k], it[s])
            weight[it[k]] += 1
    # nothing to group, which almost always means the data was never downloaded
    if not seqs:
        raise FileNotFoundError(f"{cls.name} has no items to split — "
                                "download it first, see data/download_*.sh")
    split_map = homology.deal(homology.clusters(seqs), weight)
    _write(path, split_map)
    return split_map

# read a saved split back in
def _read(path):
    with open(path, newline="") as f:
        return {r["protein"]: r["split"] for r in csv.DictReader(f)}

# save a split as a two column csv: the protein, and the split it is in
def _write(path, split_map):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("protein", "split"))
        w.writerows(sorted(split_map.items()))

from .megascale import Megascale
from .fireprot import Fireprot
from .dms import Dms
from .go import Go
from .allobench import Allobench
from .ppi import Ppi
from .shs27k import Shs27k
from .passerrank import Passerrank
from .lp_pdbbind import LpPdbbind
from .bindingdb import Bindingdb
from .atlas import Atlas


# every benchmark, looked up by name — this is what fastas.py and build_trainset loop
# over. a benchmark only turns up here if it was imported above, so keep that list current.
REGISTRY = {b.name: b for b in Benchmark.__subclasses__()}
# what `from plug.benchmarks import *` hands you: the benchmark classes, nothing else
__all__ = ["Benchmark", "REGISTRY", *(b.__name__ for b in REGISTRY.values())]