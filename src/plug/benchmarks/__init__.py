# what every benchmark shares: items, splits and leakage views
'''
Each benchmark file says only how to read its own rows and this class does the
rest, so adding a twelfth means writing one rows() and nothing else. A
benchmark declares what names the protein an item is about, where that
protein's sequence lives, and — for a benchmark whose items are about a pair or
about a protein and a ligand — both of those as tuples, after which nothing
downstream has to know which kind it is. It also declares the other ways it can
leak beyond sequence, "structure" or "ligand", which turns on that search and
puts its threshold in the split file name. Asking for a benchmark takes two
arguments: split is which split you want and how says where that split came
from, either the one the benchmark was published with, a homology split worked
out here so no two splits share similarity, or none at all, which makes the
whole benchmark one test set for a self-supervised model. A homology split is
worked out once, cached on disk under the thresholds that made it, and read
back from then on. A benchmark only appears in REGISTRY if it is imported at
the bottom of this file, and REGISTRY is what the dataset build screens
against, so an unimported benchmark is one nothing is ever held clear of.
'''

from torch.utils.data import Dataset

from .. import splitting

class Benchmark(Dataset):
    '''
    the shared half of every benchmark: which split to use, and which one to
    hand back. subclasses set the fields below and write rows().
    '''

    # what this benchmark is called, also the name of its saved split file
    name = ""
    # the field that names the protein an item is about
    # (two fields, as a tuple, when an item is about a pair of proteins)
    key = "sequence"
    # the field holding that protein's sequence (also two, for pairs)
    seq = "sequence"
    # the splits this benchmark came with, if it came with any
    published = ()
    # the other ways this benchmark can leak, on top of sequence: "structure", "ligand".
    # naming one here turns on that search and puts its threshold in the split file name.
    views = ()
    # for a benchmark whose items pair a protein with a ligand, the field holding it
    ligand = "smiles"
    # the field holding what this benchmark is SCORED on, for the stratified split. None
    # leaves the split unstratified. see strata() below for how a field becomes classes.
    label = None
    # how many bins a numeric label is cut into before it can be stratified on
    bins = 4
    # worked out the first time a homology split is asked for, then kept here
    _split = None

    @classmethod
    def rows(cls):
        '''
        each benchmark writes this one and nothing else: yield one dict per
        item. if the benchmark came with splits, give every item the split it
        belongs to.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- one item per datapoint
        '''
        raise NotImplementedError

    @classmethod
    def structures(cls):
        '''
        the structures this benchmark ships, if it ships any. a benchmark
        naming "structure" in views overrides this and its homology split then
        covers shape too.

        inputs:
            none
        outputs:
            paths: dict[str, Path] -- name -> a structure file
        '''
        return {}

    @classmethod
    def ligands(cls):
        '''
        the molecules this benchmark ships, if it ships any. a benchmark naming
        "ligand" in views overrides this and its homology split then covers
        chemistry too.

        inputs:
            none
        outputs:
            ligands: dict[str, list[str]] -- name -> its smiles
        '''
        return {}

    @classmethod
    def strata(cls, rows):
        '''
        one iterable of class names per row, in the same order. the generic
        rule reads `label`: a number is cut into `bins` quantile bins over the
        whole benchmark, a list or set puts the row in one stratum per entry,
        and anything else is a class already. a benchmark whose label is
        per-residue overrides this and summarises it first.

        inputs:
            rows: list[dict] -- the benchmark's items
        outputs:
            strata: list[tuple[str, ...]] -- the classes each row belongs to
        '''
        if not cls.label: return [() for _ in rows]
        vals = [it.get(cls.label) for it in rows]
        if all(isinstance(v, (list, set, tuple)) for v in vals):
            return [tuple(f"{cls.label}={x}" for x in sorted(set(v))) for v in vals]
        try:
            nums = [float(v) for v in vals]
        except (TypeError, ValueError):
            return [(f"{cls.label}={v}",) for v in vals]
        import numpy as np
        edges = np.quantile(nums, [i / cls.bins for i in range(1, cls.bins)])
        return [(f"{cls.label}~q{int(np.searchsorted(edges, v))}",) for v in nums]

    def sequences(self):
        '''
        every sequence in this benchmark, for the searches that keep the
        training set clean. a benchmark holding sequences its items don't carry
        — a mutant, a pdb sequence, a fasta entry with no row — adds them on
        top of these.

        inputs:
            none
        outputs:
            seqs: iterator[str] -- every sequence the benchmark holds
        '''
        for it in self.items:
            for f in self.seqfields(): yield it[f]

    def __init__(self, split="test", how=None):
        '''
        load the benchmark and keep the split asked for. how is "published" for
        the split the benchmark came with, "homology" for one worked out here
        so no two splits share similarity, or "whole" for no split at all,
        which makes everything one test set.

        inputs:
            split: str -- "train", "val", "test" or "all"
            how: str | None -- "published", "homology" or "whole"; None picks
                 published where there is one, else whole
        outputs:
            none
        '''
        how = how or ("published" if self.published else "whole")
        # read every item first — the split has to be decided over all of them
        items = list(self.rows())

        if how == "published":
            if not self.published:
                raise ValueError(f"{self.name} ships no published split — "
                                 "use how='homology' to create a split or how='whole'"
                                 " to use data as one big test set")
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

        self.items = items if split == "all" else [it for it in items if it["split"] == split]

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]

    @classmethod
    def keyfields(cls):
        '''
        key is one field for most benchmarks and two for pairs. this always
        gives back a tuple, so nothing further down has to check which it is.

        inputs:
            none
        outputs:
            fields: tuple[str, ...] -- the fields naming the proteins
        '''
        return (cls.key,) if isinstance(cls.key, str) else tuple(cls.key)

    @classmethod
    def seqfields(cls):
        '''
        the same for seq: always a tuple, however many sequences an item has.

        inputs:
            none
        outputs:
            fields: tuple[str, ...] -- the fields holding the sequences
        '''
        return (cls.seq,) if isinstance(cls.seq, str) else tuple(cls.seq)

    @classmethod
    def split_file(cls):
        '''
        where this benchmark's split is cached. the thresholds are part of the
        file name, so a split made at 20% identity is never reused after you
        change that setting to 30%.

        inputs:
            none
        outputs:
            path: Path -- the split csv
        '''
        return splitting.split_file(cls)

    @classmethod
    def split_map(cls):
        '''
        which split each entity is in. worked out once, cached on disk, and
        read back from then on — searching a benchmark against itself and
        solving takes minutes, so it should only ever happen once per setting.

        inputs:
            none
        outputs:
            split_of: dict[str, str] -- entity -> "train"/"val"/"test"
        '''
        if cls._split is None: cls._split = splitting.cached(cls)
        return cls._split

    @classmethod
    def stamp(cls, items, split_map):
        '''
        label each item with the split it belongs to, and drop the ones that
        straddle: a pair whose two proteins landed in different splits belongs
        to neither, and neither does an affinity whose protein and ligand
        disagree.

        inputs:
            items: list[dict] -- the benchmark's items
            split_map: dict[str, str] -- entity -> split
        outputs:
            items: list[dict] -- the surviving items, each carrying its split
        '''
        out = []
        for it, where in zip(items, splitting.stamp(cls, split_map, items)):
            if where:
                it["split"] = where
                out.append(it)
        return out

from .megascale import Megascale
from .fireprot import Fireprot
from .dms import Dms
from .cafa5 import Cafa5
from .allobench import Allobench
from .ppi import Ppi
from .shs27k import Shs27k
from .passerrank import Passerrank
from .lp_pdbbind import LpPdbbind
from .bindingdb import Bindingdb
from .atlas import Atlas

# every benchmark, looked up by name — this is what the dataset build loops over. a
# benchmark only turns up here if it was imported above, so keep that list current.
REGISTRY = {b.name: b for b in Benchmark.__subclasses__()}
# what `from plug.benchmarks import *` hands you: the benchmark classes, nothing else
__all__ = ["Benchmark", "REGISTRY", *(b.__name__ for b in REGISTRY.values())]
