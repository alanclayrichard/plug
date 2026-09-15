# build a benchmark's leakage-free split and prove nothing crossed
'''
This runs the method end to end: build the similarity graph, hand its groups to
the solver, and turn the answer back into a split you can use. Run it as
"python -m plug.splitting" for every benchmark, or with names for a few. The
solver is asked for 80/10/10 of whatever survives, which is the only form every
benchmark can satisfy: on ppi a pair is kept only when both ends agree, which
caps retention near 51%, so no assignment could put 80% of all pairs in train
however it was arranged. A datapoint is kept only when everything it is about
agrees on a split, so both partners of an interaction, or the protein and the
ligand of an affinity, have to land together.

Two files come out of it. data/splits/<benchmark>__<thresholds>.csv is one
entity per row, proteins and ligands side by side, because entities are what
the split is decided over, and it is what Benchmark(how="homology") reads.
Passing "--datapoints" also writes the one to hand out: one surviving datapoint
per row, indexed along Benchmark.rows(), so nobody has to join proteins back
onto ligands and rediscover which pairs agreed.

The audit at the end re-reads the raw mmseqs, foldseek and rdkit hits rather
than the groups — checking the groups would be circular, since they are what
the assignment was built from — and counts how many ended up with their two
ends in different splits. It has to be zero.
'''

import csv
import gzip
import json
import random
import sys
import time
from collections import Counter
from pathlib import Path

from . import config as c
from . import graphs
from .graphs import SPLITS
from .solver import FRAC, solve_scip

def stamp(cls, split_of, rows):
    '''
    say which split each datapoint landed in. a datapoint survives only when
    everything it is about agrees, which is one rule for every modality: both
    partners of an interaction, or the protein and the ligand of an affinity.
    an entity nothing placed is in no split, so anything about it is dropped
    for the same reason a straddler is.

    inputs:
        cls: type -- a Benchmark subclass
        split_of: dict[str, str] -- entity -> split
        rows: list[dict] -- the benchmark's items
    outputs:
        where: list[str | None] -- one split per row, None if it straddles
    '''
    out = []
    for it in rows:
        names = {split_of.get(it[k]) for k in cls.keyfields()}
        if "ligand" in cls.views and it.get(cls.ligand):
            names.add(split_of.get(it[cls.ligand]))
        where = names.pop() if len(names) == 1 else None
        out.append(where if where in SPLITS else None)
    return out

def audit(split_of, g):
    '''
    the leakage check, run against the raw search hits rather than the groups.
    checking the groups would be circular — they are what the assignment was
    built from — so these are the actual mmseqs, foldseek and rdkit pairs, and
    a mistake anywhere between the search and the assignment shows up here as a
    non-zero count.

    inputs:
        split_of: dict[str, str] -- entity -> split
        g: dict -- the graph the assignment was built from
    outputs:
        counts: dict -- per modality, how many hits crossed a split, plus the
                groups spanning splits and the total leaked edges
    '''
    out = {}
    for m, edges in g["edges"].items():
        cross, unplaced = Counter(), 0
        for a, b in edges:
            x, y = split_of.get(a), split_of.get(b)
            # an entity nothing placed is in no split, so every datapoint about it is
            # dropped and no hit touching it can leak
            if x not in SPLITS or y not in SPLITS: unplaced += 1; continue
            if x != y: cross["|".join(sorted((x, y)))] += 1
        n = sum(cross.values())
        live = len(edges) - unplaced
        out[f"{m}_edges"] = len(edges)
        out[f"{m}_edges_with_an_unplaced_end"] = unplaced
        out[f"{m}_cross_split"] = n
        out[f"{m}_cross_split_pct"] = round(100 * n / live, 4) if live else 0.0
        out[f"{m}_cross_split_by_pair"] = dict(cross)
    out["components_spanning_splits"] = sum(
        1 for gs in g["groups"].values() for gr in gs
        if len({split_of.get(k) for k in gr} & set(SPLITS)) > 1)
    out["leaked_edges"] = sum(out[f"{m}_cross_split"] for m in g["edges"])
    return out

def report(cls, split_of, g):
    '''
    retention, sizes and the audit for one assignment.

    inputs:
        cls: type -- a Benchmark subclass
        split_of: dict[str, str] -- entity -> split
        g: dict -- the graph the assignment was built from
    outputs:
        report: dict -- how many datapoints were kept and dropped, the split
                sizes both ways, and every audit count
    '''
    where = stamp(cls, split_of, g["rows"])
    n = len(where)
    counts = Counter(w for w in where if w)
    got = sum(counts.values())
    return {"datapoints": n, "retained": got, "dropped": n - got,
            "retention_pct": round(100 * got / n, 2),
            "counts": {s: counts.get(s, 0) for s in SPLITS},
            "fraction_of_original": {s: round(counts.get(s, 0) / n, 4) for s in SPLITS},
            "fraction_of_retained": {s: round(counts.get(s, 0) / max(got, 1), 4)
                                     for s in SPLITS},
            **audit(split_of, g)}

def greedy(cls, g, frac=FRAC):
    '''
    a legal assignment found greedily, to start the solver from. biggest group
    first, each going to whichever split is furthest below its share, and
    ligand groups then follow the proteins they were measured against. this is
    only ever a hint — it cannot change what counts as optimal, it just gets
    there sooner.

    inputs:
        cls: type -- a Benchmark subclass
        g: dict -- a graph from graphs.build()
        frac: tuple[float, float, float] -- the target train/val/test shares
    outputs:
        split_of: dict[str, str] -- entity -> split
    '''
    if len(cls.keyfields()) > 1:
        # a pair survives only when both ends agree, so a split holding x of the proteins
        # keeps about x*x of the pairs; dealing at the square root compensates
        r = [f ** 0.5 for f in frac]
        frac = tuple(x / sum(r) for x in r)

    groups = [list(x) for x in g["groups"]["protein"]]
    weight = g["weight"]
    # shuffle first so groups of the same size don't always fall the same way round
    random.Random(c.seed).shuffle(groups)
    groups.sort(key=lambda gr: -sum(weight[k] for k in gr))
    want = {s: f * sum(weight.values()) for s, f in zip(SPLITS, frac)}
    have = dict.fromkeys(SPLITS, 0)
    split_of = {}
    for gr in groups:
        s = max(SPLITS, key=lambda s: want[s] - have[s])
        have[s] += sum(weight[k] for k in gr)
        split_of.update(dict.fromkeys(gr, s))

    if "ligand" in cls.views:
        votes = {}
        for k, smiles in cls.ligands().items():
            for s in smiles: votes.setdefault(s, Counter())[split_of[k]] += 1
        total = sum(sum(v.values()) for v in votes.values())
        want = {s: f * total for s, f in zip(SPLITS, frac)}
        have = dict.fromkeys(SPLITS, 0)
        for gr in sorted(g["groups"]["ligand"],
                         key=lambda gr: -sum(sum(votes[k].values()) for k in gr)):
            tally = Counter()
            for k in gr: tally.update(votes[k])
            # where most of its datapoints already are, unless that split is over its share
            pick = max(SPLITS, key=lambda s: (tally[s] - max(0, have[s] - want[s]),
                                              want[s] - have[s]))
            for k in gr: split_of[k] = pick
            have[pick] += sum(tally.values())
    return split_of

def split(cls, g=None, tolerance=0.05, time_limit=600.0, warm_start=True, frac=FRAC,
          rebuild=False, verbose=True, basis="retained", balance=True, stratify=False,
          strat_tolerance=0.75):
    '''
    one benchmark's leakage-free split, as {entity: "train"/"val"/"test"}.
    the proportions are of what survives rather than of the whole benchmark,
    which is the only form every benchmark can satisfy, and the tolerance is a
    hard constraint on top of it. a second pass then pins the retention just
    proved and lands the split sizes on target, which cannot cost a datapoint.

    inputs:
        cls: type | str -- a Benchmark subclass or its name
        g: dict | None -- a prebuilt graph, or None to build one
        tolerance: float -- how far a split may sit from its share
        time_limit: float -- seconds per solver stage
        warm_start: bool -- hand the solver a greedy assignment as a bound
        frac: tuple[float, float, float] -- the target train/val/test shares
        rebuild: bool -- redo the similarity searches
        verbose: bool -- print progress
        basis: str -- "retained" (of what survives) or "original"
        balance: bool -- run the second pass that lands the sizes on target
        stratify: bool -- give each label class its share of each split too
        strat_tolerance: float -- relative, so 0.75 keeps every class between
                        a quarter and one and three quarters of its share
    outputs:
        split_of: dict[str, str] | None -- entity -> split, None if unsolved
        info: dict -- solver status, retention and runtime
    '''
    cls = graphs._resolve(cls)
    g = g if g is not None else graphs.build(cls, rebuild, verbose)
    comp_of, n_comp = graphs.index(g)
    sigs = graphs.signatures(cls, g["rows"], comp_of)
    strata = graphs.strata(cls, g["rows"], comp_of) if stratify else None
    if verbose and strata is not None:
        print(f"    stratifying on {len(strata)} classes of {cls.label!r} "
              f"at +/-{100 * strat_tolerance:.0f}% of each share", flush=True)

    hint = None
    if warm_start:
        gm = greedy(cls, g, frac)
        hint = [0] * n_comp
        for k, cid in comp_of.items(): hint[cid] = SPLITS.index(gm[k])

    t = time.time()
    if verbose:
        print(f"    basis={basis} tolerance={tolerance} {time_limit:.0f}s/stage", flush=True)
    assign, info = solve_scip(sigs, n_comp, len(g["rows"]), frac=frac, tolerance=tolerance,
                              time_limit=time_limit, hint=hint, basis=basis,
                              balance=balance, strata=strata,
                              strat_tolerance=strat_tolerance)
    info.update(requested_tolerance=tolerance, total_runtime=round(time.time() - t, 2))
    if assign is None:
        return None, info
    return {k: SPLITS[assign[cid]] for k, cid in comp_of.items()}, info

def split_file(cls):
    '''
    where a benchmark's split is kept. the thresholds are in the file name, so
    a split made at 20% identity is never reused after you change that setting.

    inputs:
        cls: type | str -- a Benchmark subclass or its name
    outputs:
        path: Path -- data/splits/<benchmark>__<thresholds>.csv
    '''
    return c.splits / f"{graphs._resolve(cls).name}__{graphs.tag(cls)}.csv"

def read(path):
    '''
    read a saved split back off disk.

    inputs:
        path: Path -- a split csv
    outputs:
        split_of: dict[str, str] -- entity -> split
    '''
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        # splits written before ligands were split too name the first column "protein"
        key = "entity" if "entity" in (r.fieldnames or ()) else "protein"
        return {row[key]: row["split"] for row in r}

def write(path, split_of):
    '''
    write a split to disk, sorted so the file is stable across runs.

    inputs:
        path: Path -- where to write it
        split_of: dict[str, str] -- entity -> split
    outputs:
        none
    '''
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("entity", "split"))
        w.writerows(sorted(split_of.items()))

def datapoint_file(cls):
    '''
    where this benchmark's per-datapoint split is kept, next to the entity one
    it was derived from.

    inputs:
        cls: type | str -- a Benchmark subclass or its name
    outputs:
        path: Path -- data/splits/<benchmark>__<thresholds>_datapoints.csv
    '''
    return c.splits / f"{graphs._resolve(cls).name}__{graphs.tag(cls)}_datapoints.csv"

def datapoint_columns(cls):
    '''
    the columns a per-datapoint file carries: the protein sequence, or both of
    them for a benchmark whose items are about a pair, then whatever names the
    structure for one that ships structures, then the ligand for one that ships
    ligands. these are exactly the things a datapoint rides on, so they are also
    exactly what had to agree for it to survive.

    inputs:
        cls: type -- a Benchmark subclass
    outputs:
        cols: list[str] -- the field names, in file order
    '''
    cols = list(cls.seqfields())
    if "structure" in cls.views:
        cols += [k for k in cls.keyfields() if k not in cols]
    if "ligand" in cls.views: cols.append(cls.ligand)
    return cols

def write_datapoints(path, cls, rows, where, compress=False):
    '''
    write the split carrying the sequence, structure and ligand a row is about
    rather than an index into something else, so a row can be read on its own.
    rows that straddle are left out, and identical rows are written once: the
    columns are exactly what the split is decided over, so a repeat carries no
    split the first one did not. that matters — proteingym names one wild type
    per assay and repeats it for every variant measured against it, which is
    2,465,767 rows over 187 sequences.

    inputs:
        path: Path -- where to write it, ".gz" appended when compressing
        cls: type -- a Benchmark subclass
        rows: list[dict] -- the benchmark's items, in rows() order
        where: list[str | None] -- one split per row, from stamp()
        compress: bool -- gzip it as well as deduplicating
    outputs:
        path: Path -- where it actually went
        n: int -- how many distinct datapoints were written
    '''
    cols = datapoint_columns(cls)
    path = Path(f"{path}.gz") if compress else Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    op = gzip.open if compress else open
    seen = set()
    with op(path, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow((*cols, "split"))
        for it, s in zip(rows, where):
            if not s: continue
            row = tuple(it[c] for c in cols)
            if row in seen: continue
            seen.add(row)
            w.writerow((*row, s))
    return path, len(seen)

def read_datapoints(path):
    '''
    read a per-datapoint file back. every row carries what it is about, so
    nothing has to be joined against the benchmark to make sense of it.

    inputs:
        path: str | Path -- a datapoints csv, gzipped or not
    outputs:
        rows: list[dict] -- one dict per surviving datapoint, "split" included
    '''
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        return list(csv.DictReader(f))

def cached(cls, **kw):
    '''
    the split for one benchmark, built and saved the first time it is asked
    for and read back from then on. searching a benchmark against itself and
    solving takes minutes, so it should only happen once per setting.

    inputs:
        cls: type | str -- a Benchmark subclass or its name
        kw: dict -- anything split() takes
    outputs:
        split_of: dict[str, str] -- entity -> split
    '''
    cls = graphs._resolve(cls)
    path = split_file(cls)
    if path.exists(): return read(path)
    split_of, info = split(cls, **kw)
    if split_of is None:
        raise RuntimeError(f"{cls.name}: no split found — {info.get('status')}")
    write(path, split_of)
    return split_of

def main(names=None, datapoints=False, compress=False, **kw):
    '''
    split the named benchmarks, or all of them, printing the problem size
    before and the retention and leaked-edge count after.

    inputs:
        names: list[str] | None -- benchmarks to do, None for all
        datapoints: bool -- also write the per-datapoint file, which is the one
                    to hand out: one surviving row per line rather than one
                    entity, so nothing has to be joined back by hand
        compress: bool -- gzip that file; a sequence repeats once per datapoint
                  about it, so the saving is large
        kw: dict -- anything split() takes
    outputs:
        none
    '''
    from .benchmarks import REGISTRY
    for name in names or list(REGISTRY):
        cls = REGISTRY[name]
        print(f"\n=== {name}", flush=True)
        g = graphs.build(cls)
        print(json.dumps({k: v for k, v in graphs.stats(cls, g).items()
                          if k in ("datapoints", "entities", "edges", "components",
                                   "largest_component_pct", "anything_to_solve")}, indent=1))
        split_of, info = split(cls, g, **kw)
        if split_of is None:
            print(f"  no split: {info.get('status')}", flush=True)
            continue
        m = report(cls, split_of, g)
        write(split_file(cls), split_of)
        if datapoints:
            out, n = write_datapoints(datapoint_file(cls), cls, g["rows"],
                                      stamp(cls, split_of, g["rows"]), compress)
            print(f"  {n:,} datapoints -> {out}", flush=True)
        print(f"  retention {m['retention_pct']:.2f}%  "
              f"{m['counts']['train']:,}/{m['counts']['val']:,}/{m['counts']['test']:,}"
              f"  leaked edges {m['leaked_edges']:,}"
              f"  ({info.get('total_runtime', 0):.0f}s"
              f", {info.get('peak_memory_mb', 0):,.0f} MB)  -> {split_file(cls)}", flush=True)

if __name__ == "__main__":
    # bare benchmark names, plus --key value for anything split() takes
    args, kw, it = [], {}, iter(sys.argv[1:])
    for a in it:
        if a in ("--datapoints", "--compress"):
            kw[a[2:]] = True
        elif a.startswith("--"):
            v = next(it)
            kw[a[2:].replace("-", "_")] = (None if v == "none" else
                                           float(v) if v.replace(".", "").isdigit() else v)
        else: args.append(a)
    main(args or None, **kw)
