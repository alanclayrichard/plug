# group a benchmark's entries that are too alike to separate
'''
This is the first half of the method and it only decides what is linked to
what; splitting.py decides where the groups go. Entries can look alike in more
than one way and the rule is the same every time: search the benchmark against
itself, link anything over the bar in either direction, and keep whole groups
of linked things together. Proteins are compared on sequence with mmseqs in
both coverage modes, and also on shape with foldseek for a benchmark that ships
structures; both describe the same proteins, so they share one graph and either
kind of match counts. Ligands are compared on ecfp4 tanimoto and get a graph of
their own, because putting them in with the proteins lets one shared molecule
chain protein to protein until the whole benchmark is a single group. Linking
carries over, so something alike something alike lands in the same group, and
it is the groups rather than the individual hits that the split treats as
indivisible. The finished graph is cached under the work directory under a name
carrying every threshold that shaped it, so a graph built at one identity is
never reused after that setting changes.
'''

import pickle
import re
import subprocess
import tempfile
import time
from collections import Counter
from pathlib import Path

from . import config as c

SPLITS = ("train", "val", "test")

# no e-value cutoff here. an e-value moves with how big the searched set is, so the same
# match scores far worse in a whole benchmark search than in the smaller split against
# split check anyone runs afterwards. any cutoff quietly drops borderline matches that
# then turn up as shared sequences across splits. identity and coverage don't move with
# the size of the search, so they decide on their own.
graph_evalue = "inf"

def sequence_links(seqs, cov_modes=("1", "2")):
    '''
    search the benchmark against itself with mmseqs and yield the pairs it says
    are the same family. mode 1 asks for the found protein to be covered and
    mode 2 the searched one; either way the two are related, so both get taken.
    a name can carry more than one sequence — allobench lists several chains
    under one target — and every one of them is searched, because a sequence
    that is never searched is a match that is never found.

    inputs:
        seqs: dict[str, list[str]] -- name -> the sequences it holds
        cov_modes: tuple[str, ...] -- mmseqs coverage modes to run
    outputs:
        pairs: iterator[tuple[str, str]] -- names that look alike
    '''
    pairs = [(k, s) for k, v in seqs.items() for s in v]
    c.work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=c.work) as tmp:
        tmp = Path(tmp)
        fa, m8 = tmp / "in.fasta", tmp / "hits.m8"
        # number the proteins: their names aren't always safe to put in a fasta
        with open(fa, "w") as o:
            for i, (_, seq) in enumerate(pairs): o.write(f">s{i}\n{seq}\n")
        for m in cov_modes:
            # no cap on matches per protein (--max-seqs). in a benchmark with many near
            # copies a cap fills up with those and hides the distant ones, and a hidden
            # match is a leak. mmseqs sizes the search by proteins x cap, so this only
            # stays affordable while the benchmarks are small.
            subprocess.run([str(c.mmseqs), "easy-search", str(fa), str(fa), str(m8),
                            str(tmp / f"s{m}"), "-s", c.sens, "-e", graph_evalue, "-c", c.cov,
                            "--cov-mode", m, "--max-seqs", str(len(pairs)), "--threads", c.threads,
                            "-v", "1", "--format-output", "query,target,fident"], check=True)
            for ln in open(m8):
                q, t, ident = ln.split()[:3]
                # too different to count as alike
                if float(ident) <= c.min_id: continue
                yield pairs[int(q[1:])][0], pairs[int(t[1:])][0]

def structure_links(paths):
    '''
    the same idea on shape, so two proteins that fold alike can't be split
    apart however far their sequences have drifted. the tm score is read both
    ways round, which is the shape version of the two coverage modes: qtmscore
    measures the match over the searched chain and ttmscore over the found one.
    foldseek's own alntmscore measures it over the matched stretch instead,
    which scores one domain buried in a much longer chain as a whole-chain
    match — six times too many hits.

    inputs:
        paths: dict[str, str | Path] -- name -> a structure file
    outputs:
        pairs: iterator[tuple[str, str]] -- names that fold the same way
    '''
    keys = list(paths)
    c.work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=c.work) as tmp:
        tmp = Path(tmp)
        # foldseek names its hits after the files, so hand it numbered copies
        pool = tmp / "pool"; pool.mkdir()
        for i, k in enumerate(keys):
            (pool / f"s{i}.pdb").write_bytes(Path(paths[k]).read_bytes())
        db, res, m8 = tmp / "db", tmp / "res", tmp / "hits.m8"
        run = lambda *a: subprocess.run([str(c.foldseek), *a, "-v", "1"], check=True)
        run("createdb", str(pool), str(db))
        # tmalign with backtrace, so the score is a real tm score
        run("search", str(db), str(db), str(res), str(tmp / "ss"), "-e", graph_evalue,
            "--alignment-type", "1", "-a", "--max-seqs", str(len(keys)), "--threads", c.threads,
            *(("--gpu", "1") if c.gpu else ()))
        run("convertalis", str(db), str(db), str(res), str(m8),
            "--format-output", "query,target,qtmscore,ttmscore", "--threads", c.threads)
        for ln in open(m8):
            q, t, qs, ts = ln.split("\t")[:4]
            # not the same shape either way round
            if max(float(qs), float(ts)) <= c.split_tm: continue
            # foldseek hands names back as "s12.pdb" or "s12.pdb_A"
            yield keys[int(re.match(r"s(\d+)", q).group(1))], keys[int(re.match(r"s(\d+)", t).group(1))]

def ligand_links(ligands):
    '''
    the same idea on chemistry. two entries sharing a molecule are obviously
    linked, and so are two whose molecules score over the tanimoto bar, which
    is what stops a test ligand reappearing in train as a near copy of itself.
    the tanimoto is done in blocks — the shared bits of every pair at once —
    because that is the only way it finishes on the big sets.

    inputs:
        ligands: dict[str, list[str]] -- name -> its smiles
    outputs:
        pairs: iterator[tuple[str, str]] -- names with the same chemistry
    '''
    try:
        import numpy as np
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
    except ImportError as e:
        raise ImportError("splitting on ligands needs rdkit — `uv sync --extra chem`") from e
    RDLogger.DisableLog("rdApp.*")

    # ecfp4 — morgan radius 2, 2048 bits. stereochemistry is left out, so two
    # enantiomers score 1.0 and always land in the same split: the safe way round.
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    owners, rows = [], []                       # one entry per distinct molecule
    seen, unreadable = {}, {}
    for k, smiles in ligands.items():
        for s in smiles:
            if s in seen: owners[seen[s]].append(k); continue
            if s in unreadable: unreadable[s].append(k); continue
            m = Chem.MolFromSmiles(s)
            # a few rows aren't readable chemistry. they still can't straddle a split,
            # so keep them by name — they just can't be compared to anything else.
            if m is None: unreadable[s] = [k]; continue
            seen[s] = len(rows); owners.append([k])
            rows.append(gen.GetFingerprintAsNumPy(m).astype(np.float32))

    # the same molecule used by two entries already links them
    for own in (*owners, *unreadable.values()):
        for o in own[1:]: yield own[0], o
    if not rows: return

    X = np.asarray(rows, dtype=np.float32)
    bits = X.sum(1)
    for lo in range(0, len(X), 1500):
        hi = min(lo + 1500, len(X))
        shared = X[lo:hi] @ X.T
        score = shared / (bits[lo:hi, None] + bits[None, :] - shared + 1e-9)
        # only look at each pair once
        for r in range(hi - lo): score[r, : lo + r + 1] = 0
        for a, b in zip(*np.nonzero(score >= c.tanimoto)):
            yield owners[lo + int(a)][0], owners[int(b)][0]

def merge(names, *links):
    '''
    union-find over every name and any number of link sources. transitivity is
    the whole point: something alike something alike lands in the same group
    and so cannot end up in a different split.

    inputs:
        names: list[str] -- every entity of one modality
        links: iterable[tuple[str, str]] -- one or more streams of linked pairs
    outputs:
        groups: list[list[str]] -- the groups that have to stay together
    '''
    at = {k: i for i, k in enumerate(names)}
    # everything starts in its own group, then they join up as links come in
    up = list(range(len(names)))

    # follow a name to its group, flattening the path on the way
    def find(i):
        while up[i] != i: up[i] = up[up[i]]; i = up[i]
        return i

    for pairs in links:
        for a, b in pairs:
            x, y = find(at[a]), find(at[b])
            if x != y: up[x] = y
    out = {}
    for k in names: out.setdefault(find(at[k]), []).append(k)
    return list(out.values())

def tag(cls):
    '''
    every threshold that shaped this benchmark's graph, for its file name.

    inputs:
        cls: type -- a Benchmark subclass
    outputs:
        tag: str -- e.g. "id20_cov80_tan60"
    '''
    t = f"id{int(c.min_id * 100)}_cov{int(float(c.cov) * 100)}"
    if "structure" in cls.views: t += f"_tm{int(c.split_tm * 100)}"
    if "ligand" in cls.views: t += f"_tan{int(c.tanimoto * 100)}"
    return t

def thresholds(cls):
    '''
    the same settings again, spelled out, so a saved result records what
    produced it.

    inputs:
        cls: type -- a Benchmark subclass
    outputs:
        thresholds: dict[str, float | str] -- the settings in force
    '''
    out = {"min_id": c.min_id, "coverage": float(c.cov), "cov_modes": "1 union 2",
           "sensitivity": c.sens}
    if "structure" in cls.views: out["split_tm"] = c.split_tm
    if "ligand" in cls.views: out["tanimoto"] = c.tanimoto
    return out

def _resolve(cls):
    '''
    take either a benchmark class or its name and give back the class.

    inputs:
        cls: type | str -- a Benchmark subclass or its registered name
    outputs:
        cls: type -- the Benchmark subclass
    '''
    # imported here rather than at the top: benchmarks imports splitting, which imports
    # this, so a module-level import would be a cycle
    if isinstance(cls, str):
        from .benchmarks import REGISTRY
        return REGISTRY[cls]
    return cls

def build(cls, rebuild=False, verbose=True):
    '''
    the nodes, edges and groups for one benchmark. shape joins the protein
    graph rather than starting its own, since it describes the same entities;
    ligands get their own. searching a benchmark against itself takes minutes,
    so the answer is cached under a name carrying every threshold.

    inputs:
        cls: type | str -- a Benchmark subclass or its name
        rebuild: bool -- redo the searches even if a cache is there
        verbose: bool -- print progress
    outputs:
        g: dict -- rows, weight, nodes, edges, groups and thresholds
    '''
    cls = _resolve(cls)
    cache = c.work / f"graph_{cls.name}__{tag(cls)}.pkl"
    if cache.exists() and not rebuild:
        return pickle.loads(cache.read_bytes())

    say = print if verbose else (lambda *a, **k: None)
    rows = list(cls.rows())
    seqs, weight = {}, Counter()
    for it in rows:
        for k, s in zip(cls.keyfields(), cls.seqfields()):
            v = seqs.setdefault(it[k], [])
            if it[s] not in v: v.append(it[s])
            weight[it[k]] += 1
    if not seqs:
        raise FileNotFoundError(f"{cls.name} has no items to split — "
                                "download it first, see data/download_*.sh")

    say(f"  sequence search over {sum(len(v) for v in seqs.values()):,} sequences "
        f"from {len(seqs):,} proteins ...", flush=True)
    t = time.time()
    edges = {"protein": list(sequence_links(seqs))}
    say(f"    {len(edges['protein']):,} hits in {time.time() - t:.0f}s", flush=True)

    if "structure" in cls.views:
        # no fallback if foldseek is missing — a graph without the structure half is a
        # different graph and must not be cached under this name
        t = time.time()
        edges["protein"] += list(structure_links(cls.structures()))
        say(f"    + structure, {len(edges['protein']):,} total in {time.time()-t:.0f}s",
            flush=True)

    nodes = {"protein": list(seqs)}
    if "ligand" in cls.views:
        nodes["ligand"] = list(dict.fromkeys(s for v in cls.ligands().values() for s in v))
        say(f"  tanimoto over {len(nodes['ligand']):,} molecules ...", flush=True)
        t = time.time()
        edges["ligand"] = list(ligand_links({s: [s] for s in nodes["ligand"]}))
        say(f"    {len(edges['ligand']):,} hits in {time.time() - t:.0f}s", flush=True)

    groups = {}
    for m in nodes:
        # one undirected pair per row: a hit is reported under each coverage mode, and
        # again in the other direction when the match holds both ways round
        edges[m] = sorted({tuple(sorted(e)) for e in edges[m]})
        groups[m] = merge(list(nodes[m]), edges[m])

    g = dict(benchmark=cls.name, rows=rows, weight=weight, nodes=nodes, edges=edges,
             groups=groups, thresholds=thresholds(cls))
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_bytes(pickle.dumps(g))
    return g

def index(g):
    '''
    number every group of every modality once, so the solver can talk about
    them as integers.

    inputs:
        g: dict -- a graph from build()
    outputs:
        comp_of: dict[str, int] -- entity -> group number
        n: int -- how many groups there are
    '''
    comp_of, at = {}, 0
    for m in g["groups"]:
        for gr in g["groups"][m]:
            for k in gr: comp_of[k] = at
            at += 1
    return comp_of, at

def signatures(cls, rows, comp_of):
    '''
    each datapoint as the set of groups it rides on, and how many datapoints
    share that set. one group for a single-protein measurement, two for an
    interaction pair, a protein and a ligand for an affinity — identical from
    here on, which is what lets one optimiser handle every modality with no
    special case per benchmark.

    inputs:
        cls: type -- a Benchmark subclass
        rows: list[dict] -- the benchmark's items
        comp_of: dict[str, int] -- entity -> group number
    outputs:
        sigs: Counter[frozenset[int], int] -- group set -> datapoints on it
    '''
    out = Counter()
    for it in rows:
        sig = {comp_of[it[k]] for k in cls.keyfields()}
        if "ligand" in cls.views and it.get(cls.ligand): sig.add(comp_of[it[cls.ligand]])
        out[frozenset(sig)] += 1
    return out

def strata(cls, rows, comp_of):
    '''
    the same signatures, counted one label class at a time. a split can be
    leakage-free, the right size, and still ask a different question of the
    test set than it trained on — proteingym's test set holds four of five
    assay types — and the fix is the size constraint written once more per
    class over these same signatures.

    inputs:
        cls: type -- a Benchmark subclass
        rows: list[dict] -- the benchmark's items
        comp_of: dict[str, int] -- entity -> group number
    outputs:
        strata: dict[str, Counter[frozenset[int], int]] -- class -> its signatures
    '''
    out = {}
    for it, names in zip(rows, cls.strata(rows)):
        sig = {comp_of[it[k]] for k in cls.keyfields()}
        if "ligand" in cls.views and it.get(cls.ligand): sig.add(comp_of[it[cls.ligand]])
        sig = frozenset(sig)
        for n in names: out.setdefault(n, Counter())[sig] += 1
    return out

def stats(cls, g):
    '''
    how big the problem is, before anyone tries to solve it. the largest
    component is the number to read first: a benchmark whose biggest group
    carries most of its datapoints has no legal 80/10/10 at all.

    inputs:
        cls: type | str -- a Benchmark subclass or its name
        g: dict -- a graph from build()
    outputs:
        stats: dict -- counts, component sizes and the thresholds used
    '''
    cls = _resolve(cls)
    comp_of, n_comp = index(g)
    sigs = signatures(cls, g["rows"], comp_of)
    multi = [k for k in sigs if len(k) > 1]
    n = len(g["rows"])
    w = Counter()                       # datapoints riding on each group
    for k, v in sigs.items():
        for cid in k: w[cid] += v
    return {"benchmark": cls.name, "datapoints": n,
            "entities": {m: len(g["nodes"][m]) for m in g["nodes"]},
            "edges": {m: len(g["edges"][m]) for m in g["edges"]},
            "components": {m: len(g["groups"][m]) for m in g["groups"]},
            "components_total": n_comp,
            "largest_component_pct": {
                m: round(100 * max(w[comp_of[gr[0]]] for gr in g["groups"][m]) / n, 2)
                for m in g["groups"]},
            "signatures": len(sigs), "signatures_multi_component": len(multi),
            "max_signature_size": max(len(k) for k in sigs),
            "anything_to_solve": bool(multi),
            "thresholds": g["thresholds"]}
