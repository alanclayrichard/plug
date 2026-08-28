# split a benchmark by similarity: link its own entries to each other, then deal whole
# groups into train/val/test so no two splits hold entries that look alike.
# different job from build_unlabeled_trainset.py — that one throws out reservoir sequences
# that look like a benchmark, this one carves up one benchmark into splits.
#
# entries can look alike in more than one way, and the rule is the same every time: search
# the benchmark against itself, link anything over the bar in either direction, and keep
# whole groups of linked things together.
#   proteins  — sequence (mmseqs, both coverage modes), plus shape (foldseek, tm read both
#               ways round) for a benchmark that ships structures. both describe the same
#               proteins, so they share one graph and either kind of match counts.
#   ligands   — chemistry (ecfp4 tanimoto). a different kind of thing from a protein, so
#               its own graph: put them together and one shared molecule chains protein to
#               protein until the whole benchmark is one group.
# proteins are dealt first, then each ligand cluster follows the proteins it was measured
# against. an item is kept only when its protein and its ligand agree on a split, which is
# the same rule a ppi pair already obeys — so nothing similar spans a split, either way.
import random, re, subprocess, tempfile
from collections import Counter
from pathlib import Path

from . import config as c

# no e-value cutoff here. an e-value moves with how big the searched set is, so the same
# match scores far worse in a whole benchmark search than in the smaller split against
# split check anyone runs afterwards. any cutoff quietly drops borderline matches that
# then turn up as shared sequences across splits. identity and coverage don't move with
# the size of the search, so they decide on their own.
graph_evalue = "inf"


# ---------- the ways two entries can look alike; each yields pairs of names ----------

# {name: sequence} -> pairs mmseqs says are the same family.
# search the benchmark against itself and link any two that clear the identity and
# coverage bar. both coverage modes: mode 1 asks for the found protein to be covered,
# mode 2 the searched one. either way the two are related, so both get taken.
def sequence_links(seqs, cov_modes=("1", "2")):
    keys = list(seqs)
    c.work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=c.work) as tmp:
        tmp = Path(tmp)
        fa, m8 = tmp / "in.fasta", tmp / "hits.m8"
        # number the proteins: their names aren't always safe to put in a fasta
        with open(fa, "w") as o:
            for i, k in enumerate(keys): o.write(f">s{i}\n{seqs[k]}\n")
        for m in cov_modes:
            # no cap on matches per protein (--max-seqs). in a benchmark with many near
            # copies a cap fills up with those and hides the distant ones, and a hidden
            # match is a leak. mmseqs sizes the search by proteins x cap, so this only
            # stays affordable while the benchmarks are small.
            subprocess.run([str(c.mmseqs), "easy-search", str(fa), str(fa), str(m8),
                            str(tmp / f"s{m}"), "-s", c.sens, "-e", graph_evalue, "-c", c.cov,
                            "--cov-mode", m, "--max-seqs", str(len(keys)), "--threads", c.threads,
                            "-v", "1", "--format-output", "query,target,fident"], check=True)
            for ln in open(m8):
                q, t, ident = ln.split()[:3]
                # too different to count as alike
                if float(ident) <= c.min_id: continue
                yield keys[int(q[1:])], keys[int(t[1:])]


# {name: path to a structure} -> pairs foldseek says fold the same way.
# same idea as above but on shape, so two proteins that fold alike can't be split apart
# even when their sequences have drifted past the identity bar.
# the tm score is read both ways round, which is the shape version of the two coverage
# modes on sequences: qtmscore measures the match over the searched chain, ttmscore over
# the found one, and either clearing the bar means they are the same shape. foldseek's
# own alntmscore measures it over the matched stretch instead, which scores a single
# domain buried in a much longer chain as a whole-chain match — six times too many hits.
def structure_links(paths):
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


# {name: [smiles, ...]} -> pairs whose ligands are the same chemistry.
# two entries sharing a molecule are obviously linked; so are two whose molecules score
# over the tanimoto bar, which is what stops a test ligand reappearing in train as a
# near copy of itself.
def ligand_links(ligands):
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

    # ...and so does a molecule close enough to another one. done in blocks: the shared
    # bits of every pair at once, which is the only way this finishes on the big sets.
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


# ---------- the groups, and dealing them out ----------

# every name, plus any number of link sources -> the groups that have to stay together.
# links carry over, so something alike something alike lands in the same group and can't
# end up in a different split.
def merge(names, *links):
    at = {k: i for i, k in enumerate(names)}
    # everything starts in its own group, then they join up as links come in
    up = list(range(len(names)))

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


# a benchmark whose items pair a protein with a ligand has two different kinds of thing in
# it, and they must not go in the same graph: one shared ligand would chain protein to
# protein until the whole benchmark is a single group. so the ligands get their own graph
# and are dealt afterwards, each cluster going where its own items' proteins already went.
# the promise still holds both ways — a protein family lands in one split and a chemical
# series lands in one split — but far fewer items have to be dropped for disagreeing.
def follow(groups, votes, frac=(0.8, 0.1, 0.1)):
    names = ("train", "val", "test")
    total = sum(sum(v.values()) for v in votes.values())
    want = {s: f * total for s, f in zip(names, frac)}
    have = dict.fromkeys(names, 0)
    out = {}
    # biggest first, same as deal
    for g in sorted(groups, key=lambda g: -sum(sum(votes[k].values()) for k in g)):
        tally = Counter()
        for k in g: tally.update(votes[k])
        # where most of its items already are, unless that part is over its share
        pick = max(names, key=lambda s: (tally[s] - max(0, have[s] - want[s]),
                                         want[s] - have[s]))
        for k in g: out[k] = pick
        have[pick] += sum(tally.values())
    return out


# deal the groups out 80/10/10 -> {name: "train"/"val"/"test"}.
# weight is how many items each name has, so the shares come out even in items and not
# just in names. biggest group first, each going to whichever split is furthest below
# its share, and a whole group always lands in one split.
def deal(groups, weight, frac=(0.8, 0.1, 0.1)):
    # shuffle first so groups of the same size don't always fall the same way round
    random.Random(c.seed).shuffle(groups)
    groups.sort(key=lambda g: -sum(weight[k] for k in g))
    names = ("train", "val", "test")
    want = {s: f * sum(weight.values()) for s, f in zip(names, frac)}
    have = dict.fromkeys(names, 0)
    out = {}
    for g in groups:
        s = max(names, key=lambda s: want[s] - have[s])
        have[s] += sum(weight[k] for k in g)
        out.update(dict.fromkeys(g, s))
    return out
