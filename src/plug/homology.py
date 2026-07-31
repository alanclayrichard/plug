# split a benchmark by similarity: group its own proteins with mmseqs, then deal whole
# groups into train/val/test so no two splits hold sequences that look alike.
# different job from build_unlabeled_trainset.py — that one throws out reservoir sequences
# that look like a benchmark, this one carves up one benchmark into splits.
import random, subprocess, tempfile
from pathlib import Path

from . import config as c

# no e-value cutoff here. an e-value moves with how big the searched set is, so the same
# match scores far worse in a whole benchmark search than in the smaller split against
# split check anyone runs afterwards. any cutoff quietly drops borderline matches that
# then turn up as shared sequences across splits. identity and coverage don't move with
# the size of the search, so they decide on their own.
graph_evalue = "inf"


# {name: sequence} -> [[name, ...]], the groups of proteins that look alike.
# search the benchmark against itself and join any two that clear the identity and
# coverage bar. joins carry over, so something similar to something similar lands in the
# same group and can't end up in a different split.
# both coverage modes: mode 1 asks for the found protein to be covered, mode 2 the
# searched one. either way the two are related, so the groups take both.
def clusters(seqs, cov_modes=("1", "2")):
    keys = list(seqs)
    c.work.mkdir(parents=True, exist_ok=True)
    # every protein starts in its own group, then they merge as matches come in
    up = list(range(len(keys)))

    def find(i):
        while up[i] != i: up[i] = up[up[i]]; i = up[i]
        return i

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
                a, b = find(int(q[1:])), find(int(t[1:]))
                if a != b: up[a] = b
    out = {}
    for i, k in enumerate(keys): out.setdefault(find(i), []).append(k)
    return list(out.values())


# deal the groups out 80/10/10 -> {protein: "train"/"val"/"test"}.
# weight is how many items each protein has, so the shares come out even in items and not
# just in proteins. biggest group first, each going to whichever split is furthest below
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
