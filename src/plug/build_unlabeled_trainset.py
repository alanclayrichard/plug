# build a training set with no leakage: sample from a reservoir, throw out anything that
# looks like a benchmark protein, and keep sampling until there are enough.
#   seq     sample uniref90, drop sequences matching a benchmark    (mmseqs identity)
#   struct  sample the pdb, drop structures matching a benchmark    (foldseek tm score)
#   both    sample the pdb and drop on either kind of match
# the reservoir is never held in memory — it is read once, in a stream, and only the
# benchmark fastas are searched against.
# run: python -m plug.build_unlabeled_trainset
import random, shutil, subprocess, tempfile
from pathlib import Path

from . import config as c

# structure file endings we know how to read
STRUCTS = (".pdb", ".cif", ".ent", ".mmcif", ".pdb.gz", ".cif.gz", ".ent.gz", ".mmcif.gz")

# no e-value cutoff, same reason as in homology.py: identity and coverage decide, and a
# cutoff moves with how big the searched set is
search_evalue = "inf"

# the sequences that survived, keyed by uniprot code. loads and indexes like a benchmark
# does, so it can be handed straight to a torch DataLoader for pretraining.
class Trainset:
    def __init__(self, seqs):
        self.seqs = seqs
        self.items = [{"sequence": s, "id": name} for name, s in seqs.items()]

    # read a training set written earlier, so it doesn't have to be built again
    @classmethod
    def from_fasta(cls, path=None):
        return cls(dict(c.iter_fasta(path or c.trainset)))

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]

    # write them out as a fasta
    def write_fasta(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as o:
            for name, s in self.seqs.items(): o.write(f">{name}\n{s}\n")
        return path


# run an mmseqs command quietly
def mmseqs(*args):
    subprocess.run([str(c.mmseqs), *args, "--threads", c.threads, "-v", "1"], check=True)


# run a foldseek command quietly
def foldseek(*args):
    subprocess.run([str(c.foldseek), *args, "-v", "1"], check=True)


# every benchmark fasta as one searchable database (padded if searching on a gpu)
def benchmark_target(tmp):
    fa = tmp / "benchmarks.fasta"
    fa.write_text("".join(p.read_text() for p in sorted(c.formatted.glob("*.fasta"))))
    db = tmp / "bench_db"
    mmseqs("createdb", str(fa), str(db), "--dbtype", "1")
    if c.gpu:
        mmseqs("makepaddedseqdb", str(db), str(tmp / "bench_pad"))
        return tmp / "bench_pad"
    return db


# the benchmark structures as one searchable database
def structure_target(tmp):
    db = tmp / "sbench_db"
    foldseek("createdb", str(c.struct_tests), str(db))
    if c.gpu:
        foldseek("makepaddedseqdb", str(db), str(tmp / "sbench_pad"))
        return tmp / "sbench_pad"
    return db


# candidates that match a benchmark sequence closely enough to count, as pairs.
# COV_MODE=union runs both coverage modes and takes everything either one finds, the same
# way the benchmark splits are built: mode 1 catches a candidate holding a whole benchmark
# protein inside it, mode 2 a candidate that sits inside one.
def sequence_matches(fa, target, tmp):
    qdb, out = tmp / "q_db", tmp / "q_edges"
    modes = ("1", "2") if c.cov_mode == "union" else (c.cov_mode,)
    mmseqs("createdb", str(fa), str(qdb), "--dbtype", "1")
    with open(out, "w") as o:
        for m in modes:
            res, m8 = tmp / f"q_res{m}", tmp / f"q{m}.m8"
            mmseqs("search", str(qdb), str(target), str(res), str(tmp / f"qs{m}"),
                   "-s", c.sens, "-e", search_evalue, "-c", c.cov, "--cov-mode", m,
                   # one match is enough to drop a candidate, so don't collect them all
                   "--max-seqs", "50", *(("--gpu", "1") if c.gpu else ()))
            mmseqs("convertalis", str(qdb), str(target), str(res), str(m8),
                   "--format-output", "query,target,fident")
            for ln in open(m8):
                q, t, ident = ln.split("\t")[:3]
                if float(ident) > c.min_id: o.write(f"{q}\t{t}\n")
    return out


# candidates whose shape matches a benchmark structure closely enough to count
def structure_matches(sub, target, tmp):
    res, m8, out = tmp / "s_res", tmp / "s.m8", tmp / "s_edges"
    # tmalign with backtrace, so the score is a real tm score
    foldseek("search", str(sub), str(target), str(res), str(tmp / "ss"),
             "-e", search_evalue, "--alignment-type", "1", "-a", "--threads", c.threads,
             *(("--gpu", "1") if c.gpu else ()))
    foldseek("convertalis", str(sub), str(target), str(res), str(m8),
             "--format-output", "query,target,alntmscore", "--threads", c.threads)
    with open(out, "w") as o:
        for ln in open(m8):
            q, t, score = ln.split("\t")[:3]
            if float(score) > c.tm: o.write(f"{q}\t{t}\n")
    return out


# join the matches into groups and return the candidates sitting in a group with a
# benchmark protein. done over both kinds of match at once, so a candidate related to a
# benchmark protein by either one is thrown out.
def leakers(match_files):
    parent, benchmarks, nodes = {}, set(), set()

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for path in match_files:
        for ln in open(path):
            q, t = ln.split()[:2]
            # keep the two apart by name, so a candidate sharing an id with a benchmark
            # protein doesn't get merged into it by accident
            q, t = "q:" + q, "t:" + t
            nodes |= {q, t}
            benchmarks.add(t)
            parent[find(q)] = find(t)
    roots = {find(t) for t in benchmarks}
    return {n[2:] for n in nodes if n not in benchmarks and find(n) in roots}


# take n items from a stream without loading it, each one equally likely
def sample_stream(items, n, rng):
    picked, seen = [], 0
    for it in items:
        seen += 1
        if len(picked) < n:
            picked.append(it)
        else:
            j = rng.randrange(seen)
            if j < n: picked[j] = it
        if seen % 20_000_000 == 0: print(f"    ...read {seen:,}", flush=True)
    return picked


# n uniref sequences we haven't drawn before, within the length bounds
def sequence_sample(n, used, rng):
    items = ((a, s) for a, s in c.iter_fasta(c.uniref_fa)
             if a not in used and c.min_len <= len(s) <= c.max_len)
    pool = sample_stream(items, n, rng)
    used.update(a for a, _ in pool)
    return pool


# n structure files we haven't drawn before, plus their sequences, and a foldseek
# database of just this round left at tmp/sub
def structure_sample(n, used, rng, tmp):
    files = (p for p in c.pdb_dir.rglob("*") if p.name.endswith(STRUCTS) and p.name not in used)
    picks = sample_stream(files, n, rng)
    if not picks: return []
    used.update(p.name for p in picks)
    pool = tmp / "pool_structs"
    shutil.rmtree(pool, ignore_errors=True)
    pool.mkdir()
    for p in picks: (pool / p.name).symlink_to(p.resolve())
    for f in tmp.glob("sub*"): f.unlink()
    foldseek("createdb", str(pool), str(tmp / "sub"))
    foldseek("convert2fasta", str(tmp / "sub"), str(tmp / "sub.fasta"))
    return list(c.iter_fasta(tmp / "sub.fasta"))


# a uniref cluster name -> the uniprot code it stands for
def code(name):
    return name.split("_", 1)[1] if name.startswith("UniRef90_") else name


# sample until there are `quota` sequences that match no benchmark protein.
# everything defaults to the setting in config, pass one to override it.
def build_trainset(*, reservoir=None, quota=None, oversample=None, seed=None, verbose=True):
    reservoir = reservoir or c.reservoir_kind
    quota = quota or c.quota
    oversample = oversample or c.oversample
    rng = random.Random(c.seed if seed is None else seed)
    by_seq, by_struct = reservoir in ("seq", "both"), reservoir in ("struct", "both")
    c.work.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=c.work) as tmp:
        tmp = Path(tmp)
        seq_target = benchmark_target(tmp) if by_seq else None
        struct_target_db = structure_target(tmp) if by_struct else None
        kept, used, rnd = {}, set(), 0

        while len(kept) < quota:
            rnd += 1
            n = int((quota - len(kept)) * oversample)
            if verbose:
                print(f"round {rnd}: drawing {n:,} candidates (have {len(kept):,}/{quota:,})",
                      flush=True)
            pool = (structure_sample(n, used, rng, tmp) if by_struct
                    else sequence_sample(n, used, rng))
            if not pool:
                print("reservoir is used up", flush=True)
                break

            matches = []
            if by_seq:
                fa = tmp / "pool.fasta"
                with open(fa, "w") as o:
                    for a, s in pool: o.write(f">{a}\n{s}\n")
                matches.append(sequence_matches(fa, seq_target, tmp))
            if by_struct:
                matches.append(structure_matches(tmp / "sub", struct_target_db, tmp))
            bad = leakers(matches)
            if verbose:
                print(f"  threw out {len(bad):,} of {len(pool):,}", flush=True)

            for a, s in pool:
                if a not in bad and len(kept) < quota: kept[code(a)] = s

    return Trainset(kept)


if __name__ == "__main__":
    trainset = build_trainset()
    path = trainset.write_fasta(c.trainset)
    print(f"wrote {path} — {len(trainset):,} sequences with no leakage", flush=True)
