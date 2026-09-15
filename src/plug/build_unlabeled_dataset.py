# sample a reservoir into a leakage-free unlabeled dataset
'''
This is the same leakage rule as splitting.py and the same thresholds, pointed
the other way: that one carves one benchmark into splits, this one holds a
whole reservoir clear of every benchmark at once. Run it as "python -m
plug.build_unlabeled_dataset". The reservoir is sampled, everything matching a
benchmark protein is thrown out, and sampling repeats until the quota is met;
RESERVOIR=seq draws from uniref90 and drops on mmseqs identity, struct draws
from the pdb and drops on foldseek tm score, and both draws structures and
drops on either kind of match. The reservoir is never held in memory — it is
read once, in a stream, and only the benchmark fastas are searched against, so
a hundred million sequences cost one pass and not one database. A candidate is
thrown out when it lands in a group with a benchmark protein, which is the
transitive rule again, and the two are kept apart by name so a candidate
sharing an id with a benchmark protein is not merged into it by accident. The
formatted fasta folder is the whole screen and nothing else checks it, so a
benchmark missing from it is rebuilt here rather than silently skipped.
'''

import random, shutil, subprocess, tempfile
from pathlib import Path

from . import config as c

# structure file endings we know how to read
STRUCTS = (".pdb", ".cif", ".ent", ".mmcif", ".pdb.gz", ".cif.gz", ".ent.gz", ".mmcif.gz")

# no e-value cutoff, same reason as in graphs.py: identity and coverage decide, and a
# cutoff moves with how big the searched set is
search_evalue = "inf"

class Unlabeled:
    '''
    the sequences that survived, keyed by uniprot code. loads and indexes like
    a benchmark does, so it can be handed straight to a torch DataLoader for
    pretraining.
    '''

    def __init__(self, seqs):
        '''
        inputs:
            seqs: dict[str, str] -- "<code> <what it is>" -> sequence
        outputs:
            none
        '''
        self.seqs = seqs
        # the key is "<code> <what it is>", so the id an item carries is the code alone
        self.items = [{"sequence": s, "id": name.split(None, 1)[0], "name": name}
                      for name, s in seqs.items()]

    @classmethod
    def from_fasta(cls, path=None):
        '''
        read a dataset written earlier, so it doesn't have to be built again.

        inputs:
            path: Path | None -- the fasta, defaulting to config's trainset
        outputs:
            data: Unlabeled -- the loaded dataset
        '''
        return cls({f"{a} {d}".rstrip(): s
                    for a, d, s in c.iter_fasta(path or c.trainset, desc=True)})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]

    def write_fasta(self, path):
        '''
        write the dataset out as a fasta.

        inputs:
            path: str | Path -- where to write it
        outputs:
            path: Path -- where it went
        '''
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as o:
            for name, s in self.seqs.items(): o.write(f">{name}\n{s}\n")
        return path

def mmseqs(*args):
    '''
    run an mmseqs command quietly.

    inputs:
        args: str -- the command and its arguments
    outputs:
        none
    '''
    subprocess.run([str(c.mmseqs), *args, "--threads", c.threads, "-v", "1"], check=True)

def foldseek(*args):
    '''
    run a foldseek command quietly.

    inputs:
        args: str -- the command and its arguments
    outputs:
        none
    '''
    subprocess.run([str(c.foldseek), *args, "-v", "1"], check=True)

def benchmark_fasta(cls):
    '''
    write one benchmark's sequences to data/formatted/<name>.fasta, skipping
    repeats. this is what the reservoir is held clear of, and a benchmark's own
    train split counts too — a model being tested should never have seen any of
    it, whichever split it sits in.

    inputs:
        cls: type -- a Benchmark subclass
    outputs:
        path: Path -- where the fasta went
        n: int -- how many distinct sequences went in
    '''
    c.formatted.mkdir(parents=True, exist_ok=True)
    out = c.formatted / f"{cls.name}.fasta"
    seen = set()
    with open(out, "w") as o:
        # the whole benchmark, no split
        for s in cls(split="all", how="whole").sequences():
            s = (s or "").strip().upper()
            if not s or s in seen: continue
            o.write(f">{cls.name}_{len(seen)}\n{s}\n")
            seen.add(s)
    return out, len(seen)

def benchmark_target(tmp):
    '''
    every benchmark fasta as one searchable database, padded if searching on a
    gpu. a benchmark with no fasta here is one the dataset was never held clear
    of, and it fails silently — the build still finishes and the count still
    looks right, so the only sign is inflated scores later. so whatever is out
    of step gets rebuilt, and said out loud.

    inputs:
        tmp: Path -- scratch directory
    outputs:
        db: Path -- the mmseqs database to search against
    '''
    # imported here, not at the top: the struct reservoir never calls this and shouldn't
    # have to have torch installed to build a trainset
    from .benchmarks import REGISTRY
    found = {p.stem: p for p in sorted(c.formatted.glob("*.fasta"))}
    for name in sorted(set(REGISTRY) - set(found)):
        out, n = benchmark_fasta(REGISTRY[name])
        print(f"  wrote {n:,} sequences for {name} -> {out}", flush=True)
        found[name] = out
    for name in sorted(set(found) - set(REGISTRY)):
        print(f"  ignoring {found.pop(name).name}: no benchmark by that name", flush=True)
    fa = tmp / "benchmarks.fasta"
    fa.write_text("".join(p.read_text() for p in found.values()))
    db = tmp / "bench_db"
    mmseqs("createdb", str(fa), str(db), "--dbtype", "1")
    if c.gpu:
        mmseqs("makepaddedseqdb", str(db), str(tmp / "bench_pad"))
        return tmp / "bench_pad"
    return db

def structure_target(tmp):
    '''
    the benchmark structures as one searchable database.

    inputs:
        tmp: Path -- scratch directory
    outputs:
        db: Path -- the foldseek database to search against
    '''
    db = tmp / "sbench_db"
    foldseek("createdb", str(c.struct_tests), str(db))
    if c.gpu:
        foldseek("makepaddedseqdb", str(db), str(tmp / "sbench_pad"))
        return tmp / "sbench_pad"
    return db

def sequence_matches(fa, target, tmp):
    '''
    candidates that match a benchmark sequence closely enough to count, as
    pairs. COV_MODE=union runs both coverage modes and takes everything either
    one finds, the same way the benchmark splits are built: mode 1 catches a
    candidate holding a whole benchmark protein inside it, mode 2 a candidate
    that sits inside one.

    inputs:
        fa: Path -- the candidate fasta
        target: Path -- the benchmark database
        tmp: Path -- scratch directory
    outputs:
        path: Path -- a tsv of (candidate, benchmark protein) pairs
    '''
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

def structure_matches(sub, target, tmp):
    '''
    candidates whose shape matches a benchmark structure closely enough to
    count, as pairs.

    inputs:
        sub: Path -- this round's candidate database
        target: Path -- the benchmark structure database
        tmp: Path -- scratch directory
    outputs:
        path: Path -- a tsv of (candidate, benchmark protein) pairs
    '''
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

def leakers(match_files):
    '''
    join the matches into groups and give back the candidates sitting in a
    group with a benchmark protein. done over both kinds of match at once, so
    a candidate related to a benchmark protein by either one is thrown out.

    inputs:
        match_files: list[Path] -- the tsvs of matched pairs
    outputs:
        bad: set[str] -- candidate ids to drop
    '''
    parent, benchmarks, nodes = {}, set(), set()

    # follow a name to its group, flattening the path on the way
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

def sample_stream(items, n, rng):
    '''
    reservoir sampling: take n items from a stream without loading it, each one
    equally likely.

    inputs:
        items: iterable -- the stream to draw from
        n: int -- how many to take
        rng: random.Random -- the seeded generator
    outputs:
        picked: list -- up to n items
    '''
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

def sequence_sample(n, used, rng):
    '''
    n uniref sequences we haven't drawn before, within the length bounds. the
    header's description rides along so the written fasta can say what each
    sequence is.

    inputs:
        n: int -- how many to draw
        used: set[str] -- accessions already drawn, updated in place
        rng: random.Random -- the seeded generator
    outputs:
        pool: list[tuple[str, str, str]] -- (accession, description, sequence)
    '''
    items = ((a, d, s) for a, d, s in c.iter_fasta(c.uniref_fa, desc=True)
             if a not in used and c.min_len <= len(s) <= c.max_len)
    pool = sample_stream(items, n, rng)
    used.update(a for a, _, _ in pool)
    return pool

def structure_sample(n, used, rng, tmp):
    '''
    n structure files we haven't drawn before, plus their sequences, leaving a
    foldseek database of just this round at tmp/sub.

    inputs:
        n: int -- how many to draw
        used: set[str] -- file names already drawn, updated in place
        rng: random.Random -- the seeded generator
        tmp: Path -- scratch directory
    outputs:
        pool: list[tuple[str, str, str]] -- (name, "", sequence)
    '''
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
    # same (id, description, sequence) shape the sequence reservoir returns; a structure
    # file's header is just its name, so there is nothing to describe
    return [(a, "", s) for a, s in c.iter_fasta(tmp / "sub.fasta")]

def code(name):
    '''
    a uniref cluster name -> the uniprot code it stands for.

    inputs:
        name: str -- e.g. "UniRef90_A0ACG4HM60"
    outputs:
        code: str -- e.g. "A0ACG4HM60"
    '''
    return name.split("_", 1)[1] if name.startswith("UniRef90_") else name

def build_unlabeled_dataset(*, reservoir=None, quota=None, oversample=None, seed=None,
                           verbose=True):
    '''
    sample until there are `quota` sequences that match no benchmark protein.
    everything defaults to the setting in config, pass one to override it.

    inputs:
        reservoir: str | None -- "seq", "struct" or "both"
        quota: int | None -- how many clean sequences to end up with
        oversample: float | None -- how far over the shortfall to draw each round
        seed: int | None -- the sampling seed
        verbose: bool -- print progress
    outputs:
        data: Unlabeled -- the sequences that survived
    '''
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
                # the accession alone is the query name: mmseqs cuts a header at the first
                # space, and the names it reports back have to match what `bad` holds
                with open(fa, "w") as o:
                    for a, _, s in pool: o.write(f">{a}\n{s}\n")
                matches.append(sequence_matches(fa, seq_target, tmp))
            if by_struct:
                matches.append(structure_matches(tmp / "sub", struct_target_db, tmp))
            bad = leakers(matches)
            if verbose:
                print(f"  threw out {len(bad):,} of {len(pool):,}", flush=True)

            # the uniprot code names the sequence and the description says what it is, so
            # the written header reads ">A0ACG4HM60 LacI family ... Tax=Nocardia ...".
            # reading the fasta back keys on the first word, so the code still names it.
            for a, d, s in pool:
                if a not in bad and len(kept) < quota:
                    kept[f"{code(a)} {d}".rstrip()] = s

    return Unlabeled(kept)

if __name__ == "__main__":
    data = build_unlabeled_dataset()
    path = data.write_fasta(c.trainset)
    print(f"wrote {path} — {len(data):,} sequences with no leakage", flush=True)
