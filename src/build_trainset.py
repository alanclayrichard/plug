# build a leakage-free training set: sample a reservoir, drop anything that leaks a functional
# test set, and resample until the quota is met.
#   RESERVOIR=seq     sample uniref90, drop seqs that hit a test seq       (mmseqs identity)
#   RESERVOIR=struct  sample the pdb, drop structures that hit a test one  (foldseek TM-score)
#   RESERVOIR=both    sample the pdb, drop via a union graph over the sequence + structure hits
# the reservoirs are never loaded whole; only the small test sets are searched. writes the kept
# set as a fasta keyed by uniprot code.
import random, shutil, subprocess, sys, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c

STRUCTS = (".pdb", ".cif", ".ent", ".mmcif", ".pdb.gz", ".cif.gz", ".ent.gz", ".mmcif.gz")  # structure file suffixes


def mmseqs(*args):  # run an mmseqs subcommand quietly
    subprocess.run([str(c.mmseqs), *args, "--threads", c.threads, "-v", "1"], check=True)


def foldseek(*args):  # run a foldseek subcommand quietly
    subprocess.run([str(c.foldseek), *args, "-v", "1"], check=True)


def test_target(tmp):  # the combined sequence test-set db (gpu-padded if GPU=1)
    fa = tmp / "tests.fasta"
    fa.write_text("".join(p.read_text() for p in sorted(c.formatted.glob("*.fasta"))))
    db = tmp / "tests_db"
    mmseqs("createdb", str(fa), str(db), "--dbtype", "1")
    if c.gpu:
        mmseqs("makepaddedseqdb", str(db), str(tmp / "tests_pad"))
        return tmp / "tests_pad"
    return db


def struct_target(tmp):  # foldseek db of the test-set structures (gpu-padded if GPU=1)
    db = tmp / "stests_db"
    foldseek("createdb", str(c.struct_tests), str(db))
    if c.gpu:
        foldseek("makepaddedseqdb", str(db), str(tmp / "stests_pad"))
        return tmp / "stests_pad"
    return db


def seq_edges(fa, target, tmp):  # mmseqs identity hits: candidate<->test edges above MIN_ID
    qdb, res, m8, out = tmp / "q_db", tmp / "q_res", tmp / "q.m8", tmp / "q_edges"
    mmseqs("createdb", str(fa), str(qdb), "--dbtype", "1")
    mmseqs("search", str(qdb), str(target), str(res), str(tmp / "qs"),
           "-s", c.sens, "-e", c.evalue, "-c", c.cov, "--cov-mode", c.cov_mode,
           "--max-seqs", "50", *(("--gpu", "1") if c.gpu else ()))
    mmseqs("convertalis", str(qdb), str(target), str(res), str(m8), "--format-output", "query,target,fident")
    with open(out, "w") as o:
        for ln in open(m8):
            q, t, ident = ln.split("\t")[:3]
            if float(ident) > c.min_id: o.write(f"{q}\t{t}\n")
    return out


def struct_edges(sub, target, tmp):  # foldseek TM-score hits: candidate<->test edges above TM
    res, m8, out = tmp / "s_res", tmp / "s.m8", tmp / "s_edges"
    foldseek("search", str(sub), str(target), str(res), str(tmp / "ss"),
             "-e", c.evalue, "--alignment-type", "1", "-a", "--threads", c.threads,  # tmalign + backtrace -> real TM-score
             *(("--gpu", "1") if c.gpu else ()))
    foldseek("convertalis", str(sub), str(target), str(res), str(m8),
             "--format-output", "query,target,alntmscore", "--threads", c.threads)
    with open(out, "w") as o:
        for ln in open(m8):
            q, t, score = ln.split("\t")[:3]
            if float(score) > c.tm: o.write(f"{q}\t{t}\n")
    return out


def union(edges):  # union graph over candidate<->test edge files; drop candidates whose component holds a test
    parent, tests, nodes = {}, set(), set()
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for e in edges:
        for ln in open(e):
            q, t = ln.split()[:2]
            q, t = "q:" + q, "t:" + t          # namespace: a candidate whose id equals a test's id must not merge into the test
            nodes |= {q, t}; tests.add(t)
            parent[find(q)] = find(t)
    roots = {find(t) for t in tests}
    return {n[2:] for n in nodes if n not in tests and find(n) in roots}


def sample_stream(items, n, rng):  # reservoir-sample n items from a stream in one pass (no full load)
    res, seen = [], 0
    for it in items:
        seen += 1
        if len(res) < n:
            res.append(it)
        else:
            j = rng.randrange(seen)
            if j < n: res[j] = it
        if seen % 20_000_000 == 0: print(f"    ...scanned {seen:,}", flush=True)
    return res


def sequence_sample(n, used, rng):  # reservoir-sample n unused, in-range uniref seqs
    items = ((a, s) for a, s in c.iter_fasta(c.uniref_fa)
             if a not in used and c.min_len <= len(s) <= c.max_len)
    pool = sample_stream(items, n, rng)
    used.update(a for a, _ in pool)
    return pool


def structure_sample(n, used, rng, tmp):  # reservoir-sample n unused structure files -> foldseek db + [(name, seq)]
    files = (p for p in c.pdb_dir.rglob("*") if p.name.endswith(STRUCTS) and p.name not in used)
    picks = sample_stream(files, n, rng)
    if not picks: return []
    used.update(p.name for p in picks)
    pool = tmp / "pool_structs"                    # a fresh db of just this round's structures
    shutil.rmtree(pool, ignore_errors=True); pool.mkdir()
    for p in picks: (pool / p.name).symlink_to(p.resolve())
    for f in tmp.glob("sub*"): f.unlink()
    foldseek("createdb", str(pool), str(tmp / "sub"))
    foldseek("convert2fasta", str(tmp / "sub"), str(tmp / "sub.fasta"))
    return list(c.iter_fasta(tmp / "sub.fasta"))


def code(uid):  # uniref90 cluster id -> representative uniprot accession
    return uid.split("_", 1)[1] if uid.startswith("UniRef90_") else uid


if __name__ == "__main__":
    rng = random.Random(c.seed)
    c.work.mkdir(parents=True, exist_ok=True)
    c.trainset.parent.mkdir(parents=True, exist_ok=True)
    seq, struct = c.reservoir_kind in ("seq", "both"), c.reservoir_kind in ("struct", "both")
    with tempfile.TemporaryDirectory(dir=c.work) as tmp:
        tmp = Path(tmp)
        qtarget = test_target(tmp) if seq else None       # sequence test db
        starget = struct_target(tmp) if struct else None  # structure test db

        def draw(n, used):  # sample a candidate pool -> [(name, seq)] (struct also leaves a foldseek db at tmp/sub)
            return structure_sample(n, used, rng, tmp) if struct else sequence_sample(n, used, rng)

        def leakers(pool):  # names in pool that leak, via the union of the sequence + structure graphs
            edges = []
            if seq:
                fa = tmp / "pool.fasta"
                with open(fa, "w") as o:
                    for a, s in pool: o.write(f">{a}\n{s}\n")
                edges.append(seq_edges(fa, qtarget, tmp))
            if struct:
                edges.append(struct_edges(tmp / "sub", starget, tmp))
            return union(edges)

        kept, used, rnd = {}, set(), 0          # code->seq, sampled ids, round counter
        while len(kept) < c.quota:
            rnd += 1
            n = int((c.quota - len(kept)) * c.oversample)
            print(f"round {rnd}: sampling {n:,} candidates (have {len(kept):,}/{c.quota:,})", flush=True)
            pool = draw(n, used)
            if not pool:
                print("reservoir exhausted", flush=True); break
            bad = leakers(pool)
            print(f"  dropped {len(bad):,} leakers of {len(pool):,}", flush=True)
            for a, s in pool:
                if a not in bad and len(kept) < c.quota: kept[code(a)] = s
        with open(c.trainset, "w") as o:
            for cd, s in kept.items(): o.write(f">{cd}\n{s}\n")
    print(f"wrote {c.trainset} — {len(kept):,} leakage-free seqs", flush=True)
