# format functional-benchmark test sets into uniform fastas for leakage checking.
# each benchmark = one ADAPTERS entry (name -> fn yielding its test seqs); write() dedups +
# uppercases to <formatted>/<name>.fasta. add a benchmark: add an adapter (+ a
# data/download_<name>.sh), then `python src/format_benchmark.py <name>`.
import csv, sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c

csv.field_size_limit(sys.maxsize)


def _col(path, *names):  # non-empty upper values from the given csv columns
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            for n in names:
                v = (row.get(n) or "").strip().upper()
                if v: yield v


def ddg():  # thermompnn megascale + fireprot test splits
    d = c.data / "ThermoMPNN" / "data_all" / "testing"
    yield from _col(d / "mega_test.csv", "aa_seq", "wt_seq")
    yield from _col(d / "fireprot_test.csv", "pdb_sequence")


def dms():  # proteingym target seqs (eval-only benchmark)
    yield from _col(c.data / "proteingym" / "DMS_substitutions.csv", "target_seq")


def go():  # bioreason-pro temporal-holdout test seqs
    p = c.data / "bioreason_go" / "test" / "data" / "test-00000-of-00001.parquet"
    for s in pd.read_parquet(p, columns=["sequence"])["sequence"]:
        if s: yield str(s).strip().upper()


def allobench():  # allobench allosteric/active-site proteins (one seq per row)
    yield from _col(c.data / "allobench" / "AlloBench.csv", "sequence")


def ppi():  # human ppi gold-standard: sequences of every protein in the test pairs
    d = c.data / "ppi"
    accs = {a for f in d.glob("Intra*_rr.txt") for line in open(f) for a in line.split()}
    for a, s in c.iter_fasta(d / "human_swissprot_oneliner.fasta"):
        if a in accs and s: yield s.upper()


def passerrank():  # passerrank allosteric proteins (uniprot seqs fetched from asd accessions)
    for _a, s in c.iter_fasta(c.data / "passerrank" / "passerrank.fasta"):
        if s: yield s.upper()


def pdbbind():  # LP-PDBBind protein-ligand binding affinity, leak-proof TEST split only
    with open(c.data / "pdbbind" / "LP_PDBBind.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r.get("new_split") != "test": continue
            s = (r.get("seq") or "").strip().upper()
            if s: yield s


def atlas():  # atlas MD sim proteins: each chain's sequence (the seq column of its per-residue RMSF.tsv)
    for p in sorted((c.data / "atlas" / "rmsf").glob("*.tsv")):
        with open(p, newline="") as f:
            seq = "".join((r.get("seq") or "").strip() for r in csv.DictReader(f, delimiter="\t"))
        if seq: yield seq.upper()


def bindingdb():  # BindingDB articles, single-chain target sequences (multichain dropped)
    chains = "Number of Protein Chains in Target (>1 implies a multichain complex)"
    with open(c.data / "bindingdb" / "BindingDB_BindingDB_Articles.tsv", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
            if (r.get(chains) or "").strip() != "1": continue
            s = (r.get("BindingDB Target Chain Sequence 1") or "").strip().upper()
            if s: yield s


# name -> test-seq iterator
ADAPTERS = {"ddg": ddg, "dms": dms, "go": go, "allobench": allobench, "ppi": ppi,
            "passerrank": passerrank, "pdbbind": pdbbind, "bindingdb": bindingdb, "atlas": atlas}


def write(name):  # dedup + write <formatted>/<name>.fasta, return (path, count)
    c.formatted.mkdir(parents=True, exist_ok=True)
    out = c.formatted / f"{name}.fasta"
    n = 0
    with open(out, "w") as o:
        for i, s in enumerate(dict.fromkeys(ADAPTERS[name]())):
            o.write(f">{name}_{i}\n{s}\n"); n += 1
    return out, n


if __name__ == "__main__":
    for name in sys.argv[1:] or list(ADAPTERS):
        out, n = write(name)
        print(f"{name}: {n} test seqs -> {out}")
