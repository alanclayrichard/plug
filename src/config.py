# paths + knobs, all env-overridable (see .env.example). 
import gzip, os
from pathlib import Path

repo = Path(__file__).resolve().parent.parent
data = repo / "data"

# load repo/.env (KEY=VALUE, $REFS and inline #comments allowed) so scripts just run
_env = repo / ".env"
if _env.exists():
    for line in _env.read_text().splitlines():
        line = line.split(" #", 1)[0].strip()
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            os.environ.setdefault(k.strip(), os.path.expandvars(v.strip().strip("\"'")))

# reservoirs to sample + the search binaries (point at gpu builds to use GPU=1)
uniref_fa = Path(os.environ.get("UNIREF_FASTA", repo / "uniref" / "uniref90.fasta.gz"))  # sequence reservoir (gzipped fasta)
pdb_dir   = Path(os.environ.get("PDB_DIR", data / "pdb"))               # structure reservoir (directory of pdb/mmcif files, see download_pdb.sh)
mmseqs    = os.environ.get("MMSEQS", "mmseqs")
foldseek  = os.environ.get("FOLDSEEK", "foldseek")
gpu       = os.environ.get("GPU", "")             # GPU=1 (+ gpu mmseqs/foldseek builds) runs the check on cuda
reservoir_kind = os.environ.get("RESERVOIR", "seq")   # seq (uniref90) | struct (pdb) | both (union graph)

# repo inputs / outputs
formatted    = Path(os.environ.get("FORMATTED_DIR", data / "formatted"))   # formatted test-set fastas
struct_tests = Path(os.environ.get("STRUCT_TESTS", data / "structs"))      # test-set structures for the foldseek check
trainset     = Path(os.environ.get("TRAINSET", repo / "trainset.fasta"))   # final leakage-free set
work         = Path(os.environ.get("WORK", repo / "work"))                 # mmseqs/foldseek scratch (ephemeral)

# leakage rule + mmseqs search
threads  = os.environ.get("THREADS", "8")
sens     = os.environ.get("SENS", "7.5")            # mmseqs sensitivity
evalue   = os.environ.get("EVALUE", "1e-3")         # mmseqs e-value cutoff
cov      = os.environ.get("COV", "0.8")             # min coverage to count as leakage
cov_mode = os.environ.get("COV_MODE", "1")          # 1 = cover the test seq (catches test-as-subdomain), 0 = both, 2 = train seq
min_id   = float(os.environ.get("MIN_ID", "0.2"))   # >this identity to a test seq counts as sequence leakage
tm       = float(os.environ.get("TM", "0.5"))       # >this foldseek TM-score to a test structure counts as structural leakage

# sampling: collect this many leakage-free seqs, oversampling each round to cover dropped leakers
quota      = int(os.environ.get("QUOTA", "1000000"))
oversample = float(os.environ.get("OVERSAMPLE", "1.5"))
min_len    = int(os.environ.get("MIN_LEN", "20"))      # only sample proteins with this many residues..
max_len    = int(os.environ.get("MAX_LEN", "1000"))    # ..up to this many (drops fragments + giant seqs)
seed       = int(os.environ.get("SEED", "0"))


def iter_fasta(path):  # (id, sequence) pairs from an optionally-gzipped fasta
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        acc, seq = None, []
        for line in f:
            if line.startswith(">"):
                if acc is not None: yield acc, "".join(seq)
                acc, seq = line[1:].split(None, 1)[0], []
            else: seq.append(line.strip())
        if acc is not None: yield acc, "".join(seq)
