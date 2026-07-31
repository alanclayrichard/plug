# config file to load env variables for paths, thresholds, parameters, etc
import gzip, os
from pathlib import Path

# the repo parent path 
repo = Path(__file__).resolve().parent.parent.parent
# the data path where raw data was downloaded and formatted
data = repo / "data"

# load repo/.env variables and store as key value pairs to be used when config imported
_env = repo / ".env"
if _env.exists():
    for line in _env.read_text().splitlines():
        line = line.split(" #", 1)[0].strip()
        if "=" in line and not line.startswith("#"):
            k, v = line.split("=", 1)
            # set the key value pairs as os.environment dict
            os.environ.setdefault(k.strip(), os.path.expandvars(v.strip().strip("\"'")))

# path to the revervoirs to sample observations for training un/self supervised models
uniref_fa = Path(os.environ.get("UNIREF_FASTA", repo / "uniref" / "uniref90.fasta.gz"))
pdb_dir   = Path(os.environ.get("PDB_DIR", data / "pdb"))
reservoir_kind = os.environ.get("RESERVOIR", "seq")

# binary path for similarity searches
mmseqs    = os.environ.get("MMSEQS", "mmseqs")
foldseek  = os.environ.get("FOLDSEEK", "foldseek")

# variables guiding the similarity search
gpu       = os.environ.get("GPU", "")
threads  = os.environ.get("THREADS", "8")

# directory with formatted fastas or structures 
formatted    = Path(os.environ.get("FORMATTED_DIR", data / "formatted")) 
struct_tests = Path(os.environ.get("STRUCT_TESTS", data / "structs"))

# output file for leakage free unlabeled training set
trainset     = Path(os.environ.get("TRAINSET", repo / "trainset.fasta"))

# tmp working directory to do mmseqs work in (store databases and tmp files)
work         = Path(os.environ.get("WORK", repo / "work"))  

# leakage thresholds (adjust these to characterize similarity)
# maximum sequence similarity, coverage frac, mode and sensitivity
min_id   = float(os.environ.get("MIN_ID", "0.2"))
cov      = os.environ.get("COV", "0.8")
cov_mode = os.environ.get("COV_MODE", "union")
sens     = os.environ.get("SENS", "7.5")
# evalue   = os.environ.get("EVALUE", "1e-3")
# structural similarity via TM score
tm       = float(os.environ.get("TM", "0.5"))

# sampling parameters
# number of observations for unlabeled training set sample
quota      = int(os.environ.get("QUOTA", "1000000"))
# oversample fraction to do rounds of screening for compute 
oversample = float(os.environ.get("OVERSAMPLE", "1.5"))
# dont sample proteins shorter than this
min_len    = int(os.environ.get("MIN_LEN", "20"))
# dont sample proteins longer than this
max_len    = int(os.environ.get("MAX_LEN", "1000"))
# seed for all random processes
seed       = int(os.environ.get("SEED", "0"))

# where homology splits are saved
splits = Path(os.environ.get("SPLITS_DIR", data / "splits"))

# create a tag for thresholds for caching homology splits 
def split_tag(): return f"id{int(min_id * 100)}_cov{int(float(cov) * 100)}"

# (id, sequence) pairs from a fasta, gzipped or not — ppi and passerrank use this
def iter_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        acc, seq = None, []
        for line in f:
            if line.startswith(">"):
                if acc is not None: yield acc, "".join(seq)
                acc, seq = line[1:].split(None, 1)[0], []
            else: seq.append(line.strip())
        if acc is not None: yield acc, "".join(seq)
