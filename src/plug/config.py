# paths, thresholds and sampling parameters, read from the environment
'''
Every setting plug uses lives here, so nothing else has to know where a file
sits or what counts as too similar. The repo's .env is read into os.environ on
import and each setting falls back to a default, which means a fresh clone runs
without being configured first. The thresholds are the interesting part: they
are what "leakage" means in this repo, and they are read by the reservoir
screen and by the benchmark splits alike, so the two can never drift apart.
Changing one changes the name of every file built from it, so a split made at
one identity is never silently reused at another.
'''

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

# path to the reservoirs to sample observations for training un/self supervised models
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
min_id   = float(os.environ.get("MIN_ID", "0.3"))
cov      = os.environ.get("COV", "0.8")
cov_mode = os.environ.get("COV_MODE", "union")
sens     = os.environ.get("SENS", "7.5")
# structural similarity via TM score
tm       = float(os.environ.get("TM", "0.5"))
# the bar the homology split uses for structures, deliberately higher than the one above.
# a split joins up matches that carry over (alike to alike to alike), so a loose bar drags
# every fold into one group — at 0.5 atlas comes out as a single cluster of 99% of itself.
split_tm = float(os.environ.get("SPLIT_TM", "0.6"))
# ligand similarity via tanimoto on ecfp4 (morgan radius 2, 2048 bits, stereochemistry
# ignored). the fingerprint has to be named with the number: the same two molecules score
# 0.43 on ecfp4 and 0.79 on the rdkit fingerprint, so a bare "tanimoto 0.5" means nothing.
# random pairs of these ligands top out around 0.23 (bindingdb) to 0.30 (lp-pdbbind), so
# 0.5 is well clear of noise. it is deliberately stricter than the ~0.4 a chemist would
# call similar, because the split closes links up transitively: at 0.4 bindingdb's biggest
# ligand cluster is 87% of the set, too big to fit a train split at all.
tanimoto = float(os.environ.get("TANIMOTO", "0.5"))

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

# sifts table for turning pdb residue numbers into uniprot ones (see download_sifts.sh)
sifts = Path(os.environ.get("SIFTS", data / "sifts" / "uniprot_segments_observed.csv.gz"))

def iter_fasta(path, desc=False):
    '''
    read a fasta, gzipped or not, one record at a time. the id is the first
    whitespace-delimited word of the header; desc=True also yields the rest of
    it, which is how the training set keeps uniref's protein name and organism.

    inputs:
        path: str | Path -- the fasta to read
        desc: bool -- yield the header's description as well as the id
    outputs:
        records: iterator[tuple[str, str]] -- (id, sequence), or
                 (id, description, sequence) when desc is True
    '''
    op = gzip.open if str(path).endswith(".gz") else open

    # one record, shaped by whether the caller wanted the description
    def out(acc, note, seq):
        return (acc, note, "".join(seq)) if desc else (acc, "".join(seq))

    with op(path, "rt") as f:
        acc, note, seq = None, "", []
        for line in f:
            if line.startswith(">"):
                if acc is not None: yield out(acc, note, seq)
                head = line[1:].strip().split(None, 1)
                acc, note, seq = head[0], (head[1] if len(head) > 1 else ""), []
            else: seq.append(line.strip())
        if acc is not None: yield out(acc, note, seq)
