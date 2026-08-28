#!/bin/sh
# download atlas MD per-residue RMSF (conformational dynamics) into data/atlas.
# atlas ships each protein chain as one ~70MB analysis zip that's almost all trajectory —
# we http-range-extract just the ~5KB <pdb>_RMSF.tsv and the ~200KB <pdb>.pdb out of each
# remote zip, no trajectories pulled.
# each RMSF.tsv is one row per residue: seq (1-letter aa), RMSF_R1, RMSF_R2, RMSF_R3 (Å, 3 MD replicas).
# the .pdb is the starting structure, used to split atlas on shape as well as sequence.
set -e
d="$(dirname "$0")/atlas"
mkdir -p "$d/rmsf" "$d/pdb"
# the list of PDB chains in ATLAS v2 (Nov 2024)
curl -fsSL https://www.dsimb.inserm.fr/ATLAS/data/download/distributions/2024_11_18_ATLAS_pdb.txt -o "$d/atlas_pdb.txt"
python3 - "$d" <<'PY'
import io, os, sys, urllib.request, zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

d = Path(sys.argv[1]); out, structs = d / "rmsf", d / "pdb"
api = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{}"
chains = [ln.strip() for ln in open(d / "atlas_pdb.txt") if ln.strip()]


class RangeFile(io.RawIOBase):  # seekable file backed by http range requests (zipfile only reads the bytes it needs)
    def __init__(self, url):
        self.url, self.pos = url, 0
        r = urllib.request.urlopen(urllib.request.Request(url, headers={"Range": "bytes=0-0"}))
        self.size = int(r.headers["Content-Range"].split("/")[1])

    def seek(self, off, whence=0): self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]; return self.pos
    def tell(self): return self.pos
    def seekable(self): return True

    def read(self, n=-1):
        n = self.size - self.pos if n < 0 else n
        if not n: return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        data = urllib.request.urlopen(req).read(); self.pos += len(data); return data


def pull(chain):  # extract only the rmsf table and the structure from the remote analysis zip
    tsv, pdb = out / f"{chain}.tsv", structs / f"{chain}.pdb"
    if tsv.exists() and pdb.exists(): return "skip"
    try:
        zf = zipfile.ZipFile(RangeFile(api.format(chain)))
        for dst, end in ((tsv, "_RMSF.tsv"), (pdb, ".pdb")):
            if dst.exists(): continue
            name = next(n for n in zf.namelist() if n.endswith(end))
            dst.write_bytes(zf.read(name))
        return "ok"
    except Exception as e:
        print(f"  {chain}: {e}", flush=True); return "fail"


# each file is a handful of small http range requests, so this is all waiting on the
# network — more workers finish sooner. ATLAS_JOBS=32 turns it up on a good connection.
with ThreadPoolExecutor(max_workers=int(os.environ.get("ATLAS_JOBS", "8"))) as ex:
    for i, _ in enumerate(ex.map(pull, chains), 1):
        if i % 100 == 0: print(f"  {i}/{len(chains)}", flush=True)
print(f"atlas: {len(list(out.glob('*.tsv')))} RMSF profiles -> {out}")
print(f"atlas: {len(list(structs.glob('*.pdb')))} structures -> {structs}")
PY
