# atlas molecular dynamics flexibility: one simulated pdb chain per file, with a
# per residue rmsf. ships no split of its own.
import csv

from . import Benchmark
from .. import config as c


class Atlas(Benchmark):
    name = "atlas"
    # every item is about one pdb chain
    key = "pdb_chain"
    published = ()
    # shape leaks too, not just sequence
    views = ("structure",)

    @classmethod
    def rows(cls):
        for path in sorted((c.data / "atlas" / "rmsf").glob("*.tsv")):
            seq, replicas = [], ([], [], [])
            with open(path, newline="") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    aa = (r.get("seq") or "").strip()
                    if not aa: continue
                    seq.append(aa)
                    for k in range(3): replicas[k].append(float(r[f"RMSF_R{k + 1}"]))
            if not seq: continue
            yield {"sequence": "".join(seq).upper(),
                   "pdb_chain": path.stem,
                   # the mean rmsf over the three runs, and the three runs themselves
                   "rmsf": [sum(v) / 3 for v in zip(*replicas)],
                   "rmsf_replicas": [list(v) for v in replicas]}

    # atlas is about how a shape moves, so two chains that fold the same way can't sit in
    # different splits however far their sequences have drifted. the split checks both.
    @classmethod
    def structures(cls):
        d = c.data / "atlas" / "pdb"
        found = {p.stem: p for p in sorted(d.glob("*.pdb"))} if d.exists() else {}
        if not found:
            raise FileNotFoundError(f"atlas has no structures in {d} — "
                                    "rerun data/download_atlas.sh to fetch them")
        return found
