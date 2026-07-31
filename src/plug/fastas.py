# turn every benchmark into one fasta of the sequences it holds, which is what the
# training set is kept clear of. a benchmark's own train split counts too — a model being
# tested should never have seen any of it, whichever split it sits in.
# run: python -m plug.fastas [name ...]
import sys

from . import config as c
from .benchmarks import REGISTRY


# write one benchmark's sequences to data/formatted/<name>.fasta, skipping repeats.
# returns where it went and how many went in.
def write(cls):
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


# one benchmark at a time, so only one is ever held in memory
def main(names=None):
    for name in names or REGISTRY:
        out, n = write(REGISTRY[name])
        print(f"{name}: {n} sequences -> {out}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or None)
