# bindingdb measured binding affinities, the articles subset, single chain targets only.
# one item per measurement. ships no split of its own.
import csv
import sys

from . import Benchmark
from .. import config as c

# this file has very long fields
csv.field_size_limit(sys.maxsize)

# how many chains the target has; anything above one is a complex, which is skipped
CHAINS = "Number of Protein Chains in Target (>1 implies a multichain complex)"
SEQUENCE = "BindingDB Target Chain Sequence 1"
# one affinity per row, taken in this order
AFFINITIES = (("IC50", "IC50 (nM)"), ("Ki", "Ki (nM)"), ("Kd", "Kd (nM)"), ("EC50", "EC50 (nM)"))


# an affinity cell like ">40" -> (">", 40.0), or None if it can't be read
def affinity(v):
    v = (v or "").strip()
    relation = "="
    while v and v[0] in "<>=~":
        relation, v = v[0], v[1:].strip()
    try: return relation, float(v)
    except ValueError: return None


class Bindingdb(Benchmark):
    name = "bindingdb"
    # thousands of rows carry no uniprot id, so the sequence itself names the protein
    key = "sequence"
    published = ()

    @classmethod
    def rows(cls):
        for r in cls.singlechain():
            hit = None
            for measure, col in AFFINITIES:
                pair = affinity(r.get(col))
                if pair: hit = (measure, pair); break
            # skip a row with no affinity that can be read
            if not hit: continue
            measure, (relation, value) = hit
            yield {"sequence": (r.get(SEQUENCE) or "").strip().upper(),
                   "uniprot": (r.get("UniProt (SwissProt) Primary ID of Target Chain 1") or "").strip(),
                   "target_name": (r.get("Target Name") or "").strip(),
                   "organism": (r.get("Target Source Organism According to Curator or DataSource") or "").strip(),
                   "smiles": (r.get("Ligand SMILES") or "").strip(),
                   "measure": measure, "value": value, "relation": relation}

    # every single chain row, affinity readable or not
    @classmethod
    def singlechain(cls):
        path = c.data / "bindingdb" / "BindingDB_BindingDB_Articles.tsv"
        with open(path, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                if (r.get(CHAINS) or "").strip() != "1": continue
                if not (r.get(SEQUENCE) or "").strip(): continue
                yield r

    # rows whose affinity couldn't be read are still in the benchmark's fasta, so their
    # sequences are held out as well
    def sequences(self):
        yield from super().sequences()
        for r in self.singlechain():
            yield (r.get(SEQUENCE) or "").strip().upper()
