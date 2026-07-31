# proteingym deep mutational scanning. one item per measured variant, grouped by assay.
# every assay targets one protein, so the assay id is what names the protein here.
# ships no split of its own.
import csv
import sys

from . import Benchmark
from .. import config as c

# target sequences are long enough to trip the default csv field limit
csv.field_size_limit(sys.maxsize)


class Dms(Benchmark):
    name = "dms"
    # every item is about the protein its assay targets
    key = "dms_id"
    published = ()

    # pass assay to keep one kind of measurement (Activity, Binding, Expression,
    # OrganismalFitness, Stability), or dms_id to keep a single assay
    def __init__(self, assay=None, dms_id=None, split="test", how=None):
        super().__init__(split, how)
        if assay: self.items = [it for it in self.items if it["assay"] == assay]
        if dms_id: self.items = [it for it in self.items if it["dms_id"] == dms_id]

    @classmethod
    def rows(cls):
        folder = c.data / "proteingym"
        # one row per assay: its target sequence, what it measured, and the protein
        with open(folder / "DMS_substitutions.csv", newline="") as f:
            meta = {r["DMS_id"]: (r["target_seq"].strip().upper(),
                                  r["coarse_selection_type"], r["UniProt_ID"])
                    for r in csv.DictReader(f)}
        for dms_id, (seq, assay, uniprot) in meta.items():
            path = folder / "DMS_ProteinGym_substitutions" / f"{dms_id}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    if not (r.get("DMS_score") or "").strip(): continue
                    yield {"sequence": seq, "mutation": r["mutant"],
                           "score": float(r["DMS_score"]), "assay": assay,
                           "dms_id": dms_id, "uniprot": uniprot}
