# bindingdb measured affinities, single chain targets, one item per measurement
'''
The articles subset of BindingDB, restricted to single chain targets — anything
above one chain is a complex and is skipped. It ships no split of its own.
Thousands of rows carry no uniprot id, so the sequence itself is what names the
protein here. One affinity is taken per row in the order IC50, Ki, Kd, EC50,
and the cell is parsed for a relation as well as a number since plenty read
">40". A smiles never has a space in it, but a couple hundred rows tack an
annotation block on the end describing stereochemistry, which is not part of
the molecule and makes the whole string unreadable, so only the molecule is
kept. The label is about a protein and a ligand together, so this benchmark
declares the ligand view and the split keeps the chemistry apart too. Rows
whose affinity couldn't be read are not items but their sequences are still in
the benchmark, so sequences() holds those out as well.
'''

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

# cells that say "no molecule here" rather than naming one. left alone they would read as
# a ligand that lots of rows share, which would wrongly tie those rows together.
BLANK = {"", "null", "none", "nan", "n/a"}

def smiles(v):
    '''
    keep only the molecule from a smiles cell, dropping the annotation block a
    couple hundred rows tack on the end ("CCO |TLB:19:20:24|"). a cell that
    says "no molecule here" comes back empty rather than as a shared ligand.

    inputs:
        v: str -- the raw cell
    outputs:
        smiles: str -- the molecule, or "" if there isn't one
    '''
    v = (v or "").strip()
    v = v.split(maxsplit=1)[0] if v else ""
    return "" if v.lower() in BLANK else v

def affinity(v):
    '''
    read an affinity cell, which carries a relation as well as a number.

    inputs:
        v: str -- the raw cell, e.g. ">40"
    outputs:
        pair: tuple[str, float] | None -- (relation, value), or None if it
              can't be read
    '''
    v = (v or "").strip()
    relation = "="
    while v and v[0] in "<>=~":
        relation, v = v[0], v[1:].strip()
    try: return relation, float(v)
    except ValueError: return None

class Bindingdb(Benchmark):
    '''measured binding affinities for single chain targets.'''

    name = "bindingdb"
    # thousands of rows carry no uniprot id, so the sequence itself names the protein
    key = "sequence"
    # what this benchmark is scored on, for the stratified split
    label = "value"
    published = ()
    # the ligand leaks too, not just the protein
    views = ("ligand",)

    @classmethod
    def rows(cls):
        '''
        one item per measurement whose affinity can be read, taking the first
        of IC50, Ki, Kd, EC50 that is there.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, ids, smiles, measure, value and
                  relation
        '''
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
                   "smiles": smiles(r.get("Ligand SMILES")),
                   "measure": measure, "value": value, "relation": relation}

    @classmethod
    def ligands(cls):
        '''
        the molecules each protein was measured against, so the split can keep
        the chemistry apart as well as the proteins.

        inputs:
            none
        outputs:
            ligands: dict[str, list[str]] -- sequence -> its smiles
        '''
        out = {}
        for it in cls.rows():
            if it["smiles"]: out.setdefault(it[cls.key], []).append(it["smiles"])
        return out

    @classmethod
    def singlechain(cls):
        '''
        every single chain row of the articles subset, affinity readable or
        not.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- the raw tsv rows
        '''
        path = c.data / "bindingdb" / "BindingDB_BindingDB_Articles.tsv"
        with open(path, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                if (r.get(CHAINS) or "").strip() != "1": continue
                if not (r.get(SEQUENCE) or "").strip(): continue
                yield r

    def sequences(self):
        '''
        rows whose affinity couldn't be read are still in the benchmark, so
        their sequences are held out as well.

        inputs:
            none
        outputs:
            seqs: iterator[str] -- every single chain target sequence
        '''
        yield from super().sequences()
        for r in self.singlechain():
            yield (r.get(SEQUENCE) or "").strip().upper()
