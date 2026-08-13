# the allosteric benchmarks number their site residues the way the pdb file does.
# everything else in plug counts along the uniprot sequence. sifts is the table that
# lines the two up, so this turns the first kind of number into the second.
import csv
import gzip
import re
from functools import lru_cache

from . import config as c

THREE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
         "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
         "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


# a pdb number can carry a letter on the end ("27A"), we only want the number
def number(v):
    m = re.match(r"^(-?\d+)", (v or "").strip())
    return int(m.group(1)) if m else None


# (pdb, chain, accession) -> [(first pdb number, last pdb number, shift), ...]
# one stretch per run of residues, so gaps in the structure stay lined up
@lru_cache(maxsize=1)
def segments():
    if not c.sifts.exists():
        raise FileNotFoundError(f"{c.sifts} is missing — run data/download_sifts.sh")
    out = {}
    with gzip.open(c.sifts, "rt") as f:
        next(f)                                  # a dated comment sits above the header
        for r in csv.DictReader(f):
            beg, end, sp = number(r["PDB_BEG"]), number(r["PDB_END"]), number(r["SP_BEG"])
            if beg is None or end is None or sp is None: continue
            out.setdefault((r["PDB"].lower(), r["CHAIN"], r["SP_PRIMARY"]), []).append(
                (beg, end, sp - beg))
    return out


# one pdb residue number -> where it sits in the uniprot sequence, or None if sifts
# doesn't cover that chain
def to_uniprot(pdb, chain, position, accession):
    for beg, end, shift in segments().get((pdb.lower(), chain, accession), ()):
        if beg <= position <= end: return position + shift
    return None


# labels are (chain, residue letter, pdb number). keep the ones that convert and then
# land on the residue they name, so every position we ship is checked against the
# sequence we ship — the same promise the mutation benchmarks make.
def sites(labels, pdb, accession, sequence):
    out = set()
    for chain, residue, position in labels:
        u = to_uniprot(pdb, chain, position, accession)
        if u and 1 <= u <= len(sequence) and sequence[u - 1] == residue: out.add(u)
    return sorted(out)
