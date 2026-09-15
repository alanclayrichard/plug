# renumber pdb residue numbers into uniprot sequence positions
'''
The allosteric benchmarks number their site residues the way the pdb file does
and everything else in plug counts along the uniprot sequence, so the two have
to be lined up before a site can be trusted. Sifts publishes that mapping as
one row per observed stretch of chain, which is what segments() reads: keeping
the stretches separate rather than storing a single offset is what makes gaps
in a structure stay aligned. A position that converts is then checked against
the sequence it lands on, so every residue plug ships is one whose identity
matches the residue the benchmark named.
'''

import csv
import gzip
import re
from functools import lru_cache

from . import config as c

THREE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
         "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
         "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def number(v):
    '''
    pull the number off a pdb residue number, which can carry an insertion
    letter on the end ("27A").

    inputs:
        v: str -- the raw cell
    outputs:
        n: int | None -- the number, or None if there isn't one
    '''
    m = re.match(r"^(-?\d+)", (v or "").strip())
    return int(m.group(1)) if m else None

@lru_cache(maxsize=1)
def segments():
    '''
    read the sifts table into the stretches it maps. one entry per run of
    residues rather than one offset per chain, so gaps in the structure stay
    lined up. cached, since the table is large and never changes.

    inputs:
        none
    outputs:
        segs: dict[tuple[str, str, str], list[tuple[int, int, int]]] --
              (pdb, chain, accession) -> [(first pdb number, last, shift), ...]
    '''
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

def to_uniprot(pdb, chain, position, accession):
    '''
    move one pdb residue number onto the uniprot sequence.

    inputs:
        pdb: str -- pdb code
        chain: str -- chain id
        position: int -- the residue number as the pdb file writes it
        accession: str -- uniprot accession the chain maps to
    outputs:
        pos: int | None -- the uniprot position, or None if sifts doesn't
             cover that chain
    '''
    for beg, end, shift in segments().get((pdb.lower(), chain, accession), ()):
        if beg <= position <= end: return position + shift
    return None

def sites(labels, pdb, accession, sequence):
    '''
    renumber a set of site residues and keep only the ones that convert and
    then land on the residue they name. every position plug ships is checked
    against the sequence plug ships, which is the same promise the mutation
    benchmarks make.

    inputs:
        labels: list[tuple[str, str, int]] -- (chain, residue letter, pdb number)
        pdb: str -- pdb code
        accession: str -- uniprot accession
        sequence: str -- the uniprot sequence being numbered against
    outputs:
        positions: list[int] -- sorted uniprot positions that checked out
    '''
    out = set()
    for chain, residue, position in labels:
        u = to_uniprot(pdb, chain, position, accession)
        if u and 1 <= u <= len(sequence) and sequence[u - 1] == residue: out.add(u)
    return sorted(out)
