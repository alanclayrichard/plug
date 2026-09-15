# allobench allosteric and active sites, one item per entry
'''
Proteins with a known allosteric site and a known active site. It ships no
split of its own. The site residues arrive numbered the way the pdb file writes
them, as labels like "B-THR-7", so they are renumbered against the uniprot
sequence through sifts and then checked against the residue they name; a
position that does not convert or does not match is dropped rather than shipped
wrong. Scoring is one auroc per protein over its residues, so what tells one
protein from another for the stratified split is how much of it is site, which
is why strata() summarises the residue list as a fraction before binning.
'''

import ast
import csv
import re

from . import Benchmark
from .. import config as c
from ..sifts import THREE, sites

def terms(v):
    '''
    a residue column arrives either as a list or as a string that looks like
    one.

    inputs:
        v: str | list -- the raw cell
    outputs:
        items: list -- what it holds
    '''
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []

def labels(v):
    '''
    parse site residue labels into the shape sifts renumbers, so "B-THR-7"
    becomes ("B", "T", 7).

    inputs:
        v: str | list -- the raw site residue cell
    outputs:
        labels: list[tuple[str, str, int]] -- (chain, residue letter, number)
    '''
    out = []
    for lab in terms(v):
        m = re.match(r"^([A-Za-z0-9]+)-([A-Z]{3})-(-?\d+)", str(lab).strip())
        if m and m.group(2) in THREE: out.append((m.group(1), THREE[m.group(2)], int(m.group(3))))
    return out

class Allobench(Benchmark):
    '''allosteric and active site residues, numbered against uniprot.'''

    name = "allobench"
    # every item is about one target protein
    key = "target_id"
    published = ()

    # the label is a set of site residues, and the score is one auroc per protein, so
    # what distinguishes one protein from another is how much of it is site.
    label = "allosteric_site"

    @classmethod
    def strata(cls, rows):
        '''
        bin proteins by what fraction of them is allosteric site, since the
        label is per-residue and the score is per-protein.

        inputs:
            rows: list[dict] -- the benchmark's items
        outputs:
            strata: list[tuple[str, ...]] -- the bin each protein falls in
        '''
        rows = [{**it, "allosteric_site": len(it["allosteric_site"]) / max(len(it["sequence"]), 1)}
                for it in rows]
        return Benchmark.strata.__func__(cls, rows)

    @classmethod
    def rows(cls):
        '''
        one item per entry, with both sites counted along the sequence shipped
        with it.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, ids, and the two site residue lists
        '''
        with open(c.data / "allobench" / "AlloBench.csv", newline="") as f:
            for r in csv.DictReader(f):
                s = (r.get("sequence") or "").strip().upper()
                if not s: continue
                pdb, acc = r.get("allosteric_pdb", ""), r.get("pdb_uniprot", "")
                yield {"sequence": s,
                       "target_id": r.get("target_id", ""),
                       "gene": r.get("target_gene", ""),
                       "organism": r.get("organism", ""),
                       "uniprot": acc,
                       "pdb": pdb,
                       # the residues making up each site, counted along the sequence above
                       "allosteric_site": sites(labels(r.get("allosteric_site_residue", "")),
                                                pdb, acc, s),
                       "active_site": [p for p in terms(r.get("active_site_residue", ""))
                                       if isinstance(p, int) and 1 <= p <= len(s)]}
