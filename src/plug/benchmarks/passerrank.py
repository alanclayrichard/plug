# passerrank allosteric proteins, one item per allosteric site entry
'''
The allosteric proteins from the Allosteric Database release PASSerRank used.
What gets downloaded is the raw ASD table rather than the paper's own
train/test, so it ships no split of its own. Site residues arrive as strings
like "Chain A:PRO150,GLN151; Chain B:ASP6" and are renumbered against the
uniprot sequence through sifts, then checked against the residue they name, the
same as allobench. Sequences come from a fasta fetched from uniprot and keyed
by accession; some of those have no entry in the ASD table, so they are in the
benchmark's fasta but in none of its items, and sequences() holds them out too.
'''

import csv
import re

from . import Benchmark
from .. import config as c
from ..sifts import THREE, sites

def site(s):
    '''
    parse a site residue string into the shape sifts renumbers, so "Chain
    A:PRO150,GLN151" becomes [("A","P",150), ("A","Q",151)].

    inputs:
        s: str -- the raw site residue cell
    outputs:
        labels: list[tuple[str, str, int]] -- (chain, residue letter, number)
    '''
    out = []
    for part in (s or "").split(";"):
        if ":" not in part: continue
        chain, residues = part.split(":", 1)
        ch = chain.replace("Chain", "").strip()
        for res in residues.split(","):
            m = re.match(r"^([A-Z]{3})(-?\d+)", res.strip())
            if m and m.group(1) in THREE: out.append((ch, THREE[m.group(1)], int(m.group(2))))
    return out

class Passerrank(Benchmark):
    '''allosteric site residues from the asd, numbered against uniprot.'''

    name = "passerrank"
    # every item is about one uniprot protein
    key = "uniprot"
    published = ()

    # as allobench: a set of site residues per protein, scored one protein at a time
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
        one item per allosteric site entry whose sequence was fetched from
        uniprot, with the site counted along that sequence.

        inputs:
            none
        outputs:
            rows: iterator[dict] -- sequence, uniprot, gene, organism, pdb and
                  the site residues
        '''
        folder = c.data / "passerrank"
        seqs = cls.fasta()
        with open(folder / "ASD_Release_201909_AS.txt", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                acc = (r.get("pdb_uniprot") or "").strip()
                # skip an entry unless its sequence was fetched from uniprot
                if acc not in seqs: continue
                pdb = r.get("allosteric_pdb", "")
                yield {"sequence": seqs[acc], "uniprot": acc,
                       "gene": r.get("target_gene", ""),
                       "organism": r.get("organism", ""),
                       "pdb": pdb,
                       # the residues making up the site, counted along the sequence above
                       "allosteric_site": sites(site(r.get("allosteric_site_residue", "")),
                                                pdb, acc, seqs[acc])}

    @classmethod
    def fasta(cls):
        '''
        the fetched sequences, keyed by the accession in the uniprot header.

        inputs:
            none
        outputs:
            seqs: dict[str, str] -- accession -> sequence
        '''
        return {h.split("|")[1] if "|" in h else h: s
                for h, s in c.iter_fasta(c.data / "passerrank" / "passerrank.fasta")}

    def sequences(self):
        '''
        some fetched sequences have no entry in the asd table, so they are in
        the benchmark's fasta but in none of its items. hold those out too.

        inputs:
            none
        outputs:
            seqs: iterator[str] -- every sequence in the benchmark's fasta
        '''
        yield from super().sequences()
        yield from self.fasta().values()
