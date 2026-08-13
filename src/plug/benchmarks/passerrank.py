# passerrank allosteric proteins, taken from the allosteric database (asd) release the
# paper used. one item per allosteric site entry. ships no split of its own — the raw
# asd table is what gets downloaded, not the paper's own train/test.
import csv
import re

from . import Benchmark
from .. import config as c
from ..sifts import THREE, sites


# "Chain A:PRO150,GLN151; Chain B:ASP6" -> [("A","P",150), ("A","Q",151), ("B","D",6)]
def site(s):
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
    name = "passerrank"
    # every item is about one uniprot protein
    key = "uniprot"
    published = ()

    @classmethod
    def rows(cls):
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

    # the fasta is keyed by accession, pulled out of the uniprot header
    @classmethod
    def fasta(cls):
        return {h.split("|")[1] if "|" in h else h: s
                for h, s in c.iter_fasta(c.data / "passerrank" / "passerrank.fasta")}

    # some fetched sequences have no entry in the asd table, so they are in the
    # benchmark's fasta but in none of its items. hold those out too.
    def sequences(self):
        yield from super().sequences()
        yield from self.fasta().values()
