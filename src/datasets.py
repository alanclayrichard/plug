# torch datasets: the leakage-free training set + the functional eval benchmarks
# (ddg, dms, go, allobench, ppi, passerrank, pdbbind, bindingdb, atlas).
import ast, csv, sys
from pathlib import Path
import pandas as pd
from torch.utils.data import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as c

csv.field_size_limit(sys.maxsize)
repo = c.repo


class TrainSet(Dataset):
    # the generated leakage-free training set (make it with build_trainset.py).
    # item: sequence, id. it's already sampled to size, so we load it in memory.
    def __init__(self, path=c.trainset):
        self.items = [{"sequence": s, "id": a} for a, s in c.iter_fasta(path)]

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class DdgDataset(Dataset):
    # thermompnn ddG test set. item: sequence, mutation, ddg, pdb, dataset.
    # dataset is megascale|fireprot (pass to filter); sign follows thermompnn.
    files = {
        "megascale": "testing/mega_test.csv",
        "fireprot":  "testing/fireprot_test.csv",
    }

    def __init__(self, dataset=None, data_dir=repo / "data" / "ThermoMPNN" / "data_all"):
        self.items = []
        for ds, rel in self.files.items():
            if dataset and ds != dataset: continue
            self.items += list(self._load(Path(data_dir) / rel, ds))

    def _load(self, path, ds):
        if not path.exists(): return
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                it = self._megascale(row) if ds == "megascale" else self._fireprot(row)
                if it: it["dataset"] = ds; yield it

    @staticmethod
    def _megascale(r):
        mut, ddg = r.get("mut_type", "").strip(), r.get("ddG_ML", "").strip()
        if not mut or mut == "wt" or ddg in ("", "-"): return None
        if any(x in mut for x in ("ins", "del", ":")): return None  # single subs only
        return {"sequence": r["wt_seq"].strip().upper(), "mutation": mut,
                "ddg": -float(ddg), "pdb": r.get("WT_name", "").replace(".pdb", "")}

    @staticmethod
    def _fireprot(r):
        if not (r.get("ddG") or "").strip(): return None
        return {"sequence": r["pdb_sequence"].strip().upper(),
                "mutation": f"{r['wild_type']}{r['position']}{r['mutation']}",
                "ddg": float(r["ddG"]), "pdb": r.get("pdb_id_corrected", "")}

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class GoDataset(Dataset):
    # bioreason-pro GO function test set (cafa5, temporal holdout).
    # item: sequence, protein_id, organism, go_bp, go_mf, go_cc (lists of GO ids).
    def __init__(self, data_dir=repo / "data" / "bioreason_go"):
        df = pd.read_parquet(Path(data_dir) / "test" / "data" / "test-00000-of-00001.parquet")
        self.items = [{"sequence": str(r.sequence).strip().upper(),
                       "protein_id": r.protein_id, "organism": r.organism,
                       "go_bp": _terms(r.go_bp), "go_mf": _terms(r.go_mf),
                       "go_cc": _terms(r.go_cc)}
                      for r in df.itertuples()]

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


def _terms(v):  # go column: array or list-repr string, nan if absent
    if isinstance(v, str): return ast.literal_eval(v) if v.startswith("[") else []
    return list(v) if hasattr(v, "__len__") else []


def _asd_site(s):  # "Chain A:PRO150,GLN151; Chain B:ASP6" -> ["A:PRO150", "A:GLN151", "B:ASP6"]
    out = []
    for part in (s or "").split(";"):
        if ":" not in part: continue
        chain, resis = part.split(":", 1)
        ch = chain.replace("Chain", "").strip()
        out += [f"{ch}:{r.strip()}" for r in resis.split(",") if r.strip()]
    return out


class DmsDataset(Dataset):
    # proteingym dms data. item: sequence, mutation, score, assay, dms_id, uniprot.
    # pass assay (Activity|Binding|Expression|OrganismalFitness|Stability) or dms_id to filter.
    def __init__(self, assay=None, dms_id=None, data_dir=repo / "data" / "proteingym"):
        data_dir = Path(data_dir)
        self.items = []
        with open(data_dir / "DMS_substitutions.csv", newline="") as f:
            meta = {r["DMS_id"]: (r["target_seq"].strip().upper(),
                                  r["coarse_selection_type"], r["UniProt_ID"])
                    for r in csv.DictReader(f)}
        for did, (seq, atype, uni) in meta.items():
            if (assay and atype != assay) or (dms_id and did != dms_id): continue
            path = data_dir / "DMS_ProteinGym_substitutions" / f"{did}.csv"
            if not path.exists(): continue
            with open(path, newline="") as f:
                for r in csv.DictReader(f):
                    if not (r.get("DMS_score") or "").strip(): continue
                    self.items.append({"sequence": seq, "mutation": r["mutant"],
                        "score": float(r["DMS_score"]), "assay": atype,
                        "dms_id": did, "uniprot": uni})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class AllobenchDataset(Dataset):
    # allobench allosteric/active-site proteins. item: sequence, target_id, gene, organism,
    # uniprot, allosteric_site, active_site (residue lists).
    def __init__(self, data_dir=repo / "data" / "allobench"):
        self.items = []
        with open(Path(data_dir) / "AlloBench.csv", newline="") as f:
            for r in csv.DictReader(f):
                s = (r.get("sequence") or "").strip().upper()
                if not s: continue
                self.items.append({"sequence": s, "target_id": r.get("target_id", ""),
                    "gene": r.get("target_gene", ""), "organism": r.get("organism", ""),
                    "uniprot": r.get("pdb_uniprot", ""),
                    "allosteric_site": _terms(r.get("allosteric_site_residue", "")),
                    "active_site": _terms(r.get("active_site_residue", ""))})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class PpiDataset(Dataset):
    # human ppi gold standard. item: sequence_a, sequence_b, id_a, id_b, label, level.
    # label 1=interacting 0=non; level 0/1/2 = redundancy-reduction stringency (pass level to filter).
    def __init__(self, level=None, data_dir=repo / "data" / "ppi"):
        data_dir = Path(data_dir)
        seqs = dict(c.iter_fasta(data_dir / "human_swissprot_oneliner.fasta"))
        self.items = []
        for lvl in (0, 1, 2):
            if level is not None and lvl != level: continue
            for label, pol in ((1, "pos"), (0, "neg")):
                path = data_dir / f"Intra{lvl}_{pol}_rr.txt"
                if not path.exists(): continue
                for line in open(path):
                    p = line.split()
                    if len(p) >= 2 and p[0] in seqs and p[1] in seqs:
                        self.items.append({"sequence_a": seqs[p[0]], "sequence_b": seqs[p[1]],
                            "id_a": p[0], "id_b": p[1], "label": label, "level": lvl})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class PasserrankDataset(Dataset):
    # passerrank allosteric proteins (asd). item: sequence, uniprot, gene, organism, pdb, allosteric_site.
    def __init__(self, data_dir=repo / "data" / "passerrank"):
        data_dir = Path(data_dir)
        seqs = {h.split("|")[1] if "|" in h else h: s
                for h, s in c.iter_fasta(data_dir / "passerrank.fasta")}
        self.items = []
        with open(data_dir / "ASD_Release_201909_AS.txt", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                acc = (r.get("pdb_uniprot") or "").strip()
                if acc not in seqs: continue
                self.items.append({"sequence": seqs[acc], "uniprot": acc,
                    "gene": r.get("target_gene", ""), "organism": r.get("organism", ""),
                    "pdb": r.get("allosteric_pdb", ""),
                    "allosteric_site": _asd_site(r.get("allosteric_site_residue", ""))})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class PdbbindDataset(Dataset):
    # LP-PDBBind protein-ligand binding affinity (leak-proof split).
    # item: sequence, pdb, value (pKd/pKi affinity), smiles, split.
    # pass split (train|val|test) to filter; default test (the eval split).
    def __init__(self, split="test", data_dir=repo / "data" / "pdbbind"):
        self.items = []
        with open(Path(data_dir) / "LP_PDBBind.csv", newline="") as f:
            for r in csv.DictReader(f):
                if split and r.get("new_split") != split: continue
                s = (r.get("seq") or "").strip().upper()
                v = (r.get("value") or "").strip()
                if not s or not v: continue
                self.items.append({"sequence": s, "pdb": r.get("", ""),
                    "value": float(v), "smiles": (r.get("smiles") or "").strip(),
                    "split": r.get("new_split", "")})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


class AtlasDataset(Dataset):
    # atlas MD per-residue flexibility (conformational dynamics), one PDB chain per simulation.
    # item: sequence, pdb_chain, rmsf (mean over 3 replicas, Å), rmsf_replicas (the 3 per-residue lists).
    def __init__(self, data_dir=repo / "data" / "atlas"):
        self.items = []
        for p in sorted((Path(data_dir) / "rmsf").glob("*.tsv")):
            seq, reps = [], ([], [], [])
            with open(p, newline="") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    aa = (r.get("seq") or "").strip()
                    if not aa: continue
                    seq.append(aa)
                    for k in range(3): reps[k].append(float(r[f"RMSF_R{k+1}"]))
            if not seq: continue
            self.items.append({"sequence": "".join(seq).upper(), "pdb_chain": p.stem,
                "rmsf": [sum(v) / 3 for v in zip(*reps)], "rmsf_replicas": [list(v) for v in reps]})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]


def _affinity(v):  # BindingDB affinity cell -> (relation, nM float); None if unparseable
    v = (v or "").strip()
    rel = "="
    while v and v[0] in "<>=~":
        rel, v = v[0], v[1:].strip()
    try: return rel, float(v)
    except ValueError: return None


class BindingdbDataset(Dataset):
    # BindingDB curated protein-ligand binding affinities (articles subset), single-chain only.
    # item: sequence, uniprot, target_name, organism, smiles, measure, value (nM), relation.
    # one affinity per row by priority IC50 > Ki > Kd > EC50; rows w/o a parseable affinity skipped.
    _aff = [("IC50", "IC50 (nM)"), ("Ki", "Ki (nM)"), ("Kd", "Kd (nM)"), ("EC50", "EC50 (nM)")]

    def __init__(self, data_dir=repo / "data" / "bindingdb"):
        self.items = []
        chains = "Number of Protein Chains in Target (>1 implies a multichain complex)"
        with open(Path(data_dir) / "BindingDB_BindingDB_Articles.tsv", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                if (r.get(chains) or "").strip() != "1": continue
                s = (r.get("BindingDB Target Chain Sequence 1") or "").strip().upper()
                if not s: continue
                hit = None
                for m, col in self._aff:
                    p = _affinity(r.get(col))
                    if p: hit = (m, p); break
                if not hit: continue
                measure, (relation, value) = hit
                self.items.append({"sequence": s,
                    "uniprot": (r.get("UniProt (SwissProt) Primary ID of Target Chain 1") or "").strip(),
                    "target_name": (r.get("Target Name") or "").strip(),
                    "organism": (r.get("Target Source Organism According to Curator or DataSource") or "").strip(),
                    "smiles": (r.get("Ligand SMILES") or "").strip(),
                    "measure": measure, "value": value, "relation": relation})

    def __len__(self): return len(self.items)

    def __getitem__(self, i): return self.items[i]
