# PLUG

**Protein Leakage-free evaluation for Unbiased Generators of function.**

![PLUG overview](plug.png)

build a **leakage-free** training subset — sampled from a sequence reservoir (uniref90), a structure reservoir (pdb), or both — to align/train a protein model on, then test whether alignment improves prediction on mutliple functional benchmarks (megascale + fireprot ddG, proteingym dms fitness, bioreason-pro GO function, allobench + passerrank allosteric sites, human PPI + string, lp-pdbbind/bindingdb affinity, atlas MD per-residue RMSF, and whatever else you want).

default leakage rule: a sampled candidate is dropped if it looks like any **test** protein —
- **sequence** (`RESERVOIR=seq`, uniref90): mmseqs2 alignment over >80% of a test seq's length at >20% identity (`COV`/`MIN_ID`), so the test protein's content is present in the train seq. Run bidirectionally (union of mode 1 and 2) so that subdomains are found.
- **structure** (`RESERVOIR=struct`, pdb): foldseek TMalign at >0.5 TM-score to a test structure (`TM`), i.e. the same fold.
- **both** (`RESERVOIR=both`): a union graph over the sequence and structure hits — a candidate is dropped if its component touches a test protein in *either* view, so a seq+struct model is safe against leakage through either lens.

thresholds/params are deliberately user-defined so you can characterize generalizability as a function of how much leakage you allow.

## leakage coverage modes (`COV_MODE`)

mmseqs `--cov-mode` decides which sequence the `COV` (0.8) threshold applies to, which matters for **multidomain** proteins. for example: a 600 aa uniref protein whose 200 aa domain B equals a 200 aa test protein — the alignment covers 100% of the test but only 33% of the train seq, so `COV_MODE=1` flags it (removed) while `COV_MODE=0` keeps it.

- `COV_MODE=union` (default) — run mode 1 **and** mode 2 and drop a candidate that either one flags. a subdomain match is leakage whichever protein is the longer one, so the honest rule takes both. on a 1M uniref sample this flags 25.2% of candidates, where mode 0 alone misses 44% of that.
- `COV_MODE=1` — cover the **test** seq. drops a train seq if it covers >80% of a test protein *regardless of the train protein's length*, so it catches a long multidomain uniref protein that contains a whole test protein as one domain.
- `COV_MODE=2` — cover the **train** seq. catches the reverse (a short train protein that is a sub-domain of a longer test protein).
- `COV_MODE=0` — cover **both** seqs (>80% of each). only removes near-duplicates of similar length; it would **keep** the multidomain protein above (the test domain is <80% of the long train seq), letting a test protein leak in as a sub-domain.

this setting applies to the **trainset scan**. the homology splits always use the union, so a benchmark's splits don't change meaning when you retune the reservoir filter.

## usage pipeline
to create your own leakage-free unlabeled training set:

```
data/download_*.sh          download benchmark(s) (all of its splits) -> data/<bench>/
plug/benchmarks/<name>.py   one file per benchmark: how to read it    -> load as a torch dataset
plug/fastas.py              each benchmark -> uniform fasta           -> data/formatted/<bench>.fasta
plug/homology.py            re-split a benchmark by similarity        -> data/splits/<bench>__<thresholds>.csv
plug/build_unlabeled_trainset.py   sample a reservoir, drop leakers   -> trainset.fasta
```

instead of scanning a whole reservoir (~121M uniref90 seqs, or the whole pdb), the trainset builder **reservoir-samples** a batch of candidates in one streaming pass so that the reservoir is never held in memory. it checks *every* sampled candidate against the combined benchmarks (~310k sequences, all splits), keeps the clean ones, and resamples until the quota is met:

- **sequence scan** (`RESERVOIR=seq`): stream the uniref90 fasta, sample a batch, `mmseqs` them against the combined test fastas, and drop any that clear the identity + coverage cutoff.
- **structure scan** (`RESERVOIR=struct`): stream the pdb directory, sample a batch of structures, `foldseek` TMalign them against the test structures, and drop any that clear the TM-score cutoff.
- **both** (`RESERVOIR=both`): sample structures, run *both* searches on them (their sequences come straight out of the structures), and drop via a union graph — any candidate whose sequence *or* structure neighborhood reaches a test protein is removed. the result is leakage-proof for a model that sees sequence and structure together.

checking a small sample against the benchmark fastas is cheap and exact; the reservoir is never copied, only the sampled training set is written. sampling is uniform (to cover lengths and families) but could — in theory — be stratified on a feature (composition, family, fold, organism, …).

## setup

the reservoirs are downloaded separately — point `.env` at them (and at the search binaries). copy `.env.example` to `.env` and edit; `config.py` auto-loads it. then:

```bash
uv sync --extra torch                    # create .venv (torch is only needed to load the benchmarks)
uv pip install -e .                      # so `import plug` works from anywhere
source .venv/bin/activate                # then plain `python ...` uses the venv
for s in data/download_*.sh; do sh "$s"; done   # fetch every benchmark
```

- **sequence reservoir**: uniref90 (https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/), point `UNIREF_FASTA` at it. needs mmseqs2 — set `MMSEQS` (or `brew install mmseqs2`).
- **structure reservoir**: a directory of pdb/mmcif files — `sh data/download_pdb.sh` mirrors the whole pdb into `data/pdb`, point `PDB_DIR` at it. needs foldseek — grab the static build for your platform from https://github.com/steineggerlab/foldseek/releases (`foldseek-osx-universal.tar.gz`, `foldseek-linux-avx2.tar.gz`, `foldseek-linux-arm64.tar.gz`, or `foldseek-linux-gpu.tar.gz`), untar it, and point `FOLDSEEK` at the extracted `bin/foldseek`. structural leakage also needs the test-set structures in `STRUCT_TESTS` (a dir of pdb/mmcif files).
- for the gpu check, point `MMSEQS`/`FOLDSEEK` at gpu builds and set `GPU=1`.

## run

```bash
python -m plug.fastas                    # all benchmarks -> data/formatted/*.fasta
RESERVOIR=seq   python -m plug.build_unlabeled_trainset   # sample uniref90, drop sequence leakers
RESERVOIR=struct python -m plug.build_unlabeled_trainset  # sample the pdb, drop TM-score leakers
RESERVOIR=both  python -m plug.build_unlabeled_trainset   # sample the pdb, drop union (seq ∪ struct) leakers
```

`.env` controls it: `RESERVOIR` picks seq/struct/both, `QUOTA` is how many leakage-free seqs to collect, `OVERSAMPLE` how much extra to draw each round to cover dropped leakers, `MIN_LEN`/`MAX_LEN` bound the sampled protein length (drops fragments + giant non-physiological seqs), `MIN_ID`/`TM` are the sequence/structure leakage thresholds, `GPU=1` runs the check on a cuda gpu.

## use

needs torch (`uv sync --extra torch`), with the venv activated:

```python
from plug.benchmarks import *
from plug.build_unlabeled_trainset import Trainset

train = Trainset.from_fasta()         # {sequence, id} — the leakage-free align set (id = uniprot code)
mega  = Megascale()                   # {sequence, mutation, ddg, pdb, split} — megascale ddG
fire  = Fireprot()                    # {sequence, mutation, ddg, pdb, split} — fireprot ddG
dms   = Dms(assay="Binding")          # {sequence, mutation, score, assay, dms_id, uniprot, split}
go    = Go()                          # {sequence, protein_id, organism, go_bp, go_mf, go_cc, split}
allo  = Allobench()                   # {sequence, target_id, gene, organism, uniprot, allosteric_site, active_site, split}
ppi   = Ppi()                         # {sequence_a, sequence_b, id_a, id_b, label, split} — human ppi
shs   = Shs27k(mode="binding")        # {sequence_a, sequence_b, id_a, id_b, types, score, split} — STRING PPI subset
pr    = Passerrank()                  # {sequence, uniprot, gene, organism, pdb, allosteric_site, split}
lppdb = LpPdbbind()                   # {sequence, pdb, value, smiles, split} — LP-PDBBind affinity
bdb   = Bindingdb()                   # {sequence, uniprot, target_name, organism, smiles, measure, value, relation, split}
atlas = Atlas()                       # {sequence, pdb_chain, rmsf, rmsf_replicas, split} — per-residue MD flexibility (Å)
```

every one of them takes `split=` and `how=`. the four that ship a split default to its **test** split; the other seven default to the **whole** benchmark (there is nothing to hold back until you split it yourself). see [benchmark splits](#benchmark-splits). `REGISTRY` maps each name to its class, for looping over all of them.

## benchmark splits

for **supervised models**, benchmarks can be broken into train/val/test splits. note that there is frequently more data available to train on in the literature. **every** benchmark takes the same two arguments. `split` picks the split you want, `how` says where it came from — published (from literature) or homology (generated here). this allows for direct comparison to literature for some benchmarks or a consistent homology splitting strategy for all benchmarks.

```python
Megascale()                     # published split, test (the default)
Megascale(split="train")        # "val", or "all" for everything
Megascale(how="homology")       # split rebuilt here, test
Atlas()                         # ships no split: the whole benchmark, as one test set
Atlas(how="published")          # ValueError — atlas never published one
Atlas(how="homology")           # so build one: test of a fresh 80/10/10
```

- `how="published"` — the split the benchmark itself shipped. four of them do (megascale, fireprot, lp-pdbbind, ppi); asking the other seven for one is a `ValueError` rather than a silently made-up split. asking for a split it doesn't publish is an error too, so `split="val"` never returns a quiet empty list.
- `how="homology"` — a fresh split built here (below). available for all eleven.
- `how="whole"` — no split at all: every item is `test`, because the whole benchmark is the test set. the default for the seven that publish no split, so a plain `Atlas()` is the entire benchmark.

`how="homology"` runs [`homology.py`](src/plug/homology.py): search the benchmark's **own** proteins against themselves with mmseqs, link any two that clear the leakage rule's bar (`MIN_ID` identity over `COV` coverage) under **either cov-mode** — mode 1 catches a hit that only covers the target, mode 2 one that only covers the query, so a protein that is a subdomain of another gets caught whichever way round it is — take the union of both graphs, then deal whole connected components 80/10/10 into train/val/test, biggest first, each to whichever split is furthest below its share. no two splits share a homolog (or a homolog of a homolog — the components are transitively closed), and every item of a protein moves together.

pair benchmarks (ppi, shs27k) pass both proteins — `("id_a", "id_b")` with `("sequence_a", "sequence_b")` — so a pair contributes both to the clustering. the fractions then apply to the *proteins*, and a pair whose two proteins land in different splits is dropped (it belongs to neither), so retained pair counts skew hard: shs27k 7,624 → 4,953/102/93, ppi 274,500 → 176,832/3,051/3,055. the published ppi split is the paper's own version of the same idea (divide proteins into three blocks, keep only within-block pairs) and is far better balanced — 163,192/59,260/52,048, 50/50 positives in every split — so prefer it unless you want the split rebuilt at this repo's stricter identity cutoff than that.


## datasets

unique sequences per benchmark — the whole thing, every published split included, which is what
`data/formatted/<bench>.fasta` holds and what the trainset is kept clear of:

| benchmark | source | unique seqs | published split |
|---|---|---|---|
| megascale | megascale ddG (tsuboyama 2023, thermompnn split) | 271,526 | train/val/test |
| fireprot | fireprot ddG (fireprotdb, thermompnn split) | 193 | train/val/test |
| dms | proteingym DMS_substitutions | 187 | — |
| go  | bioreason-pro (cafa5 temporal holdout) | 8,528 | — |
| allobench | allobench allosteric/active sites | 425 | — |
| ppi | human ppi gold standard (figshare) | 11,019 | train/val/test |
| shs27k | STRING human PPI subset (GNN-PPI / PIPR) | 1,690 | — |
| passerrank | passerrank allosteric set (ASD) | 333 | — |
| lp-pdbbind | LP-PDBBind protein–ligand affinity | 12,718 | train/val/test |
| bindingdb | bindingdb articles, protein–ligand affinity (single-chain targets) | 2,157 | — |
| atlas | atlas MD per-residue RMSF (conformational dynamics) | 1,693 | — |

**310,469 sequences** all told (308,121 unique — ~2.3k proteins appear in more than one benchmark). the three that grew when their train splits were added: megascale 28,199 → 271,526, fireprot 53 → 193, lp-pdbbind 2,644 → 12,718. ppi always shipped all three of its blocks and the other seven publish no split, so those were already whole. megascale dominates because it counts every mutant sequence as well as every wild type — 215,731 train substitutions of 239 wild-type proteins — so it is 271k sequences of only a few hundred distinct proteins, and a candidate homologous to one of them is homologous to all of its mutants anyway.

the **published split** column is what `how="published"` will hand you. a `—` means the benchmark ships no train/val/test to select from, so `how="published"` raises. several of those seven are still distributed as evaluation-only data — go is a cafa5 temporal holdout, dms and bindingdb are assay collections — and for those the whole benchmark *is* the test data, which is exactly what the default `how="whole"` gives you.

reservoirs: uniref90 (~121M sequences) and the pdb (~all deposited structures, as a directory of mmcif files).

## licenses

while this repo doesn't distribute actual data, it is important to review [`licenses/`](licenses/) for the full texts and an attribution index. all sources are freely redistributable with attribution **except atlas (CC BY-NC) and lp-pdbbind (research/non-profit only)**, which are non-commercial. the generated trainset (uniref90 CC BY + pdb CC0) is clean for any use, with attribution.

## adding a benchmark

1. add `data/download_<name>.sh` to fetch it — all of its splits, if it ships splits.
2. add `src/plug/benchmarks/<name>.py`: one class subclassing `Benchmark`, four attributes and one
   method. that's the whole contract.

```python
class Atlas(Benchmark):
    name, key = "atlas", "pdb_chain"   # key names the field holding an item's protein
    published = ()                     # ("train","val","test") if it ships a split, () if not

    @classmethod
    def rows(cls):                     # yield one dict per item, unfiltered
        yield {"sequence": ..., "pdb_chain": ...}
```

pair benchmarks set `key = ("id_a", "id_b")` and `seq = ("sequence_a", "sequence_b")`. if the benchmark ships a split, give every item its `"split"`. if it holds sequences its items don't carry (mutants, extra chains), override `sequences()` — everything it yields is held out of the trainset.
3. import it in [`benchmarks/__init__.py`](src/plug/benchmarks/__init__.py) — a benchmark missing from that list silently disappears from `REGISTRY`.
4. `python -m plug.fastas <name>`, then `python -m plug.build_unlabeled_trainset` — the sample is checked against every `data/formatted/*.fasta` (and, for `RESERVOIR=struct|both`, every structure in `STRUCT_TESTS`), so the new benchmark is protected automatically. `split=`/`how=` already work.


## TODO:
- add other sampling strategies for train set (stratify on length/family/fold/organism/…)
- add more benchmarks (more conformational/functional signals)
- fetch each benchmark's test structures automatically (if available) for the foldseek check
