# dataset licenses

plug's git repo ships **download scripts, not data** — `data/*` is gitignored except the `.sh` files.
so publishing this repo redistributes nothing but pointers to public sources, which is unrestricted.

these licenses matter the moment you redistribute the *fetched data itself* — a cached csv/fasta under
`data/`, the `data/formatted/*.fasta` bundle, or a released benchmark tarball. each source keeps its own
terms; the full texts are in this directory.

## sources

| dataset | script | source | license | text | commercial? |
|---|---|---|---|---|---|
| uniref90 (seq reservoir) | — (`UNIREF_FASTA`) | UniProt | CC BY 4.0 | [CC-BY-4.0.txt](CC-BY-4.0.txt) | yes |
| pdb (struct reservoir) | `download_pdb.sh` | RCSB / wwPDB | CC0 1.0 | [CC0-1.0.txt](CC0-1.0.txt) | yes |
| megascale (thermompnn) | `download_megascale.sh` | github Kuhlman-Lab/ThermoMPNN | MIT | [MIT-ThermoMPNN.txt](MIT-ThermoMPNN.txt) | yes |
| fireprot (thermompnn) | `download_fireprot.sh` | github Kuhlman-Lab/ThermoMPNN | MIT | [MIT-ThermoMPNN.txt](MIT-ThermoMPNN.txt) | yes |
| dms (proteingym) | `download_proteingym.sh` | github OATML-Markslab/ProteinGym | MIT (data: cite source assays) | [MIT-ProteinGym.txt](MIT-ProteinGym.txt) | yes |
| cafa5 (via bioreason-pro) | `download_cafa5.sh` | hf wanglab/bioreason-pro-test-data | Apache 2.0 | [Apache-2.0.txt](Apache-2.0.txt) | yes |
| allobench | `download_allobench.sh` | github djmaity/allobench | MIT (data ex-ASD, CC BY) | [MIT-allobench.txt](MIT-allobench.txt) | yes |
| ppi | `download_ppi.sh` | figshare 21591618 (Bernett) | CC BY 4.0 | [CC-BY-4.0.txt](CC-BY-4.0.txt) | yes |
| shs27k | `download_shs27k.sh` | github lvguofeng/GNN_PPI (STRING v10.5 subset) | CC BY 4.0 (see note) | [CC-BY-4.0.txt](CC-BY-4.0.txt) | yes |
| passerrank | `download_passerrank.sh` | ASD table + UniProt seqs | CC BY 4.0 (see note) | [CC-BY-4.0.txt](CC-BY-4.0.txt) | yes |
| bindingdb | `download_bindb.sh` | bindingdb.org (Articles subset) | CC BY 4.0 | [CC-BY-4.0.txt](CC-BY-4.0.txt) | yes |
| **lp-pdbbind** | `download_lp_pdbbind.sh` | github THGLab/LP-PDBBind | **UC Regents — research / non-profit only** | [LP-PDBBind-LICENSE.txt](LP-PDBBind-LICENSE.txt) | no |
| **atlas** | `download_atlas.sh` | dsimb.inserm.fr | **CC BY-NC 4.0** | [CC-BY-NC-4.0.txt](CC-BY-NC-4.0.txt) | no |

## watch-outs

**two sources are non-commercial** — everything else is freely redistributable with attribution:

- **atlas** — CC BY-NC 4.0. commercial use needs the authors' permission.
- **lp-pdbbind** — UC Regents license grants use/copy/modify/**distribute** for *educational, research, and
  not-for-profit* purposes only. the affinity values ultimately come from **PDBbind**, whose own license
  forbids redistribution; lp-pdbbind is the republication vehicle, so stay within its non-commercial bound
  and cite it. this is the softest spot in the set.

if you ever bundle *all* the formatted test fastas into one release, these two NC terms cover the whole
bundle. keeping them download-on-demand (as the scripts already do) avoids that.

**the generated trainset is clean.** it derives only from uniref90 (CC BY 4.0) + pdb (CC0), so it's freely
redistributable — commercial included — with a UniProt attribution line.

**minor:**
- **passerrank** — the github repo has no LICENSE file, but the script only pulls the ASD table (ASD v3+ is
  CC BY) and fetches sequences straight from UniProt (CC BY). the *content* is CC BY on both sides.
- **shs27k** — redistributed via the GNN-PPI repo, but the *content* is a subset of STRING (interactions +
  sequences), which is CC BY 4.0. cite GNN-PPI and PIPR for the subset, STRING for the underlying data.
- **bindingdb** — the *Articles* subset used here is bindingdb's own literature curation (CC BY 4.0). the
  more restrictive ChEMBL terms apply only to the separate ChEMBL download, which this repo doesn't touch.

## cite the sources

- **uniprot / uniref90** — The UniProt Consortium, *UniProt: the universal protein knowledgebase in 2023*, NAR.
- **pdb** — Berman et al., *The Protein Data Bank*, NAR 2000.
- **thermompnn** — Dieckhaus et al., *ThermoMPNN*, PNAS 2024 (megascale: Tsuboyama et al., Nature 2023; FireProtDB).
- **proteingym** — Notin et al., *ProteinGym*, NeurIPS 2023 (+ each assay's source paper).
- **cafa5 / bioreason-pro** — wanglab bioreason-pro (arXiv:2505.23579); GO annotations from CAFA5 / UniProt-GOA.
- **allobench** — Maity et al., *AlloBench*, ACS Omega 2025 (data from ASD).
- **ppi** — Bernett et al., figshare 10.6084/m9.figshare.21591618; interactions from HIPPIE v2.3.
- **shs27k** — Lv et al., *GNN-PPI*, IJCAI 2021 (subset recipe: Chen et al., *PIPR*, Bioinformatics 2019); interactions + sequences from STRING v10.5 (Szklarczyk et al., NAR).
- **passerrank** — Xiao et al., *PASSerRank*; allosteric sites from ASD.
- **bindingdb** — Gilson et al., *BindingDB*, NAR.
- **lp-pdbbind** — Li et al., *Leak Proof PDBBind*, 2023 (PDBbind: Wang et al.).
- **atlas** — Vander Meersche et al., *ATLAS*, NAR 2024.
