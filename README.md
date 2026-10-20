# Liverwort genomes through economic complexity: code and data

Code and data for Vergara A, Vergara M, Gutiérrez RA, *Fewer copies, not fewer functions? Liverwort genomes through economic complexity* (manuscript in preparation).

The analysis compares the orthogroup content of eleven green plant proteomes with Revealed Comparative Advantage (RCA), decomposes gene number into orthogroup coverage and genes per occupied orthogroup, and tests the resulting ranking with a set of controls (bootstrap, species removal, rarefaction, a degree-preserving null model, phylogenetic regression and simulation). The repository mirrors the working layout, so every script runs unchanged from its root.

In the scripts and result files the RCA score is called **MCi-RCA** (or `MCi_rca`) and the presence score **MCi-raw** (`MCi_raw`).

## Contents

| Path | Content |
|---|---|
| `01`–`35` scripts | Pipeline, in the order of execution (table below) |
| `plots/` | One script per figure (`fig01`–`fig06`, `figS1`–`figS5`) and a plain style module (see `plots/README.md`) |
| `Matias/alexander_metodologia/scripts/06_go_filter_multispecies.py` | Gene Ontology filter on the native annotation of each species |
| `data/orthofinder_v3/Results_Oct05/` | OrthoFinder 3.1.4 output read by the analyses: orthogroups, gene counts, statistics, gene duplication events, species tree and resolved gene trees (large files gzip-compressed) |
| `data/proteomes/gene_ids/` | Identifiers of the genes analysed for each species (one longest protein per gene, after the corrections below) |
| `data/proteomes/corrections/` | Gene lists used to correct the *Picea abies* and *Selaginella moellendorffii* proteomes |
| `data/busco/` | BUSCO short summaries |
| `data/filters/`, `Matias/alexander_metodologia/GO_files/` | Orthogroups assigned to each functional context under KEGG and Gene Ontology |
| `data/results_*` | Every result table behind the figures and the supplementary tables |
| `supplementary_data/` | Dataset S1: complete result tables cited in the manuscript and Supplementary Note S1, with a README describing each file |
| `gene_lists/` | Filter B orthogroups and their genes |

## Not included

These files are public and are downloaded by the scripts, or cannot be redistributed:

- Proteome FASTA files: retrieved by `01_download_proteomes.sh` (sources below).
- KEGG data: `06_kegg_filter.py` queries the KEGG REST API and writes `data/kegg_cache/`; the KOfam profiles and module definitions of KEGG release 2025-02-04 are at https://figshare.com/articles/dataset/KEGG_build_2025-02-04/28346942 and go in `data/kofam/`.
- PlantTFDB v5 transcription factor list of *Arabidopsis thaliana* (`data/tf_families/Ath_TF_list.txt`), from https://planttfdb.gao-lab.org.
- Gene Ontology annotation files of each species and `go-basic.obo` (sources in `Matias/alexander_metodologia/GO_files/download_info.txt`).
- BUSCO lineage datasets (`embryophyta_odb10`, `chlorophyta_odb10`) and OrthoFinder intermediate files (alignments, orthogroup sequences, working directory), which OrthoFinder regenerates.

## Proteomes

| Species | Source | Annotation |
|---|---|---|
| *Arabidopsis thaliana* | Ensembl Plants release 62 | TAIR10 (Araport11) |
| *Chlamydomonas reinhardtii* | Ensembl Plants release 62 | v5.5 |
| *Marchantia polymorpha* | Ensembl Plants release 62 | MpTak1 v7.1 |
| *Physcomitrium patens* | Ensembl Plants release 62 | Phypa V3 |
| *Selaginella moellendorffii* | Ensembl Plants release 62 | v1.0, allelic copies collapsed (below) |
| *Solanum lycopersicum* | Ensembl Plants release 62 | SL3.0 |
| *Azolla filiculoides* | FernBase | v1.1 high-confidence gene set |
| *Salvinia cucullata* | FernBase | v1.2 high-confidence gene set |
| *Lunularia cruciata* | CoGe, genome id 64629 | manual download, see `01_download_proteomes.sh` |
| *Persea americana* | NCBI GCA_029852735.1 | |
| *Picea abies* | V2 annotation (Picab02), provided by N. R. Street | transposable element genes removed (below) |

For proteomes with several isoforms the longest protein per gene was kept (`02_filter_primary_transcripts.py`). Three proteomes were corrected before orthogroup inference, and the resulting gene sets are listed in `data/proteomes/gene_ids/`:

- *Picea abies*: the 11,100 gene models annotated as transposable elements were removed, leaving the 43,410 protein-coding genes in `corrections/Pabies_noTE_genes.txt`.
- *Azolla filiculoides*: the FernBase v1.1 high-confidence set (20,203 genes) replaces the v1.0 set downloaded by `01_download_proteomes.sh`.
- *Selaginella moellendorffii*: the Ensembl annotation holds the gene models of both haplotypes. An all-against-all DIAMOND search grouped proteins sharing at least 95% identity over at least 80% of both lengths, and the longest protein of each group was kept (34,799 to 24,228 genes); the removed models are in `corrections/Smoellendorffii_dropped_alleles.txt`.

## Running

```bash
conda env create -f environment.yaml      # or: pip install -r requirements.txt for the Python steps only
conda activate liverwort-rca
```

Steps 01–04 rebuild the proteomes, BUSCO summaries and orthogroups and take most of the computing time. With the OrthoFinder output included here, the analyses can start at step 05.

| Step | Script | Output | Used in |
|---|---|---|---|
| 1 | `01_download_proteomes.sh`, `02_filter_primary_transcripts.py` | `data/proteomes/` | Table 1 |
| 2 | `03_busco_qc.sh`, `03b_busco_summary.py` | `data/busco/`, `data/busco/busco_summary.tsv` | Table 1, Table S9 |
| 3 | `04_run_orthofinder.sh`, then `38_fixed_species_tree.py` after an OrthoFinder rerun with the accepted species tree (`data/species_tree/species_tree_fixed.nwk`; command in the script) | `data/orthofinder_v3/Results_Oct05/`, `data/orthofinder_v3/Results_Oct05_fixedtree/` | all; duplications and hierarchical orthogroups from the second |
| 4 | `07_run_all_analyses.sh` (runs `06_kegg_filter.py` and `05_build_matrix.py` for the whole genome and each KEGG context) | `data/filters/`, `data/results_v3/`, `data/results_*_v3/` | Fig. 2, Table 2 |
| 5 | `Matias/.../06_go_filter_multispecies.py`, `31_go_rankings.py` | GO filters, `data/results_go/` | Fig. 2c, Table 2 |
| 6 | `08_extract_functional_genes.py`, `29_filterB_pathways.py`, `32_filterB_enrichment.py` | `gene_lists/`, `data/results_filterB/` | Fig. 6, Table S7 |
| 7 | `09a_bootstrap_mci.py`, `09b_fitness_complexity.py`, `09c_expansion_contraction.py` | `data/results_bootstrap*/`, `data/results_fitness_complexity/`, `data/results_expansion_contraction/` | Fig. 2b, Note S1 |
| 8 | `10_tf_family_redundancy.py`, `11_regulatory_subclasses.py`, `27_regulatory_contrast.py` | `data/results_tf_redundancy/`, `data/results_regulatory_*/` | Figs 4 and 5 |
| 9 | `12`–`21`, `24` | `data/results_controls/`, `data/results_simulation/` | Note S1, Figs S1–S5, Tables S2 and S11–S16 |
| 10 | `25_presence_coverage.py`, `26_coverage_redundancy.py`, `35_bootstrap_ratios.py` | `data/results_controls/` | Fig. 3, Table 3 |
| 11 | `28_isoforms_per_gene.py` | `data/results_isoforms/` | Methods (isoforms per gene) |
| 12 | `33_kegg_module_completeness.py`, `34_eci_kegg_modules.py` | `data/results_kegg_modules/` | Results 2.3, Table S5, Note S1 |
| 13 | `37_overlap_and_sensitivity.py` | `data/results_controls/` | Results 2.2 and 2.5, Tables S3–S5 |
| 14 | `18_make_table_s7.py`, `22_build_supplementary.py` | supplementary tables, `supplementary_data/` | Tables S1–S9, Dataset S1 |
| 15 | `plots/fig*.py` | `plots/figures/` | Figures 1–6 and S1–S5, in a plain style; the published versions were refined afterwards for presentation (`plots/README.md`) |

Every stochastic step uses a fixed random seed, 42 by default (`--seed`). The simulation parameters (60 capabilities, orthogroup requirements of size 1 + Poisson(2), gene copies 1 + Poisson(d) with a species-specific duplication rate d, 11 and 40 species) are the defaults of `simulate()` in `16_simulation_mci.py`, and the four relaxations are defined in `17_simulation_extensions.py`.

Every external tool (busco, orthofinder, diamond, famsa, fasttree, hmmsearch, hmmfetch) is called from `PATH`, which the conda environment provides once activated; the HMMER tools used by `33_kegg_module_completeness.py` can also be pointed to explicitly with the `HMMER_BIN` environment variable.

## Citation and licence

Please cite the article once published; the archived DOI (Zenodo) will be added at release. The scripts are released under the MIT licence (`LICENSE`) and the data files under CC BY 4.0 (`LICENSE-DATA`).
