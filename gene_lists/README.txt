# gene_lists/ — Functional Gene Lists

Generated: 2026-05-05  |  Script: 08_extract_functional_genes.py
Reproduce: python3 08_extract_functional_genes.py

## Biological question
Which orthogroups — within each metabolic context — are over-represented
(RCA ≥ 1) in liverworts but NOT over-represented in any angiosperm?
These OGs represent the part of the liverwort metabolic repertoire that
cannot be explained by genome size and is absent from angiosperm genomes.

## Three filters

### Filter A — liverwort-specific  (primary gene lists)
  Criteria: RCA=1 in ≥1 liverwort (M. polymorpha or L. cruciata)
            AND RCA=0 in ALL angiosperms (A. thaliana, S. lycopersicum, P. americana)
  Sorted: ubiquity ascending (most lineage-restricted OGs first),
          then |GCi_rca| descending.

### Filter B — both liverworts  (conservative, shared by genus pair)
  Criteria: RCA=1 in BOTH M. polymorpha AND L. cruciata
            AND RCA=0 in ALL angiosperms
  Interpretation: shared by two liverwort genera; more likely ancestral
  liverwort character than species-specific expansion.

### Filter C — bryophyte-conserved  (pre-angiosperm)
  Criteria: RCA=1 in ≥1 liverwort AND RCA=1 in P. patens (moss)
            AND RCA=0 in ALL angiosperms
  Interpretation: present in liverworts AND mosses but lost in angiosperms.
  Suggests the over-representation predates the liverwort/moss split
  and was selectively lost in vascular plant lineages.

## Note on GCi sign convention
  In functional subsets the GCi sign can invert relative to the whole-genome
  analysis. OGs sorted by ubiquity (ascending) are the most biologically
  informative first-pass list. GCi_rca is reported as auxiliary information.
  See Methods note in script docstring for full explanation.

## Directory structure

  gene_lists/
    00_cross_context_summary.tsv         OG and gene counts per context and filter
    C_metabolism/
      01_OGs_liverwort_specific.tsv      Filter A OG table
      02_OGs_both_liverworts.tsv         Filter B OG table
      03_OGs_bryophyte_conserved.tsv     Filter C OG table
      04_genes_by_species.tsv            Flat gene table (filter A, all species)
      05_genes_Mpolymorpha.txt           Marchantia gene IDs
      05_genes_Lcruciata.txt             Lunularia gene IDs
      05_genes_Ppatens.txt               Physcomitrium gene IDs
      05_genes_Creinhardtii.txt          Chlamydomonas gene IDs
      05_genes_Afiliculoides.txt         Azolla gene IDs
      05_genes_Scucullata.txt            Salvinia gene IDs
      06_summary.txt                     Summary statistics
    N_metabolism/    (same structure)
    primary_metabolism/    (same structure)

## Gene IDs
  IDs match the source proteome files in data/proteomes/filtered/.
  For KEGG/GO enrichment of Marchantia genes, use:
    gene_lists/<context>/05_genes_Mpolymorpha.txt
  These IDs follow the Ensembl Plants r62 naming (MpXXXXXX format).
  Back-reference to orthogroups via 04_genes_by_species.tsv (column: orthogroup).
