
Gene-tree outputs reconciled with the accepted species tree.

The first OrthoFinder run inferred its own species tree, which placed
Selaginella as sister to all other land plants. OrthoFinder was rerun from the
same DIAMOND results with the accepted rooted species tree
(data/species_tree/species_tree_fixed.nwk: bryophytes sister to tracheophytes,
Chlamydomonas as outgroup):

    orthofinder -b data/orthofinder_v3/Results_Oct05/WorkingDirectory \
        -s data/species_tree/species_tree_fixed.nwk -M msa -T fasttree -t 12 -a 12 -n fixedtree

The rerun reproduces the orthogroups (18,879 of 18,902 with identical gene
sets; the other 23 differ by 19 genes in total) but numbers them differently.
The analyses keep the orthogroups of the first run, so this script copies the
outputs that depend on the species tree (gene duplication events and
hierarchical orthogroups) into data/orthofinder_v3/Results_Oct05_fixedtree/
with every orthogroup identifier translated to the first run, by identical gene
set or, for the 23 that differ, by the largest overlap.

Node labels of the fixed tree: N1 land plants, N2 bryophytes, N3 tracheophytes,
N4 liverworts, N5 euphyllophytes, N6 ferns, N7 seed plants, N8 angiosperms.

Usage:
    python3 38_fixed_species_tree.py
