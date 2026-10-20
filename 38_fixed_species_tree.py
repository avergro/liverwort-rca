#!/usr/bin/env python3
"""
Gene-tree outputs reconciled with the accepted species tree.

The first OrthoFinder run inferred its own species tree, which placed
Selaginella as sister to all other land plants. OrthoFinder was rerun from the
same DIAMOND results with the accepted rooted species tree
(data/species_tree/species_tree_fixed.nwk: bryophytes sister to tracheophytes,
Chlamydomonas as outgroup):

    orthofinder -b data/orthofinder_v3/Results_Oct05/WorkingDirectory \\
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
"""

import csv
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OLD = ROOT / "data/orthofinder_v3/Results_Oct05"
NEW = OLD / "WorkingDirectory/OrthoFinder/Results_fixedtree"
OUT = ROOT / "data/orthofinder_v3/Results_Oct05_fixedtree"


def load(path):
    out = {}
    with open(path) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if r[0] == "Orthogroup":
                continue
            out[r[0]] = frozenset(g.strip() for cell in r[1:] for g in cell.split(",") if g.strip())
    return out


def mapping():
    old, new = load(OLD / "Orthogroups/Orthogroups.tsv"), load(NEW / "Orthogroups/Orthogroups.tsv")
    by_set = {v: k for k, v in old.items()}
    gene_old = {g: k for k, v in old.items() for g in v}
    m, inexact = {}, 0
    for k, v in new.items():
        if v in by_set:
            m[k] = by_set[v]
            continue
        hits = {}
        for g in v:
            if g in gene_old:
                hits[gene_old[g]] = hits.get(gene_old[g], 0) + 1
        m[k] = max(hits, key=hits.get)
        inexact += 1
    assert len(set(m.values())) == len(m), "orthogroup mapping is not one to one"
    return m, inexact


def main():
    if not (NEW / "Orthogroups/Orthogroups.tsv").exists():
        raise SystemExit("the OrthoFinder run with the fixed species tree is not present; its remapped outputs are "
                         f"provided in {OUT.relative_to(ROOT)}, and the command that produced it is given above")
    m, inexact = mapping()
    print(f"{len(m)} orthogroups mapped, {inexact} by largest overlap")
    pat = re.compile(r"\bOG\d{7}\b")
    sub = lambda text: pat.sub(lambda x: m[x.group(0)], text)
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "Gene_Duplication_Events").mkdir(parents=True)
    (OUT / "Phylogenetic_Hierarchical_Orthogroups").mkdir()
    shutil.copytree(NEW / "Species_Tree", OUT / "Species_Tree")
    for f in ["Duplications.tsv"]:
        (OUT / "Gene_Duplication_Events" / f).write_text(sub((NEW / "Gene_Duplication_Events" / f).read_text()))
    for f in sorted((NEW / "Phylogenetic_Hierarchical_Orthogroups").glob("N*.tsv")):
        (OUT / "Phylogenetic_Hierarchical_Orthogroups" / f.name).write_text(sub(f.read_text()))
    (OUT / "README.txt").write_text(__doc__)
    print("written", OUT)


if __name__ == "__main__":
    main()
