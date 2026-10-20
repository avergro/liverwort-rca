#!/usr/bin/env python3
"""KEGG pathway summary of the Filter B orthogroups (Fig. 6a, 6b).

Reconstructs gene_lists/<ctx>/pathway_summary_FilterB.tsv, which in the original
project was produced interactively: each Filter B orthogroup takes the union of the
KEGG pathways of its Arabidopsis genes; the table counts orthogroups per pathway and
their mean ubiquity.
Usage: python3 29_filterB_pathways.py [--root DIR] [--out-root DIR]
"""
import argparse
import re
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
ap.add_argument("--out-root", type=Path, default=None)
a = ap.parse_args()
root, out_root = a.root, (a.out_root or a.root)

tsv = sorted((root / "data/orthofinder_v3").glob("Results_*"))[-1] / "Orthogroups" / "Orthogroups.tsv"
og = pd.read_csv(tsv, sep="\t", index_col=0, usecols=["Orthogroup", "Athaliana"])
ath = {o: {re.sub(r"\.\d+$", "", g.strip()) for g in str(c).split(",")} for o, c in og.Athaliana.dropna().items()}
gp = pd.read_csv(root / "data/kegg_cache/ath_gene_pathway.tsv", sep="\t", header=None, names=["gene", "path"])
gp["gene"] = gp.gene.str.replace("ath:", "", regex=False)
gp["path"] = gp.path.str.replace("path:", "", regex=False)
gene2path = gp.groupby("gene").path.apply(set).to_dict()
names = pd.read_csv(Path(__file__).resolve().parent / "data/kegg_cache/ath_pathway_names.tsv", sep="\t",
                    header=None, names=["path", "name"])
names["name"] = names.name.str.replace(r" - Arabidopsis thaliana \(thale cress\)$", "", regex=True)
name = dict(zip(names.path, names.name))

for ctx in ("primary_metabolism", "C_metabolism", "N_metabolism"):
    b = pd.read_csv(root / "gene_lists" / ctx / "02_OGs_both_liverworts.tsv", sep="\t", index_col=0)
    rows = []
    for o, r in b.iterrows():
        paths = set().union(*[gene2path.get(g, set()) for g in ath.get(o, ())])
        rows += [(p, r.ubiquity) for p in paths]
    d = pd.DataFrame(rows, columns=["pathway", "ubiquity"])
    d = d[~d.pathway.isin(["ath01100", "ath01110"])]   # global maps, excluded as in the original tables
    s = (d.groupby("pathway").ubiquity.agg(n_OGs="size", mean_ubiquity="mean").reset_index()
         .sort_values(["n_OGs", "pathway"], ascending=[False, True]))
    s.insert(1, "pathway_name", s.pathway.map(name))
    out = out_root / "gene_lists" / ctx / "pathway_summary_FilterB.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(out, sep="\t", index=False)
    print(ctx, len(b), "Filter B OGs;", len(s), "pathways ->", out)
