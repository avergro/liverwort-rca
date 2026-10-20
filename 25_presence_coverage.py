#!/usr/bin/env python3
"""
Presence-based repertoire breadth, requested by review.

MCi-RCA and the diversity statistic both count orthogroups above an RCA
threshold, which is a share and therefore falls as the proteome grows. The
claim that liverworts hold a broad repertoire is a presence claim, so it needs
a presence statistic: how many orthogroups of each functional context does a
species occupy at all, with at least one gene.

Reported per species for the whole genome, for the four KEGG contexts and for
the four Gene Ontology contexts, as counts and as the fraction of the context
occupied.

Usage: python3 25_presence_coverage.py
Output: data/results_controls/K_presence_coverage.tsv
"""

import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"
GO_DIR = ROOT / "Matias/alexander_metodologia/GO_files"
KEGG_DIR = ROOT / "data/filters"
CONTEXTS = {
    "whole_genome": None,
    "KEGG primary": KEGG_DIR / "primary_metabolism_ogs.txt",
    "KEGG carbon": KEGG_DIR / "C_metabolism_ogs.txt",
    "KEGG nitrogen": KEGG_DIR / "N_metabolism_ogs.txt",
    "KEGG regulatory": KEGG_DIR / "regulatory_networks_ogs.txt",
    "GO primary": GO_DIR / "primary_metabolism_ogs.txt",
    "GO carbon": GO_DIR / "C_metabolism_ogs.txt",
    "GO nitrogen": GO_DIR / "N_metabolism_ogs.txt",
    "GO regulatory": GO_DIR / "regulatory_networks_ogs.txt",
}
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    bm = load_bm()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    counts = counts[(counts > 0).sum(axis=1) >= 2]          # same min_species filter
    genes = counts.sum(axis=0)

    rows = []
    for ctx, filt in CONTEXTS.items():
        sub = counts if filt is None else counts[counts.index.isin(
            set(filt.read_text().split()))]
        present = (sub > 0).sum(axis=0)
        for sp in counts.columns:
            rows.append({"context": ctx, "species": sp,
                         "orthogroups_in_context": sub.shape[0],
                         "occupied": int(present[sp]),
                         "fraction_occupied": round(present[sp] / sub.shape[0], 4),
                         "genes": int(genes[sp])})
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "K_presence_coverage.tsv", sep="\t", index=False)

    wide = df.pivot(index="species", columns="context", values="fraction_occupied")
    order = ["whole_genome", "KEGG primary", "KEGG carbon", "KEGG nitrogen",
             "KEGG regulatory", "GO primary", "GO carbon", "GO nitrogen",
             "GO regulatory"]
    wide = (wide[order] * 100).round(1)
    occ = df[df.context == "whole_genome"].set_index("species").occupied
    wide["genes"] = genes
    wide["genes_per_occupied_OG"] = (genes / occ).round(2)
    wide = wide.sort_values("genes")
    print("\nfraction of each context occupied (%, presence, >=1 gene)\n")
    print(wide.to_string())

    print("\nliverwort mean / angiosperm mean, by context\n")
    for ctx in order:
        lw = df[(df.context == ctx) & (df.species.isin(LIVERWORTS))].occupied.mean()
        an = df[(df.context == ctx) & (df.species.isin(ANGIOSPERMS))].occupied.mean()
        n = df[df.context == ctx].orthogroups_in_context.iloc[0]
        print(f"  {ctx:<18} liverworts {lw:8.0f}   angiosperms {an:8.0f}   "
              f"ratio {lw/an:5.2f}   (context has {n:,} orthogroups)")


if __name__ == "__main__":
    main()
