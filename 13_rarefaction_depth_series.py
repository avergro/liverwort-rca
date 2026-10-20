#!/usr/bin/env python3
"""
Rarefaction depth series for MCi-RCA (follow-up to 12_gene_number_controls.py).

Full equalization to the smallest proteome (Chlamydomonas, 12,470 genes assigned
to orthogroups) removes so many genes that functional contexts lose power. Here
each depth T caps every proteome at T genes (species below T keep all their
genes), so the gene-number disparity is reduced progressively instead of being
collapsed in one step. A second series repeats this with Chlamydomonas dropped,
where full equalization is 14,137 genes (Marchantia).

Usage:
    caffeinate -i python3 13_rarefaction_depth_series.py --reps 100
Outputs: data/results_controls/E_depth_series*.tsv
"""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"
CONTEXTS = {
    "whole_genome": None,
    "primary_metabolism": ROOT / "data/filters/primary_metabolism_ogs.txt",
    "C_metabolism": ROOT / "data/filters/C_metabolism_ogs.txt",
    "N_metabolism": ROOT / "data/filters/N_metabolism_ogs.txt",
    "regulatory_networks": ROOT / "data/filters/regulatory_networks_ogs.txt",
}
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
DEPTHS = [12702, 14160, 17455, 24542, 34211, None]  # None = no capping


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cap_proteomes(counts, target, rng):
    """Subsample every species with more than `target` genes down to `target`."""
    out = {}
    for sp in counts.columns:
        col = counts[sp].values.astype(int)
        total = col.sum()
        out[sp] = rng.multivariate_hypergeometric(col, target) if total > target else col
    return pd.DataFrame(out, index=counts.index)


def mci_rca(counts, bm, og):
    if og is not None:
        counts = counts[counts.index.isin(og)]
    M = bm.filter_matrix(bm.build_rca_matrix(counts), 2)
    mci, _ = bm.compute_complexity(M)
    return mci


def run(counts, bm, reps, seed, tag):
    rng = np.random.default_rng(seed)
    rows = []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        for depth in DEPTHS:
            n = 1 if depth is None else reps
            for r in range(n):
                sub = counts if depth is None else cap_proteomes(counts, depth, rng)
                mci = mci_rca(sub, bm, og)
                rank = mci.rank(ascending=False)
                for sp in counts.columns:
                    rows.append({"set": tag, "context": ctx,
                                 "depth": depth if depth else "none",
                                 "rep": r, "species": sp,
                                 "MCi": mci[sp], "rank": rank[sp]})
            print(f"    {tag} {ctx} depth={depth}: done")
    return pd.DataFrame(rows)


def summarize(df):
    s = (df.groupby(["set", "context", "depth", "species"])
           .agg(mean_rank=("rank", "mean"), sd_rank=("rank", "std"),
                frac_top2=("rank", lambda r: (r <= 2).mean()))
           .reset_index())
    lw = (df[df.species.isin(LIVERWORTS)]
            .groupby(["set", "context", "depth", "rep"])["rank"]
            .apply(lambda r: (r <= 2).all())
            .groupby(["set", "context", "depth"]).mean()
            .rename("both_liverworts_top2").reset_index())
    return s, lw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    bm = load_bm()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))

    print("depth series, all 11 species")
    a = run(counts, bm, args.reps, args.seed, "all_species")
    print("depth series, Chlamydomonas excluded")
    b = run(counts.drop(columns=["Creinhardtii"]), bm, args.reps, args.seed, "no_Chlamydomonas")

    df = pd.concat([a, b], ignore_index=True)
    df.to_csv(OUT / "E_depth_series_raw.tsv", sep="\t", index=False)
    s, lw = summarize(df)
    s.to_csv(OUT / "E_depth_series_summary.tsv", sep="\t", index=False)
    lw.to_csv(OUT / "E_depth_series_liverworts.tsv", sep="\t", index=False)
    print("\nfraction of replicates with both liverworts in the top two:")
    print(lw.pivot_table(index=["set", "context"], columns="depth",
                         values="both_liverworts_top2").to_string())
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
