#!/usr/bin/env python3
"""
RCA threshold sweep and baseline comparison for MCi.

Threshold sweep: the RCA >= 1 cut is inherited convention. Here the matrix is
binarized at several thresholds and the species ranking is recomputed, so the
dependence of the result on that choice is measurable.

Baselines: MCi is compared against simpler statistics computed on the same
matrix, to show how much the eigenvector machinery adds over each one.
  diversity        number of orthogroups with RCA >= threshold
  CA axis 1        correspondence analysis of the binary matrix
  fitness          Tacchella et al. 2012 non-linear metric
  low_ubiquity     share of a species' enriched orthogroups found in <= 3 species
  MCi_raw          eigenvector on presence/absence, without RCA

Usage:
    caffeinate -i python3 15_threshold_and_baselines.py
Outputs: data/results_controls/G_threshold_sweep.tsv, G_baselines.tsv
"""

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
THRESHOLDS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def rca_values(counts):
    tot_sp = counts.sum(axis=0)
    tot_og = counts.sum(axis=1)
    grand = counts.values.sum()
    return counts.div(tot_sp, axis=1).div(tot_og / grand, axis=0)


def ca_axis1(M):
    A = M.values.astype(float)
    P = A / A.sum()
    r, c = P.sum(1), P.sum(0)
    S = (P - np.outer(r, c)) / np.sqrt(np.outer(r, c))
    _, _, Vt = np.linalg.svd(S, full_matrices=False)
    return pd.Series(Vt[0] / np.sqrt(c), index=M.columns)


def fitness_complexity(M, iters=200):
    A = M.values.astype(float)
    F = np.ones(A.shape[1]); Q = np.ones(A.shape[0])
    for _ in range(iters):
        Fn = A.T @ Q
        Qn = 1.0 / (A @ (1.0 / np.where(F == 0, 1e-12, F)))
        F = Fn / Fn.mean(); Q = Qn / Qn.mean()
    return pd.Series(F, index=M.columns)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bm = load_bm()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))

    sweep, base = [], []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[counts.index.isin(og)] if og else counts
        rca = rca_values(sub)
        ref = None
        for thr in THRESHOLDS:
            M = bm.filter_matrix((rca >= thr).astype(int), 2)
            mci, _ = bm.compute_complexity(M)
            rank = mci.rank(ascending=False)
            if thr == 1.0:
                ref = mci
            for sp in M.columns:
                sweep.append({"context": ctx, "threshold": thr, "species": sp,
                              "MCi": mci[sp], "rank": rank[sp], "n_ogs": M.shape[0]})
        # baselines at the conventional threshold
        M = bm.filter_matrix((rca >= 1).astype(int), 2)
        mci = ref
        ub = M.sum(axis=1)
        low = M[ub <= 3].sum(axis=0) / M.sum(axis=0)
        Mraw = bm.filter_matrix(bm.build_raw_matrix(sub), 2)
        mci_raw, _ = bm.compute_complexity(Mraw)
        stats = pd.DataFrame({
            "MCi": mci, "diversity": M.sum(axis=0), "CA_axis1": ca_axis1(M),
            "fitness": fitness_complexity(M), "low_ubiquity_share": low,
            "MCi_raw": mci_raw, "genes": sub.sum(axis=0),
        })
        for col in stats.columns:
            if col == "MCi":
                continue
            base.append({"context": ctx, "baseline": col,
                         "spearman_vs_MCi": abs(stats["MCi"].corr(stats[col], method="spearman")),
                         "liverworts_top2_under_baseline":
                             bool(stats[col].rank(ascending=False).reindex(LIVERWORTS).le(2).all())})
        stats.to_csv(OUT / f"G_baseline_values_{ctx}.tsv", sep="\t")

    sw = pd.DataFrame(sweep)
    sw.to_csv(OUT / "G_threshold_sweep.tsv", sep="\t", index=False)
    bl = pd.DataFrame(base)
    bl.to_csv(OUT / "G_baselines.tsv", sep="\t", index=False)

    print("=== threshold sweep: rank of each liverwort ===")
    piv = (sw[sw.species.isin(LIVERWORTS)]
             .pivot_table(index=["context", "species"], columns="threshold", values="rank"))
    print(piv.to_string())
    print("\n=== orthogroups retained per threshold ===")
    print(sw.pivot_table(index="context", columns="threshold", values="n_ogs", aggfunc="first").to_string())
    print("\n=== baselines: |Spearman| against MCi ===")
    print(bl.pivot_table(index="baseline", columns="context", values="spearman_vs_MCi").round(3).to_string())
    print("\n=== does each baseline also put both liverworts in the top two? ===")
    print(bl.pivot_table(index="baseline", columns="context",
                         values="liverworts_top2_under_baseline", aggfunc="first").to_string())
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
