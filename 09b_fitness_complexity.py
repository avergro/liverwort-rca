#!/usr/bin/env python3
"""
Fitness-Complexity Index (Tacchella et al. 2012) for metabolic complexity.

The Fitness-Complexity (FC) algorithm uses coupled non-linear equations
that are less sensitive to ubiquitous products than the linear eigenvector
method of Hidalgo & Hausmann.

Algorithm:
  F_s^(0) = 1         (fitness — species complexity)
  Q_g^(0) = 1         (complexity — orthogroup complexity)

  F_s^(n+1) = Σ_g M_{g,s} · Q_g^(n)
  Q_g^(n+1) = 1 / Σ_s M_{g,s} · (1 / F_s^(n))

  Normalize after each iteration: mean(F) = 1, mean(Q) = 1
  Iterate until convergence (||F^(n+1) - F^(n)|| < tol)

Usage:
    python3 09b_fitness_complexity.py [--max-iter 100] [--tol 1e-6]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib.machinery import SourceFileLoader
bm = SourceFileLoader("05_build_matrix", "05_build_matrix.py").load_module()


def fitness_complexity(M: pd.DataFrame, max_iter: int = 100, tol: float = 1e-6):
    """
    Compute Fitness (species) and Complexity (orthogroups) à la Tacchella et al. 2012.

    M: binary matrix (n_ogs × n_species)
    Returns: (F_s Series, Q_g Series, n_iter)
    """
    mat = M.values.astype(float).T  # shape: (n_species, n_ogs)
    n_spp, n_ogs = mat.shape

    F = np.ones(n_spp)
    Q = np.ones(n_ogs)

    for iteration in range(1, max_iter + 1):
        Q_new = 1.0 / (mat.T @ (1.0 / (F + 1e-12)) + 1e-12)
        F_new = mat @ Q_new

        # Normalize
        Q_new /= Q_new.mean()
        F_new /= F_new.mean()

        delta = np.max(np.abs(F_new - F))
        F, Q = F_new, Q_new

        if delta < tol:
            break

    return (pd.Series(F, index=M.columns, name="Fitness_FC"),
            pd.Series(Q, index=M.index, name="Complexity_FC"),
            iteration)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    parser.add_argument("--out-dir", default="data/results_fitness_complexity", type=Path)
    parser.add_argument("--max-iter", default=100, type=int)
    parser.add_argument("--tol", default=1e-6, type=float)
    parser.add_argument("--og-filter", default=None, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    tsv = bm.find_orthogroups_tsv(args.orthofinder_dir)
    counts = bm.parse_orthogroups(tsv)

    if args.og_filter:
        og_ids = set(args.og_filter.read_text().split())
        counts = counts[counts.index.isin(og_ids)]
        filter_label = args.og_filter.stem
    else:
        filter_label = "whole_genome"

    # Build RCA matrix
    M_rca = bm.build_rca_matrix(counts)
    M_rca_f = bm.filter_matrix(M_rca)

    # ── Hidalgo MCi-RCA (for comparison) ──
    MCi_rca, GCi_rca = bm.compute_complexity(M_rca_f)

    # ── Fitness-Complexity ──
    F_fc, Q_fc, n_iter = fitness_complexity(M_rca_f, args.max_iter, args.tol)

    print(f"\nFitness-Complexity converged in {n_iter} iterations")
    print(f"Matrix: {M_rca_f.shape[0]:,} OGs × {M_rca_f.shape[1]} species\n")

    # Comparison table
    results = pd.DataFrame({
        "MCi_Hidalgo": MCi_rca,
        "Fitness_FC":  F_fc,
    })
    results["Rank_Hidalgo"] = results["MCi_Hidalgo"].rank(ascending=False)
    results["Rank_FC"]     = results["Fitness_FC"].rank(ascending=False)
    results["Rank_delta"]  = results["Rank_Hidalgo"] - results["Rank_FC"]
    results = results.sort_values("Fitness_FC", ascending=False)

    print(f"  {'Species':<25} {'Fitness_FC':>10} {'Rank_FC':>6}  "
          f"{'MCi_Hidalgo':>10} {'Rank_H':>6}  {'ΔRank':>6}")
    print(f"  {'─'*75}")
    for sp, row in results.iterrows():
        print(f"  {sp:<25} {row['Fitness_FC']:>+10.6f} {int(row['Rank_FC']):>6}  "
              f"{row['MCi_Hidalgo']:>+10.6f} {int(row['Rank_Hidalgo']):>6}  "
              f"{int(row['Rank_delta']):>+6}")

    # Spearman correlation between methods (using numpy rank-based Pearson)
    def spearman_r(x, y):
        rx = x.rank()
        ry = y.rank()
        return rx.corr(ry)
    rho = spearman_r(MCi_rca, F_fc)
    print(f"\n  Spearman ρ(Hidalgo, FC) = {rho:.4f}")

    # Agreement on top species
    top3_hidalgo = set(results.sort_values("MCi_Hidalgo", ascending=False).head(3).index)
    top3_fc = set(results.head(3).index)
    print(f"  Top-3 agreement: {len(top3_hidalgo & top3_fc)}/3 species in common")
    print(f"  Top-3 Hidalgo: {top3_hidalgo}")
    print(f"  Top-3 FC:      {top3_fc}")

    # Save
    results.to_csv(args.out_dir / f"species_complexity_fc_{filter_label}.tsv", sep="\t")

    og_results = pd.DataFrame({
        "GCi_Hidalgo":    GCi_rca,
        "Complexity_FC":  Q_fc,
    })
    og_results.to_csv(args.out_dir / f"og_complexity_fc_{filter_label}.tsv", sep="\t")

    print(f"\nOutputs saved to: {args.out_dir}/")


if __name__ == "__main__":
    main()
