#!/usr/bin/env python3
"""
Bootstrap resampling of MCi-RCA to compute:
  1. Confidence intervals for MCi-RCA per species
  2. Ranking stability matrix (how often each species occupies each rank)
  3. Probability of liverworts occupying top-2 ranks (hypothesis test)

Method:
  - Resample orthogroups with replacement from the raw copy-number matrix
  - Recompute RCA and MCi on each bootstrap sample
  - B = 2000 iterations (default)

Usage:
    python3 09_bootstrap_mci.py [--n-bootstrap 2000]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse helper from 05_build_matrix.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib.machinery import SourceFileLoader
bm = SourceFileLoader("05_build_matrix", "05_build_matrix.py").load_module()


def bootstrap_mci(counts: pd.DataFrame, n_bootstrap: int = 2000, seed: int = 42):
    """
    Returns:
      ci_df  — pd.DataFrame with species rows, columns: mean, lower_95, upper_95, std
      rank_matrix — pd.DataFrame: rows=species, cols=rank, values=count (how many
                     times each species occupies each rank position)
      raw_vals — list of pd.Series (MCi-RCA per bootstrap iteration)
    """
    rng = np.random.default_rng(seed)
    n_ogs, n_spp = counts.shape
    species = counts.columns.tolist()
    n_ranks = n_spp

    mci_samples = np.zeros((n_bootstrap, n_spp))
    rank_counts = np.zeros((n_spp, n_ranks), dtype=int)

    for b in range(n_bootstrap):
        # Resample orthogroups with replacement
        idx = rng.integers(0, n_ogs, size=n_ogs)
        boot_counts = counts.iloc[idx]

        # Recompute RCA matrix
        M_rca = bm.build_rca_matrix(boot_counts)
        M_rca_f = bm.filter_matrix(M_rca, min_species=2)

        try:
            MCi_rca, _ = bm.compute_complexity(M_rca_f)
        except Exception:
            # Degenerate matrix — use NaN
            MCi_rca = pd.Series(np.nan, index=M_rca_f.columns)

        # Reindex to full species set
        MCi_aligned = MCi_rca.reindex(species)
        mci_samples[b, :] = MCi_aligned.values

        # Ranking (NaN-safe: put NaN at bottom)
        ranks = np.argsort(np.argsort(-MCi_aligned.values))  # 0 = highest
        for i, r in enumerate(ranks):
            if not np.isnan(mci_samples[b, i]):
                rank_counts[i, r] += 1

    # Confidence intervals (ignoring NaN iterations)
    mci_mean = np.nanmean(mci_samples, axis=0)
    mci_lower = np.nanpercentile(mci_samples, 2.5, axis=0)
    mci_upper = np.nanpercentile(mci_samples, 97.5, axis=0)
    mci_std = np.nanstd(mci_samples, axis=0)

    ci_df = pd.DataFrame({
        "species":             species,
        "MCi_RCA_mean":        mci_mean,
        "MCi_RCA_lower95":     mci_lower,
        "MCi_RCA_upper95":     mci_upper,
        "MCi_RCA_std":         mci_std,
        "MCi_RCA_observed":    np.nan,  # filled below
    }).set_index("species")

    rank_df = pd.DataFrame(
        rank_counts,
        index=species,
        columns=[f"rank_{r+1}" for r in range(n_ranks)]
    )
    rank_df.index.name = "species"

    return ci_df, rank_df, mci_samples, species


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    parser.add_argument("--n-bootstrap", default=2000, type=int)
    parser.add_argument("--out-dir", default="data/results_bootstrap", type=Path)
    parser.add_argument("--og-filter", default=None, type=Path,
                        help="KEGG filter file (OG IDs, one per line)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load copy-number matrix
    tsv = bm.find_orthogroups_tsv(args.orthofinder_dir)
    counts = bm.parse_orthogroups(tsv)
    print(f"Loaded: {counts.shape[0]:,} OGs × {counts.shape[1]} species")

    if args.og_filter:
        og_ids = set(args.og_filter.read_text().split())
        counts = counts[counts.index.isin(og_ids)]
        filter_label = args.og_filter.stem
        print(f"OG filter: {counts.shape[0]:,} OGs ({filter_label})")
    else:
        filter_label = "whole_genome"

    # ── Observed MCi-RCA ──
    M_rca = bm.build_rca_matrix(counts)
    M_rca_f = bm.filter_matrix(M_rca)
    MCi_obs, _ = bm.compute_complexity(M_rca_f)

    print(f"\n{'═'*50}")
    print(f"  Observed MCi-RCA [{filter_label}]")
    print(f"{'═'*50}")
    for sp in MCi_obs.sort_values(ascending=False).index:
        print(f"  {sp:<25} {MCi_obs[sp]:>+10.4f}")

    # ── Bootstrap ──
    print(f"\nRunning bootstrap ({args.n_bootstrap} iterations)...")
    ci_df, rank_df, samples, sp_order = bootstrap_mci(counts, args.n_bootstrap, seed=42)

    ci_df["MCi_RCA_observed"] = MCi_obs
    ci_df = ci_df.sort_values("MCi_RCA_mean", ascending=False)

    print(f"\n{'═'*50}")
    print(f"  Bootstrap 95% CI [{filter_label}]")
    print(f"{'═'*50}")
    print(f"  {'Species':<25} {'Mean':>8}  {'Lower':>8}  {'Upper':>8}  {'Observed':>8}")
    print(f"  {'─'*65}")
    for sp, row in ci_df.iterrows():
        flag = "*" if (row["MCi_RCA_observed"] < row["MCi_RCA_lower95"] or
                       row["MCi_RCA_observed"] > row["MCi_RCA_upper95"]) else " "
        print(f"  {sp:<25} {row['MCi_RCA_mean']:>+8.4f}  {row['MCi_RCA_lower95']:>+8.4f}  "
              f"{row['MCi_RCA_upper95']:>+8.4f}  {row['MCi_RCA_observed']:>+8.4f} {flag}")

    # ── Ranking stability ──
    rank_pct = rank_df.div(rank_df.sum(axis=1), axis=0) * 100

    print(f"\n{'═'*60}")
    print(f"  Ranking stability — probability of each rank (%)")
    print(f"{'═'*60}")
    species_order = MCi_obs.sort_values(ascending=False).index.tolist()
    rank_labels_top = [f"#{i+1}" for i in range(min(5, len(species_order)))]
    cols = [f"rank_{i+1}" for i in range(min(5, len(species_order)))]
    print(f"  {'Species':<25}", end="")
    for rl in rank_labels_top:
        print(f"  {rl:>7}", end="")
    print()
    print(f"  {'─'*25}", end="")
    for _ in rank_labels_top:
        print(f"  {'─'*7}", end="")
    print()
    for sp in species_order:
        print(f"  {sp:<25}", end="")
        for c in cols:
            print(f"  {rank_pct.loc[sp, c]:>6.1f}%", end="")
        print()

    # ── Hypothesis test: liverwort in top-2 probability ──
    liverworts = ["Mpolymorpha", "Lcruciata"]
    angiosperms = ["Athaliana", "Slycopersicum", "Pamericana"]
    # P(both liverworts in top 2 simultaneously)
    top2_both = 0
    top2_any_lw = 0
    for b in range(args.n_bootstrap):
        vals = pd.Series(samples[b], index=sp_order).dropna().sort_values(ascending=False)
        if all(lw in vals.index[:2].tolist() for lw in liverworts):
            top2_both += 1
        if any(lw in vals.index[:2].tolist() for lw in liverworts):
            top2_any_lw += 1

    # P(any liverwort > all angiosperms)
    lw_above_all_angio = 0
    for b in range(args.n_bootstrap):
        vals = pd.Series(samples[b], index=sp_order).dropna()
        lw_vals = vals[liverworts].dropna()
        ang_vals = vals[angiosperms].dropna()
        if len(lw_vals) > 0 and len(ang_vals) > 0:
            if lw_vals.max() > ang_vals.max():
                lw_above_all_angio += 1

    print(f"\n{'═'*60}")
    print(f"  Hypothesis tests (n = {args.n_bootstrap} bootstrap replicates)")
    print(f"{'═'*60}")
    print(f"  P(both liverworts in top 2 simultaneously):  {top2_both/args.n_bootstrap:.4f} ({top2_both}/{args.n_bootstrap})")
    print(f"  P(at least one liverwort in top 2):           {top2_any_lw/args.n_bootstrap:.4f} ({top2_any_lw}/{args.n_bootstrap})")
    print(f"  P(max(liverwort MCi) > max(angiosperm MCi)):   {lw_above_all_angio/args.n_bootstrap:.4f} ({lw_above_all_angio}/{args.n_bootstrap})")

    # ── Save ──
    ci_df.to_csv(args.out_dir / f"bootstrap_ci_{filter_label}.tsv", sep="\t")
    rank_pct.to_csv(args.out_dir / f"bootstrap_rankprob_{filter_label}.tsv", sep="\t")

    # Save all bootstrap MCi values as a matrix
    boot_df = pd.DataFrame(
        samples, columns=sp_order,
        index=range(1, args.n_bootstrap + 1)
    )
    boot_df.index.name = "iteration"
    boot_df.to_csv(args.out_dir / f"bootstrap_mci_values_{filter_label}.tsv", sep="\t")

    print(f"\nOutputs saved to: {args.out_dir}/")


if __name__ == "__main__":
    main()
