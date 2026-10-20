#!/usr/bin/env python3
"""
Controls for the gene-number confound in MCi-RCA (panel review, 2026-09-17).

Four analyses, all on the OrthoFinder v3 copy-number matrix:

  A. Correlation between MCi-RCA and gene number / diversity / BUSCO, per context.
  B. Rarefaction: every proteome subsampled to the smallest gene total
     (multivariate hypergeometric over each species' orthogroup counts),
     MCi-RCA recomputed per replicate.
  C. Species jackknife: leave-one-species-out, including the Chlamydomonas re-run.
  D. Degree-preserving null model (curveball) on the binarized RCA matrix,
     giving a z-score per species against matrices with identical row and
     column sums.

MCi eigenvector entries have arbitrary scale, so every replicate is compared by
rank and by z-scored MCi within its own run.

Usage:
    caffeinate -i python3 12_gene_number_controls.py --reps 100 --null-reps 500
Outputs: data/results_controls/
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


def load_module():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mci_rca(counts: pd.DataFrame, bm, og_filter=None) -> pd.Series:
    """MCi-RCA for a copy-number matrix, mirroring 05_build_matrix.py."""
    if og_filter is not None:
        counts = counts[counts.index.isin(og_filter)]
    M = bm.build_rca_matrix(counts)
    M = bm.filter_matrix(M, 2)
    if M.shape[0] < 10:
        return pd.Series(np.nan, index=counts.columns)
    mci, _ = bm.compute_complexity(M)
    return mci


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=0)


# ── A. Correlation with gene number ───────────────────────────────────────

def analysis_a(counts, bm, busco):
    rows = []
    genes = counts.sum(axis=0)
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[counts.index.isin(og)] if og else counts
        M = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        mci = mci_rca(counts, bm, og)
        div = M.sum(axis=0)
        rows.append({
            "context": ctx,
            "spearman_MCi_vs_genes": mci.corr(genes, method="spearman"),
            "pearson_MCi_vs_log_genes": mci.corr(np.log(genes), method="pearson"),
            "spearman_MCi_vs_diversity": mci.corr(div, method="spearman"),
            "spearman_MCi_vs_busco": mci.corr(busco, method="spearman") if busco is not None else np.nan,
        })
    return pd.DataFrame(rows)


# ── B. Rarefaction ────────────────────────────────────────────────────────

def analysis_b(counts, bm, reps, seed):
    rng = np.random.default_rng(seed)
    target = int(counts.sum(axis=0).min())
    print(f"  rarefying every proteome to {target:,} genes, {reps} replicates")
    records = []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        for r in range(reps):
            rare = pd.DataFrame(
                {sp: rng.multivariate_hypergeometric(counts[sp].values.astype(int), target)
                 for sp in counts.columns},
                index=counts.index,
            )
            mci = mci_rca(rare, bm, og)
            z = zscore(mci)
            rank = mci.rank(ascending=False)
            for sp in counts.columns:
                records.append({"context": ctx, "rep": r, "species": sp,
                                "MCi": mci[sp], "MCi_z": z[sp], "rank": rank[sp]})
        print(f"    {ctx}: done")
    return pd.DataFrame(records)


# ── C. Species jackknife ──────────────────────────────────────────────────

def analysis_c(counts, bm):
    records = []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        for dropped in ["none"] + list(counts.columns):
            sub = counts if dropped == "none" else counts.drop(columns=[dropped])
            mci = mci_rca(sub, bm, og)
            rank = mci.rank(ascending=False)
            for sp in sub.columns:
                records.append({"context": ctx, "dropped": dropped, "species": sp,
                                "MCi": mci[sp], "rank": rank[sp],
                                "n_species": sub.shape[1]})
    return pd.DataFrame(records)


# ── D. Degree-preserving null (curveball) ─────────────────────────────────

def curveball(mat: np.ndarray, rng, n_swaps: int) -> np.ndarray:
    """Strona et al. 2014. Preserves both row and column sums exactly."""
    presences = [set(np.flatnonzero(row)) for row in mat]
    n = len(presences)
    for _ in range(n_swaps):
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        a, b = presences[i], presences[j]
        shared = a & b
        only_a = list(a - shared)
        only_b = list(b - shared)
        if not only_a or not only_b:
            continue
        pool = only_a + only_b
        rng.shuffle(pool)
        k = len(only_a)
        presences[i] = shared | set(pool[:k])
        presences[j] = shared | set(pool[k:])
    out = np.zeros_like(mat)
    for i, cols in enumerate(presences):
        if cols:
            out[i, list(cols)] = 1
    return out


def analysis_d(counts, bm, null_reps, seed):
    rng = np.random.default_rng(seed)
    records = []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[counts.index.isin(og)] if og else counts
        M = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        obs, _ = bm.compute_complexity(M)
        mat = M.values.astype(int)
        n_swaps = 5 * mat.shape[0]
        null = []
        for _ in range(null_reps):
            rand = pd.DataFrame(curveball(mat, rng, n_swaps),
                                index=M.index, columns=M.columns)
            rand = rand[rand.sum(axis=1) >= 2]
            mci, _ = bm.compute_complexity(rand)
            null.append(zscore(mci))
        null = pd.DataFrame(null)
        obs_z = zscore(obs)
        for sp in M.columns:
            greater = (null[sp] >= obs_z[sp]).sum()
            records.append({
                "context": ctx, "species": sp,
                "MCi_z_observed": obs_z[sp],
                "null_mean": null[sp].mean(), "null_sd": null[sp].std(ddof=0),
                "z_vs_null": (obs_z[sp] - null[sp].mean()) / null[sp].std(ddof=0),
                "p_one_sided": (greater + 1) / (null_reps + 1),
            })
        print(f"    {ctx}: done ({null_reps} null matrices)")
    return pd.DataFrame(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthofinder-dir", default=ROOT / "data/orthofinder_v3", type=Path)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--null-reps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    bm = load_module()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(args.orthofinder_dir))
    print(f"counts matrix: {counts.shape[0]:,} OGs x {counts.shape[1]} species")

    busco_file = ROOT / "data/busco/busco_summary.tsv"
    busco = None
    if busco_file.exists():
        b = pd.read_csv(busco_file, sep="\t", index_col=0)
        col = [c for c in b.columns if c.lower().startswith("c")]
        if col:
            busco = b[col[0]].reindex(counts.columns)

    print("\nA. correlation of MCi-RCA with gene number")
    a = analysis_a(counts, bm, busco)
    a.to_csv(OUT / "A_correlations.tsv", sep="\t", index=False)
    print(a.to_string(index=False))

    print("\nB. rarefaction")
    b = analysis_b(counts, bm, args.reps, args.seed)
    b.to_csv(OUT / "B_rarefaction_raw.tsv", sep="\t", index=False)
    summ = (b.groupby(["context", "species"])
              .agg(mean_rank=("rank", "mean"), sd_rank=("rank", "std"),
                   best_rank=("rank", "min"), worst_rank=("rank", "max"),
                   mean_MCi_z=("MCi_z", "mean"),
                   frac_top2=("rank", lambda r: (r <= 2).mean()))
              .reset_index().sort_values(["context", "mean_rank"]))
    summ.to_csv(OUT / "B_rarefaction_summary.tsv", sep="\t", index=False)
    both = (b[b.species.isin(LIVERWORTS)].groupby(["context", "rep"])["rank"]
              .apply(lambda r: (r <= 2).all()).groupby("context").mean())
    both.to_csv(OUT / "B_rarefaction_liverworts_top2.tsv", sep="\t")
    print(summ.to_string(index=False))
    print("\nfraction of replicates with both liverworts in the top two:")
    print(both.to_string())

    print("\nC. species jackknife")
    c = analysis_c(counts, bm)
    c.to_csv(OUT / "C_jackknife.tsv", sep="\t", index=False)
    piv = (c[c.species.isin(LIVERWORTS)]
             .pivot_table(index=["context", "dropped"], columns="species", values="rank"))
    print(piv.to_string())

    print("\nD. degree-preserving null model (curveball)")
    d = analysis_d(counts, bm, args.null_reps, args.seed)
    d.to_csv(OUT / "D_null_model.tsv", sep="\t", index=False)
    print(d[d.species.isin(LIVERWORTS)].to_string(index=False))

    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
