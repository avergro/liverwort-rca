#!/usr/bin/env python3
"""
Two controls requested by review (follow-up to 13_rarefaction_depth_series.py).

A. Sign rule. The index is the second eigenvector of the species transition
   matrix, whose sign is arbitrary. 05_build_matrix.py orients it by its
   correlation with diversity, which is one of the estimators the index is
   later compared against. This script recomputes every real-data index under
   both that rule and the largest-loading-positive rule, across contexts and
   across the RCA threshold sweep, and reports where the two disagree.

B. Symmetric rarefaction. 13_rarefaction_depth_series.py caps every proteome
   at T genes and leaves species below T untouched, so at T = 17,413 the two
   liverworts carry no resampling noise while every competitor does. Here every
   depth is below the smallest proteome, so all eleven species are subsampled
   at every depth and the noise is symmetric.

Usage:
    python3 24_sign_rule_and_symmetric_rarefaction.py --reps 100
Outputs: data/results_controls/J_sign_rule*.tsv, J_symmetric_rarefaction*.tsv
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
THRESHOLDS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def second_eigenvector(M: pd.DataFrame):
    """Unoriented second eigenvector of the species transition matrix."""
    mat = M.values.astype(float)
    k_sp = mat.sum(axis=0)
    k_og = mat.sum(axis=1)
    k_sp_safe = np.where(k_sp == 0, 1, k_sp)
    k_og_safe = np.where(k_og == 0, 1, k_og)
    T = (mat.T / k_sp_safe[:, None]) @ (mat / k_og_safe[:, None])
    vals, vecs = np.linalg.eig(T)
    idx = np.argsort(vals.real)[::-1]
    return vecs[:, idx[1]].real, k_sp


def orient(vec, k_sp, rule):
    if rule == "diversity":
        return -vec if np.corrcoef(k_sp, vec)[0, 1] < 0 else vec
    if rule == "max_loading":
        return -vec if vec[np.argmax(np.abs(vec))] < 0 else vec
    raise ValueError(rule)


def rca_binary(counts, bm, og, threshold=1.0):
    if og is not None:
        counts = counts[counts.index.isin(og)]
    shares = counts / counts.sum(axis=0)
    global_share = counts.sum(axis=1) / counts.values.sum()
    rca = shares.div(global_share, axis=0)
    M = (rca >= threshold).astype(int)
    return M.loc[M.sum(axis=1) >= 2]


# ── A. sign rule ───────────────────────────────────────────────────────────

def sign_rule(counts, bm):
    rows, disagree = [], []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        for thr in THRESHOLDS:
            M = rca_binary(counts, bm, og, thr)
            vec, k_sp = second_eigenvector(M)
            per_rule = {}
            for rule in ("diversity", "max_loading"):
                v = orient(vec, k_sp, rule)
                s = pd.Series(v, index=M.columns)
                per_rule[rule] = s
                r = s.rank(ascending=False, method="min")
                for sp in M.columns:
                    rows.append({"context": ctx, "threshold": thr, "rule": rule,
                                 "species": sp, "MCi": s[sp], "rank": int(r[sp]),
                                 "n_orthogroups": M.shape[0]})
            flipped = bool(np.sign(per_rule["diversity"]).ne(
                np.sign(per_rule["max_loading"])).all())
            rho = per_rule["diversity"].corr(per_rule["max_loading"], method="spearman")
            disagree.append({"context": ctx, "threshold": thr,
                             "n_orthogroups": M.shape[0],
                             "spearman_between_rules": round(rho, 4),
                             "sign_flip": flipped,
                             "leader_diversity_rule":
                                 per_rule["diversity"].idxmax(),
                             "leader_max_loading_rule":
                                 per_rule["max_loading"].idxmax()})
    return pd.DataFrame(rows), pd.DataFrame(disagree)


# ── B. symmetric rarefaction ───────────────────────────────────────────────

def subsample_all(counts, target, rng):
    """Subsample EVERY species to `target` genes. Requires target <= min total."""
    out = {}
    for sp in counts.columns:
        col = counts[sp].values.astype(int)
        out[sp] = rng.multivariate_hypergeometric(col, target)
    return pd.DataFrame(out, index=counts.index)


def symmetric_rarefaction(counts, bm, depths, reps, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for depth in depths:
        for rep in range(reps):
            sub = subsample_all(counts, depth, rng)
            for ctx, filt in CONTEXTS.items():
                og = set(filt.read_text().split()) if filt else None
                M = rca_binary(sub, bm, og)
                vec, k_sp = second_eigenvector(M)
                s = pd.Series(orient(vec, k_sp, "diversity"), index=M.columns)
                r = s.rank(ascending=False, method="min")
                rows.append({"depth": depth, "rep": rep, "context": ctx,
                             "liverworts_top2": bool(
                                 set(LIVERWORTS) <= set(r.nsmallest(2).index)),
                             **{f"rank_{sp}": int(r[sp]) for sp in LIVERWORTS}})
        print(f"  depth {depth}: {reps} replicates done", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    bm = load_bm()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    totals = counts.sum(axis=0)
    print("genes per species:\n", totals.sort_values().to_string())

    print("\n[A] sign rule, both orientations x 5 contexts x 6 thresholds")
    per_species, summary = sign_rule(counts, bm)
    per_species.to_csv(OUT / "J_sign_rule_per_species.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "J_sign_rule_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    smallest = int(totals.min())
    depths = [d for d in (6000, 8000, 10000, 12000) if d < smallest]
    print(f"\n[B] symmetric rarefaction, every species subsampled "
          f"(smallest proteome = {smallest:,}); depths {depths}")
    raw = symmetric_rarefaction(counts, bm, depths, args.reps, args.seed)
    raw.to_csv(OUT / "J_symmetric_rarefaction_raw.tsv", sep="\t", index=False)
    summ = (raw.groupby(["context", "depth"])
               .agg(liverworts_top2=("liverworts_top2", "mean"),
                    mean_rank_Mpolymorpha=("rank_Mpolymorpha", "mean"),
                    mean_rank_Lcruciata=("rank_Lcruciata", "mean"))
               .round(3).reset_index())
    summ.to_csv(OUT / "J_symmetric_rarefaction_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
