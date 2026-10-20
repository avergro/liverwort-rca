#!/usr/bin/env python3
"""
Simulation study: when does MCi-RCA recover a planted capability signal, and
when does it recover proteome size instead?

Generative model of a species x orthogroup count matrix with two independent
species-level parameters:

  endowment a_s   number of capabilities a species holds, under the capability
                  model of Hausmann and Hidalgo (2011). An orthogroup is
                  occupied only by species holding every capability it requires,
                  so a_s is the quantity the index is supposed to recover.
  duplication d_s mean extra copies per occupied orthogroup, which inflates the
                  gene total N_s without adding new functions. This is the
                  nuisance variable the RCA normalisation is supposed to remove.

endowment and duplication are drawn independently, so any correlation between the index and N_s
is an artefact of the estimator rather than of the data-generating process.

Estimators compared on each simulated matrix:
  MCi_rca      the published index
  MCi_raw      eigenvector on presence/absence
  diversity    number of orthogroups with RCA >= 1
  fitness      Tacchella et al. 2012
  MCi_rarefied MCi-RCA after subsampling every species to the smallest gene total

Reported per cell of the parameter grid:
  rho(estimator, r)   recovery of the planted signal
  rho(estimator, N)   contamination by proteome size

The eigenvector is defined only up to sign, so recovery is scored by the
absolute partial correlation: the published convention orients it by its
correlation with diversity, which is one of the estimators under comparison,
and the alternative rule implemented here (largest loading positive) simply
flips it. Taking absolute values removes the convention from the comparison and
leaves the ordering of estimators unchanged.

Usage:
    caffeinate -i python3 16_simulation_mci.py --reps 20
Outputs: data/results_simulation/
"""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_simulation"


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def simulate(rng, S, G, cap_prob, d, n_cap=60, req_lambda=2.0):
    """
    Capability model of Hausmann and Hidalgo (2011), the generative model that
    motivates the complexity index.

      - there are `n_cap` capabilities;
      - species s carries each capability with probability cap_prob[s], so its
        endowment a_s is the ground truth the index should recover;
      - orthogroup g requires a random subset of capabilities (size ~ 1 +
        Poisson(req_lambda)); a species occupies g only if it holds all of them;
      - occupied entries get 1 + Poisson(d_s) gene copies, so d_s inflates the
        gene total without adding any function.

    This produces the nested species x orthogroup structure the eigenvector
    method was designed to read, and separates capability breadth (signal) from
    gene number (nuisance) by construction.
    """
    has_cap = rng.random((S, n_cap)) < cap_prob[:, None]
    endowment = has_cap.sum(axis=1)
    counts = np.zeros((G, S), dtype=int)
    for g in range(G):
        k = 1 + rng.poisson(req_lambda)
        req = rng.choice(n_cap, size=min(k, n_cap), replace=False)
        occ = np.flatnonzero(has_cap[:, req].all(axis=1))
        if occ.size == 0:
            continue
        counts[g, occ] = 1 + rng.poisson(d[occ])
    return counts, endowment


def ca_axis1(M):
    A = M.values.astype(float)
    P = A / A.sum()
    rr, cc = P.sum(1), P.sum(0)
    Sm = (P - np.outer(rr, cc)) / np.sqrt(np.outer(rr, cc))
    _, _, Vt = np.linalg.svd(Sm, full_matrices=False)
    return pd.Series(Vt[0] / np.sqrt(cc), index=M.columns)


def fitness_complexity(M, iters=200):
    A = M.values.astype(float)
    F = np.ones(A.shape[1]); Q = np.ones(A.shape[0])
    for _ in range(iters):
        Fn = A.T @ Q
        Qn = 1.0 / (A @ (1.0 / np.where(F == 0, 1e-12, F)))
        F = Fn / Fn.mean(); Q = Qn / Qn.mean()
    return pd.Series(F, index=M.columns)


def oriented(series, reference):
    """Eigenvector sign is arbitrary; orient it the way 05_build_matrix.py does."""
    return series if series.corr(reference) >= 0 else -series


def orient_max_loading(series):
    """
    Diversity-free sign rule: make the entry of largest magnitude positive.

    05_build_matrix.py orients the eigenvector by its correlation with diversity,
    which is one of the estimators being compared here, so using that rule in the
    simulation would let the competitor set the sign of the estimator under test.
    """
    return series if series[series.abs().idxmax()] >= 0 else -series


def partial_spearman(x, y, z):
    """Spearman correlation of x and y after removing the linear effect of z on ranks."""
    rx, ry, rz = (pd.Series(v).rank().values for v in (x, y, z))
    Z = np.column_stack([np.ones_like(rz), rz])
    def resid(v):
        beta = np.linalg.lstsq(Z, v, rcond=None)[0]
        return v - Z @ beta
    ex, ey = resid(rx), resid(ry)
    denom = np.sqrt((ex ** 2).sum() * (ey ** 2).sum())
    return float((ex * ey).sum() / denom) if denom else np.nan


def estimators(counts_df, bm, rng, rarefy_reps=5, sign_rule="diversity"):
    M_rca = bm.filter_matrix(bm.build_rca_matrix(counts_df), 2)
    mci_rca, _ = bm.compute_complexity(M_rca)
    M_raw = bm.filter_matrix(bm.build_raw_matrix(counts_df), 2)
    mci_raw, _ = bm.compute_complexity(M_raw)
    div = M_rca.sum(axis=0)
    if sign_rule == "max_loading":
        mci_rca = orient_max_loading(mci_rca)
        mci_raw = orient_max_loading(mci_raw)

    target = int(counts_df.sum(axis=0).min())
    rare = []
    for _ in range(rarefy_reps):
        sub = pd.DataFrame(
            {sp: rng.multivariate_hypergeometric(counts_df[sp].values.astype(int), target)
             for sp in counts_df.columns}, index=counts_df.index)
        Mr = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        m, _ = bm.compute_complexity(Mr)
        if sign_rule == "max_loading":
            m = orient_max_loading(m)
        # the eigenvector sign is arbitrary per replicate: align to the first one
        rare.append(m if not rare else oriented(m, rare[0]))
    mci_rarefied = pd.concat(rare, axis=1).mean(axis=1)

    return {
        "MCi_rca": mci_rca,
        "MCi_raw": mci_raw,
        "diversity": div,
        "CA_axis1": oriented(ca_axis1(M_rca), div),
        "fitness": fitness_complexity(M_rca),
        "MCi_rarefied": mci_rarefied,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--species", type=int, default=11)
    ap.add_argument("--orthogroups", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--sign-rule", choices=("diversity", "max_loading"),
                    default="diversity",
                    help="how to orient the eigenvector; 'diversity' reproduces "
                         "05_build_matrix.py, 'max_loading' is independent of the "
                         "competing estimator")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    bm = load_bm()
    rng = np.random.default_rng(args.seed)
    S, G = args.species, args.orthogroups

    # grid: how strong the planted rarity signal is, and how unequal duplication is
    cap_spreads = [0.0, 0.10, 0.20]           # heterogeneity of capability endowment
    dup_spreads = [0.0, 0.25, 0.5, 1.0, 2.0]  # heterogeneity of duplication across species

    rows = []
    for rs in cap_spreads:
        for ds in dup_spreads:
            for rep in range(args.reps):
                cap_prob = np.clip(0.55 + rs * rng.standard_normal(S), 0.15, 0.95)
                d = np.clip(0.5 + ds * rng.standard_normal(S), 0.0, None)
                counts, endowment = simulate(rng, S, G, cap_prob, d)
                species = [f"sp{j:02d}" for j in range(S)]
                cdf = pd.DataFrame(counts, columns=species,
                                   index=[f"OG{g:05d}" for g in range(G)])
                cdf = cdf[cdf.sum(axis=1) > 0]
                N = cdf.sum(axis=0)
                truth = pd.Series(endowment, index=species)
                for name, est in estimators(cdf, bm, rng,
                                            sign_rule=args.sign_rule).items():
                    rows.append({
                        "cap_spread": rs, "dup_spread": ds, "rep": rep,
                        "estimator": name,
                        "rho_signal": est.corr(truth, method="spearman"),
                        "rho_size": est.corr(N, method="spearman"),
                        "partial_signal_given_size": partial_spearman(est.values, truth.values, N.values),
                        "partial_size_given_signal": partial_spearman(est.values, N.values, truth.values),
                    })
            print(f"  cap_spread={rs} dup={ds}: done")

    df = pd.DataFrame(rows)
    tag = "" if args.sign_rule == "diversity" else "_maxload"
    df.to_csv(OUT / f"S1_grid_raw{tag}.tsv", sep="\t", index=False)
    df["abs_partial_signal"] = df.partial_signal_given_size.abs()
    summ = (df.groupby(["cap_spread", "dup_spread", "estimator"])
              [["rho_signal", "rho_size", "partial_signal_given_size",
                "abs_partial_signal", "partial_size_given_signal"]].mean().reset_index())
    summ.to_csv(OUT / f"S1_grid_summary{tag}.tsv", sep="\t", index=False)

    print("\n=== recovery of the planted signal, rho(estimator, r) ===")
    print(summ.pivot_table(index=["cap_spread", "dup_spread"],
                           columns="estimator", values="rho_signal").round(2).to_string())
    print("\n=== contamination by proteome size, rho(estimator, N) ===")
    print(summ.pivot_table(index=["cap_spread", "dup_spread"],
                           columns="estimator", values="rho_size").round(2).to_string())
    print("\n=== signal recovered after removing proteome size, |partial| ===")
    print(summ.pivot_table(index=["cap_spread", "dup_spread"], columns="estimator",
                           values="abs_partial_signal").round(2).to_string())
    print("\n=== size dependence after removing the signal (partial) ===")
    print(summ.pivot_table(index=["cap_spread", "dup_spread"], columns="estimator",
                           values="partial_size_given_signal").round(2).to_string())
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
