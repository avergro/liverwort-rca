#!/usr/bin/env python3
"""
Three extensions of the capability simulation (16_simulation_mci.py), each one
relaxing an idealisation that the first version made:

  phylogeny     capabilities are inherited along a tree instead of being drawn
                independently per species, so related species share repertoire
                and the columns of the matrix are not independent;
  completeness  a species-specific fraction of genes is dropped after the matrix
                is generated, mimicking incomplete proteome annotation, with
                completeness drawn independently of the capability endowment;
  neofunction   a fraction of duplicated copies gains a new capability, so
                duplication is no longer purely redundant and gene number becomes
                causally linked to functional breadth.

Each scenario is run alone and in combination with the others, and the same
estimators as in 16_simulation_mci.py are scored against the true endowment.

Usage:
    caffeinate -i python3 17_simulation_extensions.py --reps 20
Outputs: data/results_simulation/S2_extensions*.tsv
"""

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_simulation"


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def random_tree(S, rng, height=1.0):
    """Balanced-ish ultrametric tree: returns per-species list of branch segments."""
    clades = [([i], height) for i in range(S)]
    segments = {i: [] for i in range(S)}   # list of (segment_id, length)
    seg_id = 0
    while len(clades) > 1:
        rng.shuffle(clades)
        (a, ha), (b, hb) = clades.pop(), clades.pop()
        join = min(ha, hb) * rng.uniform(0.3, 0.8)
        for members, h in ((a, ha), (b, hb)):
            for m in members:
                segments[m].append((seg_id, h - join))
            seg_id += 1
        clades.append((a + b, join))
    root = clades[0][1]
    for m in range(S):
        segments[m].append((seg_id, root))
    return segments


def capabilities_on_tree(S, n_cap, rng, rate=2.5):
    """Simulate presence/absence of each capability along the tree."""
    segments = random_tree(S, rng)
    has = np.zeros((S, n_cap), dtype=bool)
    for c in range(n_cap):
        root_state = rng.random() < 0.6
        seg_flip = {}
        for sp in range(S):
            state = root_state
            for seg_id, length in reversed(segments[sp]):
                if seg_id not in seg_flip:
                    seg_flip[seg_id] = rng.random() < (1 - np.exp(-rate * length))
                if seg_flip[seg_id]:
                    state = not state
            has[sp, c] = state
    return has


def build_matrix(rng, S, G, has_cap, d, req_lambda=2.0):
    n_cap = has_cap.shape[1]
    counts = np.zeros((G, S), dtype=int)
    for g in range(G):
        k = 1 + rng.poisson(req_lambda)
        req = rng.choice(n_cap, size=min(k, n_cap), replace=False)
        occ = np.flatnonzero(has_cap[:, req].all(axis=1))
        if occ.size:
            counts[g, occ] = 1 + rng.poisson(d[occ])
    return counts


def apply_completeness(counts, completeness, rng):
    out = counts.copy()
    for s, c in enumerate(completeness):
        col = out[:, s]
        out[:, s] = rng.binomial(col, c)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--species", type=int, default=11)
    ap.add_argument("--orthogroups", type=int, default=5000)
    ap.add_argument("--n-cap", type=int, default=60)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    bm = load("bm", "05_build_matrix.py")
    sim = load("sim", "16_simulation_mci.py")
    rng = np.random.default_rng(args.seed)
    S, G, n_cap = args.species, args.orthogroups, args.n_cap

    scenarios = {
        "baseline": dict(phylo=False, incomplete=False, neofunc=0.0),
        "phylogeny": dict(phylo=True, incomplete=False, neofunc=0.0),
        "incomplete": dict(phylo=False, incomplete=True, neofunc=0.0),
        "neofunction": dict(phylo=False, incomplete=False, neofunc=0.5),
        "all_three": dict(phylo=True, incomplete=True, neofunc=0.5),
    }

    rows = []
    for name, cfg in scenarios.items():
        for rep in range(args.reps):
            d = np.clip(0.5 + 1.0 * rng.standard_normal(S), 0.0, None)
            if cfg["phylo"]:
                has_cap = capabilities_on_tree(S, n_cap, rng)
            else:
                p = np.clip(0.55 + 0.2 * rng.standard_normal(S), 0.15, 0.95)
                has_cap = rng.random((S, n_cap)) < p[:, None]
            if cfg["neofunc"]:
                # duplication-rich species gain extra capabilities
                extra = (cfg["neofunc"] * d * 10).astype(int)
                for s in range(S):
                    off = np.flatnonzero(~has_cap[s])
                    if off.size and extra[s]:
                        pick = rng.choice(off, size=min(extra[s], off.size), replace=False)
                        has_cap[s, pick] = True
            endowment = has_cap.sum(axis=1)
            counts = build_matrix(rng, S, G, has_cap, d)
            completeness = np.ones(S)
            if cfg["incomplete"]:
                completeness = rng.uniform(0.70, 1.0, S)
                counts = apply_completeness(counts, completeness, rng)

            species = [f"sp{j:02d}" for j in range(S)]
            cdf = pd.DataFrame(counts, columns=species,
                               index=[f"OG{g:05d}" for g in range(G)])
            cdf = cdf[cdf.sum(axis=1) > 0]
            if cdf.shape[0] < 100 or (cdf.sum(axis=0) == 0).any():
                continue
            N = cdf.sum(axis=0)
            truth = pd.Series(endowment, index=species)
            comp = pd.Series(completeness, index=species)
            for est_name, est in sim.estimators(cdf, bm, rng).items():
                rows.append({
                    "scenario": name, "rep": rep, "estimator": est_name,
                    "rho_signal": est.corr(truth, method="spearman"),
                    "partial_signal_given_size":
                        sim.partial_spearman(est.values, truth.values, N.values),
                    "rho_size": est.corr(N, method="spearman"),
                    "rho_completeness": est.corr(comp, method="spearman"),
                })
        print(f"  {name}: done")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "S2_extensions_raw.tsv", sep="\t", index=False)
    summ = (df.groupby(["scenario", "estimator"])
              [["rho_signal", "partial_signal_given_size", "rho_size", "rho_completeness"]]
              .mean().reset_index())
    summ.to_csv(OUT / "S2_extensions_summary.tsv", sep="\t", index=False)

    for col, label in [("partial_signal_given_size", "signal recovered, controlling for size"),
                       ("rho_size", "correlation with gene number"),
                       ("rho_completeness", "correlation with annotation completeness")]:
        print(f"\n=== {label} ===")
        print(summ.pivot_table(index="scenario", columns="estimator", values=col)
                  .round(2).to_string())
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
