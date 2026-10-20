#!/usr/bin/env python3
"""
Every estimator, every species, every context.

The paper recommends reporting diversity and Fitness-Complexity rather than the
eigenvector, so those rankings have to be on the table for the real genomes and
not only as correlations with the index. This script produces that table.

Usage:
    caffeinate -i python3 21_estimator_table.py
Outputs: data/results_controls/J_estimators_values.tsv, J_estimators_ranks.tsv
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
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
SEED_PLANTS = ["Athaliana", "Slycopersicum", "Pamericana", "Pabies"]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    bm = load("bm", "05_build_matrix.py")
    sim = load("sim", "16_simulation_mci.py")
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))

    values, ranks = [], []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[counts.index.isin(og)] if og else counts
        M = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        mci, _ = bm.compute_complexity(M)
        est = {
            "MCi_RCA": mci,
            "MCi_RCA_maxload": sim.orient_max_loading(mci),
            "diversity": M.sum(axis=0).astype(float),
            "CA_axis1": sim.oriented(sim.ca_axis1(M), M.sum(axis=0)),
            "fitness": sim.fitness_complexity(M),
        }
        for name, series in est.items():
            r = series.rank(ascending=False)
            for sp in M.columns:
                values.append({"context": ctx, "estimator": name, "species": sp,
                               "value": float(series[sp]), "rank": float(r[sp])})
                ranks.append({"context": ctx, "estimator": name, "species": sp,
                              "rank": float(r[sp])})

    v = pd.DataFrame(values); r = pd.DataFrame(ranks)
    OUT.mkdir(parents=True, exist_ok=True)
    v.to_csv(OUT / "J_estimators_values.tsv", sep="\t", index=False)
    r.to_csv(OUT / "J_estimators_ranks.tsv", sep="\t", index=False)

    print("=== ranks, whole genome ===")
    print(r[r.context == "whole_genome"].pivot_table(index="species", columns="estimator",
                                                     values="rank").to_string())
    print("\n=== do both liverworts hold the first two ranks? ===")
    lw = (r[r.species.isin(LIVERWORTS)].groupby(["context", "estimator"])["rank"]
          .apply(lambda x: bool((x <= 2).all())).unstack())
    print(lw.to_string())
    print("\n=== does every non-seed species outrank every seed plant? ===")
    def separates(g):
        worst_nonseed = g[~g.species.isin(SEED_PLANTS)]["rank"].max()
        best_seed = g[g.species.isin(SEED_PLANTS)]["rank"].min()
        return bool(worst_nonseed < best_seed)
    sep = (r.groupby(["context", "estimator"]).apply(separates, include_groups=False)
            .unstack())
    print(sep.to_string())
    print("\nwritten to", OUT)


if __name__ == "__main__":
    main()
