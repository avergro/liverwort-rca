#!/usr/bin/env python3
"""
Eigenvalue spectrum, correspondence-analysis inertia and nestedness.

Two questions raised in review: is the index a graded axis or a two-block split
of the species set, and how nested is each matrix. Nestedness matters because a
perfectly nested matrix makes any eigenvector ranking a function of the
marginals, which would explain why the index tracks gene number.

NODF follows Almeida-Neto et al. 2008 and is compared against the same
degree-preserving (curveball) null used elsewhere. It is computed exactly, not
by subsampling: with eleven species every orthogroup is one of at most 2^11
binary patterns, so rows are grouped by pattern and pairs are weighted by the
multiplicity of each pattern. That turns roughly 121 million row pairs into at
most a few million pattern pairs and gives the same number.

Usage:
    caffeinate -i python3 20_spectrum_nestedness.py --null-reps 200
Outputs: data/results_controls/I_spectrum.tsv, I_nestedness.tsv
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


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def spectrum(M):
    A = M.values.astype(float)
    ks = A.sum(0); kg = A.sum(1)
    ks_s = np.where(ks == 0, 1, ks); kg_s = np.where(kg == 0, 1, kg)
    T = (A.T / ks_s[:, None]) @ (A / kg_s[:, None])
    ev = np.sort(np.linalg.eigvals(T).real)[::-1]
    P = A / A.sum()
    r, c = P.sum(1), P.sum(0)
    S = (P - np.outer(r, c)) / np.sqrt(np.outer(r, c))
    sv = np.linalg.svd(S, compute_uv=False)
    inertia = sv ** 2 / (sv ** 2).sum()
    return ev, inertia


def _paired_small(mat):
    """Exact NODF contribution for a matrix with few rows."""
    n = mat.shape[0]
    fills = mat.sum(1)
    total, pairs = 0.0, 0
    for i in range(n):
        for j in range(i + 1, n):
            if fills[i] == fills[j] or min(fills[i], fills[j]) == 0:
                continue
            hi, lo = (i, j) if fills[i] > fills[j] else (j, i)
            total += 100.0 * np.logical_and(mat[hi], mat[lo]).sum() / fills[lo]
            pairs += 1
    return total, pairs


def _paired_patterns(mat):
    """
    Exact NODF contribution over rows, grouped by binary pattern.

    Rows with the same pattern contribute nothing between themselves (equal
    fill) and contribute identically against every other pattern, so the sum
    over row pairs equals the sum over pattern pairs weighted by multiplicity.
    """
    patterns, counts = np.unique(mat, axis=0, return_counts=True)
    fills = patterns.sum(1)
    total, pairs = 0.0, 0
    for i in range(len(patterns)):
        for j in range(i + 1, len(patterns)):
            if fills[i] == fills[j] or min(fills[i], fills[j]) == 0:
                continue
            hi, lo = (i, j) if fills[i] > fills[j] else (j, i)
            w = counts[i] * counts[j]
            total += w * 100.0 * np.logical_and(patterns[hi], patterns[lo]).sum() / fills[lo]
            pairs += w
    return total, pairs


def nodf(A):
    """NODF for a binary matrix (Almeida-Neto et al. 2008), rows and columns."""
    tr, pr = _paired_patterns(A)
    tc, pc = _paired_small(A.T)
    return (tr + tc) / (pr + pc) if (pr + pc) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    bm = load("bm", "05_build_matrix.py")
    ctrl = load("ctrl", "12_gene_number_controls.py")
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    rng = np.random.default_rng(args.seed)

    spec_rows, nest_rows = [], []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[counts.index.isin(og)] if og else counts
        M = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        ev, inertia = spectrum(M)
        spec_rows.append({"context": ctx, "n_orthogroups": M.shape[0],
                          "lambda2": ev[1], "lambda3": ev[2], "lambda4": ev[3],
                          "gap_2_3": ev[1] - ev[2],
                          "ratio_3_2": ev[2] / ev[1],
                          "CA_inertia_axis1": inertia[0],
                          "CA_inertia_axis2": inertia[1]})

        A = M.values.astype(int)
        obs = nodf(A)
        null = []
        for _ in range(args.null_reps):
            rand = ctrl.curveball(A, rng, 5 * A.shape[0])
            rand = rand[rand.sum(1) >= 2]
            null.append(nodf(rand))
        null = np.array(null)
        nest_rows.append({"context": ctx, "NODF_observed": obs,
                          "NODF_null_mean": null.mean(), "NODF_null_sd": null.std(ddof=0),
                          "z": (obs - null.mean()) / null.std(ddof=0),
                          "p_one_sided": ((null >= obs).sum() + 1) / (args.null_reps + 1)})
        print(f"  {ctx}: done")

    sp = pd.DataFrame(spec_rows); ne = pd.DataFrame(nest_rows)
    sp.to_csv(OUT / "I_spectrum.tsv", sep="\t", index=False)
    ne.to_csv(OUT / "I_nestedness.tsv", sep="\t", index=False)
    print("\n=== eigenvalue spectrum and CA inertia ===")
    print(sp.round(3).to_string(index=False))
    print("\n=== nestedness (NODF) against the degree-preserving null ===")
    print(ne.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
