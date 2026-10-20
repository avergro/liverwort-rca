#!/usr/bin/env python3
"""
Phylogenetic generalised least squares of MCi-RCA against proteome size.

Answers the panel's request for a test that separates "liverworts rank first"
from "bryophytes are an outgroup to the seed plants that dominate the sample".

Model:  MCi ~ log(genes assigned to orthogroups)  [+ BUSCO completeness]
        errors ~ MVN(0, sigma^2 * V(lambda))
V is the phylogenetic variance-covariance matrix of a dated topology; Pagel's
lambda is estimated by maximum likelihood on a grid (lambda = 0 is ordinary
least squares, lambda = 1 is Brownian motion). Residuals per species say which
species exceed the expectation given their proteome size and their position in
the tree.

Divergence times (Ma) are the midpoints of published 95% HPD intervals: Harris et
al. 2022 (Nat Ecol Evol) for the deep nodes, Morris et al. 2018 (PNAS) for the seed
plants, Magallon et al. 2015 (New Phytol) for the angiosperm nodes, Villarreal et al.
2016 and Flores et al. 2023 (New Phytol) for the liverwort split and the Fern Tree of
Life (Nitta et al. 2022) for Azolla and Salvinia (Table S1). The tree is hard-coded
below and printed with the results. Bryophyte monophyly is assumed.

Usage:
    caffeinate -i python3 14_pgls_mci.py
Outputs: data/results_controls/F_pgls*.tsv
"""

import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"

# Dated topology, branch lengths in Ma. Bryophyte monophyly assumed.
NEWICK = (
    "(Creinhardtii:1000,"
    "(((Mpolymorpha:280,Lcruciata:280):185,Ppatens:465):39,"
    "(Smoellendorffii:450,((Scucullata:137,Afiliculoides:137):286,"
    "(Pabies:348,(Pamericana:136,(Athaliana:124,Slycopersicum:124):12):212):75):27):54):496);"
)

CONTEXTS = {
    "whole_genome": None,
    "primary_metabolism": ROOT / "data/filters/primary_metabolism_ogs.txt",
    "C_metabolism": ROOT / "data/filters/C_metabolism_ogs.txt",
    "N_metabolism": ROOT / "data/filters/N_metabolism_ogs.txt",
    "regulatory_networks": ROOT / "data/filters/regulatory_networks_ogs.txt",
}
BUSCO = {  # complete %; Afiliculoides updated for the v1.1 proteome (RG_corrected)
    "Athaliana": 99.6, "Creinhardtii": 99.0, "Pamericana": 92.7, "Slycopersicum": 91.3,
    "Ppatens": 88.5, "Pabies": 88.2, "Mpolymorpha": 87.4, "Smoellendorffii": 76.9,
    "Afiliculoides": 73.7, "Lcruciata": 74.0, "Scucullata": 72.3,
}
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]


# ── phylogeny ─────────────────────────────────────────────────────────────

def parse_newick(s: str):
    """Return nested (children, length) tree; tips are (name, length)."""
    tokens = re.findall(r"\(|\)|,|[^(),:]+|:[0-9.]+", s.replace(";", ""))
    pos = 0

    def node():
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = [node()]
            while pos < len(tokens) and tokens[pos] == ",":
                pos += 1
                children.append(node())
            pos += 1  # ')'
            name = None
            if pos < len(tokens) and not tokens[pos].startswith(":") \
                    and tokens[pos] not in (",", ")"):
                name = tokens[pos]; pos += 1
        else:
            children, name = None, tokens[pos]
            pos += 1
        length = 0.0
        if pos < len(tokens) and tokens[pos].startswith(":"):
            length = float(tokens[pos][1:]); pos += 1
        return {"children": children, "name": name, "length": length}

    return node()


def vcv(tree):
    """Phylogenetic variance-covariance: shared root-to-tip path length."""
    tips = []

    def collect(nd, depth):
        if nd["children"] is None:
            tips.append((nd["name"], depth + nd["length"]))
            return [nd["name"]]
        names = []
        for ch in nd["children"]:
            names += collect(ch, depth + nd["length"])
        nd["_tips"] = names
        nd["_depth"] = depth
        return names

    collect(tree, 0.0)
    names = [t[0] for t in tips]
    height = max(t[1] for t in tips)
    V = pd.DataFrame(0.0, index=names, columns=names)

    def fill(nd, depth):
        if nd["children"] is None:
            return
        shared = depth + nd["length"]
        for a in nd["_tips"]:
            for b in nd["_tips"]:
                V.loc[a, b] = max(V.loc[a, b], shared)
        for ch in nd["children"]:
            fill(ch, depth + nd["length"])

    fill(tree, 0.0)
    for n, d in tips:
        V.loc[n, n] = d
    return V / height


# ── GLS with Pagel's lambda ───────────────────────────────────────────────

def gls(y, X, V):
    Vi = np.linalg.inv(V)
    XtVi = X.T @ Vi
    beta = np.linalg.solve(XtVi @ X, XtVi @ y)
    resid = y - X @ beta
    n, p = X.shape
    s2 = float(resid.T @ Vi @ resid) / (n - p)
    cov = s2 * np.linalg.inv(XtVi @ X)
    sign, logdet = np.linalg.slogdet(V)
    ll = -0.5 * (n * np.log(2 * np.pi * s2) + logdet + float(resid.T @ Vi @ resid) / s2)
    return beta, np.sqrt(np.diag(cov)), resid, ll, s2


def lambda_transform(V, lam):
    Vl = V * lam
    np.fill_diagonal(Vl, np.diag(V))
    return Vl


def pgls(y, X, V):
    grid = np.linspace(0, 1, 101)
    best = None
    for lam in grid:
        try:
            out = gls(y, X, lambda_transform(V.copy(), lam))
        except np.linalg.LinAlgError:
            continue
        if best is None or out[3] > best[1][3]:
            best = (lam, out)
    return best


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    bm = importlib.util.module_from_spec(spec); spec.loader.exec_module(bm)
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))

    V_full = vcv(parse_newick(NEWICK))
    species = list(counts.columns)
    missing = set(species) - set(V_full.index)
    if missing:
        raise SystemExit(f"species missing from the tree: {missing}")
    V = V_full.loc[species, species].values
    genes = counts.sum(axis=0).reindex(species)
    busco = pd.Series(BUSCO).reindex(species)

    fits, residuals = [], []
    for ctx, filt in CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[counts.index.isin(og)] if og else counts
        M = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        mci, _ = bm.compute_complexity(M)
        y = mci.reindex(species).values
        for label, cols in [("log_genes", [np.log(genes.values)]),
                            ("log_genes+busco", [np.log(genes.values), busco.values])]:
            X = np.column_stack([np.ones(len(species))] + cols)
            lam, (beta, se, resid, ll, s2) = pgls(y, X, V)
            t = beta / se
            fits.append({
                "context": ctx, "model": label, "lambda_ML": lam,
                "beta_log_genes": beta[1], "se_log_genes": se[1], "t_log_genes": t[1],
                "beta_busco": beta[2] if len(beta) > 2 else np.nan,
                "t_busco": t[2] if len(beta) > 2 else np.nan,
                "logLik": ll,
            })
            if label == "log_genes":
                for sp, r in zip(species, resid):
                    residuals.append({"context": ctx, "species": sp, "residual": r})

    fits = pd.DataFrame(fits)
    res = pd.DataFrame(residuals)
    res["rank_residual"] = res.groupby("context")["residual"].rank(ascending=False)
    fits.to_csv(OUT / "F_pgls_fits.tsv", sep="\t", index=False)
    res.to_csv(OUT / "F_pgls_residuals.tsv", sep="\t", index=False)

    print("Tree (Ma, bryophyte monophyly assumed):")
    print(NEWICK)
    print("\nPGLS fits:")
    print(fits.round(4).to_string(index=False))
    print("\nResiduals (MCi minus the phylogenetic expectation given proteome size):")
    piv = res.pivot_table(index="species", columns="context", values="residual")
    rnk = res.pivot_table(index="species", columns="context", values="rank_residual")
    print(piv.round(3).to_string())
    print("\nResidual ranks (1 = most in excess of expectation):")
    print(rnk.to_string())
    print("\nLiverwort residual ranks:")
    print(rnk.loc[LIVERWORTS].to_string())
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
