#!/usr/bin/env python3
"""Pathway enrichment of the Filter B orthogroups against three backgrounds.

Answers the reviewer objection that Filter B (RCA >= 1 in both liverworts, < 1 in
all three angiosperms) recovers conserved low-copy orthogroups mechanically, so
that its pathway composition needs a null before it is interpreted.

For each KEGG context (primary, C, N):
  1. Rebuild the context RCA matrix from Orthogroups.tsv and check that it equals
     data/results_<ctx>_v3/matrix_rca.tsv and that Filter B equals
     gene_lists/<ctx>/02_OGs_both_liverworts.tsv.
  2. Hypergeometric over-representation of each pathway in Filter B against the
     orthogroups of the context (BH-FDR within context).
  3. Copy-permutation null: within every species, the copy numbers of its occupied
     orthogroups are shuffled among those orthogroups (gene total, occupancy and
     copy-number distribution of each species kept), RCA and Filter B are
     recomputed. The size of Filter B is compared with the null, and the share of
     Filter B taken by each pathway is compared with its share in the null sets.
  4. Stratified test: among orthogroups occupied by both liverworts, the
     association between pathway membership and Filter B is tested within strata
     of mean angiosperm copy number (Cochran-Mantel-Haenszel), which asks whether
     a pathway is in Filter B beyond its low copy number in angiosperms.
Also writes the pathway overlap of Filter B orthogroups (shared KEGG maps) and a
per-orthogroup table (supplementary).

Pathways = union of the KEGG maps of the Arabidopsis genes of each orthogroup;
global and overview maps (ath011xx, ath012xx) are excluded from testing.

Usage: python3 32_filterB_enrichment.py [--reps 2000]
Outputs: data/results_filterB/
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_filterB"
LIVER = ["Mpolymorpha", "Lcruciata"]
ANGIO = ["Athaliana", "Slycopersicum", "Pamericana"]
CTX = ["primary_metabolism", "C_metabolism", "N_metabolism"]
MIN_B = 3          # pathways with fewer Filter B orthogroups are not tested
STRATA = [0, 1.0, 2.0, 4.0, np.inf]   # mean angiosperm copies: (0,1], (1,2], (2,4], >4


def rca_binary(c):
    rca = c.div(c.sum(axis=0), axis=1).div(c.sum(axis=1) / c.values.sum(), axis=0)
    return (rca >= 1).astype(int)


def filter_b(M):
    return (M[LIVER].min(axis=1) == 1) & (M[ANGIO].sum(axis=1) == 0)


def load():
    tsv = sorted((ROOT / "data/orthofinder_v3").glob("Results_*"))[-1] / "Orthogroups" / "Orthogroups.tsv"
    raw = pd.read_csv(tsv, sep="\t", index_col=0)
    counts = raw.map(lambda x: 0 if pd.isna(x) or not str(x).strip() else len(str(x).split(",")))
    ath = {o: {re.sub(r"\.\d+$", "", g.strip()) for g in str(c).split(",")}
           for o, c in raw.Athaliana.dropna().items()}
    gp = pd.read_csv(ROOT / "data/kegg_cache/ath_gene_pathway.tsv", sep="\t", header=None, names=["gene", "path"])
    gp["gene"] = gp.gene.str.replace("ath:", "", regex=False)
    gp["path"] = gp.path.str.replace("path:", "", regex=False)
    gene2path = gp.groupby("gene").path.apply(set).to_dict()
    og2path = {o: {p for p in set().union(*[gene2path.get(g, set()) for g in gs])
                   if not re.match(r"ath01[12]\d\d$", p)} for o, gs in ath.items()}
    names = pd.read_csv(ROOT / "data/kegg_cache/ath_pathway_names.tsv", sep="\t", header=None, names=["path", "name"])
    names["name"] = names.name.str.replace(r" - Arabidopsis thaliana \(thale cress\)$", "", regex=True)
    return raw, counts, og2path, dict(zip(names.path, names.name))


def permute_counts(c, rng):
    """Shuffle each species' copy numbers among the orthogroups it occupies."""
    a = c.values.copy()
    for j in range(a.shape[1]):
        nz = np.flatnonzero(a[:, j])
        a[nz, j] = rng.permutation(a[nz, j])
    return pd.DataFrame(a, index=c.index, columns=c.columns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    raw, counts, og2path, name = load()
    rng = np.random.default_rng(args.seed)
    all_rows, overlap_rows, og_rows, null_size = [], [], [], []

    for ctx in CTX:
        ids = set((ROOT / "data/filters" / f"{ctx}_ogs.txt").read_text().split())
        c = counts[counts.index.isin(ids)]
        M = rca_binary(c)
        M = M[M.sum(axis=1) >= 2]
        ref = pd.read_csv(ROOT / f"data/results_{ctx}_v3/matrix_rca.tsv", sep="\t", index_col=0)
        assert M.index.equals(ref.index) and (M[ref.columns].values == ref.values).all(), f"{ctx}: RCA matrix differs"
        B = set(M.index[filter_b(M)])
        refB = set(pd.read_csv(ROOT / "gene_lists" / ctx / "02_OGs_both_liverworts.tsv", sep="\t", index_col=0).index)
        assert B == refB, f"{ctx}: Filter B differs ({len(B)} vs {len(refB)})"
        U = list(M.index)
        print(f"{ctx}: {len(U)} orthogroups, Filter B {len(B)} (matrix and Filter B reproduced)")

        # pathway membership over the context background
        paths = sorted({p for o in U for p in og2path.get(o, ())})
        memb = pd.DataFrame({p: [p in og2path.get(o, ()) for o in U] for p in paths}, index=U)
        inB = memb.index.isin(B)
        kB = memb[inB].sum()
        tested = [p for p in paths if kB[p] >= MIN_B]

        # copy-permutation null
        nullk = np.zeros((args.reps, len(tested)), dtype=int)
        nulln = np.zeros(args.reps, dtype=int)
        tcol = memb[tested]
        for r in range(args.reps):
            Mr = rca_binary(permute_counts(c, rng))
            Mr = Mr[Mr.sum(axis=1) >= 2]
            Br = Mr.index[filter_b(Mr)]
            nullk[r] = tcol.reindex(Br).fillna(False).astype(int).sum().values
            nulln[r] = len(Br)
            null_size.append((ctx, r, len(Br)))

        # strata of mean angiosperm copy number, among orthogroups occupied by both liverworts
        occ = c.loc[U, LIVER].min(axis=1) > 0
        stratum = pd.cut(c.loc[U, ANGIO].mean(axis=1), STRATA, labels=False, include_lowest=True)

        for i, p in enumerate(tested):
            K, N, n, k = int(memb[p].sum()), len(U), len(B), int(kB[p])
            ph = hypergeom.sf(k - 1, N, K, n)
            # composition: share of Filter B in pathway p, observed vs null
            share_null = nullk[:, i] / np.maximum(nulln, 1)
            null_share = share_null.mean()
            p_perm = (1 + (share_null >= k / n).sum()) / (args.reps + 1)
            tables = []
            for s in sorted(stratum[occ].dropna().unique()):
                sel = occ & (stratum == s)
                a = int((sel & memb[p] & inB).sum()); b = int((sel & memb[p] & ~inB).sum())
                cc = int((sel & ~memb[p] & inB).sum()); d = int((sel & ~memb[p] & ~inB).sum())
                if min(a + b, cc + d, a + cc, b + d) > 0:
                    tables.append(np.array([[a, b], [cc, d]]) + 0.5)
            if tables:
                st = StratifiedTable(tables)
                or_mh, p_mh = st.oddsratio_pooled, st.test_null_odds(correction=True).pvalue
            else:
                or_mh, p_mh = np.nan, np.nan
            all_rows.append(dict(context=ctx, pathway=p, pathway_name=name.get(p, p),
                                 n_filterB=k, n_context=K, filterB_size=n, context_size=N,
                                 fold_vs_context=(k / n) / (K / N), p_hypergeom=ph,
                                 share_filterB=k / n, share_null=null_share,
                                 fold_vs_null=(k / n) / null_share if null_share else np.nan,
                                 p_permutation=p_perm, OR_MH=or_mh, p_MH=p_mh))

        # overlap of KEGG maps among Filter B orthogroups (top pathways)
        top = kB[kB >= MIN_B].sort_values(ascending=False).index[:12]
        for x in top:
            for y in top:
                if x < y:
                    sh = int((memb.loc[inB, x] & memb.loc[inB, y]).sum())
                    if sh:
                        overlap_rows.append(dict(context=ctx, pathway_a=name.get(x, x), pathway_b=name.get(y, y),
                                                 n_a=int(kB[x]), n_b=int(kB[y]), shared=sh))
        for o in sorted(B):
            og_rows.append(dict(context=ctx, orthogroup=o,
                                **{f"copies_{s}": int(counts.at[o, s]) for s in counts.columns},
                                kegg_pathways="; ".join(sorted(name.get(p, p) for p in og2path.get(o, ()))),
                                Athaliana=raw.at[o, "Athaliana"], Mpolymorpha=raw.at[o, "Mpolymorpha"],
                                Lcruciata=raw.at[o, "Lcruciata"]))

    res = pd.DataFrame(all_rows)
    for col in ["p_hypergeom", "p_permutation", "p_MH"]:
        q = np.full(len(res), np.nan)
        for ctx, g in res.groupby("context"):
            ok = g[col].notna()
            q[g.index[ok]] = multipletests(g.loc[ok, col], method="fdr_bh")[1]
        res[col.replace("p_", "q_")] = q
    res = res.sort_values(["context", "n_filterB"], ascending=[True, False])
    res.to_csv(OUT / "filterB_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(overlap_rows).to_csv(OUT / "filterB_pathway_overlap.tsv", sep="\t", index=False)
    pd.DataFrame(og_rows).to_csv(OUT / "filterB_orthogroups.tsv", sep="\t", index=False)
    ns = pd.DataFrame(null_size, columns=["context", "rep", "filterB_size"])
    ns.to_csv(OUT / "filterB_null_size.tsv", sep="\t", index=False)
    print(ns.groupby("context").filterB_size.describe(percentiles=[.025, .975]).round(1))
    show = ["context", "pathway_name", "n_filterB", "n_context", "fold_vs_context", "q_hypergeom",
            "share_null", "fold_vs_null", "q_permutation", "OR_MH", "q_MH"]
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_colwidth", 42):
        print(res[show].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
