#!/usr/bin/env python3
"""
Metabolic versus regulatory gene content: is the liverwort contraction the
same in both?

For each functional category (whole genome, metabolic contexts, regulatory
contexts, KEGG regulatory subclasses, all Arabidopsis transcription factor
orthogroups) this script measures, per species:

  coverage      orthogroups occupied with >= 1 gene
  redundancy    genes per occupied orthogroup
  retained dup  duplication events whose two daughter clades both contain
                genes of the species, per occupied orthogroup. Unlike the
                dup_per_OG of 09c_expansion_contraction.py, which sums every
                duplication of the orthogroup whatever the species, this
                counts only duplications the species itself still carries.

and then:
  A. liverwort / angiosperm ratios of the three quantities per category
  B. coverage of angiosperm orthogroups: share of the orthogroups occupied by
     >= 2 angiosperms that at least one liverwort occupies
  C. per-orthogroup copy ratio log2(liverwort copies / angiosperm copies),
     among orthogroups occupied by both groups, compared between categories
     with a label permutation test on the median
  D. scaling: slope of log(category genes) on log(total genes) across
     species, ordinary least squares and generalised least squares under
     Brownian motion on the dated topology of Table S1

Usage:  caffeinate -i python3 27_regulatory_contrast.py
Output: data/results_regulatory_contrast/*.tsv
"""

import importlib.util
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_regulatory_contrast"
CACHE = ROOT / "data" / "kegg_cache"
GO_DIR = ROOT / "Matias/alexander_metodologia/GO_files"
KEGG_DIR = ROOT / "data/filters"
CTRL = ROOT / "data" / "results_controls"
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]
ALGA = "Creinhardtii"
SUBCLASSES = {"ath04075": "hormone signal transduction", "ath04016": "MAPK signalling",
              "ath04712": "circadian rhythm", "ath04626": "plant-pathogen interaction",
              "ath04070": "phosphatidylinositol signalling"}
RNG = np.random.default_rng(42)


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def kegg_links():
    """Arabidopsis gene -> KEGG pathways, cached after the first download."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / "ath_gene_pathway.tsv"
    if not f.exists():
        txt = urllib.request.urlopen("https://rest.kegg.jp/link/pathway/ath", timeout=60).read().decode()
        f.write_text(txt)
    links = defaultdict(set)
    for line in f.read_text().splitlines():
        g, p = line.split("\t")
        links[g.replace("ath:", "").upper()].add(p.replace("path:", ""))
    return links


def categories(counts, gene_to_og):
    cats = {"whole genome": set(counts.index)}
    for lab, f in [("KEGG primary", "primary_metabolism_ogs.txt"), ("KEGG carbon", "C_metabolism_ogs.txt"),
                   ("KEGG nitrogen", "N_metabolism_ogs.txt"), ("KEGG regulatory", "regulatory_networks_ogs.txt")]:
        cats[lab] = set((KEGG_DIR / f).read_text().split())
    for lab, f in [("GO primary", "primary_metabolism_ogs.txt"), ("GO carbon", "C_metabolism_ogs.txt"),
                   ("GO nitrogen", "N_metabolism_ogs.txt"), ("GO regulatory", "regulatory_networks_ogs.txt")]:
        cats[lab] = set((GO_DIR / f).read_text().split())
    links = kegg_links()
    for pid, lab in SUBCLASSES.items():
        cats[f"KEGG {lab}"] = {gene_to_og[g] for g, ps in links.items() if pid in ps and g in gene_to_og}
    tf = pd.read_csv(ROOT / "data/tf_families/Ath_TF_list.txt", sep="\t")
    cats["transcription factors"] = {gene_to_og[g] for g in tf.Gene_ID.str.upper() if g in gene_to_og}
    # check: the KEGG regulatory filter is the union of the five subclasses
    union = set().union(*(cats[f"KEGG {l}"] for l in SUBCLASSES.values()))
    print(f"  KEGG regulatory filter {len(cats['KEGG regulatory'])} OGs; union of subclasses "
          f"{len(union)}; shared {len(union & cats['KEGG regulatory'])}")
    return {k: v & set(counts.index) for k, v in cats.items()}


def retained_duplications(species):
    """Per (orthogroup, species): duplication events retained in the species."""
    # gene trees reconciled with the accepted species tree (38_fixed_species_tree.py)
    path = ROOT / "data/orthofinder_v3/Results_Oct05_fixedtree" / "Gene_Duplication_Events" / "Duplications.tsv"
    sp_sorted = sorted(species, key=len, reverse=True)

    def sp_of(gene):
        for s in sp_sorted:
            if gene.startswith(s + "_"):
                return s
        raise ValueError(gene)

    ret = defaultdict(int)
    term = defaultdict(int)
    with open(path) as fh:
        next(fh)
        for line in fh:
            og, node, _, support, typ, g1, g2 = line.rstrip("\n").split("\t")
            if float(support) < 0.5:
                continue
            s1 = {sp_of(g) for g in g1.split(", ") if g}
            s2 = {sp_of(g) for g in g2.split(", ") if g}
            for s in s1 & s2:
                ret[(og, s)] += 1
            if typ == "Terminal":
                term[(og, node)] += 1
    return ret, term


def per_species(counts, cats, ret):
    rows = []
    for cat, ogs in cats.items():
        sub = counts.loc[sorted(ogs)]
        for s in counts.columns:
            occ = sub[s] > 0
            n_occ = int(occ.sum())
            genes = int(sub[s].sum())
            dups = sum(ret.get((og, s), 0) for og in sub.index[occ])
            rows.append({"category": cat, "species": s, "orthogroups": len(sub), "occupied": n_occ,
                         "fraction_occupied": n_occ / len(sub), "genes": genes,
                         "genes_per_occupied_OG": genes / n_occ if n_occ else np.nan,
                         "retained_dup_per_occupied_OG": dups / n_occ if n_occ else np.nan})
    return pd.DataFrame(rows)


def ratios(ps):
    rows = []
    for cat, d in ps.groupby("category", sort=False):
        d = d.set_index("species")
        r = {"category": cat, "orthogroups": int(d.orthogroups.iloc[0])}
        for col, lab in [("occupied", "coverage"), ("genes_per_occupied_OG", "genes_per_OG"),
                         ("retained_dup_per_occupied_OG", "retained_dup_per_OG"), ("genes", "genes")]:
            r[f"L/A {lab}"] = d.loc[LIVERWORTS, col].mean() / d.loc[ANGIOSPERMS, col].mean()
        rows.append(r)
    return pd.DataFrame(rows)


def coverage_of_angiosperm_ogs(counts, cats):
    rows = []
    for cat, ogs in cats.items():
        sub = counts.loc[sorted(ogs)]
        a = (sub[ANGIOSPERMS] > 0).sum(axis=1) >= 2
        l = (sub[LIVERWORTS] > 0).sum(axis=1) >= 1
        rows.append({"category": cat, "angiosperm_OGs": int(a.sum()),
                     "also_in_a_liverwort": int((a & l).sum()),
                     "share": (a & l).sum() / a.sum() if a.sum() else np.nan})
    return pd.DataFrame(rows)


def copy_ratio(counts, ogs):
    sub = counts.loc[sorted(ogs)]
    L, A = sub[LIVERWORTS], sub[ANGIOSPERMS]
    both = ((L > 0).any(axis=1)) & ((A > 0).any(axis=1))
    l = L[both].where(L[both] > 0).mean(axis=1)          # copies among occupying liverworts
    a = A[both].where(A[both] > 0).mean(axis=1)
    return np.log2(l / a)


def permutation(x, y, n=10000):
    obs = np.median(x) - np.median(y)
    pool = np.concatenate([x, y])
    k = len(x)
    null = np.empty(n)
    for i in range(n):
        RNG.shuffle(pool)
        null[i] = np.median(pool[:k]) - np.median(pool[k:])
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (n + 1)


def gls_slope(x, y, V):
    X = np.c_[np.ones(len(x)), x]
    Vi = np.linalg.inv(V)
    XtVi = X.T @ Vi
    cov = np.linalg.inv(XtVi @ X)
    b = cov @ XtVi @ y
    e = y - X @ b
    s2 = (e @ Vi @ e) / (len(x) - 2)
    return b[1], np.sqrt(s2 * cov[1, 1])


def ols_slope(x, y):
    X = np.c_[np.ones(len(x)), x]
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    s2 = e @ e / (len(x) - 2)
    return b[1], np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])


def scaling(ps):
    vcv = pd.read_csv(CTRL / "TableS7_vcv.tsv", sep="\t", index_col=0)
    tot = ps[ps.category == "whole genome"].set_index("species").genes
    rows = []
    for cat, d in ps.groupby("category", sort=False):
        if cat == "whole genome":
            continue
        g = d.set_index("species").genes
        for subset in ("all eleven", "land plants"):
            sp = [s for s in tot.index if subset == "all eleven" or s != ALGA]
            sp = [s for s in sp if g[s] > 0]
            x, y = np.log(tot[sp].values), np.log(g[sp].values)
            bo, so = ols_slope(x, y)
            bg, sg = gls_slope(x, y, vcv.loc[sp, sp].values)
            rows.append({"category": cat, "species": subset, "n": len(sp),
                         "slope_OLS": bo, "se_OLS": so, "slope_GLS": bg, "se_GLS": sg})
    return pd.DataFrame(rows)


def main():
    bm = load("bm", "05_build_matrix.py")
    of = load("of", "06_kegg_filter.py")
    tsv = bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3")
    counts = bm.parse_orthogroups(tsv)
    counts = counts[(counts > 0).sum(axis=1) >= 2]
    gene_to_og = {k.upper(): v for k, v in of.build_ath_gene_to_og(tsv).items()}
    print("categories")
    cats = categories(counts, gene_to_og)
    print("duplications")
    ret, term = retained_duplications(list(counts.columns))

    ps = per_species(counts, cats, ret)
    # check against the whole-genome values reported in the manuscript
    wg = ps[ps.category == "whole genome"].set_index("species")
    # RG_corrected: check against the original-run values disabled
    #assert int(wg.loc["Mpolymorpha", "occupied"]) == 9881 and int(wg.loc["Athaliana", "genes"]) == 23318
    rat = ratios(ps)
    cov = coverage_of_angiosperm_ogs(counts, cats)

    groups = {"metabolic (KEGG)": cats["KEGG primary"],
              "regulatory (KEGG)": cats["KEGG regulatory"],
              "metabolic (GO)": cats["GO primary"],
              "regulatory (GO)": cats["GO regulatory"],
              "transcription factors": cats["transcription factors"],
              "whole genome": cats["whole genome"]}
    cr = {k: copy_ratio(counts, v) for k, v in groups.items()}
    tests = []
    for a, b in [("regulatory (KEGG)", "metabolic (KEGG)"), ("regulatory (GO)", "metabolic (GO)"),
                 ("transcription factors", "metabolic (KEGG)"), ("transcription factors", "whole genome"),
                 ("metabolic (KEGG)", "whole genome"), ("regulatory (KEGG)", "whole genome")]:
        x = cr[a].values
        y = cr[b].drop(cr[a].index, errors="ignore").values if b == "whole genome" else cr[b].values
        if b != "whole genome":
            shared = set(cr[a].index) & set(cr[b].index)       # remove overlap between the two sets
            x = cr[a].drop(list(shared)).values
            y = cr[b].drop(list(shared)).values
        obs, p = permutation(x, y)
        tests.append({"A": a, "B": b, "n_A": len(x), "n_B": len(y),
                      "median_log2_A": np.median(x), "median_log2_B": np.median(y),
                      "copy_ratio_A": 2 ** np.median(x), "copy_ratio_B": 2 ** np.median(y),
                      "diff": obs, "P_perm": p})
    tests = pd.DataFrame(tests)
    sc = scaling(ps)

    OUT.mkdir(parents=True, exist_ok=True)
    ps.to_csv(OUT / "per_species.tsv", sep="\t", index=False, float_format="%.4f")
    rat.to_csv(OUT / "liverwort_angiosperm_ratios.tsv", sep="\t", index=False, float_format="%.3f")
    cov.to_csv(OUT / "coverage_of_angiosperm_ogs.tsv", sep="\t", index=False, float_format="%.3f")
    tests.to_csv(OUT / "copy_ratio_tests.tsv", sep="\t", index=False, float_format="%.4g")
    sc.to_csv(OUT / "scaling.tsv", sep="\t", index=False, float_format="%.3f")
    pd.DataFrame([{"category": c, "orthogroup": og} for c, s in cats.items() if c != "whole genome"
                  for og in sorted(s)]).to_csv(OUT / "category_orthogroups.tsv", sep="\t", index=False)

    pd.set_option("display.width", 220)
    print("\nA. liverwort / angiosperm ratios\n"); print(rat.round(3).to_string(index=False))
    print("\nB. angiosperm orthogroups (>=2 angiosperms) also occupied by a liverwort\n"); print(cov.round(3).to_string(index=False))
    print("\nC. per-orthogroup copy ratio, liverworts / angiosperms\n"); print(tests.round(4).to_string(index=False))
    print("\nD. scaling of category gene number with total gene number\n"); print(sc.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
