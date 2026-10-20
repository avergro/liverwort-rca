#!/usr/bin/env python3
"""
Sensitivity of the phylogenetic regression to the outgroup and to the root age.

Chlamydomonas is the only species that uses the root calibration, the least
supported age in the tree, and the only one that is not a land plant. Here the
regression is repeated (i) with the alga removed, where the tree reduces to the
embryophyte clade and the root age no longer enters, and (ii) with the root age
set to 750, 1000 and 1500 Ma with the alga kept, and (iii) with every land plant
node at the younger or at the older bound of the intervals given in Table S1.

Usage:
    caffeinate -i python3 19_pgls_sensitivity.py
Outputs: data/results_controls/H_pgls_sensitivity_fits.tsv, H_pgls_sensitivity_ranks.tsv
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"

EMBRYO = ("(((Mpolymorpha:280,Lcruciata:280):185,Ppatens:465):39,"
          "(Smoellendorffii:450,((Scucullata:137,Afiliculoides:137):286,"
          "(Pabies:348,(Pamericana:136,(Athaliana:124,Slycopersicum:124):12):212):75):27):54);")


def dated(a):
    """Newick with the alga at 1000 Ma and the land plant node ages of dict a (Ma)."""
    return (f"(Creinhardtii:1000,(((Mpolymorpha:{a['mlu']},Lcruciata:{a['mlu']}):{a['bry'] - a['mlu']},"
            f"Ppatens:{a['bry']}):{a['emb'] - a['bry']},(Smoellendorffii:{a['tra']},"
            f"((Scucullata:{a['salv']},Afiliculoides:{a['salv']}):{a['euph'] - a['salv']},"
            f"(Pabies:{a['sperm']},(Pamericana:{a['meso']},(Athaliana:{a['penta']},Slycopersicum:{a['penta']})"
            f":{a['meso'] - a['penta']}):{a['sperm'] - a['meso']}):{a['euph'] - a['sperm']}):{a['tra'] - a['euph']})"
            f":{a['emb'] - a['tra']}):{1000 - a['emb']});")


# every node at the younger or the older bound of the intervals of Table S1
AGE_BOUNDS = {
    "younger": dict(emb=493, bry=450, mlu=228, tra=447, euph=414, salv=89, sperm=331, meso=135, penta=121),
    "older": dict(emb=515, bry=479, mlu=311, tra=452, euph=432, salv=153, sperm=365, meso=137, penta=127),
}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fit(pg, bm, counts, newick, species, label):
    V_full = pg.vcv(pg.parse_newick(newick))
    V = V_full.loc[species, species].values
    genes = counts[species].sum(axis=0)
    busco = pd.Series(pg.BUSCO).reindex(species)
    rows, ranks = [], []
    for ctx, filt in pg.CONTEXTS.items():
        og = set(filt.read_text().split()) if filt else None
        sub = counts[species]
        sub = sub[sub.index.isin(og)] if og else sub
        M = bm.filter_matrix(bm.build_rca_matrix(sub), 2)
        mci, _ = bm.compute_complexity(M)
        y = mci.reindex(species).values
        for model, cols in [("log_genes", [np.log(genes.values)]),
                            ("log_genes+busco", [np.log(genes.values), busco.values])]:
            X = np.column_stack([np.ones(len(species))] + cols)
            lam, (beta, se, resid, ll, s2) = pg.pgls(y, X, V)
            rows.append({"setting": label, "context": ctx, "model": model,
                         "lambda_ML": lam, "beta_log_genes": beta[1],
                         "t_log_genes": beta[1] / se[1],
                         "beta_busco": beta[2] if len(beta) > 2 else np.nan,
                         "t_busco": (beta[2] / se[2]) if len(beta) > 2 else np.nan})
            if model == "log_genes":
                r = pd.Series(resid, index=species).rank(ascending=False)
                for sp in species:
                    ranks.append({"setting": label, "context": ctx, "species": sp,
                                  "residual": float(pd.Series(resid, index=species)[sp]),
                                  "rank_residual": float(r[sp])})
    return rows, ranks


def main():
    pg = load("pg", "14_pgls_mci.py")
    bm = load("bm", "05_build_matrix.py")
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    all_species = list(counts.columns)
    land = [s for s in all_species if s != "Creinhardtii"]

    fits, ranks = [], []
    for root_age in (750, 1000, 1500):
        nw = pg.NEWICK.replace("Creinhardtii:1000", f"Creinhardtii:{root_age}") \
                      .replace("):496);", f"):{root_age - 504});")
        f, r = fit(pg, bm, counts, nw, all_species, f"all species, root {root_age} Ma")
        fits += f; ranks += r
    for bound, ages in AGE_BOUNDS.items():
        f, r = fit(pg, bm, counts, dated(ages), all_species, f"all species, node ages at the {bound} bounds")
        fits += f; ranks += r
    f, r = fit(pg, bm, counts, EMBRYO, land, "Chlamydomonas removed")
    fits += f; ranks += r

    fits = pd.DataFrame(fits); ranks = pd.DataFrame(ranks)
    OUT.mkdir(parents=True, exist_ok=True)
    fits.to_csv(OUT / "H_pgls_sensitivity_fits.tsv", sep="\t", index=False)
    ranks.to_csv(OUT / "H_pgls_sensitivity_ranks.tsv", sep="\t", index=False)

    print("=== lambda and gene-number effect, model with gene number only ===")
    sub = fits[fits.model == "log_genes"]
    for col in ("lambda_ML", "beta_log_genes", "t_log_genes"):
        print(f"\n{col}")
        print(sub.pivot_table(index="setting", columns="context", values=col)
                 .round(2).to_string())
    print("\n=== residual ranks (model with gene number only) ===")
    piv = ranks.pivot_table(index=["setting", "species"], columns="context",
                            values="rank_residual")
    print(piv.to_string())
    print("\nwritten to", OUT)


if __name__ == "__main__":
    main()
