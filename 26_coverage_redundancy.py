#!/usr/bin/env python3
"""
Coverage and redundancy: the two components of gene number.

Gene number = orthogroups occupied x genes per occupied orthogroup. This script
measures both components per species and context, asks how much of the
variation in gene number among land plants each component carries, and how
closely the RCA ranking follows each of them. It also counts how many Filter B
orthogroups (RCA >= 1 in both liverworts, below threshold in all angiosperms)
are enriched in the moss and in the lycophyte.

Unlike K_presence_coverage.tsv, whose `genes` column is the whole-genome total
in every context, genes here are counted inside each context.

Usage: python3 26_coverage_redundancy.py
Output: data/results_controls/L_coverage_redundancy.tsv
        data/results_controls/L_decomposition.tsv
        data/results_controls/L_rank_vs_redundancy.tsv
        data/results_controls/L_filterB_sharing.tsv
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"
GO_DIR = ROOT / "Matias/alexander_metodologia/GO_files"
KEGG_DIR = ROOT / "data/filters"
CONTEXTS = {
    "whole_genome": None,
    "KEGG primary": KEGG_DIR / "primary_metabolism_ogs.txt",
    "KEGG carbon": KEGG_DIR / "C_metabolism_ogs.txt",
    "KEGG nitrogen": KEGG_DIR / "N_metabolism_ogs.txt",
    "KEGG regulatory": KEGG_DIR / "regulatory_networks_ogs.txt",
    "GO primary": GO_DIR / "primary_metabolism_ogs.txt",
    "GO carbon": GO_DIR / "C_metabolism_ogs.txt",
    "GO nitrogen": GO_DIR / "N_metabolism_ogs.txt",
    "GO regulatory": GO_DIR / "regulatory_networks_ogs.txt",
}
# estimator table uses the KEGG context names of 05_build_matrix.py
EST_CONTEXT = {"whole_genome": "whole_genome", "KEGG primary": "primary_metabolism",
               "KEGG carbon": "C_metabolism", "KEGG nitrogen": "N_metabolism",
               "KEGG regulatory": "regulatory_networks"}
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]
ALGA = "Creinhardtii"


def load_bm():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def coverage_table(counts):
    rows = []
    for ctx, filt in CONTEXTS.items():
        sub = counts if filt is None else counts[counts.index.isin(set(filt.read_text().split()))]
        occ = (sub > 0).sum(axis=0)
        genes = sub.sum(axis=0)
        for sp in counts.columns:
            rows.append({"context": ctx, "species": sp, "orthogroups_in_context": sub.shape[0],
                         "occupied": int(occ[sp]), "fraction_occupied": occ[sp] / sub.shape[0],
                         "genes_in_context": int(genes[sp]),
                         "genes_per_occupied_OG": genes[sp] / occ[sp]})
    return pd.DataFrame(rows)


def decomposition(cov):
    """Share of var(log genes) carried by log coverage and by log genes per OG.

    log G = log occupied + log (G / occupied), so the two covariances with log G
    sum to var(log G) exactly; each share is cov(component, log G) / var(log G).
    """
    rows = []
    for ctx, d in cov.groupby("context", sort=False):
        d = d[d.species != ALGA]
        lg = np.log(d.genes_in_context)
        lo = np.log(d.occupied)
        lr = np.log(d.genes_per_occupied_OG)
        v = lg.var()
        s_occ, s_red = np.cov(lo, lg)[0, 1] / v, np.cov(lr, lg)[0, 1] / v
        assert abs(s_occ + s_red - 1) < 1e-9, (ctx, s_occ, s_red)
        cv = lambda x: x.std() / x.mean()
        rows.append({"context": ctx, "n_species": len(d),
                     "cv_occupied": cv(d.occupied), "cv_genes_per_OG": cv(d.genes_per_occupied_OG),
                     "share_var_log_genes_coverage": s_occ,
                     "share_var_log_genes_redundancy": s_red})
    return pd.DataFrame(rows)


def rank_vs_redundancy(cov):
    est = pd.read_csv(OUT / "J_estimators_values.tsv", sep="\t")
    wg_genes = cov[cov.context == "whole_genome"].set_index("species").genes_in_context
    rows = []
    for ctx, ectx in EST_CONTEXT.items():
        c = cov[cov.context == ctx].set_index("species")
        for name in ("MCi_RCA", "diversity", "fitness"):
            v = est[(est.context == ectx) & (est.estimator == name)].set_index("species").value
            for subset, sp in (("all eleven", list(c.index)),
                               ("land plants", [s for s in c.index if s != ALGA])):
                d = pd.DataFrame({"v": v, "red": c.genes_per_occupied_OG,
                                  "occ": c.occupied, "genes": wg_genes}).loc[sp]
                r = d.corr(method="spearman")["v"]
                rows.append({"context": ctx, "estimator": name, "species": subset,
                             "rho_genes_per_OG": r.red, "rho_occupied": r.occ,
                             "rho_genome_gene_number": r.genes})
    return pd.DataFrame(rows)


def filterB_sharing():
    rows = []
    for ctx in ("primary_metabolism", "C_metabolism", "N_metabolism"):
        b = pd.read_csv(ROOT / "gene_lists" / ctx / "02_OGs_both_liverworts.tsv",
                        sep="\t", index_col=0)
        assert (b[LIVERWORTS] == 1).all().all() and (b[ANGIOSPERMS] == 0).all().all()
        rows.append({"context": ctx, "filterB": len(b),
                     "Ppatens": int(b.Ppatens.sum()), "Smoellendorffii": int(b.Smoellendorffii.sum()),
                     "Ppatens_or_Smoellendorffii": int(((b.Ppatens + b.Smoellendorffii) > 0).sum()),
                     "Ppatens_and_Smoellendorffii": int(((b.Ppatens + b.Smoellendorffii) == 2).sum()),
                     "any_other_non_seed": int((b[["Ppatens", "Smoellendorffii", "Afiliculoides",
                                                   "Scucullata"]].sum(axis=1) > 0).sum())})
    return pd.DataFrame(rows)


def main():
    bm = load_bm()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    counts = counts[(counts > 0).sum(axis=1) >= 2]          # same min_species filter as 05

    cov = coverage_table(counts)
    # check against the whole-genome values already reported (K table, manuscript 3.9)
    wg = cov[cov.context == "whole_genome"].set_index("species")
    k = pd.read_csv(OUT / "K_presence_coverage.tsv", sep="\t")
    kwg = k[k.context == "whole_genome"].set_index("species")
    assert (wg.occupied == kwg.occupied).all() and (wg.genes_in_context == kwg.genes).all()
    # RG_corrected: check against the original-run values disabled
    #assert round(wg.loc["Mpolymorpha", "genes_per_occupied_OG"], 2) == 1.37

    dec = decomposition(cov)
    rvr = rank_vs_redundancy(cov)
    share = filterB_sharing()

    OUT.mkdir(parents=True, exist_ok=True)
    cov.to_csv(OUT / "L_coverage_redundancy.tsv", sep="\t", index=False, float_format="%.4f")
    dec.to_csv(OUT / "L_decomposition.tsv", sep="\t", index=False, float_format="%.4f")
    rvr.to_csv(OUT / "L_rank_vs_redundancy.tsv", sep="\t", index=False, float_format="%.3f")
    share.to_csv(OUT / "L_filterB_sharing.tsv", sep="\t", index=False)

    pd.set_option("display.width", 200)
    print("\ncoverage (fraction occupied) / genes per occupied orthogroup\n")
    print(cov.pivot(index="species", columns="context", values="fraction_occupied")
          [list(CONTEXTS)].round(2).to_string())
    print()
    print(cov.pivot(index="species", columns="context", values="genes_per_occupied_OG")
          [list(CONTEXTS)].round(2).to_string())
    print("\nliverwort / angiosperm ratios\n")
    for ctx in CONTEXTS:
        c = cov[cov.context == ctx].set_index("species")
        print(f"  {ctx:<16} coverage {c.loc[LIVERWORTS].occupied.mean() / c.loc[ANGIOSPERMS].occupied.mean():.2f}"
              f"   genes per OG {c.loc[LIVERWORTS].genes_per_occupied_OG.mean() / c.loc[ANGIOSPERMS].genes_per_occupied_OG.mean():.2f}")
    print("\nvariance decomposition of log gene number, ten land plants\n")
    print(dec.round(3).to_string(index=False))
    print("\nSpearman of each estimator with redundancy, coverage and gene number\n")
    print(rvr.round(2).to_string(index=False))
    print("\nFilter B orthogroups also enriched in the moss and the lycophyte\n")
    print(share.to_string(index=False))


if __name__ == "__main__":
    main()
