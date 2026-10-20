#!/usr/bin/env python3
"""
Bootstrap intervals for the three liverwort/angiosperm ratios of the main text.

  coverage    mean orthogroups occupied by the two liverworts / mean of the three angiosperms (1.02)
  redundancy  mean genes per occupied orthogroup, same groups (0.53)
  genes/KO    genes per KEGG ortholog of Marchantia / mean of the three angiosperms (0.62)

Coverage and redundancy are resampled over the 11,383 orthogroups present in at least
two species (as in 26_coverage_redundancy.py), whole genome and each functional context.
Genes per KO is resampled over the 825 KOs of the 118 modules complete in at least one
land plant (as in 33_kegg_module_completeness.py), and is also given for both liverworts.
2,000 replicates, seed 42, percentile intervals.

Usage: python3 35_bootstrap_ratios.py
Output: data/results_controls/M_bootstrap_ratios.tsv
"""

import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"
KM = ROOT / "data" / "results_kegg_modules"
REPS, SEED = 2000, 42
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / file)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def og_ratios(c):
    """c: orthogroups x species gene counts (numpy, columns liverworts then angiosperms)."""
    occ = (c > 0).sum(axis=0)
    red = c.sum(axis=0) / occ
    return occ[:2].mean() / occ[2:].mean(), red[:2].mean() / red[2:].mean()


def ko_ratio(c, who):
    """c: KOs x species gene counts; genes per KO = mean over the KOs a species carries."""
    with np.errstate(invalid="ignore"):
        gpk = np.nanmean(np.where(c > 0, c, np.nan), axis=0)
    return gpk[who].mean() / gpk[2:].mean()


def interval(fun, mat, rng):
    n = mat.shape[0]
    b = np.array([fun(mat[rng.integers(0, n, n)]) for _ in range(REPS)])
    return np.percentile(b, [2.5, 97.5], axis=0)


def main():
    rng = np.random.default_rng(SEED)
    cr = load("cr", "26_coverage_redundancy.py")
    bm = cr.load_bm()
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    counts = counts[(counts > 0).sum(axis=1) >= 2]
    cols = LIVERWORTS + ANGIOSPERMS
    rows = []
    for ctx, filt in cr.CONTEXTS.items():
        sub = counts if filt is None else counts[counts.index.isin(set(filt.read_text().split()))]
        m = sub[cols].to_numpy()
        est = og_ratios(m)
        lo, hi = interval(lambda x: np.array(og_ratios(x)), m, rng)
        for i, name in enumerate(("coverage", "redundancy")):
            rows.append(dict(measure=name, context=ctx, units=len(sub), estimate=est[i],
                             ci_low=lo[i], ci_high=hi[i]))

    # genes per KO over the KOs of the plant-relevant modules
    ev = pd.read_csv(KM / "module_completeness.tsv", sep="\t")
    mods = pd.read_csv(KM / "modules.tsv", sep="\t")
    land = [s for s in ev.species.unique() if s != "Creinhardtii"]
    wide = ev.pivot(index="module", columns="species", values="complete")
    relevant = wide.index[wide[land].any(axis=1)]
    kos = sorted({k for d in mods[mods.module.isin(relevant)].definition for k in re.findall(r"K\d{5}", d)})
    a = pd.read_csv(KM / "ko_assignments.tsv", sep="\t")
    kc = (a[a.ko.isin(kos)].groupby(["ko", "species"]).gene.nunique().unstack(fill_value=0)
          .reindex(index=kos, columns=cols, fill_value=0))
    summ = pd.read_csv(KM / "species_summary.tsv", sep="\t", index_col=0)
    m = kc.to_numpy()
    for label, who in (("genes_per_KO Marchantia", [0]), ("genes_per_KO liverworts", [0, 1])):
        est = ko_ratio(m, who)
        check = summ.genes_per_ko[[cols[i] for i in who]].mean() / summ.genes_per_ko[ANGIOSPERMS].mean()
        assert abs(est - check) < 1e-3, (label, est, check)
        lo, hi = interval(lambda x: ko_ratio(x, who), m, rng)
        rows.append(dict(measure=label, context="KEGG modules", units=len(kos), estimate=est,
                         ci_low=lo, ci_high=hi))

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "M_bootstrap_ratios.tsv", sep="\t", index=False, float_format="%.4f")
    with pd.option_context("display.width", 160):
        print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
