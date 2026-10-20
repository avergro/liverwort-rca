#!/usr/bin/env python3
"""Economic complexity on the species x KEGG module matrix.

The KEGG module is the closest biological equivalent of a product in the
capability model: a module is complete only when every step (KO) has a gene,
as a product is made only when every required capability is held. This script
applies the Economic Complexity Index (ECI, second eigenvector, 05_build_matrix)
and the Fitness-Complexity metric to that matrix, from the uniform KOfam
annotation of 33_kegg_module_completeness.py.

Variants
  complete      M = 1 when the module is complete (presence, no RCA: a module has no volume)
  near          M = 1 when at most one step is missing
  rca_genes     genes assigned to the KOs of each module, RCA >= 1 (the orthogroup analogue,
                where copy number enters)
each on all 11 species and on the 10 land plants. Modules complete (or with RCA >= 1)
in no species are dropped. Ranks are compared with diversity, BUSCO completeness and
gene number, and the rank of every species is bootstrapped over modules (2,000, seed 42).

Usage: python3 34_eci_kegg_modules.py
Outputs: data/results_kegg_modules/eci_*.tsv
"""
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
KM = ROOT / "data" / "results_kegg_modules"
REPS, SEED = 2000, 42


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, ROOT / file)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


bm = load("bm", "05_build_matrix.py")
fc = load("fc", "15_threshold_and_baselines.py")


def scores(M):
    eci, pci = bm.compute_complexity(M)
    return pd.DataFrame({"diversity": M.sum(axis=0), "ECI": eci,
                         "fitness": fc.fitness_complexity(M)}), pci


def main():
    ev = pd.read_csv(KM / "module_completeness.tsv", sep="\t")
    mods = pd.read_csv(KM / "modules.tsv", sep="\t").set_index("module")
    ko = pd.read_csv(KM / "ko_assignments.tsv", sep="\t")
    s1 = pd.read_csv(ROOT / "data/busco/busco_summary.tsv", sep="\t")      # 03b_busco_summary.py
    short = {"Mpolymorpha": "Marchantia polymorpha", "Lcruciata": "Lunularia cruciata",
             "Smoellendorffii": "Selaginella moellendorffii", "Ppatens": "Physcomitrium patens",
             "Afiliculoides": "Azolla filiculoides", "Scucullata": "Salvinia cucullata",
             "Creinhardtii": "Chlamydomonas reinhardtii", "Athaliana": "Arabidopsis thaliana",
             "Slycopersicum": "Solanum lycopersicum", "Pamericana": "Persea americana",
             "Pabies": "Picea abies"}
    s1 = s1.set_index(s1.Species.map({v: k for k, v in short.items()}))
    busco = s1["C (%)"]
    genes = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3")).sum(axis=0)

    complete = ev.pivot(index="module", columns="species", values="complete").astype(int)
    near = ev.pivot(index="module", columns="species", values="near_complete").astype(int)
    # genes per module: genes carrying any KO of the module definition
    kos = {m: set(re.findall(r"K\d{5}", d)) for m, d in mods.definition.items()}
    k2m = {}
    for m, ks in kos.items():
        for k in ks:
            k2m.setdefault(k, []).append(m)
    rows = [(r.species, r.gene, m) for r in ko.itertuples() for m in k2m.get(r.ko, ())]
    gm = pd.DataFrame(rows, columns=["species", "gene", "module"]).drop_duplicates()
    counts = gm.groupby(["module", "species"]).gene.nunique().unstack(fill_value=0)
    counts = counts.reindex(index=complete.index, columns=complete.columns, fill_value=0)

    variants = {"complete": complete, "near": near}
    out, pcis, boots = [], [], []
    rng = np.random.default_rng(SEED)
    for sp_set, cols in [("all", list(complete.columns)),
                         ("land", [c for c in complete.columns if c != "Creinhardtii"])]:
        mats = {k: v[cols] for k, v in variants.items()}
        c = counts[cols]
        c = c[c.sum(axis=1) > 0]
        mats["rca_genes"] = bm.build_rca_matrix(c)
        for name, M in mats.items():
            M = M[M.sum(axis=1) > 0]
            sc, pci = scores(M)
            sc["rank_ECI"] = sc.ECI.rank(ascending=False).astype(int)
            sc["rank_fitness"] = sc.fitness.rank(ascending=False).astype(int)
            sc["BUSCO"] = busco.reindex(sc.index)
            sc["genes"] = genes.reindex(sc.index)
            # bootstrap over modules
            br = []
            for _ in range(REPS):
                Mb = M.iloc[rng.integers(0, len(M), len(M))]
                Mb = Mb[(Mb.sum(axis=1) > 0)]
                if (Mb.sum(axis=0) == 0).any():
                    continue
                br.append(bm.compute_complexity(Mb)[0].rank(ascending=False))
            br = pd.DataFrame(br)
            sc["rank_ECI_lo"] = br.quantile(0.025).reindex(sc.index)
            sc["rank_ECI_hi"] = br.quantile(0.975).reindex(sc.index)
            sc.insert(0, "variant", name); sc.insert(0, "species_set", sp_set)
            sc["n_modules"] = len(M)
            out.append(sc.reset_index(names="species"))
            land = sc.index != "Creinhardtii"
            r = lambda a, b: spearmanr(sc.loc[land, a], sc.loc[land, b])[0]
            print(f"\n== {sp_set} / {name}: {len(M)} modules ==")
            print(sc.sort_values("ECI", ascending=False)[
                ["diversity", "ECI", "rank_ECI", "rank_ECI_lo", "rank_ECI_hi", "fitness", "rank_fitness",
                 "BUSCO"]].round(3).to_string())
            print(f"land plants, Spearman: ECI~diversity {r('ECI', 'diversity'):.2f}, "
                  f"ECI~BUSCO {r('ECI', 'BUSCO'):.2f}, ECI~genes {r('ECI', 'genes'):.2f}, "
                  f"fitness~BUSCO {r('fitness', 'BUSCO'):.2f}, fitness~genes {r('fitness', 'genes'):.2f}")
            p = pci.rename("PCI").to_frame().join(mods["name"])
            p["ubiquity"] = M.sum(axis=1)
            p.insert(0, "variant", name); p.insert(0, "species_set", sp_set)
            pcis.append(p.reset_index(names="module"))
    pd.concat(out).to_csv(KM / "eci_species.tsv", sep="\t", index=False, float_format="%.4g")
    pd.concat(pcis).to_csv(KM / "eci_modules.tsv", sep="\t", index=False, float_format="%.4g")
    top = pd.concat(pcis).query("species_set == 'land' and variant == 'complete'").sort_values("PCI")
    with pd.option_context("display.width", 200, "display.max_colwidth", 70):
        print("\nland/complete, most and least complex modules:")
        print(top.tail(12)[["module", "name", "ubiquity", "PCI"]].to_string(index=False))
        print(top.head(6)[["module", "name", "ubiquity", "PCI"]].to_string(index=False))


if __name__ == "__main__":
    main()
