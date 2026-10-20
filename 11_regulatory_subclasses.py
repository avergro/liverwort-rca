#!/usr/bin/env python3
"""
MCi-RCA por clase funcional regulatoria individual (no agregada).

06_kegg_filter.py agrupa 5 pathways KEGG en un solo filtro "regulatory_networks".
Este script calcula el MCi-RCA por separado para cada uno de los 5, igual que
ya se hace con C_metabolism vs N_metabolism (en vez de un solo "metabolism").

Usage:
    python3 11_regulatory_subclasses.py [--orthofinder-dir data/orthofinder_v3]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib.machinery import SourceFileLoader
bm = SourceFileLoader("05_build_matrix", "05_build_matrix.py").load_module()
kf = SourceFileLoader("06_kegg_filter", "06_kegg_filter.py").load_module()

REGULATORY_SUBCLASSES = {
    "ath04075": "Plant hormone signal transduction",
    "ath04016": "MAPK signaling - plant",
    "ath04712": "Circadian rhythm - plant",
    "ath04626": "Plant-pathogen interaction",
    "ath04070": "Phosphatidylinositol signaling",
}

LIVERWORTS  = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    parser.add_argument("--out-dir", default="data/results_regulatory_subclasses", type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("[1] Parsing OrthoFinder results ...")
    tsv = bm.find_orthogroups_tsv(args.orthofinder_dir)
    counts = bm.parse_orthogroups(tsv)
    gene_to_og = kf.build_ath_gene_to_og(tsv)

    print("\n[2] Downloading KEGG gene-pathway links ...")
    gene_to_pathways = kf.download_kegg_gene_pathway_links("ath")

    print("\n[3] Computing MCi-RCA per regulatory subclass ...")
    all_mci = {}
    og_counts = {}

    for pathway_id, label in REGULATORY_SUBCLASSES.items():
        og_set = kf.og_set_for_pathways(gene_to_og, gene_to_pathways, {pathway_id})
        og_counts[label] = len(og_set)

        sub_counts = counts.loc[counts.index.isin(og_set)]
        M_rca = bm.build_rca_matrix(sub_counts)
        M_rca = bm.filter_matrix(M_rca, min_species=2)

        if M_rca.shape[0] < 3:
            print(f"  SKIP {label}: solo {M_rca.shape[0]} OGs informativos tras filtro")
            continue

        MCi, _ = bm.compute_complexity(M_rca)
        all_mci[label] = MCi

        out_path = args.out_dir / f"species_complexity_{pathway_id}.tsv"
        MCi.sort_values(ascending=False).to_csv(out_path, sep="\t", header=["MCi_rca"])

    result_df = pd.DataFrame(all_mci)
    result_df.to_csv(args.out_dir / "mci_by_regulatory_subclass.tsv", sep="\t")

    print("\n" + "═" * 90)
    print("  MCi-RCA por clase funcional regulatoria (columnas = pathway)")
    print("═" * 90)
    print(f"  OGs por clase: " + "  ".join(f"{l}={n}" for l, n in og_counts.items()))
    print()
    print(result_df.round(3).to_string())

    print("\n" + "═" * 90)
    print("  Hepáticas vs. Angiospermas — MCi-RCA medio por clase regulatoria")
    print("═" * 90)
    print(f"  {'Clase funcional':<35} {'Hepáticas':>12} {'Angiospermas':>14} {'Hep > Ang?':>12}")
    print("  " + "─" * 78)
    for label in result_df.columns:
        lw_mean = result_df.loc[LIVERWORTS, label].mean()
        ang_mean = result_df.loc[ANGIOSPERMS, label].mean()
        print(f"  {label:<35} {lw_mean:>12.3f} {ang_mean:>14.3f} {'sí' if lw_mean > ang_mean else 'no':>12}")

    print(f"\nOutputs guardados en: {args.out_dir}/")


if __name__ == "__main__":
    main()
