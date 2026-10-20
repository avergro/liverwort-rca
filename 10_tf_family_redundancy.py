#!/usr/bin/env python3
"""
Redundancia de paralogs en familias de factores de transcripción (TF) reguladoras.

A diferencia de 06_kegg_filter.py + 05_build_matrix.py (que miden DIVERSIDAD de
funciones distintas presentes, vía MCi-RCA), este script mide EXPANSIÓN DENTRO
de cada familia: cuántas copias (paralogs) tiene cada especie en orthogroups que
pertenecen a una familia de TF dada. Es el eje de "redundancia regulatoria" que
describe la literatura (p.ej. Bowman et al. 2017, genoma de Marchantia): pocas
copias por familia = baja redundancia regulatoria.

Fuente de clasificación de familias: PlantTFDB v5 (Ath_TF_list.txt),
descargado de planttfdb.gao-lab.org y mapeado a orthogroups vía la columna
'Athaliana' de Orthogroups.tsv (mismo método que 06_kegg_filter.py).

Usage:
    python3 10_tf_family_redundancy.py [--orthofinder-dir data/orthofinder_v3]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib.machinery import SourceFileLoader
bm = SourceFileLoader("05_build_matrix", "05_build_matrix.py").load_module()
kf = SourceFileLoader("06_kegg_filter", "06_kegg_filter.py").load_module()

# Familias de TF de desarrollo/hormonales clásicamente citadas en la literatura
# de expansión génica entre briófitas y angiospermas (Lang et al. 2010;
# de Mendoza et al. 2013; Bowman et al. 2017).
FAMILIES_OF_INTEREST = [
    "ARF",         # Auxin Response Factor — señalización de auxinas
    "GRAS",        # vías de giberelinas / desarrollo
    "MIKC_MADS",   # identidad de órganos reproductivos
    "NAC",         # desarrollo / estrés
    "MYB",         # regulación amplia, secundaria incl.
    "WRKY",        # defensa / desarrollo
    "bHLH",        # red regulatoria amplia
    "AP2",         # identidad floral / estrés
    "ERF",         # señalización de etileno
    "HD-ZIP",      # desarrollo, polaridad
    "Dof",         # luz / desarrollo de semillas
    "TCP",         # arquitectura de la planta
]

LIVERWORTS  = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]


def load_tf_families(tf_list_path: Path) -> dict:
    """Returns {AT_locus: family_name}, deduplicated by gene (not isoform)."""
    df = pd.read_csv(tf_list_path, sep="\t")
    gene_to_family = {}
    for _, row in df.iterrows():
        gene_to_family[row["Gene_ID"]] = row["Family"]
    return gene_to_family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    parser.add_argument("--tf-list", default="data/tf_families/Ath_TF_list.txt", type=Path)
    parser.add_argument("--out-dir", default="data/results_tf_redundancy", type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("[1] Parsing OrthoFinder results ...")
    tsv = bm.find_orthogroups_tsv(args.orthofinder_dir)
    counts = bm.parse_orthogroups(tsv)          # OG × species, gene copy number
    gene_to_og = kf.build_ath_gene_to_og(tsv)    # AT locus -> OG
    species = counts.columns.tolist()

    print("\n[2] Loading PlantTFDB family classification ...")
    gene_to_family = load_tf_families(args.tf_list)
    print(f"  {len(gene_to_family):,} Arabidopsis TF genes classified into "
          f"{len(set(gene_to_family.values())):,} families")

    print("\n[3] Mapping TF families to orthogroups and counting paralogs per species ...")
    rows = []
    family_ogs = {}
    for family in FAMILIES_OF_INTEREST:
        genes = [g for g, f in gene_to_family.items() if f == family]
        ogs = sorted({gene_to_og[g] for g in genes if g in gene_to_og})
        family_ogs[family] = ogs
        if not ogs:
            print(f"  WARNING: no orthogroups found for family {family}")
            continue

        for sp in species:
            sp_counts = counts.loc[ogs, sp]
            n_ogs_present = int((sp_counts > 0).sum())
            total_copies = int(sp_counts.sum())
            rows.append({
                "family": family,
                "species": sp,
                "n_ath_genes": len(genes),
                "n_OGs_total": len(ogs),
                "n_OGs_present": n_ogs_present,
                "total_copies": total_copies,
                "copies_per_OG_present": total_copies / n_ogs_present if n_ogs_present else 0.0,
            })

    results = pd.DataFrame(rows)
    results.to_csv(args.out_dir / "tf_family_paralogs.tsv", sep="\t", index=False)

    print("\n" + "═" * 78)
    print("  Copias totales (paralogs) por familia y especie")
    print("═" * 78)
    pivot_copies = results.pivot(index="family", columns="species", values="total_copies")
    pivot_copies = pivot_copies[species]
    print(pivot_copies.to_string())
    pivot_copies.to_csv(args.out_dir / "pivot_total_copies.tsv", sep="\t")

    print("\n" + "═" * 78)
    print("  Hepáticas vs. Angiospermas — copias totales por familia")
    print("═" * 78)
    print(f"  {'Familia':<14} {'Hepáticas (media)':>18} {'Angiospermas (media)':>22} {'Ratio L/A':>10}")
    print("  " + "─" * 68)
    summary_rows = []
    for family in FAMILIES_OF_INTEREST:
        if family not in pivot_copies.index:
            continue
        lw_mean = pivot_copies.loc[family, LIVERWORTS].mean()
        ang_mean = pivot_copies.loc[family, ANGIOSPERMS].mean()
        ratio = lw_mean / ang_mean if ang_mean > 0 else float("nan")
        print(f"  {family:<14} {lw_mean:>18.1f} {ang_mean:>22.1f} {ratio:>10.2f}")
        summary_rows.append({
            "family": family, "liverwort_mean_copies": lw_mean,
            "angiosperm_mean_copies": ang_mean, "ratio_L_A": ratio,
        })

    summary_df = pd.DataFrame(summary_rows).sort_values("ratio_L_A")
    summary_df.to_csv(args.out_dir / "liverwort_vs_angiosperm_summary.tsv", sep="\t", index=False)

    overall_lw = summary_df["liverwort_mean_copies"].sum()
    overall_ang = summary_df["angiosperm_mean_copies"].sum()
    print("\n" + "─" * 68)
    print(f"  TOTAL (12 familias): hepáticas = {overall_lw:.0f} copias  ·  "
          f"angiospermas = {overall_ang:.0f} copias  ·  ratio = {overall_lw/overall_ang:.2f}")

    print(f"\nOutputs guardados en: {args.out_dir}/")


if __name__ == "__main__":
    main()
