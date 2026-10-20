#!/usr/bin/env python3
"""
Gene family expansion/contraction analysis from OrthoFinder output.

Uses OrthoFinder's Gene_Duplication_Events to count duplications per
KEGG functional category, comparing liverworts vs angiosperms.

Two analyses:
  1. Per-species duplication counts, normalized by total OGs per category
  2. Duplication density: duplications per OG per species, by KEGG category

Usage:
    python3 09c_expansion_contraction.py [--orthofinder-dir data/orthofinder_v3]
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib.machinery import SourceFileLoader
bm = SourceFileLoader("05_build_matrix", "05_build_matrix.py").load_module()


def load_kegg_filter(filter_file: Path) -> set:
    """Load orthogroup IDs from KEGG filter file."""
    return set(filter_file.read_text().split())


def parse_duplications_per_node(ortho_dir: Path) -> pd.DataFrame:
    """Parse Duplications_per_Species_Tree_Node.tsv"""
    results_dir = sorted(ortho_dir.glob("Results_*"))[-1]
    path = results_dir / "Comparative_Genomics_Statistics" / "Duplications_per_Species_Tree_Node.tsv"
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df


def parse_duplications_per_og(ortho_dir: Path) -> pd.DataFrame:
    """Parse Duplications_per_Orthogroup.tsv"""
    results_dir = sorted(ortho_dir.glob("Results_*"))[-1]
    path = results_dir / "Comparative_Genomics_Statistics" / "Duplications_per_Orthogroup.tsv"
    return pd.read_csv(path, sep="\t", index_col=0)


def map_duplications_to_category(dups_per_og: pd.DataFrame,
                                  counts: pd.DataFrame,
                                  filter_ogs: dict) -> dict:
    """
    For each species and KEGG category, compute:
      - total duplications in that category
      - duplications / OG in that category
    """
    species = counts.columns.tolist()
    results = {}
    cat_ogs_map = {}

    for cat_name, og_set in filter_ogs.items():
        cat_ogs = [og for og in og_set if og in dups_per_og.index]
        if not cat_ogs:
            continue
        cat_ogs_map[cat_name] = cat_ogs

        cat_dups = dups_per_og.loc[cat_ogs, "Duplications (all)"]

        for sp in species:
            # Count how many of these OGs actually have genes in this species
            sp_counts = counts[sp]
            cat_present = [og for og in cat_ogs
                          if og in sp_counts.index and sp_counts[og] > 0]
            n_present = len(cat_present)

            # Total duplications in these OGs
            total_dup = int(cat_dups[cat_present].sum()) if cat_present else 0

            results[(cat_name, sp)] = {
                "n_OGs_present": n_present,
                "total_duplications": total_dup,
                "duplications_per_OG": total_dup / n_present if n_present > 0 else 0,
            }

    return results, cat_ogs_map


def parse_duplications_by_species(ortho_dir: Path, species_list: list) -> pd.DataFrame:
    """
    From Duplications.tsv, count species-specific (Terminal) duplications
    per species by scanning the gene IDs in each duplication event.
    """
    results_dir = sorted(ortho_dir.glob("Results_*"))[-1]
    path = results_dir / "Gene_Duplication_Events" / "Duplications.tsv"
    sp_dup_counts = {sp: 0 for sp in species_list}

    # Species prefixes for identification in gene IDs
    sp_prefixes = {
        "Afiliculoides": "Azfi",
        "Athaliana": "AT",
        "Creinhardtii": "Cre",
        "Lcruciata": "Lc_female",
        "Mpolymorpha": "Mp",
        "Pabies": "Pab",
        "Pamericana": "KAJ",
        "Ppatens": "Pp",
        "Scucullata": "Sacu",
        "Slycopersicum": "Solyc",
        "Smoellendorffii": "Smo",
    }

    with open(path) as f:
        for line in f:
            if line.startswith("Orthogroup"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            og, species_node, gene_node, support, dup_type = parts[0], parts[1], parts[2], parts[3], parts[4]

            if dup_type != "Terminal":
                continue

            # Check if the node is a species (not an internal N_ node)
            if species_node not in species_list:
                continue

            sp_dup_counts[species_node] += 1

    return pd.Series(sp_dup_counts, name="terminal_duplications")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    parser.add_argument("--out-dir", default="data/results_expansion_contraction", type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ortho_dir = args.orthofinder_dir

    # Load data
    print("Loading OrthoFinder data...")
    tsv = bm.find_orthogroups_tsv(ortho_dir)
    counts = bm.parse_orthogroups(tsv)
    species = counts.columns.tolist()

    # Duplications per orthogroup (from Comparative Genomics)
    dups_per_og = parse_duplications_per_og(ortho_dir)
    print(f"  {len(dups_per_og):,} OGs with duplication data")

    # Duplications per species tree node
    dups_per_node = parse_duplications_per_node(ortho_dir)

    # Terminal duplications per species
    terminal_dups = parse_duplications_by_species(ortho_dir, species)

    # ── Load KEGG filter sets ──
    filter_dir = Path("data/filters")
    filter_ogs = {}
    filter_labels = {
        "primary_metabolism_ogs": "Primary metabolism",
        "C_metabolism_ogs": "C metabolism",
        "N_metabolism_ogs": "N metabolism",
        "regulatory_networks_ogs": "Regulatory networks",
    }

    for fname, label in filter_labels.items():
        fpath = filter_dir / f"{fname}.txt"
        if fpath.exists():
            filter_ogs[label] = load_kegg_filter(fpath)
            print(f"  {label}: {len(filter_ogs[label]):,} OGs")

    # ── 1. Total duplications per species × category ──
    print("\n\n" + "═" * 70)
    print("  Gene Family Expansion: Duplications per KEGG Category")
    print("═" * 70)

    results, cat_ogs_map = map_duplications_to_category(dups_per_og, counts, filter_ogs)

    # Build pivot table
    for cat_name in filter_labels.values():
        if cat_name not in cat_ogs_map:
            continue

        rows = []
        for sp in species:
            if (cat_name, sp) in results:
                rows.append({
                    "species": sp,
                    "n_OGs_present": results[(cat_name, sp)]["n_OGs_present"],
                    "total_dup": results[(cat_name, sp)]["total_duplications"],
                    "dup_per_OG": results[(cat_name, sp)]["duplications_per_OG"],
                })

        df = pd.DataFrame(rows).sort_values("dup_per_OG", ascending=False)

        print(f"\n  ── {cat_name} ──")
        print(f"  {'Species':<25} {'OGs present':>12} {'Total dups':>12} {'Dups/OG':>8}")
        print(f"  {'─'*60}")
        for _, row in df.iterrows():
            print(f"  {row['species']:<25} {int(row['n_OGs_present']):>12,} "
                  f"{int(row['total_dup']):>12,} {row['dup_per_OG']:>8.2f}")

        df.to_csv(args.out_dir / f"duplications_{cat_name.replace(' ', '_')}.tsv",
                  sep="\t", index=False)

    # ── 2. Duplications per species tree node ──
    print("\n\n" + "═" * 70)
    print("  Duplications per Species Tree Node (from OrthoFinder)")
    print("═" * 70)
    print(dups_per_node.to_string())

    # ── 3. Terminal duplications per species ──
    print("\n\n" + "═" * 70)
    print("  Terminal Duplications per Species")
    print("═" * 70)
    td_df = terminal_dups.sort_values(ascending=False)
    for sp, val in td_df.items():
        print(f"  {sp:<25} {val:>8,}")

    # ── 4. Liverwort vs Angiosperm comparison ──
    print("\n\n" + "═" * 70)
    print("  Liverwort vs Angiosperm Duplication Density Comparison")
    print("═" * 70)

    liverworts = ["Mpolymorpha", "Lcruciata"]
    angiosperms = ["Athaliana", "Slycopersicum", "Pamericana"]

    for cat_name in filter_labels.values():
        if cat_name not in cat_ogs_map:
            continue

        lw_vals = []
        ang_vals = []
        for sp in liverworts:
            if (cat_name, sp) in results:
                lw_vals.append(results[(cat_name, sp)]["duplications_per_OG"])
        for sp in angiosperms:
            if (cat_name, sp) in results:
                ang_vals.append(results[(cat_name, sp)]["duplications_per_OG"])

        if lw_vals and ang_vals:
            print(f"\n  {cat_name}:")
            print(f"    Liverworts:   {np.mean(lw_vals):.2f} dups/OG (values: {[f'{v:.2f}' for v in lw_vals]})")
            print(f"    Angiosperms:  {np.mean(ang_vals):.2f} dups/OG (values: {[f'{v:.2f}' for v in ang_vals]})")
            print(f"    Ratio L/A:    {np.mean(lw_vals)/np.mean(ang_vals):.2f}")

    # ── 5. Duplications in liverwort-specific OGs ──
    print("\n\n" + "═" * 70)
    print("  Duplications in Filter B OGs (liverwort-enriched, angiosperm-absent)")
    print("═" * 70)

    for cat_abbr, cat_name in [("N_metabolism", "N_metabolism_ogs"),
                                ("C_metabolism", "C_metabolism_ogs")]:
        filter_b_path = Path(f"gene_lists/{cat_abbr}/02_OGs_both_liverworts.tsv")
        if not filter_b_path.exists():
            continue

        filter_b_ogs = set(pd.read_csv(filter_b_path, sep="\t")["Orthogroup"])
        fb_with_dup = [og for og in filter_b_ogs if og in dups_per_og.index]

        if fb_with_dup:
            fb_dups = dups_per_og.loc[fb_with_dup, "Duplications (all)"]
            total = int(fb_dups.sum())
            mean_dup = float(fb_dups.mean()) if len(fb_dups) > 0 else 0

            # Compare with all OGs in same category
            cat_ogs_for_comparison = filter_ogs.get(
                {"N_metabolism_ogs": "N metabolism",
                 "C_metabolism_ogs": "C metabolism"}.get(cat_name, ""),
                set()
            )
            cat_ogs_with_dup = [og for og in cat_ogs_for_comparison if og in dups_per_og.index]
            if cat_ogs_with_dup:
                cat_dups = dups_per_og.loc[cat_ogs_with_dup, "Duplications (all)"]
                cat_mean = float(cat_dups.mean())
            else:
                cat_mean = 0

            print(f"\n  {cat_abbr} Filter B OGs:")
            print(f"    OGs:              {len(fb_with_dup)}")
            print(f"    Total dups:       {total:,}")
            print(f"    Mean dups/OG:     {mean_dup:.2f}")
            print(f"    Category mean:    {cat_mean:.2f}")
            print(f"    Filter B / all:   {mean_dup/cat_mean:.2f}" if cat_mean > 0 else "")

    # Save summary
    terminal_dups.to_csv(args.out_dir / "terminal_duplications_per_species.tsv", sep="\t")
    dups_per_node.to_csv(args.out_dir / "duplications_per_node.tsv", sep="\t")

    print(f"\nOutputs saved to: {args.out_dir}/")


if __name__ == "__main__":
    import numpy as np
    main()
