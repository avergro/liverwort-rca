#!/usr/bin/env python3
"""
Download KEGG pathway-gene associations for Arabidopsis thaliana (ath)
and map them to OrthoFinder orthogroups.

Generates three OG filter files:
  data/filters/primary_metabolism_ogs.txt  — all KEGG metabolic pathways
  data/filters/C_metabolism_ogs.txt        — carbohydrate + energy pathways
  data/filters/N_metabolism_ogs.txt        — amino acid + nitrogen assimilation

Usage:
    python3 06_kegg_filter.py [--orthofinder-dir data/orthofinder_v2]
"""

import argparse
import glob
import sys
import time
import urllib.request
from pathlib import Path
from collections import defaultdict

import pandas as pd


# ── KEGG pathway category definitions ─────────────────────────────────────
# Based on KEGG BRITE hierarchy for Arabidopsis thaliana (ath)
# Pathway IDs follow: ath00010 = Glycolysis, etc.

# Carbohydrate metabolism (KEGG category 1.1)
C_PATHWAYS_CARBOHYDRATE = {
    "ath00010",  # Glycolysis / Gluconeogenesis
    "ath00020",  # Citrate cycle (TCA cycle)
    "ath00030",  # Pentose phosphate pathway
    "ath00040",  # Pentose and glucuronate interconversions
    "ath00051",  # Fructose and mannose metabolism
    "ath00052",  # Galactose metabolism
    "ath00053",  # Ascorbate and aldarate metabolism
    "ath00500",  # Starch and sucrose metabolism
    "ath00520",  # Amino sugar and nucleotide sugar metabolism
    "ath00562",  # Inositol phosphate metabolism
    "ath00620",  # Pyruvate metabolism
    "ath00630",  # Glyoxylate and dicarboxylate metabolism
    "ath00640",  # Propanoate metabolism
    "ath00650",  # Butanoate metabolism
    "ath00660",  # C5-Branched dibasic acid metabolism
}

# Energy metabolism (KEGG category 1.2) — carbon fixation included here
C_PATHWAYS_ENERGY = {
    "ath00190",  # Oxidative phosphorylation
    "ath00195",  # Photosynthesis
    "ath00196",  # Photosynthesis - antenna proteins
    "ath00710",  # Carbon fixation in photosynthetic organisms
    "ath00720",  # Carbon fixation pathways in prokaryotes
}

C_METABOLISM_PATHWAYS = C_PATHWAYS_CARBOHYDRATE | C_PATHWAYS_ENERGY

# Amino acid metabolism (KEGG categories 1.5 + 1.6)
N_PATHWAYS_AMINO_ACIDS = {
    "ath00220",  # Arginine biosynthesis
    "ath00250",  # Alanine, aspartate and glutamate metabolism
    "ath00260",  # Glycine, serine and threonine metabolism
    "ath00270",  # Cysteine and methionine metabolism
    "ath00280",  # Valine, leucine and isoleucine degradation
    "ath00290",  # Valine, leucine and isoleucine biosynthesis
    "ath00300",  # Lysine biosynthesis
    "ath00310",  # Lysine degradation
    "ath00330",  # Arginine and proline metabolism
    "ath00340",  # Histidine metabolism
    "ath00350",  # Tyrosine metabolism
    "ath00360",  # Phenylalanine metabolism
    "ath00380",  # Tryptophan metabolism
    "ath00400",  # Phenylalanine, tyrosine and tryptophan biosynthesis
    "ath00410",  # beta-Alanine metabolism
    "ath00430",  # Taurine and hypotaurine metabolism
    "ath00460",  # Cyanoamino acid metabolism
    "ath00470",  # D-Amino acid metabolism
    "ath00480",  # Glutathione metabolism
}

# Nitrogen assimilation (GS/GOGAT, nitrate reductase)
N_PATHWAYS_ASSIMILATION = {
    "ath00910",  # Nitrogen metabolism
    "ath00910",  # (included once)
}

N_METABOLISM_PATHWAYS = N_PATHWAYS_AMINO_ACIDS | N_PATHWAYS_ASSIMILATION

# Primary metabolism = all KEGG "Metabolism" pathways
# Pathway numbers 00010–00999 are Metabolism; 01100+ are overview/global
# We include C + N + lipid + nucleotide + cofactors + secondary metabolism
LIPID_PATHWAYS = {
    "ath00061",  # Fatty acid biosynthesis
    "ath00062",  # Fatty acid elongation
    "ath00071",  # Fatty acid degradation
    "ath00072",  # Synthesis and degradation of ketone bodies
    "ath00100",  # Steroid biosynthesis
    "ath00120",  # Primary bile acid biosynthesis
    "ath00121",  # Secondary bile acid biosynthesis
    "ath00130",  # Ubiquinone and other terpenoid-quinone biosynthesis
    "ath00561",  # Glycerolipid metabolism
    "ath00564",  # Glycerophospholipid metabolism
    "ath00565",  # Ether lipid metabolism
    "ath00590",  # Arachidonic acid metabolism
    "ath00591",  # Linoleic acid metabolism
    "ath00592",  # alpha-Linolenic acid metabolism
    "ath00600",  # Sphingolipid metabolism
}

NUCLEOTIDE_PATHWAYS = {
    "ath00230",  # Purine metabolism
    "ath00240",  # Pyrimidine metabolism
}

COFACTOR_PATHWAYS = {
    "ath00730",  # Thiamine metabolism
    "ath00740",  # Riboflavin metabolism
    "ath00750",  # Vitamin B6 metabolism
    "ath00760",  # Nicotinate and nicotinamide metabolism
    "ath00770",  # Pantothenate and CoA biosynthesis
    "ath00780",  # Biotin metabolism
    "ath00785",  # Lipoic acid metabolism
    "ath00790",  # Folate biosynthesis
    "ath00860",  # Porphyrin metabolism
    "ath00130",  # Ubiquinone biosynthesis (already in lipids)
}

PRIMARY_METABOLISM_PATHWAYS = (
    C_METABOLISM_PATHWAYS
    | N_METABOLISM_PATHWAYS
    | LIPID_PATHWAYS
    | NUCLEOTIDE_PATHWAYS
    | COFACTOR_PATHWAYS
)

# Regulatory networks / signal transduction (KEGG "Environmental Information
# Processing" supercategory) — the conceptual opposite of Metabolism, used to
# test the "simple in regulation, complex in metabolism" hypothesis.
REGULATORY_NETWORK_PATHWAYS = {
    "ath04075",  # Plant hormone signal transduction
    "ath04016",  # MAPK signaling pathway - plant
    "ath04712",  # Circadian rhythm - plant
    "ath04626",  # Plant-pathogen interaction (receptor-mediated immune signaling)
    "ath04070",  # Phosphatidylinositol signaling system
}


# ── KEGG REST download helpers ─────────────────────────────────────────────

def kegg_get(url: str, retries: int = 3, delay: float = 1.0) -> str:
    """Download text from KEGG REST API with retry."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                return resp.read().decode("utf-8")
        except Exception as e:
            if attempt < retries - 1:
                print(f"  Retry {attempt+1}/{retries-1} after error: {e}")
                time.sleep(delay)
            else:
                raise


def download_kegg_gene_pathway_links(organism: str = "ath") -> dict[str, set[str]]:
    """
    Download gene→pathway links from KEGG REST API.
    Returns: {gene_locus (e.g. 'AT1G01010'): {pathway_id, ...}}
    """
    url = f"https://rest.kegg.jp/link/pathway/{organism}"
    print(f"  Downloading KEGG gene-pathway links from {url} ...")
    text = kegg_get(url)

    gene_to_pathways: dict[str, set[str]] = defaultdict(set)
    lines_parsed = 0
    for line in text.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        gene_raw, path_raw = parts
        # gene_raw: "ath:AT1G01010"  →  "AT1G01010"
        gene_id = gene_raw.split(":")[-1].strip()
        # path_raw: "path:ath00010"  →  "ath00010"
        pathway_id = path_raw.split(":")[-1].strip()
        gene_to_pathways[gene_id].add(pathway_id)
        lines_parsed += 1

    print(f"  Parsed {lines_parsed:,} gene-pathway associations for {len(gene_to_pathways):,} genes")
    return dict(gene_to_pathways)


# ── OrthoFinder parsing ────────────────────────────────────────────────────

def find_orthogroups_tsv(base_dir: Path) -> Path:
    patterns = [
        base_dir / "Results_*" / "Orthogroups" / "Orthogroups.tsv",
        base_dir / "OrthoFinder" / "Results_*" / "Orthogroups" / "Orthogroups.tsv",
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(str(pattern)))
        if matches:
            return Path(matches[-1])
    sys.exit(f"ERROR: Could not find Orthogroups.tsv under {base_dir}")


def build_ath_gene_to_og(tsv: Path) -> dict[str, str]:
    """
    Parse OrthoFinder TSV and return {arabidopsis_locus: orthogroup_id}.
    Strips isoform suffix: AT1G01010.1 → AT1G01010
    """
    raw = pd.read_csv(tsv, sep="\t", index_col=0)

    # Find the Arabidopsis column (named 'Athaliana' or contains 'thaliana')
    ath_col = None
    for col in raw.columns:
        if "thaliana" in col.lower() or col.lower().startswith("at"):
            ath_col = col
            break
    if ath_col is None:
        sys.exit("ERROR: Could not find Arabidopsis column in Orthogroups.tsv")

    print(f"  Using column '{ath_col}' for Arabidopsis gene IDs")

    gene_to_og: dict[str, str] = {}
    for og_id, row in raw.iterrows():
        cell = row[ath_col]
        if pd.isna(cell) or str(cell).strip() == "":
            continue
        for gene_full in str(cell).split(","):
            gene_full = gene_full.strip()
            if not gene_full:
                continue
            # Strip isoform suffix: AT1G01010.1 → AT1G01010
            locus = gene_full.rsplit(".", 1)[0]
            gene_to_og[locus] = og_id

    print(f"  Mapped {len(gene_to_og):,} Arabidopsis gene loci to orthogroups")
    return gene_to_og


# ── Filter logic ───────────────────────────────────────────────────────────

def og_set_for_pathways(
    gene_to_og: dict[str, str],
    gene_to_pathways: dict[str, set[str]],
    target_pathways: set[str],
) -> set[str]:
    """
    Return orthogroup IDs where at least one Arabidopsis gene belongs
    to any pathway in target_pathways.
    """
    og_ids: set[str] = set()
    for locus, og in gene_to_og.items():
        pathways = gene_to_pathways.get(locus, set())
        if pathways & target_pathways:
            og_ids.add(og)
    return og_ids


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v2", type=Path)
    parser.add_argument("--out-dir", default="data/filters", type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Parse OrthoFinder
    print("\n[1] Parsing OrthoFinder results ...")
    tsv = find_orthogroups_tsv(args.orthofinder_dir)
    print(f"  Reading: {tsv}")
    gene_to_og = build_ath_gene_to_og(tsv)

    # 2. Download KEGG
    print("\n[2] Downloading KEGG pathway annotations ...")
    gene_to_pathways = download_kegg_gene_pathway_links("ath")

    # Coverage check
    loci_in_kegg = set(gene_to_og.keys()) & set(gene_to_pathways.keys())
    print(f"  Arabidopsis loci in OrthoFinder: {len(gene_to_og):,}")
    print(f"  Arabidopsis loci with KEGG annotations: {len(gene_to_pathways):,}")
    print(f"  Overlap (usable for mapping): {len(loci_in_kegg):,}")

    # 3. Build OG filter sets
    print("\n[3] Building OG filter sets ...")

    og_primary = og_set_for_pathways(gene_to_og, gene_to_pathways, PRIMARY_METABOLISM_PATHWAYS)
    og_C       = og_set_for_pathways(gene_to_og, gene_to_pathways, C_METABOLISM_PATHWAYS)
    og_N       = og_set_for_pathways(gene_to_og, gene_to_pathways, N_METABOLISM_PATHWAYS)
    og_reg     = og_set_for_pathways(gene_to_og, gene_to_pathways, REGULATORY_NETWORK_PATHWAYS)

    total_ogs = len(set(gene_to_og.values()))
    print(f"\n  Total OGs with ≥1 Arabidopsis gene: {total_ogs:,}")
    print(f"  Primary metabolism OGs:              {len(og_primary):,}  ({100*len(og_primary)/total_ogs:.1f}%)")
    print(f"  C metabolism OGs:                    {len(og_C):,}  ({100*len(og_C)/total_ogs:.1f}%)")
    print(f"  N metabolism OGs:                    {len(og_N):,}  ({100*len(og_N)/total_ogs:.1f}%)")
    print(f"  Regulatory network OGs:               {len(og_reg):,}  ({100*len(og_reg)/total_ogs:.1f}%)")
    print(f"  C ∩ N overlap:                       {len(og_C & og_N):,} OGs")
    print(f"  Regulatory ∩ primary overlap:         {len(og_reg & og_primary):,} OGs")

    # 4. Write filter files
    print(f"\n[4] Writing filter files to {args.out_dir}/ ...")

    filters = {
        "primary_metabolism": og_primary,
        "C_metabolism": og_C,
        "N_metabolism": og_N,
        "regulatory_networks": og_reg,
    }

    for name, og_set in filters.items():
        out_path = args.out_dir / f"{name}_ogs.txt"
        sorted_ogs = sorted(og_set)
        out_path.write_text("\n".join(sorted_ogs) + "\n")
        print(f"  {out_path.name}: {len(sorted_ogs):,} orthogroups")

    # Save pathway coverage summary
    summary_rows = []
    for pathway_id in sorted(PRIMARY_METABOLISM_PATHWAYS | REGULATORY_NETWORK_PATHWAYS):
        og_count = len(og_set_for_pathways(gene_to_og, gene_to_pathways, {pathway_id}))
        gene_count = sum(1 for locus in gene_to_og if pathway_id in gene_to_pathways.get(locus, set()))
        category = (
            "C_metabolism" if pathway_id in C_METABOLISM_PATHWAYS
            else "N_metabolism" if pathway_id in N_METABOLISM_PATHWAYS
            else "regulatory_networks" if pathway_id in REGULATORY_NETWORK_PATHWAYS
            else "other_primary"
        )
        summary_rows.append({
            "pathway_id": pathway_id,
            "category": category,
            "ath_genes_in_OGs": gene_count,
            "orthogroups": og_count,
        })

    summary_df = pd.DataFrame(summary_rows).sort_values(["category", "pathway_id"])
    summary_path = args.out_dir / "pathway_coverage.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"  pathway_coverage.tsv: pathway-level summary")

    print("\nDone. Run 05_build_matrix.py with --og-filter for each subset.")
    print(f"  Example:")
    print(f"    python3 05_build_matrix.py --orthofinder-dir data/orthofinder_v2 \\")
    print(f"        --og-filter data/filters/C_metabolism_ogs.txt \\")
    print(f"        --out-dir data/results_C_metabolism")


if __name__ == "__main__":
    main()
