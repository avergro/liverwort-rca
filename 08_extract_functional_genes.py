#!/usr/bin/env python3
"""
08_extract_functional_genes.py

For each functional context (primary_metabolism, C_metabolism, N_metabolism):

  1. Load MCi-RCA matrix (presence after RCA correction) and GCi values.
  2. Identify orthogroups of interest using SPECIES PRESENCE PATTERN, not GCi sign.
     (GCi sign convention can invert in small biased subsets — see Methods note below.)

  Three filters:
    A) "Liverwort-enriched":   RCA=1 in ≥1 liverwort AND RCA=0 in ALL angiosperms
                               Sorted by ubiquity desc (broadly conserved in non-angiosperms)
    B) "Both liverworts":      RCA=1 in BOTH liverworts AND RCA=0 in ALL angiosperms
                               Most conservative; shared by Marchantia + Lunularia
    C) "Bryophyte-conserved":  RCA=1 in ≥1 liverwort AND RCA=1 in Ppatens
                               AND RCA=0 in ALL angiosperms
                               Suggests pre-angiosperm conservation across bryophytes

  3. Extract gene IDs from OrthoFinder Orthogroups.tsv.
  4. Save organised output files to gene_lists/<context>/.

METHODS NOTE — GCi sign:
  In the economic complexity framework the GCi sign convention is set so that
  lower ubiquity correlates with higher GCi. In filtered metabolic subsets the
  sign can invert because angiosperm-specific OGs (low ubiquity, very negative
  MCi) skew the ubiquity-GCi correlation positive, flipping the sign and making
  liverwort-specific OGs (low ubiquity, positive MCi) appear to have NEGATIVE
  GCi. We therefore use species presence patterns as the primary filter and
  report |GCi_rca| as an auxiliary sort key.

Output per context:
  01_OGs_liverwort_specific.tsv      Filter A: OG table, per-species RCA, ubiquity, GCi
  02_OGs_both_liverworts.tsv         Filter B
  03_OGs_bryophyte_conserved.tsv     Filter C
  04_genes_by_species.tsv            All genes in filter-A OGs, flat table
  05_genes_Mpolymorpha.txt           Marchantia gene IDs (one per line)
  05_genes_Lcruciata.txt             Lunularia gene IDs
  05_genes_Ppatens.txt               Physcomitrium gene IDs
  05_genes_Creinhardtii.txt          Chlamydomonas gene IDs
  05_genes_Afiliculoides.txt         Azolla gene IDs
  05_genes_Scucullata.txt            Salvinia gene IDs
  06_summary.txt                     Plain-text summary for records and Methods

Usage (from project root):
    python3 08_extract_functional_genes.py [--orthofinder-dir data/orthofinder_v2]
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── Constants ──────────────────────────────────────────────────────────────

LIVERWORTS  = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]
NON_ANGIO   = ["Mpolymorpha", "Lcruciata", "Ppatens",
                "Creinhardtii", "Afiliculoides", "Scucullata",
                "Smoellendorffii", "Pabies"]
ALL_SPECIES = ["Afiliculoides", "Athaliana", "Creinhardtii", "Lcruciata",
               "Mpolymorpha", "Pabies", "Pamericana", "Ppatens",
               "Scucullata", "Slycopersicum", "Smoellendorffii"]

DISPLAY = {
    "Mpolymorpha":      "Marchantia polymorpha",
    "Lcruciata":        "Lunularia cruciata",
    "Ppatens":          "Physcomitrium patens",
    "Creinhardtii":     "Chlamydomonas reinhardtii",
    "Afiliculoides":    "Azolla filiculoides",
    "Scucullata":       "Salvinia cucullata",
    "Smoellendorffii":  "Selaginella moellendorffii",
    "Pabies":           "Picea abies",
    "Athaliana":        "Arabidopsis thaliana",
    "Slycopersicum":    "Solanum lycopersicum",
    "Pamericana":       "Persea americana",
}

CONTEXTS = {
    "C_metabolism":       "data/results_C_metabolism_v3",
    "N_metabolism":       "data/results_N_metabolism_v3",
    "primary_metabolism": "data/results_primary_metabolism_v3",
}


# ── OrthoFinder helpers ────────────────────────────────────────────────────

def find_orthogroups_tsv(base_dir: Path) -> Path:
    for pattern in [base_dir / "Results_*" / "Orthogroups" / "Orthogroups.tsv",
                    base_dir / "OrthoFinder" / "Results_*" / "Orthogroups" / "Orthogroups.tsv"]:
        matches = sorted(glob.glob(str(pattern)))
        if matches:
            return Path(matches[-1])
    sys.exit(f"ERROR: Orthogroups.tsv not found under {base_dir}")


def load_gene_table(tsv: Path) -> pd.DataFrame:
    return pd.read_csv(tsv, sep="\t", index_col=0)


def genes_for_ogs(gene_table: pd.DataFrame, og_ids: list[str],
                  species_col: str) -> list[str]:
    """Flat list of unique gene IDs for given OGs and one species column."""
    valid = [og for og in og_ids if og in gene_table.index]
    genes: list[str] = []
    for cell in gene_table.loc[valid, species_col]:
        if pd.notna(cell) and str(cell).strip():
            genes.extend(g.strip() for g in str(cell).split(",") if g.strip())
    return sorted(set(genes))


# ── Core analysis ──────────────────────────────────────────────────────────

def analyse_context(context_name: str, results_dir: Path,
                    gene_table: pd.DataFrame, out_root: Path) -> dict:
    out_dir = out_root / context_name
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = results_dir / "matrix_rca.tsv"
    gci_path    = results_dir / "og_complexity.tsv"

    if not matrix_path.exists():
        print(f"  [SKIP] {matrix_path} not found")
        return {}

    mat = pd.read_csv(matrix_path, sep="\t", index_col=0)
    gci = pd.read_csv(gci_path,    sep="\t", index_col=0).reindex(mat.index)

    n_total = len(mat)

    # Presence masks
    in_any_lw   = mat[LIVERWORTS].max(axis=1) == 1
    in_both_lw  = mat[LIVERWORTS].min(axis=1) == 1
    in_ppatens  = mat["Ppatens"] == 1 if "Ppatens" in mat.columns else pd.Series(False, index=mat.index)
    absent_angio = mat[ANGIOSPERMS].sum(axis=1) == 0

    # ── Filter A ─────────────────────────────────────────────────────────
    mask_A = in_any_lw & absent_angio
    ogs_A  = _build_og_table(mat, gci, mask_A)

    # ── Filter B ─────────────────────────────────────────────────────────
    mask_B = in_both_lw & absent_angio
    ogs_B  = _build_og_table(mat, gci, mask_B)

    # ── Filter C ─────────────────────────────────────────────────────────
    mask_C = in_any_lw & in_ppatens & absent_angio
    ogs_C  = _build_og_table(mat, gci, mask_C)

    # ── Save OG tables ────────────────────────────────────────────────────
    ogs_A.to_csv(out_dir / "01_OGs_liverwort_specific.tsv", sep="\t")
    ogs_B.to_csv(out_dir / "02_OGs_both_liverworts.tsv",    sep="\t")
    ogs_C.to_csv(out_dir / "03_OGs_bryophyte_conserved.tsv", sep="\t")

    # ── Gene lists for filter A ───────────────────────────────────────────
    og_ids_A = list(ogs_A.index)

    # Flat table: one row per gene
    rows = []
    for og in og_ids_A:
        if og not in gene_table.index:
            continue
        for sp in ALL_SPECIES:
            if sp not in gene_table.columns:
                continue
            cell = gene_table.at[og, sp]
            if pd.notna(cell) and str(cell).strip():
                for g in str(cell).split(","):
                    g = g.strip()
                    if g:
                        rows.append({
                            "orthogroup":      og,
                            "species_code":    sp,
                            "species_name":    DISPLAY.get(sp, sp),
                            "gene_id":         g,
                            "ubiquity":        int(ogs_A.at[og, "ubiquity"]),
                            "GCi_rca":         ogs_A.at[og, "GCi_rca"],
                            "rca_Mpolymorpha": int(mat.at[og, "Mpolymorpha"]) if "Mpolymorpha" in mat.columns else None,
                            "rca_Lcruciata":   int(mat.at[og, "Lcruciata"])   if "Lcruciata"   in mat.columns else None,
                        })

    pd.DataFrame(rows).to_csv(out_dir / "04_genes_by_species.tsv", sep="\t", index=False)

    # Per-species flat .txt files
    gene_counts = {}
    for sp in NON_ANGIO:
        if sp not in gene_table.columns:
            continue
        genes = genes_for_ogs(gene_table, og_ids_A, sp)
        gene_counts[sp] = len(genes)
        hdr = (f"# {DISPLAY.get(sp, sp)}\n"
               f"# Context: {context_name}\n"
               f"# Filter: RCA=1 in ≥1 liverwort, RCA=0 in all angiosperms\n"
               f"# Source OGs: {len(og_ids_A)}  |  Genes in this species: {len(genes)}\n")
        (out_dir / f"05_genes_{sp}.txt").write_text(hdr + "\n".join(genes) + "\n")

    # ── Summary ───────────────────────────────────────────────────────────
    lines = [
        f"=== {context_name} — gene extraction summary (2026-05-05) ===",
        f"",
        f"Input: {n_total} orthogroups (MCi-RCA matrix, ≥2 species after singleton filter)",
        f"",
        f"--- Filter A: ≥1 liverwort present, all angiosperms absent ---",
        f"  OGs: {len(ogs_A)} / {n_total} ({100*len(ogs_A)/n_total:.1f}%)",
        f"  Ubiquity range: {ogs_A['ubiquity'].min()}–{ogs_A['ubiquity'].max()} species",
        f"  Genes per non-angiosperm species in these OGs:",
    ]
    for sp in NON_ANGIO:
        n = gene_counts.get(sp, 0)
        lines.append(f"    {DISPLAY.get(sp, sp):<35} {n:>5} genes")
    lines += [
        f"",
        f"--- Filter B: BOTH liverworts present, all angiosperms absent ---",
        f"  OGs: {len(ogs_B)} / {n_total} ({100*len(ogs_B)/n_total:.1f}%)",
        f"",
        f"--- Filter C: ≥1 liverwort + Physcomitrium present, all angiosperms absent ---",
        f"  OGs: {len(ogs_C)} / {n_total} ({100*len(ogs_C)/n_total:.1f}%)",
        f"  (Bryophyte-conserved: OGs shared by liverworts + moss, lost in angiosperms)",
        f"",
        f"--- GCi note ---",
        f"  GCi sign convention can invert in small biased subsets.",
        f"  OGs sorted by ubiquity (ascending) within filter A — lower ubiquity =",
        f"  more restricted to non-angiosperm lineages.",
        f"",
        f"--- Output files ---",
        f"  01_OGs_liverwort_specific.tsv    Filter A OG table",
        f"  02_OGs_both_liverworts.tsv       Filter B OG table",
        f"  03_OGs_bryophyte_conserved.tsv   Filter C OG table",
        f"  04_genes_by_species.tsv          All genes in filter-A OGs, all species",
        f"  05_genes_SPECIES.txt             Per-species gene ID lists (filter A)",
        f"  06_summary.txt                   This file",
    ]
    (out_dir / "06_summary.txt").write_text("\n".join(lines) + "\n")

    print(f"\n  [OK] {context_name}")
    print(f"       Filter A (≥1 liverwort, no angiosperms): {len(ogs_A)} OGs")
    print(f"       Filter B (both liverworts, no angiosperms): {len(ogs_B)} OGs")
    print(f"       Filter C (liverwort + Ppatens, no angiosperms): {len(ogs_C)} OGs")
    print(f"       Marchantia genes (filter A): {gene_counts.get('Mpolymorpha',0)}")
    print(f"       Lunularia genes  (filter A): {gene_counts.get('Lcruciata',0)}")

    return {
        "context":           context_name,
        "total_ogs":         n_total,
        "filter_A_ogs":      len(ogs_A),
        "filter_B_ogs":      len(ogs_B),
        "filter_C_ogs":      len(ogs_C),
        "genes_Mpolymorpha": gene_counts.get("Mpolymorpha", 0),
        "genes_Lcruciata":   gene_counts.get("Lcruciata",   0),
        "genes_Ppatens":     gene_counts.get("Ppatens",     0),
    }


def _build_og_table(mat: pd.DataFrame, gci: pd.DataFrame,
                    mask: pd.Series) -> pd.DataFrame:
    """Build sorted OG summary table for a given boolean mask."""
    df = mat[mask].copy()
    df["ubiquity"] = mat[mask].sum(axis=1)
    df["GCi_rca"]  = gci.loc[mask, "GCi_rca"]
    df["GCi_raw"]  = gci.loc[mask, "GCi_raw"]
    df["abs_GCi_rca"] = df["GCi_rca"].abs()
    # Primary sort: ubiquity ascending (most restricted first), then |GCi| descending
    return df.sort_values(["ubiquity", "abs_GCi_rca"], ascending=[True, False]) \
             .drop(columns=["abs_GCi_rca"])


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Extract gene lists by functional context",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    parser.add_argument("--out-root", default="gene_lists", type=Path)
    args = parser.parse_args()

    args.out_root.mkdir(parents=True, exist_ok=True)

    print("[1] Loading OrthoFinder gene table ...")
    tsv = find_orthogroups_tsv(args.orthofinder_dir)
    print(f"    {tsv}")
    gene_table = load_gene_table(tsv)
    print(f"    {gene_table.shape[0]:,} orthogroups × {gene_table.shape[1]} species")

    print("\n[2] Extracting gene lists per functional context ...")
    summaries = []
    for context_name, results_dir in CONTEXTS.items():
        s = analyse_context(context_name, Path(results_dir), gene_table, args.out_root)
        if s:
            summaries.append(s)

    # Cross-context summary
    if summaries:
        print("\n[3] Cross-context summary ...")
        df = pd.DataFrame(summaries).set_index("context")
        df.to_csv(args.out_root / "00_cross_context_summary.tsv", sep="\t")
        print(f"\n{'Context':<25} {'TotalOGs':>9} {'FilterA':>8} {'FilterB':>8} "
              f"{'FilterC':>8} {'M.poly':>8} {'L.cru':>8}")
        print("-" * 76)
        for _, r in df.iterrows():
            print(f"{r.name:<25} {int(r.total_ogs):>9} {int(r.filter_A_ogs):>8} "
                  f"{int(r.filter_B_ogs):>8} {int(r.filter_C_ogs):>8} "
                  f"{int(r.genes_Mpolymorpha):>8} {int(r.genes_Lcruciata):>8}")

    # README
    readme = """# gene_lists/ — Functional Gene Lists

Generated: 2026-05-05  |  Script: 08_extract_functional_genes.py
Reproduce: python3 08_extract_functional_genes.py

## Biological question
Which orthogroups — within each metabolic context — are over-represented
(RCA ≥ 1) in liverworts but NOT over-represented in any angiosperm?
These OGs represent the part of the liverwort metabolic repertoire that
cannot be explained by genome size and is absent from angiosperm genomes.

## Three filters

### Filter A — liverwort-specific  (primary gene lists)
  Criteria: RCA=1 in ≥1 liverwort (M. polymorpha or L. cruciata)
            AND RCA=0 in ALL angiosperms (A. thaliana, S. lycopersicum, P. americana)
  Sorted: ubiquity ascending (most lineage-restricted OGs first),
          then |GCi_rca| descending.

### Filter B — both liverworts  (conservative, shared by genus pair)
  Criteria: RCA=1 in BOTH M. polymorpha AND L. cruciata
            AND RCA=0 in ALL angiosperms
  Interpretation: shared by two liverwort genera; more likely ancestral
  liverwort character than species-specific expansion.

### Filter C — bryophyte-conserved  (pre-angiosperm)
  Criteria: RCA=1 in ≥1 liverwort AND RCA=1 in P. patens (moss)
            AND RCA=0 in ALL angiosperms
  Interpretation: present in liverworts AND mosses but lost in angiosperms.
  Suggests the over-representation predates the liverwort/moss split
  and was selectively lost in vascular plant lineages.

## Note on GCi sign convention
  In functional subsets the GCi sign can invert relative to the whole-genome
  analysis. OGs sorted by ubiquity (ascending) are the most biologically
  informative first-pass list. GCi_rca is reported as auxiliary information.
  See Methods note in script docstring for full explanation.

## Directory structure

  gene_lists/
    00_cross_context_summary.tsv         OG and gene counts per context and filter
    C_metabolism/
      01_OGs_liverwort_specific.tsv      Filter A OG table
      02_OGs_both_liverworts.tsv         Filter B OG table
      03_OGs_bryophyte_conserved.tsv     Filter C OG table
      04_genes_by_species.tsv            Flat gene table (filter A, all species)
      05_genes_Mpolymorpha.txt           Marchantia gene IDs
      05_genes_Lcruciata.txt             Lunularia gene IDs
      05_genes_Ppatens.txt               Physcomitrium gene IDs
      05_genes_Creinhardtii.txt          Chlamydomonas gene IDs
      05_genes_Afiliculoides.txt         Azolla gene IDs
      05_genes_Scucullata.txt            Salvinia gene IDs
      06_summary.txt                     Summary statistics
    N_metabolism/    (same structure)
    primary_metabolism/    (same structure)

## Gene IDs
  IDs match the source proteome files in data/proteomes/filtered/.
  For KEGG/GO enrichment of Marchantia genes, use:
    gene_lists/<context>/05_genes_Mpolymorpha.txt
  These IDs follow the Ensembl Plants r62 naming (MpXXXXXX format).
  Back-reference to orthogroups via 04_genes_by_species.tsv (column: orthogroup).
"""
    (args.out_root / "README.txt").write_text(readme)
    print(f"\nDone. All output in: {args.out_root}/")


if __name__ == "__main__":
    main()
