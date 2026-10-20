#!/usr/bin/env python3
"""
Build presence/absence matrix from OrthoFinder output and compute
Metabolic Complexity Index (MCi) and Gene Complexity Index (GCi)
following Hidalgo & Hausmann's economic complexity framework.

Two variants:
  MCi-raw : binary presence/absence matrix
  MCi-RCA : binarized by Revealed Comparative Advantage (corrects for genome size)

Usage:
    python3 05_build_matrix.py [--orthofinder-dir data/orthofinder]
    python3 05_build_matrix.py --og-filter data/filters/C_metabolism_ogs.txt \\
        --out-dir data/results_C_metabolism
"""

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── 1. Locate OrthoFinder results ─────────────────────────────────────────

def find_orthogroups_tsv(base_dir: Path) -> Path:
    patterns = [
        base_dir / "Results_*" / "Orthogroups" / "Orthogroups.tsv",
        base_dir / "OrthoFinder" / "Results_*" / "Orthogroups" / "Orthogroups.tsv",
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(str(pattern)))
        if matches:
            return Path(matches[-1])
    sys.exit(f"ERROR: Could not find Orthogroups.tsv under {base_dir}\n"
             f"Run 04_run_orthofinder.sh first.")


# ── 2. Parse OrthoFinder table → copy-number matrix ───────────────────────

def parse_orthogroups(tsv: Path) -> pd.DataFrame:
    """
    Returns DataFrame with shape (n_orthogroups × n_species).
    Values = number of genes per species per OG (copy number).
    """
    raw = pd.read_csv(tsv, sep="\t", index_col=0)

    def count_genes(cell):
        if pd.isna(cell) or str(cell).strip() == "":
            return 0
        return len(str(cell).split(","))

    counts = raw.map(count_genes)
    return counts


# ── 3. Build binary matrices ───────────────────────────────────────────────

def build_raw_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    """MCi-raw: 1 if species has ≥1 gene in OG."""
    return (counts > 0).astype(int)


def build_rca_matrix(counts: pd.DataFrame) -> pd.DataFrame:
    """
    MCi-RCA: binarized Revealed Comparative Advantage.
    RCA(species, OG) = (copies_sp / total_genes_sp) / (copies_all / total_genes_all)
    Entry = 1 if RCA ≥ 1, else 0.
    """
    total_per_species = counts.sum(axis=0)   # total genes per species
    total_per_og      = counts.sum(axis=1)   # total copies per OG across all species
    grand_total       = counts.values.sum()

    # Broadcast: RCA[og, sp] = (counts[og,sp] / total_per_species[sp])
    #                         / (total_per_og[og] / grand_total)
    rca = counts.div(total_per_species, axis=1)          # normalize by species size
    rca = rca.div(total_per_og / grand_total, axis=0)    # normalize by OG ubiquity

    return (rca >= 1).astype(int)


# ── 4. Filter low-information rows/columns ─────────────────────────────────

def filter_matrix(M: pd.DataFrame, min_species: int = 2) -> pd.DataFrame:
    """
    Remove OGs present in fewer than min_species species.
    Singleton OGs carry no comparative information.
    """
    og_breadth = M.sum(axis=1)
    return M[og_breadth >= min_species]


# ── 5. Economic complexity (Hidalgo & Hausmann 2009) ──────────────────────

def compute_complexity(M: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """
    Compute MCi (species complexity) and GCi (orthogroup complexity)
    as the second eigenvector of the transition matrices.

    M: binary matrix, shape (n_og × n_species)
    Returns: (MCi Series indexed by species, GCi Series indexed by OG)
    """
    mat = M.values.astype(float)
    species = M.columns.tolist()
    ogs     = M.index.tolist()

    k_sp = mat.sum(axis=0)   # diversity: OGs per species  (shape: n_species)
    k_og = mat.sum(axis=1)   # ubiquity:  species per OG   (shape: n_og)

    # Avoid division by zero
    k_sp_safe = np.where(k_sp == 0, 1, k_sp)
    k_og_safe = np.where(k_og == 0, 1, k_og)

    # Transition matrix for species (8×8 — tractable):
    # T_sp[s,s'] = (1/k_sp[s]) · Σ_og (M[og,s] · M[og,s'] / k_og[og])
    T_sp = (mat.T / k_sp_safe[:, None]) @ (mat / k_og_safe[:, None])

    eigenvalues, eigenvectors = np.linalg.eig(T_sp)
    idx = np.argsort(eigenvalues.real)[::-1]
    MCi = eigenvectors[:, idx[1]].real

    # Sign convention: higher diversity → higher MCi
    if np.corrcoef(k_sp, MCi)[0, 1] < 0:
        MCi = -MCi

    # GCi via method of reflections (avoids n_og × n_og matrix):
    # GCi[og] = (1/k_og[og]) · Σ_sp (M[og,sp] · MCi[sp])
    GCi = (mat @ MCi) / k_og_safe

    # Sign convention: lower ubiquity → higher GCi
    if np.corrcoef(k_og, GCi)[0, 1] > 0:
        GCi = -GCi

    return pd.Series(MCi, index=species, name="MCi"), \
           pd.Series(GCi, index=ogs,     name="GCi")


# ── 6. Summary stats ───────────────────────────────────────────────────────

def print_summary(label: str, M_bin: pd.DataFrame,
                  MCi: pd.Series, GCi: pd.Series) -> None:
    print(f"\n{'═'*60}")
    print(f"  {label}")
    print(f"{'═'*60}")
    print(f"  Matrix: {M_bin.shape[0]:,} orthogroups × {M_bin.shape[1]} species")
    print(f"  Mean OG breadth (ubiquity): {M_bin.sum(axis=1).mean():.2f} species/OG")
    print(f"  Mean species breadth (diversity): {M_bin.sum(axis=0).mean():.0f} OGs/species")

    df = pd.DataFrame({
        "diversity": M_bin.sum(axis=0),
        "MCi": MCi
    }).sort_values("MCi", ascending=False)

    print(f"\n  {'Species':<25} {'Diversity':>10} {'MCi':>10}")
    print(f"  {'-'*47}")
    for sp, row in df.iterrows():
        print(f"  {sp:<25} {int(row.diversity):>10,} {row.MCi:>10.4f}")


# ── 7. Main ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthofinder-dir", default="data/orthofinder",
                        type=Path)
    parser.add_argument("--min-species", default=2, type=int,
                        help="Min species per OG to keep (default: 2)")
    parser.add_argument("--out-dir", default="data/results", type=Path)
    parser.add_argument("--og-filter", default=None, type=Path,
                        help="File with orthogroup IDs to include (one per line). "
                             "If not given, all OGs are used.")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Locate input
    tsv = find_orthogroups_tsv(args.orthofinder_dir)
    print(f"Reading: {tsv}")

    # Parse
    counts = parse_orthogroups(tsv)
    print(f"Loaded: {counts.shape[0]:,} orthogroups × {counts.shape[1]} species")

    # Apply OG filter if provided
    if args.og_filter is not None:
        og_ids = set(args.og_filter.read_text().split())
        before = counts.shape[0]
        counts = counts[counts.index.isin(og_ids)]
        print(f"OG filter '{args.og_filter.name}': {before:,} → {counts.shape[0]:,} OGs "
              f"({counts.shape[0]/before*100:.1f}% retained)")

    # Build matrices
    M_raw = build_raw_matrix(counts)
    M_rca = build_rca_matrix(counts)

    # Filter singletons
    M_raw_f = filter_matrix(M_raw, args.min_species)
    M_rca_f = filter_matrix(M_rca, args.min_species)

    print(f"After filtering (≥{args.min_species} species): "
          f"{M_raw_f.shape[0]:,} OGs (raw), {M_rca_f.shape[0]:,} OGs (RCA)")

    # Complexity indices
    MCi_raw, GCi_raw = compute_complexity(M_raw_f)
    MCi_rca, GCi_rca = compute_complexity(M_rca_f)

    # Label for summaries
    filter_label = args.og_filter.stem if args.og_filter else "whole_genome"

    # Print summaries
    print_summary(f"MCi-raw  (presence/absence) [{filter_label}]", M_raw_f, MCi_raw, GCi_raw)
    print_summary(f"MCi-RCA  (genome-size corrected) [{filter_label}]", M_rca_f, MCi_rca, GCi_rca)

    # Save outputs
    M_raw_f.to_csv(args.out_dir / "matrix_raw.tsv", sep="\t")
    M_rca_f.to_csv(args.out_dir / "matrix_rca.tsv", sep="\t")

    results = pd.DataFrame({
        "diversity":    M_raw_f.sum(axis=0),
        "MCi_raw":      MCi_raw,
        "MCi_rca":      MCi_rca,
        "filter":       filter_label,
    }).sort_values("MCi_raw", ascending=False)
    results.to_csv(args.out_dir / "species_complexity.tsv", sep="\t")

    og_results = pd.DataFrame({
        "ubiquity":  M_raw_f.sum(axis=1),
        "GCi_raw":   GCi_raw,
        "GCi_rca":   GCi_rca,
    }).sort_values("GCi_raw", ascending=False)
    og_results.to_csv(args.out_dir / "og_complexity.tsv", sep="\t")

    print(f"\nOutputs saved to: {args.out_dir}/")
    print("  matrix_raw.tsv          — binary presence/absence matrix")
    print("  matrix_rca.tsv          — RCA-binarized matrix")
    print("  species_complexity.tsv  — MCi-raw and MCi-RCA per species")
    print("  og_complexity.tsv       — GCi-raw and GCi-RCA per orthogroup")


if __name__ == "__main__":
    main()
