#!/usr/bin/env python3
"""
06_go_filter_multispecies.py — GO functional filters using NATIVE GO
annotation from every species that has one (not just Arabidopsis).

Baseline (06_go_filter.py) tags an orthogroup by looking ONLY at Arabidopsis
genes in that row. This script tags an orthogroup if ANY gene from ANY of
the 8 species with native GO (Athaliana, Mpolymorpha, Ppatens, Slycopersicum,
Creinhardtii, Smoellendorffii, Afiliculoides, Scucullata) maps to the target
category — via GO ancestor propagation, same term sets as 06_go_filter.py.

Produces the same 4 filter files for direct comparison:
    data/filters_go_multispecies/primary_metabolism_ogs.txt
    data/filters_go_multispecies/C_metabolism_ogs.txt
    data/filters_go_multispecies/N_metabolism_ogs.txt
    data/filters_go_multispecies/regulatory_networks_ogs.txt
    data/filters_go_multispecies/comparison_vs_arabidopsis_only.tsv

Usage:
    python3 06_go_filter_multispecies.py \
        --orthofinder-dir data/orthofinder_v3 \
        --go-species-dir ../output/GO_por_especie \
        --go-obo data/go/go-basic.obo \
        --baseline-dir data/filters \
        --out-dir data/filters_go_multispecies
"""

import argparse
import glob
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ── Same GO target term sets as 06_go_filter.py ────────────────────────────

C_GO_TERMS = {
    "GO:0005975", "GO:0016051", "GO:0016052", "GO:0044262", "GO:0006091",
    "GO:0015979", "GO:0019684", "GO:0019685", "GO:0006096", "GO:0006099",
    "GO:0009052", "GO:0005982", "GO:0005986",
}
N_GO_TERMS = {
    "GO:0006520", "GO:0009308", "GO:0009309", "GO:0042126", "GO:0042128",
    "GO:0019740", "GO:0006749",
}
LIPID_GO_TERMS = {"GO:0006629", "GO:0008610", "GO:0016042", "GO:0006631"}
NUCLEOTIDE_GO_TERMS = {"GO:0009117", "GO:0006163", "GO:0006220", "GO:0006753"}
COFACTOR_GO_TERMS = {"GO:0051186", "GO:0009108", "GO:0006766"}
PRIMARY_GO_TERMS = (C_GO_TERMS | N_GO_TERMS | LIPID_GO_TERMS
                     | NUCLEOTIDE_GO_TERMS | COFACTOR_GO_TERMS)
REGULATORY_GO_TERMS = {
    "GO:0000160", "GO:0009755", "GO:0007623", "GO:0045087", "GO:0046854",
}
# GO:0006355 (regulation of DNA-templated transcription, 9,224 genes pooled)
# and GO:0007165 (signal transduction, 6,847 genes pooled) were dropped:
# on their own they cover ~2x the genes of C_metabolism and ~5x N_metabolism,
# so keeping them made "regulatory_networks" incomparable to the other
# categories (biased any C/N-vs-regulatory complexity comparison by
# methodology, not biology). Without them the category totals 7,225 genes —
# in line with C_metabolism (8,391). See go_term_popularity.py /
# regulatory_fix_options.py for the numbers behind this call.

CATEGORY_SETS = {
    "primary_metabolism": PRIMARY_GO_TERMS,
    "C_metabolism": C_GO_TERMS,
    "N_metabolism": N_GO_TERMS,
    "regulatory_networks": REGULATORY_GO_TERMS,
}

# ── GO ontology (ancestor propagation) ──────────────────────────────────────

RELATIONS = {"is_a", "part_of", "regulates",
             "positively_regulates", "negatively_regulates"}


def parse_go_obo(obo_path: Path) -> dict[str, set[str]]:
    parent_map: dict[str, set[str]] = defaultdict(set)
    current_id = None
    in_term = False
    with open(obo_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                in_term = True
                continue
            if not in_term:
                continue
            if not line:
                in_term = False
                current_id = None
                continue
            if line.startswith("id:"):
                current_id = line.split(":", 1)[1].strip()
            elif line.startswith("is_a:") and current_id:
                parent_map[current_id].add(line.split()[1])
            elif line.startswith("relationship:") and current_id:
                parts = line.split()
                if len(parts) >= 3 and parts[1] in RELATIONS:
                    parent_map[current_id].add(parts[2])
    print(f"  Parsed {len(parent_map):,} GO terms with parent relations")
    return dict(parent_map)


def parse_go_id_remap(obo_path: Path) -> tuple[set[str], dict[str, str], dict[str, str]]:
    """Secondary-ID resolution: alt_id (merged synonym IDs) and replaced_by
    (obsolete term with a single unambiguous successor). Needed because
    eggNOG-mapper's internal GO database is older than the go-basic.obo used
    here — ~24% of its emitted GO IDs don't match a current primary id
    (see go_annotation_depth_by_species.py); about 43% of those are cleanly
    recoverable through alt_id/replaced_by. Terms with only ambiguous
    `consider:` suggestions, or no successor at all, are left unresolved on
    purpose — there is no safe automatic choice for those."""
    primary_ids: set[str] = set()
    alt_to_primary: dict[str, str] = {}
    replaced_by: dict[str, str] = {}
    current_id = None
    in_term = False
    with open(obo_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line == "[Term]":
                in_term = True
                current_id = None
                continue
            if not in_term:
                continue
            if not line:
                in_term = False
                current_id = None
                continue
            if line.startswith("id:"):
                current_id = line.split(":", 1)[1].strip()
                primary_ids.add(current_id)
            elif line.startswith("alt_id:") and current_id:
                alt_to_primary[line.split(":", 1)[1].strip()] = current_id
            elif line.startswith("replaced_by:") and current_id:
                replaced_by[current_id] = line.split(":", 1)[1].strip()
    return primary_ids, alt_to_primary, replaced_by


def normalize_gene_to_go(gene_to_go: dict[str, set[str]], primary_ids: set[str],
                         alt_to_primary: dict[str, str],
                         replaced_by: dict[str, str]) -> dict[str, set[str]]:
    """Remap alt_id/obsolete-with-replacement terms to their current primary
    id; drop terms that still don't resolve (ambiguous `consider:` or no
    successor). Applied once per species right after loading, so every
    downstream script (filters, ancestor propagation, depth stats) sees only
    resolvable ids."""
    out: dict[str, set[str]] = {}
    for gid, terms in gene_to_go.items():
        resolved = set()
        for t in terms:
            if t in primary_ids:
                resolved.add(t)
            elif t in alt_to_primary:
                resolved.add(alt_to_primary[t])
            elif t in replaced_by:
                resolved.add(replaced_by[t])
        if resolved:
            out[gid] = resolved
    return out


def compute_ancestors(terms: set[str], parent_map: dict) -> dict[str, frozenset]:
    cache: dict[str, frozenset] = {}

    def ancestors(term: str) -> frozenset:
        if term in cache:
            return cache[term]
        seen: set[str] = set()
        stack = [term]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            stack.extend(parent_map.get(cur, ()))
        cache[term] = frozenset(seen)
        return cache[term]

    return {t: ancestors(t) for t in terms}


# ── Per-species native GO parsers → {normalized_gene_id: {GO_ids}} ─────────

AT_LOCUS_RE = re.compile(r"(AT[1-5MC]G\d+)", re.IGNORECASE)


def parse_athaliana_gaf(path: Path) -> dict[str, set[str]]:
    gene_to_go: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            if "NOT" in parts[3] or not parts[4].startswith("GO:"):
                continue
            m = AT_LOCUS_RE.search(parts[10]) or AT_LOCUS_RE.search(parts[1])
            if m:
                gene_to_go[m.group(1).upper()].add(parts[4])
    return dict(gene_to_go)


def parse_biomart_tsv(path: Path, id_normalizer) -> dict[str, set[str]]:
    """Generic parser for the Ensembl BioMart exports (ensembl_gene_id or
    protein_id in a fixed column, go_id in another)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    protein_cols = [c for c in df.columns if "protein" in c.lower()]
    id_col = protein_cols[0] if protein_cols else df.columns[0]
    go_col = [c for c in df.columns if "GO term accession" in c or c == "go_id"][0]
    gene_to_go: dict[str, set[str]] = defaultdict(set)
    for _, row in df.iterrows():
        raw_id = row[id_col]
        go_id = row[go_col]
        if pd.isna(raw_id) or pd.isna(go_id) or not str(go_id).startswith("GO:"):
            continue
        norm = id_normalizer(str(raw_id))
        if norm:
            gene_to_go[norm].add(go_id)
    return dict(gene_to_go)


GO_TOKEN_RE = re.compile(r"GO:\d{7}")


def parse_eggnog_annotations(path: Path) -> dict[str, set[str]]:
    """eggNOG-mapper .emapper.annotations — col 1 '#query' (gene id, matches
    Orthogroups.tsv directly since we fed it the same filtered proteome),
    col 'GOs' comma-separated GO ids or '-' if none."""
    with open(path) as fh:
        header_idx = next(i for i, l in enumerate(fh) if l.startswith("#query"))
    df = pd.read_csv(path, sep="\t", dtype=str, skiprows=header_idx)
    df = df[~df.iloc[:, 0].astype(str).str.startswith("#")]
    gene_to_go: dict[str, set[str]] = defaultdict(set)
    for _, row in df.iterrows():
        gid = row.iloc[0]
        go_field = row.get("GOs")
        if pd.isna(gid) or pd.isna(go_field) or go_field == "-" or go_field == "":
            continue
        terms = set(GO_TOKEN_RE.findall(str(go_field)))
        if terms:
            gene_to_go[str(gid).strip()] |= terms
    return dict(gene_to_go)


def parse_fernbase_annotation(path: Path) -> dict[str, set[str]]:
    """FernBase functional_annotation files: locus + 3 GO columns with
    embedded 'term name (GO:XXXXXXX) [score]' text."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    go_cols = [c for c in df.columns if c.startswith("GO_")]
    gene_to_go: dict[str, set[str]] = defaultdict(set)
    for _, row in df.iterrows():
        locus = row["locus"]
        if pd.isna(locus):
            continue
        terms = set()
        for c in go_cols:
            val = row[c]
            if isinstance(val, str):
                terms |= set(GO_TOKEN_RE.findall(val))
        if terms:
            gene_to_go[locus] = terms
    return dict(gene_to_go)


def norm_direct(x: str) -> str:
    return x.strip()


def norm_first_dot(x: str) -> str:
    return x.split(".")[0].strip()


def norm_strip_v_version(x: str) -> str:
    return re.sub(r"V\d+\.\d+$", "", x.strip())


def norm_solyc(x: str) -> str | None:
    m = re.search(r"(Solyc\d+g\d+)", x)
    return m.group(1) if m else None


SPECIES_PARSERS = {
    # column name in Orthogroups.tsv -> (file, parser, id_normalizer for Orthogroups side)
    "Athaliana": ("Athaliana/goa_arabidopsis.gaf.gz", "gaf", norm_first_dot),
    "Mpolymorpha": ("Mpolymorpha/GO_biomart_MpTak1_v7.1.tsv", "biomart", norm_first_dot),
    "Ppatens": ("Ppatens/GO_biomart_PhypaV3.tsv", "biomart", norm_strip_v_version),
    "Slycopersicum": ("Slycopersicum/GO_biomart_SL4.0.tsv", "biomart", norm_solyc),
    "Creinhardtii": ("Creinhardtii/GO_biomart_v5.5_via_protein_id.tsv", "biomart", norm_direct),
    "Smoellendorffii": ("Smoellendorffii/GO_biomart_v1.0_via_protein_id.tsv", "biomart", norm_direct),
    "Afiliculoides": ("Afiliculoides/Afiliculoides_functional_annotation_v1.1.txt", "fernbase", norm_direct),
    "Scucullata": ("Scucullata/Scucullata_functional_annotation_v2.txt", "fernbase", norm_direct),
    "Lcruciata": ("Lcruciata/Lcruciata_eggnog_annotations.tsv", "eggnog", norm_direct),
    "Pamericana": ("Pamericana/Pamericana_eggnog_annotations.tsv", "eggnog", norm_direct),
    "Pabies": ("Pabies/Pabies_eggnog_annotations.tsv", "eggnog", norm_direct),
}


def find_orthogroups_tsv(base_dir: Path) -> Path:
    patterns = [base_dir / "Results_*" / "Orthogroups" / "Orthogroups.tsv"]
    for pattern in patterns:
        matches = sorted(glob.glob(str(pattern)))
        if matches:
            return Path(matches[-1])
    sys.exit(f"ERROR: Could not find Orthogroups.tsv under {base_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orthofinder-dir", default="data/orthofinder_v3", type=Path)
    ap.add_argument("--go-species-dir", default="../output/GO_por_especie", type=Path)
    ap.add_argument("--go-obo", default="data/go/go-basic.obo", type=Path)
    ap.add_argument("--baseline-dir", default="data/filters", type=Path,
                    help="Arabidopsis-only filters, for comparison")
    ap.add_argument("--out-dir", default="data/filters_go_multispecies", type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("[1] Loading Orthogroups.tsv ...")
    tsv = find_orthogroups_tsv(args.orthofinder_dir)
    raw = pd.read_csv(tsv, sep="\t", index_col=0)
    print(f"  {raw.shape[0]:,} orthogroups x {raw.shape[1]} species")

    print("\n[2] Loading GO ontology ...")
    parent_map = parse_go_obo(args.go_obo)
    primary_ids, alt_to_primary, replaced_by = parse_go_id_remap(args.go_obo)
    print(f"  {len(alt_to_primary):,} alt_id + {len(replaced_by):,} replaced_by "
          f"mappings available for secondary-ID resolution")

    print("\n[3] Loading native GO per species ...")
    species_gene_to_go = {}
    for sp, (relpath, kind, _) in SPECIES_PARSERS.items():
        path = args.go_species_dir / relpath
        if kind == "gaf":
            raw_gene_to_go = parse_athaliana_gaf(path)
        elif kind == "biomart":
            raw_gene_to_go = parse_biomart_tsv(path, lambda x: x)
        elif kind == "fernbase":
            raw_gene_to_go = parse_fernbase_annotation(path)
        elif kind == "eggnog":
            raw_gene_to_go = parse_eggnog_annotations(path)
        species_gene_to_go[sp] = normalize_gene_to_go(
            raw_gene_to_go, primary_ids, alt_to_primary, replaced_by)
        print(f"  {sp}: {len(species_gene_to_go[sp]):,} genes with >=1 resolvable GO term")

    all_terms = set()
    for d in species_gene_to_go.values():
        for s in d.values():
            all_terms |= s
    print(f"\n  Building ancestor map for {len(all_terms):,} distinct GO terms ...")
    all_ancestors = compute_ancestors(all_terms, parent_map)

    print("\n[4] Scanning orthogroups (multi-species) ...")
    og_filters = {name: set() for name in CATEGORY_SETS}

    for og_id, row in raw.iterrows():
        og_terms = set()
        for sp, (_, _, normalizer) in SPECIES_PARSERS.items():
            cell = row.get(sp)
            if pd.isna(cell) or str(cell).strip() == "":
                continue
            gene_to_go = species_gene_to_go[sp]
            for gene_full in str(cell).split(","):
                gene_full = gene_full.strip()
                if not gene_full:
                    continue
                norm_id = normalizer(gene_full)
                if norm_id is None:
                    continue
                og_terms |= gene_to_go.get(norm_id, set())
        if not og_terms:
            continue
        og_ancestor_union = set()
        for t in og_terms:
            og_ancestor_union |= all_ancestors.get(t, {t})
        for cat_name, target_terms in CATEGORY_SETS.items():
            if og_ancestor_union & target_terms:
                og_filters[cat_name].add(og_id)

    total_ogs = raw.shape[0]
    print(f"\n  Total OGs: {total_ogs:,}")
    for name, ogs in og_filters.items():
        print(f"  {name:<22} {len(ogs):>6,} OGs  ({100*len(ogs)/total_ogs:.1f}%)")

    print(f"\n[5] Writing filter files to {args.out_dir}/ ...")
    for name, ogs in og_filters.items():
        out_path = args.out_dir / f"{name}_ogs.txt"
        out_path.write_text("\n".join(sorted(ogs)) + "\n")
        print(f"  {out_path.name}: {len(ogs):,} orthogroups")

    print("\n[6] Comparing against Arabidopsis-only baseline ...")
    comparison_rows = []
    for name in CATEGORY_SETS:
        baseline_path = args.baseline_dir / f"{name}_ogs.txt"
        if baseline_path.exists():
            baseline_ogs = set(baseline_path.read_text().split())
        else:
            baseline_ogs = set()
        multi_ogs = og_filters[name]
        comparison_rows.append({
            "category": name,
            "arabidopsis_only": len(baseline_ogs),
            "multispecies": len(multi_ogs),
            "only_in_arabidopsis": len(baseline_ogs - multi_ogs),
            "only_in_multispecies": len(multi_ogs - baseline_ogs),
            "shared": len(baseline_ogs & multi_ogs),
        })
    comp_df = pd.DataFrame(comparison_rows)
    comp_path = args.out_dir / "comparison_vs_arabidopsis_only.tsv"
    comp_df.to_csv(comp_path, sep="\t", index=False)
    print(comp_df.to_string(index=False))
    print(f"\n  {comp_path}")


if __name__ == "__main__":
    main()
