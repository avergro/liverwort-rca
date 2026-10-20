#!/usr/bin/env python3
"""
Build the supplementary tables from the result files, so that none of them is
transcribed by hand and all of them can be regenerated after any re-analysis.

Tables produced (Submission/supplementary/):
  S1  dated topology (written by 18_make_table_s7.py)
  S2  MCi-raw and MCi-RCA for every species and context
  S3-S8  shared orthogroups, missing GO terms and DYW, sensitivity, KEGG modules, Filter B, WGDs
  S9  species, source, BUSCO categories, genes assigned to orthogroups
  S10 KEGG pathway sets defining each functional context
  S12 degree-preserving null model and nestedness
  S13 phylogenetic regression: all fits plus the sensitivity analyses
  S14 simulation grid and its relaxations

Usage:
    caffeinate -i python3 22_build_supplementary.py
"""

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
CTRL = ROOT / "data" / "results_controls"
SIM = ROOT / "data" / "results_simulation"
SUPP = ROOT / "Submission" / "supplementary"
SHORT = {"Mpolymorpha": "Marchantia polymorpha", "Lcruciata": "Lunularia cruciata",
         "Smoellendorffii": "Selaginella moellendorffii", "Ppatens": "Physcomitrium patens",
         "Afiliculoides": "Azolla filiculoides", "Scucullata": "Salvinia cucullata",
         "Creinhardtii": "Chlamydomonas reinhardtii", "Athaliana": "Arabidopsis thaliana",
         "Slycopersicum": "Solanum lycopersicum", "Pamericana": "Persea americana",
         "Pabies": "Picea abies"}
SOURCE = {"Athaliana": "Ensembl Plants r62", "Creinhardtii": "Ensembl Plants r62",
          "Pamericana": "NCBI GCA_029852735.1", "Slycopersicum": "Ensembl Plants r62",
          "Ppatens": "Ensembl Plants r62", "Pabies": "Picab02 V2, TE genes removed",
          "Mpolymorpha": "Ensembl Plants r62", "Smoellendorffii": "Ensembl Plants r62, alleles collapsed",
          "Afiliculoides": "FernBase v1.1 high-confidence", "Lcruciata": "CoGe gid:64629",
          "Scucullata": "FernBase v1.2"}
CTX_FILES = {"primary_metabolism": "results_primary_metabolism_v3",
             "C_metabolism": "results_C_metabolism_v3",
             "N_metabolism": "results_N_metabolism_v3",
             "regulatory_networks": "results_regulatory_v3"}


# readable labels for the values and column names the analysis scripts write
CONTEXT = {"whole_genome": "Whole genome", "primary_metabolism": "Primary metabolism",
           "C_metabolism": "Carbon metabolism", "N_metabolism": "Nitrogen metabolism",
           "regulatory_networks": "Regulatory context"}
MODEL = {"log_genes": "log gene number", "log_genes+busco": "log gene number + BUSCO"}
# Tables moved out of the supplementary file into Dataset S1 (repository folder supplementary_data)
DATA = ROOT / "supplementary_data"
DATASET = {"Table S10. KEGG": "kegg_pathway_sets.tsv",
           "Table S12a.": "null_model.tsv", "Table S12b.": "nestedness.tsv", "Table S12c.": "eigenvalue_spectrum.tsv",
           "Table S13a.": "pgls_fits.tsv", "Table S13b.": "pgls_sensitivity.tsv", "Table S13c.": "pgls_residual_ranks.tsv",
           "Table S14a.": "simulation_grid.tsv", "Table S14b.": "simulation_relaxations.tsv",
           "Table S4b. The 25 conserved": "rna_modification_orthogroups.tsv",
           "Table S7full.": "liverwort_enriched_pathways.tsv"}
DATASET_NAMES = {"Table S12a": "null_model.tsv", "Table S12b": "nestedness.tsv", "Table S12c": "eigenvalue_spectrum.tsv",
                 "Table S13a": "pgls_fits.tsv", "Table S13b": "pgls_sensitivity.tsv", "Table S13c": "pgls_residual_ranks.tsv",
                 "Table S14a": "simulation_grid.tsv", "Table S14b": "simulation_relaxations.tsv",
                 "Table S7b": "liverwort_enriched_orthogroups.tsv",
                 "Table S10b": "go_term_sets.tsv"}
ESTIMATOR = {"MCi_rca": "RCA score", "MCi_raw": "Presence score", "MCi_rarefied": "RCA score, rarefied",
             "CA_axis1": "Correspondence axis 1", "diversity": "RCA diversity",
             "fitness": "Fitness-Complexity"}
SCENARIO = {"baseline": "Baseline", "phylogeny": "Capabilities inherited along a tree",
            "incomplete": "Incomplete annotation", "neofunction": "Duplicates gain capabilities",
            "all_three": "All three combined"}
KEGG_SET = {"c_pathways_carbohydrate": "Carbon, carbohydrates", "c_pathways_energy": "Carbon, energy",
            "n_pathways_amino_acids": "Nitrogen, amino acids", "n_pathways_assimilation": "Nitrogen, assimilation",
            "lipid_pathways": "Lipids", "nucleotide_pathways": "Nucleotides", "cofactor_pathways": "Cofactors",
            "regulatory_network_pathways": "Regulatory context"}
COLUMNS = {
    "context": "Context", "model": "Model", "setting": "Setting", "species": "Species",
    "estimator": "Estimator", "scenario": "Scenario",
    "lambda_ML": "Pagel λ", "beta_log_genes": "β, log genes", "se_log_genes": "SE, log genes",
    "t_log_genes": "t, log genes", "beta_busco": "β, BUSCO", "t_busco": "t, BUSCO", "logLik": "log-likelihood",
    "MCi_z_observed": "Observed score (standardized)", "null_mean": "Null mean", "null_sd": "Null SD",
    "z_vs_null": "z", "p_one_sided": "P (one-sided)",
    "NODF_observed": "Observed NODF", "NODF_null_mean": "Null mean", "NODF_null_sd": "Null SD",
    "n_orthogroups": "Orthogroups", "lambda2": "λ2", "lambda3": "λ3", "lambda4": "λ4",
    "gap_2_3": "λ2 − λ3", "ratio_3_2": "λ3 / λ2",
    "CA_inertia_axis1": "CA inertia, axis 1", "CA_inertia_axis2": "CA inertia, axis 2",
    "cap_spread": "Endowment heterogeneity", "dup_spread": "Duplication heterogeneity",
    "rho_signal": "ρ, endowment", "rho_size": "ρ, gene number",
    "partial_signal_given_size": "Partial ρ, endowment given gene number",
    "partial_size_given_signal": "Partial ρ, gene number given endowment",
    "rho_completeness": "ρ, completeness", "abs_partial_signal": "Mean absolute partial ρ per replicate",
}


def readable(df):
    """Full names for contexts, models, estimators, scenarios and species; readable headers."""
    df = df.copy()
    for col, mapping in (("context", CONTEXT), ("model", MODEL), ("estimator", ESTIMATOR),
                         ("scenario", SCENARIO), ("species", SHORT)):
        if col in df.columns:
            unknown = set(df[col].astype(str)) - set(mapping)
            assert not unknown, (col, unknown)
            df[col] = df[col].astype(str).map(mapping)
    return df.rename(columns=COLUMNS)


def md_table(df, floatfmt="{:.3f}"):
    def fmt(v, col=""):
        if v is None or (isinstance(v, float) and v != v):   # missing value
            return "–"
        if isinstance(v, float) and ("(%)" in col or col.startswith("BUSCO")):
            return f"{v:.1f}".replace("-", "\u2212")
        if isinstance(v, float) and "λ" in col:
            return f"{v:.2f}"
        if isinstance(v, float) and v.is_integer() and abs(v) >= 2:
            return f"{int(v):,}".replace("-", "\u2212")
        if isinstance(v, float) and v != 0 and abs(v) < 1e-3:
            m, e = f"{v:.2e}".split("e")
            return f"{float(m):.2g} × 10^{int(e)}^".replace("-", "\u2212")
        if isinstance(v, float):
            return floatfmt.format(v).replace("-", "\u2212")
        if isinstance(v, (int, np.integer)) and abs(v) >= 1000:
            return f"{v:,}"
        if isinstance(v, str) and v == "Chlamydomonas removed":
            return "*Chlamydomonas* removed"
        if col in ("Species", "species") and isinstance(v, str) and v and not v.startswith("*") and " " in v:
            return f"*{v}*"
        if isinstance(v, str) and v.strip().lower() in ("", "nan"):
            return "–"
        return str(v)
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "|" + "---|" * len(df.columns)
    rows = ["| " + " | ".join(fmt(v, str(c)) for v, c in zip(r, df.columns)) + " |"
            for r in df.itertuples(index=False)]
    return "\n".join([head, rule] + rows)


def busco_table():
    """BUSCO categories from data/busco/busco_summary.tsv (03b_busco_summary.py)."""
    b = pd.read_csv(ROOT / "data/busco/busco_summary.tsv", sep="\t")
    return pd.DataFrame({"species": b.Species, "lineage_dataset": b["BUSCO dataset"], "C": b["C (%)"],
                         "S": b["S (%)"], "D": b["D (%)"], "F": b["F (%)"], "M": b["M (%)"]})


def busco_table_from_json():
    rows = []
    for d in sorted((ROOT / "data" / "busco").iterdir()):
        js = list(d.glob("short_summary*.json"))
        if not js:
            continue
        data = json.loads(js[0].read_text())
        res = data.get("results", {})
        sp = d.name.replace("_chlorophyta", "")
        rows.append({"species": SHORT.get(sp, sp), "dir": d.name,
                     "lineage_dataset": data.get("lineage_dataset", {}).get("name", ""),
                     "C": res.get("Complete percentage"), "S": res.get("Single copy percentage"),
                     "D": res.get("Multi copy percentage"), "F": res.get("Fragmented percentage"),
                     "M": res.get("Missing percentage")})
    df = pd.DataFrame(rows)
    # Chlamydomonas was assessed twice; keep the chlorophyta run, which is the
    # applicable lineage dataset and the value reported in the manuscript
    df = df.sort_values("dir", key=lambda c: c.str.contains("chlorophyta"))
    df = df.drop_duplicates("species", keep="last").drop(columns="dir")
    return df


KEGG_NAME = {"ath00910": "Nitrogen metabolism"}


def kegg_sets():
    src = (ROOT / "06_kegg_filter.py").read_text()
    out = []
    for name in re.findall(r"^([A-Z_]+PATHWAYS?[A-Z_]*) = \{", src, re.M):
        block = src.split(f"{name} = {{", 1)[1].split("}", 1)[0]
        ids = re.findall(r'"(ath\d+)",\s*#\s*(.+)', block)
        for pid, desc in ids:
            # code comments such as "(included once)" or "(already in lipids)" are notes, not names
            desc = re.sub(r"\s*\((included once|already in [^)]*)\)", "", desc).strip()
            out.append({"set": name.lower(), "pathway": pid, "description": desc or KEGG_NAME.get(pid, pid)})
    return pd.DataFrame(out).drop_duplicates(["set", "pathway"])


def main():
    SUPP.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    bm = importlib.util.module_from_spec(spec); spec.loader.exec_module(bm)
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    genes = counts.sum(axis=0)

    # ── S1 ────────────────────────────────────────────────────────────────
    b = busco_table().set_index("species")
    s1 = pd.DataFrame({
        "Species": [SHORT[s] for s in genes.index],
        "Source": [SOURCE[s] for s in genes.index],
        "Genes in orthogroups": genes.values,
    })
    for col in ("lineage_dataset", "C", "S", "D", "F", "M"):
        s1[f"{col} (%)" if col != "lineage_dataset" else "BUSCO dataset"] = [
            b[col].get(SHORT[s], "") for s in genes.index]
    s1 = s1.sort_values("Species")

    # ── S3 ────────────────────────────────────────────────────────────────
    whole = pd.read_csv(ROOT / "data/results_v3/species_complexity.tsv", sep="\t", index_col=0)
    rca_div = pd.read_csv(ROOT / "data/results_v3/matrix_rca.tsv", sep="\t", index_col=0).sum()
    s3 = pd.DataFrame({"Species": [SHORT[s] for s in whole.index],
                       "Diversity (whole genome)": rca_div.reindex(whole.index).values,
                       "Presence score": whole.MCi_raw.values, "RCA score, whole genome": whole.MCi_rca.values})
    for ctx, d in CTX_FILES.items():
        v = pd.read_csv(ROOT / "data" / d / "species_complexity.tsv", sep="\t", index_col=0)
        s3[f"RCA score, {CONTEXT[ctx].lower()}"] = v.MCi_rca.reindex(whole.index).values
    s3 = s3.sort_values("RCA score, whole genome", ascending=False)

    # ── S4, S5, S6 ────────────────────────────────────────────────────────
    s4a = readable(pd.read_csv(CTRL / "F_pgls_fits.tsv", sep="\t"))
    s4b = readable(pd.read_csv(CTRL / "H_pgls_sensitivity_fits.tsv", sep="\t"))
    # residual ranks per context: main tree, and range over the sensitivity settings
    rk = pd.read_csv(CTRL / "F_pgls_residuals.tsv", sep="\t")
    hr = pd.read_csv(CTRL / "H_pgls_sensitivity_ranks.tsv", sep="\t")
    ctx_order = ["whole_genome", "primary_metabolism", "C_metabolism", "N_metabolism", "regulatory_networks"]
    ctx_name = {"whole_genome": "Whole genome", "primary_metabolism": "Primary metabolism",
                "C_metabolism": "Carbon metabolism", "N_metabolism": "Nitrogen metabolism",
                "regulatory_networks": "Regulatory context"}
    s4c = pd.DataFrame({"Species": [f"*{SHORT[x]}*" for x in sorted(rk.species.unique(), key=lambda x: SHORT[x])]})
    for c in ctx_order:
        main_r = rk[rk.context == c].set_index("species").rank_residual
        rng = hr[hr.context == c].groupby("species").rank_residual.agg(["min", "max"])
        col = []
        for sp in sorted(rk.species.unique(), key=lambda x: SHORT[x]):
            lo, hi = int(rng.loc[sp, "min"]), int(rng.loc[sp, "max"])
            col.append(f"{int(main_r[sp])} ({lo}–{hi})" if lo != hi else f"{int(main_r[sp])}")
        s4c[ctx_name[c]] = col
    s5a = readable(pd.read_csv(CTRL / "D_null_model.tsv", sep="\t"))
    s5b = readable(pd.read_csv(CTRL / "I_nestedness.tsv", sep="\t"))
    s5c = readable(pd.read_csv(CTRL / "I_spectrum.tsv", sep="\t"))
    s6a = readable(pd.read_csv(SIM / "S1_grid_summary.tsv", sep="\t"))
    s6b = readable(pd.read_csv(SIM / "S2_extensions_summary.tsv", sep="\t"))
    s2 = kegg_sets()
    s2["set"] = s2["set"].map(KEGG_SET)
    s2 = s2.rename(columns={"set": "Set", "pathway": "KEGG pathway", "description": "Description"})

    DATA.mkdir(exist_ok=True)
    for old in DATA.glob("*.tsv"):
        old.unlink()
    readme = []

    def dataset_file(fn, title, note, df):
        """Write one table of Dataset S1 (repository folder supplementary_data) and list it in its README."""
        df.to_csv(DATA / fn, sep="\t", index=False)
        for k, v in DATASET_NAMES.items():
            note = note.replace(k, v)
        readme.append(f"- `{fn}`: {title}." + (f" {note}" if note.strip() else ""))

    doc = ["# Supplementary tables", "",
           "All tables are generated from the analysis outputs deposited in the repository (Data and code availability).", ""]

    def section(title, note, df, fmt="{:.3f}"):
        note = re.sub(r"\s*Generated by `[^.]*\.", "", note)          # script names live in the repository
        for prefix, fn in DATASET.items():                               # detail tables go to Dataset S1
            if title.startswith(prefix):
                dataset_file(fn, title.split(". ", 1)[1], note, df)
                return
        doc.extend([f"## {title}", "", note, "", md_table(df, fmt), ""])

    section("Table S9. Species, sources and proteome completeness",
            "BUSCO categories are percentages: complete (C), complete single-copy (S), complete "
            "duplicated (D), fragmented (F) and missing (M). *Chlamydomonas reinhardtii* was assessed "
            "against the chlorophyta lineage dataset and the remaining species against embryophyta.",
            s1, "{:.1f}")
    section("Table S10. KEGG pathway sets defining each functional context",
            "Sets are those that define the KEGG functional contexts; primary metabolism is the union of all "
            "metabolic sets, and ubiquinone biosynthesis, listed under both lipids and cofactors, is counted once. "
            "The equivalent Gene Ontology term sets "
            "are listed in Table S10b.", s2)
    section("Table S2. Ranking scores for every species and context",
            "RCA diversity is the number of orthogroups with RCA at or above one in the whole genome. "
            "The presence score is computed on gene presence and the RCA score after RCA "
            "normalization. Missing values are shown as –.",
            s3)
    section("Table S13a. Phylogenetic regression, all fits",
            "Phylogenetic generalized least squares regression of the RCA score of each context on the natural "
            "logarithm of gene number, alone or with BUSCO completeness, on the dated topology of Table S1. Pagel λ "
            "is estimated by maximum likelihood between 0 and 1. The log-likelihood (logL) is computed with the "
            "residual variance estimated with n − p degrees of freedom, which shifts every profile by a constant "
            "and leaves the estimate of λ unchanged.", s4a, "{:.4f}")
    section("Table S13b. Phylogenetic regression, sensitivity to the outgroup and the root age",
            "The regression of Table S13a repeated without the alga, with the root of the topology at 750, 1,000 and "
            "1,500 million years, and with every land plant node at the younger or at the older bound of the "
            "interval given by its source in Table S1.", s4b, "{:.4f}")
    section("Table S13c. Residual rank of each species in the phylogenetic regression",
            "Rank of the residual of each species (1, largest positive residual) in the regression of the RCA score "
            "on the natural logarithm of gene number on the dated topology of Table S1. In parentheses, the range of "
            "ranks across the settings of Table S13b (without the alga, root at 750, 1,000 and 1,500 million years, "
            "and every land plant node at the younger or at the older bound of its interval); no parentheses when "
            "the rank did not change.", s4c)
    section("Table S12a. Degree-preserving null model, all species and contexts",
            "z scores are reported for reference; the standard deviations of some null "
            "distributions are very small and the corresponding z values should not be "
            "interpreted as effect sizes.", s5a, "{:.4f}")
    section("Table S12b. Nestedness against the null model",
            "NODF computed exactly by grouping orthogroups into binary patterns.", s5b, "{:.3f}")
    section("Table S12c. Eigenvalue spectrum and correspondence analysis inertia",
            "Eigenvalues of the species-by-species matrix T of equation 2 in the main text (the first, equal to one, is omitted) "
            "and the share of inertia carried by the first two axes of a correspondence analysis of the same "
            "binary matrix.", s5c)
    section("Table S14a. Simulation grid", "Recovery is the partial Spearman correlation with the "
            "planted capability endowment, controlling for gene number.", s6a)
    section("Table S14b. Relaxations of the simulation model",
            "Recovery under each relaxation of the capability model, given as the partial Spearman correlation "
            "with the planted endowment, controlling for gene number, with its sign. Scores are oriented by RCA diversity, as in the main text; "
            "Table S14a reports both the signed and the absolute value.", s6b)
    doc.extend(["## Table S1. Dated topology", "",
                "The dated topology and its variance-covariance matrix are given in the accompanying Table S1 file.", ""])

    # ── S11: pathway composition of Filter B (32_filterB_enrichment.py) ──────
    fb = pd.read_csv(ROOT / "data/results_filterB/filterB_enrichment.tsv", sep="\t")
    fb["context"] = fb.context.map(CONTEXT)
    s11a = fb.rename(columns={
        "context": "Context", "pathway_name": "KEGG pathway", "n_filterB": "Liverwort-enriched orthogroups",
        "n_context": "Context orthogroups", "fold_vs_context": "Fold, context",
        "q_hypergeom": "q, context", "fold_vs_null": "Fold, permutation null",
        "q_permutation": "q, permutation null", "OR_MH": "Odds ratio, stratified",
        "q_MH": "q, stratified"})[["Context", "KEGG pathway", "Liverwort-enriched orthogroups", "Context orthogroups",
                                   "Fold, context", "q, context", "Fold, permutation null",
                                   "q, permutation null", "Odds ratio, stratified", "q, stratified"]]
    s11b = pd.read_csv(ROOT / "data/results_filterB/filterB_orthogroups.tsv", sep="\t")
    ns = pd.read_csv(ROOT / "data/results_filterB/filterB_null_size.tsv", sep="\t")
    obs = s11b.context.value_counts()
    null_note = "; ".join(
        f"{g.filterB_size.mean():.0f} orthogroups on average in {CONTEXT[c].lower()} (95% interval "
        f"{np.percentile(g.filterB_size, 2.5):.0f} to {np.percentile(g.filterB_size, 97.5):.0f}), against "
        f"{obs[c]} observed" for c, g in ns.groupby("context", sort=False))
    s11b["context"] = s11b.context.map(CONTEXT)
    section("Table S7full. Pathway composition of the liverwort-enriched orthogroups, all pathways with at least three",
            "Columns as in Table S7 of the supplementary material.", s11a)
    top = fb.sort_values("n_filterB", ascending=False, kind="stable").groupby("context", sort=False).head(7)
    keep = top.index.union(fb.index[(fb.context == "Carbon metabolism") & (fb.pathway_name == "Pentose phosphate pathway")])
    s11a = s11a.loc[s11a.index.isin(keep)]
    section("Table S7. Pathway composition of the liverwort-enriched orthogroups",
            "Liverwort-enriched orthogroups have RCA ≥ 1 in both liverworts and RCA < 1 in all three angiosperms "
            "(section 4.5). The seven pathways with most of them in each context are listed, together with the "
            "pentose phosphate pathway of carbon metabolism, cited in the main text; all pathways with at least three "
            "are given in Dataset S1. Fold is the share of the set "
            "taken by the pathway divided by its share of the context orthogroups (hypergeometric test) or "
            "by its mean share in liverwort-enriched sets recomputed after shuffling the copy numbers of each species "
            "among the orthogroups it occupies (2,000 permutations). The stratified odds ratio is the "
            "Cochran-Mantel-Haenszel estimate among orthogroups occupied by both liverworts, stratified by "
            "mean angiosperm copy number. q values are Benjamini-Hochberg corrected within each context. "
            "The orthogroups, their copy numbers, KEGG pathways and gene identifiers are listed in "
            "Dataset S1. "
            f"Under the copy-permutation null the liverwort-enriched set had {null_note}. Generated by "
            "`32_filterB_enrichment.py`.", s11a)
    # ── S12: KEGG module completeness (33_kegg_module_completeness.py) ───────
    km = pd.read_csv(ROOT / "data/results_kegg_modules/species_summary.tsv", sep="\t")
    km["species"] = km.species.map(SHORT)
    s12a = km.merge(s1[["Species", "C (%)"]], left_on="species", right_on="Species").rename(columns={
        "C (%)": "BUSCO C (%)", "kos": "KOs", "genes": "Genes with a KO", "genes_per_ko": "Genes per KO",
        "complete": "Complete modules", "near_complete": "Near-complete modules",
        "mean_fraction": "Mean step completeness"})[
        ["Species", "BUSCO C (%)", "KOs", "Genes with a KO", "Genes per KO", "Complete modules",
         "Near-complete modules", "Mean step completeness"]].sort_values("Genes per KO")
    ev = pd.read_csv(ROOT / "data/results_kegg_modules/module_completeness.tsv", sep="\t")
    bh = pd.read_csv(CTRL / "M_kegg_best_hit.tsv", sep="\t").set_index("assignment").loc["best-scoring KO per gene"]
    best_note = (f"With each gene assigned only to its best-scoring KO, {bh.modules_complete_in_a_land_plant} modules "
                 f"were complete in at least one land plant and {bh.complete_in_all_three_angiosperms} in the three "
                 f"angiosperms, of which *Marchantia* completed {bh.of_those_complete_in_Mpolymorpha}, with the same six "
                 "gaps and one more (M00987, sulfate assimilation).")
    land = [s for s in SHORT if s != "Creinhardtii"]
    rel = ev[ev.species.isin(land) & ev.complete].module.unique()
    s12b = ev[ev.module.isin(rel)].pivot_table(index=["module", "name", "category"], columns="species",
                                               values="fraction").round(2).rename(columns=SHORT).reset_index()
    section("Table S5. KEGG module completeness",
            "Every proteome was annotated with the KOfam profiles of the KEGG pathway modules and their "
            "score thresholds. Values are computed over the 118 modules complete in at least one land "
            "plant and the KEGG orthologs (KOs) that define them. A module is complete when every step "
            "has a gene and near-complete when at most one step is missing. The fraction of steps "
            "satisfied for every module and species is listed in Dataset S1. In *Marchantia*, "
            "the six modules complete in the three angiosperms but incomplete in *Marchantia* had no gene "
            "assigned to ferulate 5-hydroxylase and caffeic acid O-methyltransferase (monolignols), allene "
            "oxide synthase and OPC-8:0 CoA ligase (jasmonic acid), the three enzymes from naringenin to "
            "pelargonidin, ACC synthase and ACC oxidase (ethylene), xanthoxin dehydrogenase and abscisic aldehyde oxidase (abscisic acid), "
            "and gibberellin 20-oxidase (gibberellins A4 and A1). "
            f"{best_note} Generated by `33_kegg_module_completeness.py`.", s12a)
    # Table S8: whole-genome duplications reported in the literature, beside redundancy
    WGD = [  # species, lineage, reported whole-genome duplications, sources
        ("Creinhardtii", "Green alga", "None reported; 317 of 798 two-gene families arranged in tandem",
         "Merchant et al. 2007"),
        ("Mpolymorpha", "Liverwort", "None", "Bowman et al. 2017"),
        ("Lcruciata", "Liverwort", "None reported in liverworts", "Linde et al. 2023; Lang et al. 2018"),
        ("Scucullata", "Fern", "One, at the base of core leptosporangiate ferns", "Li et al. 2018"),
        ("Afiliculoides", "Fern", "Two: the core leptosporangiate event and one specific to *Azolla*",
         "Li et al. 2018"),
        ("Smoellendorffii", "Lycophyte",
         "None in the genome analysis; two paleopolyploidizations (13 to 15 and 125 to 142 Ma) proposed later",
         "Banks et al. 2011; Wang et al. 2020"),
        ("Ppatens", "Moss", "Two", "Lang et al. 2018"),
        ("Athaliana", "Angiosperm",
         "Ancestral seed plant (ζ) and angiosperm (ε) events, the core eudicot triplication (γ) and two "
         "in the crucifer lineage (β, α)", "Jiao et al. 2011"),
        ("Slycopersicum", "Angiosperm", "ζ, ε, γ and a triplication in the *Solanum* lineage",
         "Jiao et al. 2011; Tomato Genome Consortium 2012"),
        ("Pamericana", "Angiosperm", "ζ, ε and two in the avocado lineage",
         "Jiao et al. 2011; Nath et al. 2022"),
        ("Pabies", "Gymnosperm", "Only the ancestral seed plant event (ζ); no recent duplication",
         "Jiao et al. 2011; Nystedt et al. 2013; Ahlgren Kalman et al. 2025"),
    ]
    cov = pd.read_csv(CTRL / "L_coverage_redundancy.tsv", sep="\t")
    gpo = cov[cov.context == "whole_genome"].set_index("species").genes_per_occupied_OG
    s13 = pd.DataFrame([{"Species": f"*{SHORT[sp]}*", "Lineage": lin,
                         "Genes per occupied orthogroup": gpo[sp],
                         "Reported whole-genome duplications": w, "Sources": src}
                        for sp, lin, w, src in WGD]).sort_values("Genes per occupied orthogroup")
    section("Table S8. Whole-genome duplications reported for each species",
            "Genes per occupied orthogroup are those of Table 3 (whole genome, 11,383 orthogroups present "
            "in at least two species). Duplications are those reported in the genome papers and in the "
            "cited analyses of ancient polyploidy; ζ and ε are the events ancestral to seed plants and to "
            "angiosperms. Sources: Ahlgren Kalman et al. 2025, Research Square, "
            "doi:10.21203/rs.3.rs-6502828/v1; Banks et al. 2011, Science 332:960; "
            "Bowman et al. 2017, Cell 171:287; Jiao et al. 2011, "
            "Nature 473:97; Lang et al. 2018, Plant J 93:515; Li et al. 2018, Nat Plants 4:460; "
            "Linde et al. 2023, Genome Biol Evol 15:evad014; Merchant et al. 2007, Science 318:245; "
            "Nath et al. 2022, Hortic Res 9:uhac157; Nystedt et al. 2013, Nature 497:579; Tomato Genome "
            "Consortium 2012, Nature 485:635; Wang et al. 2020, Genomics Proteomics Bioinformatics "
            "18:333.", s13, "{:.2f}")
    # ── S14 and S15: overlap with the angiosperms and sensitivity (37_overlap_and_sensitivity.py) ──
    ov = pd.read_csv(CTRL / "M_overlap_with_angiosperms.tsv", sep="\t")
    exp = pd.read_csv(CTRL / "M_busco_expected_loss.tsv", sep="\t")
    exp = exp[exp.unit == "orthogroups"].set_index("species")
    ref = exp.drop(index=["Mpolymorpha", "Lcruciata"])
    b, a = np.polyfit(ref.busco, ref.observed, 1)
    n, xm = len(ref), ref.busco.mean()
    sd = np.sqrt(((ref.observed - (a + b * ref.busco)) ** 2).sum() / (n - 2))
    half = stats.t.ppf(0.975, n - 2) * sd * np.sqrt(1 + 1 / n + (exp.busco - xm) ** 2 / ((ref.busco - xm) ** 2).sum())
    ov["expected_from_busco"] = ov.species.map(exp.expected)
    ov["pi_low"] = ov.species.map(exp.expected - half).where(ov.species.isin(["Mpolymorpha", "Lcruciata"]))
    ov["pi_high"] = ov.species.map(exp.expected + half).where(ov.species.isin(["Mpolymorpha", "Lcruciata"]))
    ov["species"] = ov.species.map(lambda x: f"*{SHORT[x]}*")
    s14 = ov.rename(columns={
        "species": "Species", "occupied": "Orthogroups occupied",
        "missing_of_all_three_angiosperms": "Missing of the 7,310 in all three angiosperms",
        "fraction_missing_of_all_three": "Fraction missing",
        "missing_of_conserved": "Missing of the conserved set",
        "fraction_missing_of_conserved": "Fraction missing, conserved",
        "absent_from_all_angiosperms": "Occupied, absent from all angiosperms",
        "in_all_land_plants_outside_lineage": "Reference outside the lineage",
        "missing_of_those_outside_lineage": "Missing of the reference",
        "fraction_missing_outside_lineage": "Fraction missing, reference",
        "expected_from_busco": "Expected from BUSCO",
        "pi_low": "95% prediction interval, lower", "pi_high": "95% prediction interval, upper",
        "jaccard_with_Athaliana": "Jaccard with *Arabidopsis*"})[
        ["Species", "Orthogroups occupied", "Missing of the 7,310 in all three angiosperms", "Fraction missing",
         "Missing of the conserved set", "Fraction missing, conserved", "Occupied, absent from all angiosperms",
         "Reference outside the lineage", "Missing of the reference", "Fraction missing, reference",
         "Expected from BUSCO", "95% prediction interval, lower", "95% prediction interval, upper",
         "Jaccard with *Arabidopsis*"]].sort_values("Fraction missing, reference")
    section("Table S3a. Orthogroups shared with the angiosperms, by species",
            "Computed over the 11,383 orthogroups present in at least two species. 7,310 orthogroups are present "
            "in all three angiosperms. The conserved set adds the condition that the orthogroup is also present in "
            "at least one other non-seed land plant (6,694 to 6,725 orthogroups; the species examined is excluded from "
            "the condition). The reference outside the lineage comprises the orthogroups present in every land plant "
            "outside the lineage of the species (liverworts, ferns, angiosperms, or the species alone), so that "
            "losses shared within a lineage are counted. Expected from BUSCO is the fraction of the reference missing "
            "predicted by linear regression on BUSCO completeness across the eight land plants outside the liverworts "
            "(R² = 0.75), with the 95% prediction interval of that regression for the two liverworts. "
            "Values marked – are not defined for the angiosperms.", s14)
    ll = pd.read_csv(CTRL / "M_lineage_losses.tsv", sep="\t")
    lab = {"liverworts": "Liverworts", "angiosperms": "Angiosperms", "ferns": "Ferns", "moss": "*Physcomitrium*",
           "lycophyte": "*Selaginella*", "gymnosperm": "*Picea*", "bryophytes": "Bryophytes",
           "seed_plants": "Seed plants"}
    rows = []
    for _, r in ll.iterrows():
        each = ", ".join(f"*{SHORT[c[8:]].split()[0]}* {int(r[c])}" for c in ll.columns
                         if c.startswith("missing_") and c != "missing_in_whole_lineage" and r[c] == r[c])
        rows.append({"Lineage": lab[r.lineage], "Unit": r.unit, "Present in every land plant outside the lineage":
                     int(r.reference), "Absent from the whole lineage": int(r.missing_in_whole_lineage),
                     "Absent from each species": each})
    section("Table S3b. Orthogroups, KEGG orthologs and complete modules lost by each lineage",
            "For each lineage, the reference comprises the units present (orthogroups and KEGG orthologs) or complete "
            "(modules) in every land plant outside it; the columns count those absent from all members of the "
            "lineage and from each member. Bryophytes and seed plants are added as wider lineages. KEGG orthologs "
            "are the 825 that define the modules complete in at least one land plant.", pd.DataFrame(rows))
    dd = pd.read_csv(CTRL / "M_deficit_decomposition.tsv", sep="\t")
    s14c = pd.DataFrame({"Liverwort": dd.liverwort.map(lambda x: f"*{SHORT[x]}*"),
                         "Angiosperm": dd.angiosperm.map(lambda x: f"*{SHORT[x]}*"),
                         "Gene difference": dd.gene_difference,
                         "Copies in shared orthogroups": dd.copies_in_shared_orthogroups,
                         "Orthogroups the liverwort lacks, present in another non-seed land plant":
                             dd.orthogroups_lacking_but_in_other_non_seed,
                         "Orthogroups absent from all non-seed land plants": dd.orthogroups_restricted_to_seed_plants,
                         "Orthogroups of the liverwort absent from the angiosperm":
                             dd.orthogroups_of_the_liverwort_absent_from_angiosperm})
    section("Table S3c. Decomposition of the difference in gene number between each angiosperm and each liverwort",
            "Over the 11,383 orthogroups present in at least two species. The four parts are fractions of the gene "
            "difference and sum to one; the last is negative because it counts liverwort genes in orthogroups the "
            "angiosperm lacks.", s14c)
    sens = pd.read_csv(CTRL / "M_sensitivity_ratios.tsv", sep="\t")
    sens["set"] = sens["set"].map({
        "orthogroups (reported)": "Orthogroups (main analysis)",
        "hierarchical orthogroups, land plant node (N1)": "Hierarchical orthogroups at the land plant node",
        "genes < 0.5 x orthogroup median length removed": "Gene models shorter than half the orthogroup median removed",
        "lineage-specific orthogroups removed": "Lineage-specific orthogroups removed"})
    s15 = sens.rename(columns={"set": "Orthogroup set", "coverage_ratio_L_A": "Orthogroups occupied, L / A",
                               "genes_per_OG_ratio_L_A": "Genes per occupied orthogroup, L / A"})
    section("Table S6. Sensitivity of the coverage and redundancy ratios",
            "Ratio of the mean of the two liverworts (L) to the mean of the three angiosperms (A). "
            "Hierarchical orthogroups are those inferred by OrthoFinder at the node of the land plants, with gene trees reconciled with the accepted species tree. "
            "Fragmentary genes are those shorter than half the median protein length of their orthogroup. "
            "Lineage-specific orthogroups are, for each species, those occupied only by its own lineage "
            "(liverworts, ferns or angiosperms) or by the species alone. Generated by "
            "`37_overlap_and_sensitivity.py`.", s15)
    hl = pd.read_csv(CTRL / "M_hog_lineages.tsv", sep="\t").rename(columns={
        "orthogroups": "Orthogroups", "hierarchical_orthogroups": "Hierarchical orthogroups",
        "liverwort_genes_per_angiosperm_lineage_N8": "Liverwort genes per angiosperm gene lineage",
        "liverwort_genes_per_seed_plant_lineage_N7": "Liverwort genes per seed plant gene lineage",
        "angiosperm_genes_per_angiosperm_lineage_N8": "Angiosperm genes per angiosperm gene lineage"})
    section("Table S6b. Liverwort genes per ancestral gene lineage of angiosperms and seed plants",
            "Land plant hierarchical orthogroups present in all three angiosperms and both liverworts. A gene lineage "
            "is a hierarchical orthogroup delimited by OrthoFinder at the node of the angiosperms or of the seed "
            "plants; values near one mean that the liverworts have about one gene for each gene that the common "
            "ancestor of the angiosperms (or seed plants) is inferred to have had. Orthogroups with and without an "
            "*Arabidopsis* transcription factor (PlantTFDB) are shown separately.", hl)
    st = pd.read_csv(CTRL / "M_regulatory_stratified.tsv", sep="\t")
    crt = pd.read_csv(ROOT / "data/results_regulatory_contrast/copy_ratio_tests.tsv", sep="\t").set_index(["A", "B"])
    cr_note = (f"{crt.loc[('regulatory (KEGG)', 'metabolic (KEGG)'), 'copy_ratio_A']:.2f} in the KEGG regulatory "
               f"context, {crt.loc[('regulatory (GO)', 'metabolic (GO)'), 'copy_ratio_A']:.2f} in the Gene Ontology "
               f"regulatory context and {crt.loc[('transcription factors', 'metabolic (KEGG)'), 'copy_ratio_A']:.2f} "
               "in transcription factor orthogroups, against "
               f"{crt.loc[('regulatory (KEGG)', 'metabolic (KEGG)'), 'copy_ratio_B']:.2f} in KEGG and Gene Ontology "
               "metabolism")
    for c in [c for c in st.columns if c.startswith("median")]:
        st[c] = st[c].astype(float)
    for c in [c for c in st.columns if c.startswith("P ")]:
        st[c] = st[c].map(lambda p: "–" if pd.isna(p) else ("< 0.001" if p < 0.001 else f"{p:.2f}"))
    st.columns = [c[0].upper() + c[1:] for c in st.columns.str.replace("_", " ")]
    section("Table S6c. Copy ratio of regulatory and metabolic orthogroups within strata of angiosperm copy number",
            "Orthogroups occupied by at least one liverwort and one angiosperm. The copy ratio is the mean number of "
            "genes among the liverworts that occupy an orthogroup divided by the same mean among the angiosperms; "
            "strata are defined by that angiosperm mean, as in the stratified test of section 4.5. P is the two-sided "
            "Mann-Whitney test of the regulatory against the metabolic copy ratios within each stratum. Orthogroups that "
            "belong to both categories of a comparison are removed, as in the permutation test. Over all strata, the "
            f"median copy ratio was {cr_note} (permutation test of the difference of medians, 10,000 permutations, "
            "P < 0.001 for each comparison). Generated by `37_overlap_and_sensitivity.py` and "
            "`27_regulatory_contrast.py`.", st)
    sc = pd.read_csv(ROOT / "data/results_regulatory_contrast/scaling.tsv", sep="\t")
    sc = sc[~sc.category.str.startswith("KEGG ") | sc.category.isin(
        ["KEGG primary", "KEGG carbon", "KEGG nitrogen", "KEGG regulatory"])]
    s5d = sc.rename(columns={"category": "Category", "species": "Species set", "n": "n", "slope_OLS": "Slope, OLS",
                             "se_OLS": "SE, OLS", "slope_GLS": "Slope, GLS", "se_GLS": "SE, GLS"})
    s5d["Species set"] = s5d["Species set"].map({"all eleven": "All eleven", "land plants": "Ten land plants"})
    s5d["Category"] = s5d.Category.str.replace("GO ", "Gene Ontology ").str.replace(
        "transcription factors", "Transcription factor orthogroups")
    section("Table S6d. Scaling of category gene number with total gene number",
            "Slope of the regression of log gene number in each category on the log number of genes in the 11,383 "
            "orthogroups present in at least two species, by ordinary "
            "least squares (OLS) and by generalized least squares under Brownian motion on the dated topology of "
            "Table S1 (GLS). A slope of one means that the category grows in proportion to the genome. Generated by "
            "`27_regulatory_contrast.py`.", s5d)
    go = pd.read_csv(CTRL / "M_missing_GO_enrichment.tsv", sep="\t")
    s16 = go[go.q < 0.25].rename(columns={
        "GO": "GO term", "term": "Name", "missing_orthogroups": "Missing from both liverworts",
        "background_orthogroups": "Conserved orthogroups", "fold": "Fold", "p": "P", "q": "q",
        "with_PPR_protein": "With a PPR protein"})[
        ["GO term", "Name", "Missing from both liverworts", "Conserved orthogroups", "Fold", "P", "q",
         "With a PPR protein"]]
    section("Table S4. Gene Ontology terms over-represented among the conserved orthogroups absent from both liverworts",
            "Conserved orthogroups are the 6,700 present in all three angiosperms and in at least one of "
            "*Physcomitrium*, *Selaginella*, *Azolla* and *Salvinia*; 411 of them are absent from both liverworts. "
            "Terms are the biological process annotations of their *Arabidopsis* members, propagated to ancestor "
            "terms. P values are from the hypergeometric test, and q values are Benjamini-Hochberg corrected over "
            "the terms with at least ten conserved orthogroups; terms with q < 0.25 are listed. The last column "
            "counts the orthogroups whose *Arabidopsis* members include a pentatricopeptide repeat (PPR) protein. "
            "Generated by `37_overlap_and_sensitivity.py`.", s16, "{:.3g}")
    det = pd.read_csv(CTRL / "M_missing_GO_top_term_orthogroups.tsv", sep="\t")
    spc = [c for c in det.columns if c in SHORT]
    s16b = det.rename(columns={"orthogroup": "Orthogroup", "PPR_protein": "PPR protein",
                               "evidence_IEA_or_IBA": "Annotations IEA or IBA", "evidence_other": "Other annotations"})
    s16b["PPR protein"] = s16b["PPR protein"].map({True: "yes", False: "no"})
    s16b = s16b.rename(columns={c: f"*{SHORT[c].split()[0]}*" for c in spc})
    section("Table S4b. The 25 conserved orthogroups annotated to RNA modification and absent from both liverworts",
            "Gene numbers per species, whether the *Arabidopsis* members include a pentatricopeptide repeat (PPR) "
            "protein, and the evidence codes of the direct annotations of those members to RNA modification or its "
            "descendant terms (IEA, inferred from electronic annotation; IBA, inferred from biological aspect of "
            "ancestor).", s16b)
    dyw = pd.read_csv(CTRL / "M_dyw_counts.tsv", sep="\t")
    dyw["species"] = dyw.species.map(lambda x: f"*{SHORT[x]}*")
    dog = pd.read_csv(CTRL / "M_dyw_orthogroups.tsv", sep="\t")
    moss = dog[dog.get("Ppatens", 0) > 0]
    note = "; ".join(f"{r.orthogroup}, {int(r.Ppatens)} moss protein{'s' if int(r.Ppatens) != 1 else ''}, {int(r.angiosperm_genes_in_orthogroup)} angiosperm genes"
                     for r in moss.itertuples())
    section("Table S4b. Proteins ending in the canonical DYW domain",
            "Proteins whose last 12 residues match the motif C-x-C-x(1,2)-D-Y-W, a conservative criterion that "
            "underestimates the DYW proteins described for *Arabidopsis*. Orthogroups containing the DYW proteins of "
            f"*Physcomitrium*: {note}.",
            dyw.rename(columns={"species": "Species", "DYW_proteins": "DYW proteins"}))
    # ── the rest of Dataset S1: data files, the variance-covariance matrix and the tables cited in Note S1 ──
    def md_to_df(md):
        rows = [[c.strip() for c in l.strip().strip("|").split("|")] for l in md.splitlines() if l.startswith("|")]
        rows = [r for r in rows if not all(set(c) <= set(":- ") for c in r)]
        return pd.DataFrame(rows[1:], columns=rows[0])

    def md_section(text, tid):
        block = [b for b in re.split(r"\n(?=#+ Table )", text) if re.match(r"#+ Table " + tid + r"\.", b)][0]
        lines = block.split("\n")
        title = re.sub(r"#+ Table \S+\. ", "", lines[0])
        note = " ".join(l.strip() for l in lines[1:] if l.strip() and not l.startswith("|"))
        return title, note, "\n".join(l for l in lines if l.startswith("|"))

    dataset_file("kegg_module_step_fractions.tsv", "Fraction of the steps of each KEGG module satisfied in each species",
                 "Modules complete in at least one land plant (Table S5).", s12b)
    dataset_file("liverwort_enriched_orthogroups.tsv",
                 "Liverwort-enriched orthogroups with their copy numbers, KEGG pathways and gene identifiers",
                 "Orthogroups with RCA ≥ 1 in both liverworts and RCA < 1 in all three angiosperms (Table S7).", s11b)
    vcv = pd.read_csv(CTRL / "TableS7_vcv.tsv", sep="\t")
    dataset_file("dated_topology_vcv.tsv", "Variance-covariance matrix of the dated topology of Table S1",
                 "Shared root-to-tip path length between two species divided by the total tree height.", vcv)
    go_md = (SUPP / "TableS2b_go_terms.md").read_text()
    t, nt, tb = md_section(go_md, "S10b")
    dataset_file("go_term_sets.tsv", t.rstrip("."), nt, md_to_df(tb))
    note_tables = (SUPP / "Note_S1_tables.md").read_text()
    for tid, fn in (("S11", "rarefaction.tsv"), ("S15", "simulation_recovery_summary.tsv"),
                    ("S16", "estimator_ranks.tsv")):
        t, nt, tb = md_section(note_tables, tid)
        dataset_file(fn, t.rstrip("."), nt, md_to_df(tb))
    # results cited in Note S1 that no table prints
    boot = []
    for ctx, f in (("Whole genome", "results_bootstrap/bootstrap_mci_values_whole_genome.tsv"),
                   ("Primary metabolism", "results_bootstrap_primary_metabolism/bootstrap_mci_values_primary_metabolism_ogs.tsv"),
                   ("Carbon metabolism", "results_bootstrap_C_metabolism/bootstrap_mci_values_C_metabolism_ogs.tsv"),
                   ("Nitrogen metabolism", "results_bootstrap_N_metabolism/bootstrap_mci_values_N_metabolism_ogs.tsv")):
        b = pd.read_csv(ROOT / "data" / f, sep="\t", index_col=0)
        top2 = b.apply(lambda r: set(r.nlargest(2).index) == {"Mpolymorpha", "Lcruciata"}, axis=1)
        boot.append({"Context": ctx, "Replicates": len(b), "Both liverworts first and second": round(top2.mean(), 3)})
    dataset_file("bootstrap_liverworts_top_two.tsv", "Fraction of bootstrap replicates with both liverworts first and second",
                 "Bootstrap over orthogroups, KEGG contexts.", pd.DataFrame(boot))
    dataset_file("rarefaction_symmetric.tsv", "Rarefaction below the smallest proteome",
                 "Fraction of replicates with both liverworts first and second and their mean ranks, at depths where all "
                 "eleven proteomes are subsampled.", pd.read_csv(CTRL / "J_symmetric_rarefaction_summary.tsv", sep="\t"))
    dataset_file("sign_rule.tsv", "Orientation of the axis under the two sign rules",
                 "Spearman correlation between the orderings given by the rule of the main text (positive correlation "
                 "with RCA diversity) and by the largest-component rule, for each context and RCA threshold.",
                 pd.read_csv(CTRL / "J_sign_rule_summary.tsv", sep="\t"))
    dataset_file("simulation_grid_40_species.tsv", "Simulation grid with 40 species",
                 "As simulation_grid.tsv, with 40 species.",
                 pd.read_csv(ROOT / "data/results_simulation/S1_grid_summary_S40.tsv", sep="\t"))
    ov_rows = []
    for ctx in ("primary_metabolism", "C_metabolism", "N_metabolism", "regulatory_networks"):
        kegg = set((ROOT / "data/filters" / f"{ctx}_ogs.txt").read_text().split())
        go = set((ROOT / "Matias/alexander_metodologia/GO_files" / f"{ctx}_ogs.txt").read_text().split())
        ov_rows.append({"Context": CONTEXT[ctx], "KEGG orthogroups": len(kegg), "Gene Ontology orthogroups": len(go),
                        "Shared": len(kegg & go), "Jaccard": round(len(kegg & go) / len(kegg | go), 3),
                        "Share of KEGG in Gene Ontology": round(len(kegg & go) / len(kegg), 3)})
    dataset_file("go_kegg_overlap.tsv", "Overlap between the KEGG and Gene Ontology orthogroup sets of each context",
                 "Orthogroups assigned to each context by the two annotations, before the RCA filter.", pd.DataFrame(ov_rows))
    (DATA / "README.md").write_text(
        "# Dataset S1\n\nResult tables that support the manuscript and Supplementary Note S1 but are not printed in "
        "the supplementary file. All files are tab-separated and are regenerated by `22_build_supplementary.py`.\n\n"
        + "\n".join(sorted(readme)) + "\n")
    print("Dataset S1:", len(readme), "files in", DATA)

    head, *secs = "\n".join(doc).split("\n## ")
    def order(sec):
        m = re.match(r"Table S(\d+)([a-z]?)\.", sec)
        return (int(m.group(1)), m.group(2)) if m else (99, sec)
    (SUPP / "supplementary_tables.md").write_text(head + "".join("\n## " + x for x in sorted(secs, key=order)))
    print("written", SUPP / "supplementary_tables.md")
    for name, df in [("S9", s1), ("S2", s3), ("S13a", s4a), ("S13b", s4b), ("S13c", s4c),
                     ("S12a", s5a), ("S12b", s5b), ("S12c", s5c), ("S14a", s6a), ("S14b", s6b),
                     ("S7a", s11a), ("S7b", s11b), ("S6a", s12a), ("S6b", s12b),
                     ("S3", s14), ("S5", s15), ("S5b", hl), ("S5c", st), ("S5d", s5d), ("S4", s16), ("S4b_rna", s16b), ("S4b", dyw)]:
        df.to_csv(SUPP / f"Table_{name}.tsv", sep="\t", index=False)
    print("rows:", {k: len(v) for k, v in [("S9", s1), ("S10", s2), ("S2", s3)]})


if __name__ == "__main__":
    main()
