#!/usr/bin/env python3
"""
Table S1: the dated topology used in the phylogenetic regression.

Reads the tree from 14_pgls_mci.py so the table can never drift from the tree
that the analysis actually used, and writes the node table, the Newick string
and the resulting variance-covariance matrix.

Usage:
    caffeinate -i python3 18_make_table_s7.py
Outputs: data/results_controls/TableS7_tree_nodes.tsv, TableS7_vcv.tsv,
         Submission/supplementary/TableS7_tree.md
"""

import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data" / "results_controls"
SUPP = ROOT / "Submission" / "supplementary"

SHORT = {"Mpolymorpha": "Marchantia polymorpha", "Lcruciata": "Lunularia cruciata",
         "Smoellendorffii": "Selaginella moellendorffii", "Ppatens": "Physcomitrium patens",
         "Afiliculoides": "Azolla filiculoides", "Scucullata": "Salvinia cucullata",
         "Creinhardtii": "Chlamydomonas reinhardtii", "Athaliana": "Arabidopsis thaliana",
         "Slycopersicum": "Solanum lycopersicum", "Pamericana": "Persea americana",
         "Pabies": "Picea abies"}

# Clade names and the source of each age, keyed by the set of descendant tips.
SOURCES = {
    frozenset(SHORT): ("Viridiplantae (chlorophyte and streptophyte split)",
                       "approximate outgroup calibration; 750 and 1,500 Ma tested in Dataset S1"),
    frozenset(SHORT) - {"Creinhardtii"}: ("Embryophyta, crown",
                                          "Harris et al. 2022, 515–493 Ma"),
    frozenset({"Mpolymorpha", "Lcruciata", "Ppatens"}): ("Setaphyta, crown (liverworts and mosses)",
        "Harris et al. 2022, 479–450 Ma, interval given for the divergence of Setaphyta"),
    frozenset({"Mpolymorpha", "Lcruciata"}): ("*Marchantia*–*Lunularia* split",
        "between the estimate of Villarreal et al. 2016 for this split (228 Ma; 95% HPD 188.6–294.3) and the total-evidence ages of Flores et al. 2023, in which it falls between the crowns of Marchantiales excluding *Lunularia* (248.6 Ma) and of Marchantiidae (310.5 Ma)"),
    frozenset({"Smoellendorffii", "Scucullata", "Afiliculoides", "Pabies",
               "Pamericana", "Athaliana", "Slycopersicum"}): ("Tracheophyta, crown",
                                                              "Harris et al. 2022, 452–447 Ma"),
    frozenset({"Scucullata", "Afiliculoides", "Pabies", "Pamericana", "Athaliana",
               "Slycopersicum"}): ("Euphyllophyta, crown", "Harris et al. 2022, 432–414 Ma"),
    frozenset({"Scucullata", "Afiliculoides"}): ("Salviniaceae, crown (*Azolla* and *Salvinia*)",
        "dated consensus tree of the Fern Tree of Life (Nitta et al. 2022; data release 1.9.0 of the ftolr package), 136.7 Ma, between the minimum age of the earliest *Azolla* megaspores (89 Ma; Metzgar et al. 2007) and the stem age of Salviniaceae (153.4 Ma; Testo and Sundue 2016)"),
    frozenset({"Pabies", "Pamericana", "Athaliana", "Slycopersicum"}):
        ("Spermatophyta, crown", "Morris et al. 2018, 365.0–330.9 Ma"),
    frozenset({"Pamericana", "Athaliana", "Slycopersicum"}): ("Mesangiospermae, crown (*Persea* and eudicots)",
                                                              "Magallón et al. 2015, 137–135 Ma"),
    frozenset({"Athaliana", "Slycopersicum"}): ("Pentapetalae, crown (rosids and asterids)",
                                                "Magallón et al. 2015, 126.5–120.9 Ma"),
}


def load_tree():
    spec = importlib.util.spec_from_file_location("pg", ROOT / "14_pgls_mci.py")
    pg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pg)
    return pg


def collect(node, depth, height, rows):
    """Walk the tree and record each internal node as an age before present."""
    if node["children"] is None:
        return [node["name"]], depth + node["length"]
    tips, d = [], depth + node["length"]
    for ch in node["children"]:
        t, _ = collect(ch, d, height, rows)
        tips += t
    name, source = SOURCES.get(frozenset(tips), ("unnamed node", ""))
    rows.append({"node": name, "age_Ma": round(height - d, 1),
                 "n_tips": len(tips), "source_of_age": source,
                 "descendants": ", ".join(sorted(SHORT[t].split()[0] for t in tips))})
    return tips, d


BINOMIAL = {"Creinhardtii": "Chlamydomonas reinhardtii", "Mpolymorpha": "Marchantia polymorpha",
            "Lcruciata": "Lunularia cruciata", "Ppatens": "Physcomitrium patens",
            "Smoellendorffii": "Selaginella moellendorffii", "Scucullata": "Salvinia cucullata",
            "Afiliculoides": "Azolla filiculoides", "Pabies": "Picea abies",
            "Pamericana": "Persea americana", "Athaliana": "Arabidopsis thaliana",
            "Slycopersicum": "Solanum lycopersicum"}
# first letter of the genus and the first two of the epithet
CODE = {sp: b.split()[0][0] + b.split()[1][:2] for sp, b in BINOMIAL.items()}
assert len(set(CODE.values())) == len(CODE)


def main():
    pg = load_tree()
    tree = pg.parse_newick(pg.NEWICK)
    V = pg.vcv(pg.parse_newick(pg.NEWICK))
    height = 1000.0
    rows = []
    collect(tree, 0.0, height, rows)
    nodes = pd.DataFrame(rows).sort_values("age_Ma", ascending=False)
    OUT.mkdir(parents=True, exist_ok=True)
    SUPP.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(OUT / "TableS7_tree_nodes.tsv", sep="\t", index=False)
    V.round(4).to_csv(OUT / "TableS7_vcv.tsv", sep="\t")

    md = ["# Table S1. Dated topology used in the phylogenetic regression", "",
          "Ages are node ages in millions of years before present. Bryophyte monophyly is",
          "assumed. Each age is the midpoint of the interval (95% highest posterior density) given by",
          "the source of that node, rounded to whole millions of years, except the *Marchantia*–*Lunularia*",
          "split and the crown of Salviniaceae, whose ages are explained in the table, and the root is an approximate",
          "outgroup calibration. The Pagel lambda is estimated from the data, and setting every land",
          "plant node at the younger or at the older bound of its interval changed the regression",
          "little (Dataset S1), but the topology does determine which species are treated as",
          "non-independent. Sources cited only in this table: Flores et al. 2023, New Phytol 240:2137;",
          "Magallón et al. 2015, New Phytol 207:437; Metzgar et al. 2007, Int J Plant Sci 168:1045;",
          "Nitta et al. 2022, Front Plant Sci 13:909768; Testo and Sundue 2016, Mol Phylogenet Evol",
          "105:200; Villarreal et al. 2016, New Phytol 209:1734. Harris et al. 2022 and Morris et al. 2018",
          "are listed in the main text.",
          "", "## Nodes", "",
          "| Node | Age (Ma) | Tips | Source of age |", "|---|---|---|---|"]
    for _, r in nodes.iterrows():
        md.append(f"| {r.node} | {r.age_Ma:.0f} | {r.n_tips} | {r.source_of_age} |")
    md += ["", "## Newick string", "", "```", pg.NEWICK, "```", "",
           "The variance-covariance matrix used in the regression is given in Dataset S1.", ""]
    (SUPP / "TableS7_tree.md").write_text("\n".join(md) + "\n")

    print(nodes.to_string(index=False))
    print("\nwritten to", OUT, "and", SUPP / "TableS7_tree.md")


if __name__ == "__main__":
    main()
