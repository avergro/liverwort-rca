#!/usr/bin/env python3
"""
Fig. 5 — metabolism held, regulation thinned.

Story: liverworts keep most metabolic orthogroups with about half the copies,
while regulatory orthogroups (signalling pathways and transcription factors)
are contracted about twice as strongly. Across species, regulatory families
carry as many copies per orthogroup as metabolic ones in liverworts and twice
as many in seed plants.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
RES = ROOT / "data" / "results_regulatory_contrast"
OUTDIR = ROOT / "plots" / "figures"
ORDER = ["Creinhardtii", "Mpolymorpha", "Lcruciata", "Ppatens", "Smoellendorffii",
         "Scucullata", "Afiliculoides", "Pabies", "Athaliana", "Slycopersicum", "Pamericana"]
SHORT = {"Mpolymorpha": "Marchantia", "Lcruciata": "Lunularia", "Smoellendorffii": "Selaginella",
         "Ppatens": "Physcomitrium", "Afiliculoides": "Azolla", "Scucullata": "Salvinia",
         "Creinhardtii": "Chlamydomonas", "Athaliana": "Arabidopsis", "Slycopersicum": "Solanum",
         "Pamericana": "Persea", "Pabies": "Picea"}
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]


def main():
    rat = pd.read_csv(RES / "liverwort_angiosperm_ratios.tsv", sep="\t").set_index("category")
    ps = pd.read_csv(RES / "per_species.tsv", sep="\t")
    g = ps.pivot(index="species", columns="category", values="genes_per_occupied_OG")

    # the two panels stacked, so that each spans the full width
    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[14, 11], hspace=0.42,
                          left=0.24, right=0.96, top=0.94, bottom=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    # (b) one column per measure, species in shared rows, so no marker covers another
    gb = gs[1, 0].subgridspec(1, 3, wspace=0.12)
    axes_b = [fig.add_subplot(gb[0, k]) for k in range(3)]

    # (a) liverwort / angiosperm ratios by category
    rows = [("primary (KEGG)", "KEGG primary", False), ("carbon (KEGG)", "KEGG carbon", False),
            ("nitrogen (KEGG)", "KEGG nitrogen", False), ("primary (GO)", "GO primary", False),
            ("carbon (GO)", "GO carbon", False), ("nitrogen (GO)", "GO nitrogen", False),
            ("regulatory (GO)", "GO regulatory", True), ("regulatory (KEGG)", "KEGG regulatory", True),
            ("  hormone signal transduction", "KEGG hormone signal transduction", True),
            ("  MAPK signalling", "KEGG MAPK signalling", True),
            ("  plant-pathogen interaction", "KEGG plant-pathogen interaction", True),
            ("  circadian rhythm", "KEGG circadian rhythm", True),
            ("  phosphatidylinositol signalling", "KEGG phosphatidylinositol signalling", True),
            ("transcription factors", "transcription factors", True)]
    for i, (lab, cat, reg) in enumerate(rows[::-1]):
        c = fs.RED if reg else fs.BLUE
        ax_a.scatter(rat.loc[cat, "L/A coverage"], i, s=24, facecolors="none", edgecolors=c, linewidths=1.3)
        ax_a.scatter(rat.loc[cat, "L/A genes_per_OG"], i, s=24, color=c, edgecolors="none")
        ax_a.plot([rat.loc[cat, "L/A genes_per_OG"], rat.loc[cat, "L/A coverage"]], [i, i],
                  color=c, lw=0.85, alpha=0.5)
        ax_a.text(0.125, i, lab, ha="right", va="center", fontsize=fs.LABEL_PT - 0.5,
                  color=c)
    ax_a.axvline(1.0, color=fs.INK, lw=1.1, ls=(0, (4, 3)))
    fs.shape_legend(ax_a, [("o", False, "orthogroups occupied"), ("o", True, "genes per occupied orthogroup")],
                    ncol=2, anchor=(0.0, 1.02))
    ax_a.set_xlim(0.14, 1.06)
    ax_a.set_ylim(-0.7, len(rows) - 0.3)
    ax_a.set_xticks([0.25, 0.5, 0.75, 1.0])
    ax_a.set_yticks([])
    ax_a.set_xlabel("liverworts / angiosperms")
    fs.bare(ax_a, keep=("bottom",))
    fs.panel_label(ax_a, "(a)", x=-0.30, y=1.06)

    # (b) regulatory over metabolic genes per orthogroup, per species
    series = [("KEGG", g["KEGG regulatory"] / g["KEGG primary"], "o"),
              ("GO", g["GO regulatory"] / g["GO primary"], "s"),
              ("transcription factors", g["transcription factors"] / g["KEGG primary"], "^")]
    for ax, (lab, v, m) in zip(axes_b, series):
        for i, sp in enumerate(ORDER[::-1]):
            ax.axhline(i, color=fs.GROUND, lw=0.6, alpha=0.5, zorder=1)
            ax.scatter(v[sp], i, s=22, marker=m, zorder=3, edgecolors="none",
                       color=fs.BLUE if sp in LIVERWORTS else fs.GROUND)
        ax.axvline(1.0, color=fs.INK, lw=1.1, ls=(0, (4, 3)))
        ax.set_title(lab, fontsize=fs.LABEL_PT, color=fs.INK, pad=4)
        ax.set_xlim(0.7, 4.3)
        ax.set_ylim(-0.7, len(ORDER) - 0.3)
        ax.set_xticks([1, 2, 3, 4])
        ax.set_yticks([])
        fs.bare(ax, keep=("bottom",))
    for i, sp in enumerate(ORDER[::-1]):
        axes_b[0].text(0.55, i, SHORT[sp], ha="right", va="center", fontsize=fs.LABEL_PT - 0.5,
                       color=fs.BLUE if sp in LIVERWORTS else fs.INK, fontstyle="italic",
                       transform=axes_b[0].transData, clip_on=False)
    axes_b[1].set_xlabel("genes per occupied orthogroup, regulatory / primary metabolism")
    fs.panel_label(axes_b[0], "(b)", x=-0.95, y=1.10)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_5_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_5_preliminary.pdf")


if __name__ == "__main__":
    main()
