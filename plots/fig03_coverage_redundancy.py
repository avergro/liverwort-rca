#!/usr/bin/env python3
"""
Fig. 3 — a full repertoire held with few genes.

Story: land plants occupy about the same share of orthogroups, so what
separates them is the number of genes behind each occupied orthogroup. The
liverworts sit at full coverage with the fewest genes per orthogroup, the alga
at low coverage, the seed plants at full coverage with many genes. The RCA
ranking follows that second axis, and the liverwort ratios hold in every
functional context. Liverworts carry the finding and are the coloured marks.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
CTRL = ROOT / "data" / "results_controls"
OUTDIR = ROOT / "plots" / "figures"
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]
SHORT = {"Mpolymorpha": "Marchantia", "Lcruciata": "Lunularia",
         "Smoellendorffii": "Selaginella", "Ppatens": "Physcomitrium",
         "Afiliculoides": "Azolla", "Scucullata": "Salvinia",
         "Creinhardtii": "Chlamydomonas", "Athaliana": "Arabidopsis",
         "Slycopersicum": "Solanum", "Pamericana": "Persea", "Pabies": "Picea"}
# label positions in data coordinates (x, y, ha), set by hand so that every name
# clears the marks and the other names; a thin leader joins name and point
LABEL_AT = {"Mpolymorpha": (1.50, 74.4, "right"), "Lcruciata": (1.72, 70.0, "right"),
            "Creinhardtii": (1.50, 51.0, "left"), "Scucullata": (2.05, 65.6, "right"),
            "Afiliculoides": (2.38, 65.6, "left"), "Smoellendorffii": (2.62, 61.3, "center"),
            "Ppatens": (2.90, 66.9, "left"), "Athaliana": (2.40, 73.4, "center"),
            "Slycopersicum": (3.30, 74.7, "center"), "Pamericana": (3.70, 72.4, "left"),
            "Pabies": (4.07, 70.1, "center")}


def colour(sp):
    return fs.BLUE if sp in LIVERWORTS else fs.GROUND


def main():
    cov = pd.read_csv(CTRL / "L_coverage_redundancy.tsv", sep="\t")
    est = pd.read_csv(CTRL / "J_estimators_values.tsv", sep="\t")
    rvr = pd.read_csv(CTRL / "L_rank_vs_redundancy.tsv", sep="\t")
    wg = cov[cov.context == "whole_genome"].set_index("species")
    rank = est[(est.context == "whole_genome") & (est.estimator == "MCi_RCA")].set_index("species")["rank"]

    fig = plt.figure(figsize=(7.2, 5.6))
    top = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.35,
                           left=0.09, right=0.97, top=0.93, bottom=0.57)
    bottom = fig.add_gridspec(1, 1, left=0.24, right=0.66, top=0.43, bottom=0.08)
    ax_a, ax_b = fig.add_subplot(top[0, 0]), fig.add_subplot(top[0, 1])
    ax_c = fig.add_subplot(bottom[0, 0])

    # (a) coverage against genes per occupied orthogroup
    for sp, r in wg.iterrows():
        ax_a.scatter(r.genes_per_occupied_OG, 100 * r.fraction_occupied, s=34,
                     color=colour(sp), edgecolors="none", zorder=3)
        lx, ly, ha = LABEL_AT[sp]
        ax_a.annotate(SHORT[sp], (r.genes_per_occupied_OG, 100 * r.fraction_occupied),
                      xytext=(lx, ly), textcoords="data", fontsize=fs.LABEL_PT - 1,
                      color=fs.BLUE if sp in LIVERWORTS else fs.INK,
                      ha=ha, va="center", fontstyle="italic",
                      arrowprops=dict(arrowstyle="-", color=fs.INK, lw=0.6, alpha=0.6,
                                      shrinkA=1.5, shrinkB=3.5))
    ax_a.set_xlim(0.45, 4.75)
    ax_a.set_ylim(48, 76)
    ax_a.set_xlabel("genes per occupied orthogroup")
    ax_a.set_ylabel("orthogroups occupied (%)")
    fs.bare(ax_a, keep=("left", "bottom"))
    fs.panel_label(ax_a, "(a)", x=-0.15)

    # (b) RCA rank against genes per occupied orthogroup
    for sp in wg.index:
        x = wg.loc[sp, "genes_per_occupied_OG"]
        ax_b.scatter(x, rank[sp], s=30, color=colour(sp), edgecolors="none", zorder=3)
        # one species per rank, so each name sits on its own row beside its point,
        # to the right when there is room before the axis end and to the left otherwise
        right = x < 2.4
        ax_b.text(x + (0.14 if right else -0.14), rank[sp], SHORT[sp], ha="left" if right else "right",
                  va="center", fontsize=fs.LABEL_PT - 1, fontstyle="italic",
                  color=fs.BLUE if sp in LIVERWORTS else fs.INK)
    # correlation of the plotted rank (rank 1 = highest score), so its sign matches the panel
    rho = pd.Series(rank).corr(wg.genes_per_occupied_OG, method="spearman")
    ax_b.text(1.0, 1.03, f"Spearman ρ (rank) = {rho:.2f}".replace("-", "\u2212"), transform=ax_b.transAxes, ha="right", va="bottom",
              fontsize=fs.LABEL_PT, color=fs.INK)
    ax_b.set_ylim(11.6, 0.4)
    ax_b.set_yticks([1, 6, 11])
    ax_b.set_xlim(1.0, 4.2)
    ax_b.set_xticks([1, 2, 3, 4])
    ax_b.set_xlabel("genes per occupied orthogroup")
    ax_b.set_ylabel("RCA rank, whole genome")
    fs.bare(ax_b, keep=("left", "bottom"))
    fs.panel_label(ax_b, "(b)", x=-0.2)

    # (c) liverwort / angiosperm ratios of coverage and redundancy by context
    rows = [("whole genome", "whole_genome", None), ("primary", "KEGG primary", "GO primary"),
            ("carbon", "KEGG carbon", "GO carbon"), ("nitrogen", "KEGG nitrogen", "GO nitrogen")]

    def ratio(ctx, col):
        c = cov[cov.context == ctx].set_index("species")
        return c.loc[LIVERWORTS, col].mean() / c.loc[ANGIOSPERMS, col].mean()

    for i, (lab, kctx, gctx) in enumerate(rows[::-1]):
        for j, (ctx, mark) in enumerate([(kctx, "D" if gctx is None else "o"), (gctx, "s")]):
            if ctx is None:
                continue
            dy = 0 if gctx is None else (j - 0.5) * 0.24
            ax_c.scatter(ratio(ctx, "occupied"), i + dy, s=26, marker=mark,
                         facecolors="none", edgecolors=fs.BLUE, linewidths=1.3, zorder=3)
            ax_c.scatter(ratio(ctx, "genes_per_occupied_OG"), i + dy, s=26, marker=mark,
                         color=fs.BLUE, edgecolors="none", zorder=3)
        ax_c.text(0.32, i, lab, ha="right", va="center", fontsize=fs.LABEL_PT, color=fs.INK)
    ax_c.axvline(1.0, color=fs.INK, lw=1.1, ls=(0, (4, 3)))
    # legend to the right of the panel, outside the data area
    from matplotlib.lines import Line2D
    key = dict(linestyle="none", markersize=5)
    handles = [Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor=fs.BLUE, label="orthogroups\noccupied", **key),
               Line2D([], [], marker="o", color=fs.BLUE, markeredgecolor="none", label="genes per occupied\northogroup", **key),
               Line2D([], [], marker="D", color=fs.INK, markeredgecolor="none", label="whole genome",
                      linestyle="none", markersize=4.5),
               Line2D([], [], marker="o", color=fs.INK, markeredgecolor="none", label="KEGG", **key),
               Line2D([], [], marker="s", color=fs.INK, markeredgecolor="none", label="GO", **key)]
    # framed, so the key glyphs are not read as data points
    leg = ax_c.legend(handles=handles, loc="center left", bbox_to_anchor=(1.06, 0.5), frameon=True,
                      fontsize=fs.LABEL_PT - 1, labelcolor=fs.INK, handletextpad=0.3, labelspacing=0.9,
                      borderaxespad=0, borderpad=0.6, fancybox=False)
    leg.get_frame().set_edgecolor(fs.INK)
    leg.get_frame().set_linewidth(0.6)
    leg.get_frame().set_facecolor("white")
    ax_c.set_xlim(0.34, 1.08)
    ax_c.set_ylim(-0.6, len(rows) - 0.4)
    ax_c.set_xticks([0.5, 0.75, 1.0])
    ax_c.set_yticks([])
    ax_c.set_xlabel("liverworts / angiosperms")
    fs.bare(ax_c, keep=("bottom",))
    fs.panel_label_aligned(ax_c, "(c)", ax_a, -0.15)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_3_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_3_preliminary.pdf", f"(rho {rho:.2f})")


if __name__ == "__main__":
    main()
