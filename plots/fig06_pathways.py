#!/usr/bin/env python3
"""
Fig. 6 — what the liverwort-enriched orthogroups encode, and who shares them.

Story: the Filter B orthogroups follow the pathway composition of each context
(no pathway is over-represented; 32_filterB_enrichment.py), and nearly all of
them are also enriched in the moss or the lycophyte. Panels (a, b) show the
observed Filter B orthogroups per pathway (bars, ground) against the number
expected from the share of the pathway in the context (open diamonds). KEGG
global maps (ath011xx, ath012xx) are omitted because they summarise the
specific pathways shown.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
FB = ROOT / "data" / "results_filterB" / "filterB_enrichment.tsv"
CTRL = ROOT / "data" / "results_controls"
OUTDIR = ROOT / "plots" / "figures"
TOP = 7


def pathways(ctx):
    d = pd.read_csv(FB, sep="\t")
    d = d[d.context == ctx].copy()
    d["n_OGs"] = d.n_filterB
    d["expected"] = d.n_context * d.filterB_size / d.context_size
    return d.sort_values("n_OGs", ascending=False, kind="stable").head(TOP)


def bars(ax, d, n_total, title):
    ys = list(range(len(d)))[::-1]
    ax.barh(ys, d.n_OGs, height=0.64, color=fs.GROUND)
    ax.plot(d.expected, ys, linestyle="none", marker="D", markersize=4.5,
            markerfacecolor="none", markeredgecolor=fs.BLUE, markeredgewidth=1.0, zorder=3)
    for y, (_, r) in zip(ys, d.iterrows()):
        ax.text(max(r.n_OGs + 0.6, r.expected + 1.1), y, str(r.n_OGs), va="center",
                fontsize=fs.LABEL_PT, color=fs.INK)
        ax.text(-0.6, y, r.pathway_name.replace(" / ", "/"), va="center", ha="right",
                fontsize=fs.LABEL_PT, color=fs.INK)
    ax.set_xlim(0, 24)
    ax.set_yticks([])
    ax.set_xticks([0, 10, 20])
    ax.set_title(f"{title} ({n_total} orthogroups)", loc="left", fontsize=fs.LABEL_PT + 1,
                 color=fs.INK)
    fs.bare(ax, keep=("bottom",))


def main():
    share = pd.read_csv(CTRL / "L_filterB_sharing.tsv", sep="\t").set_index("context")
    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.75], height_ratios=[1, 1],
                          wspace=0.45, hspace=0.55, left=0.36, right=0.98, top=0.86, bottom=0.17)
    ax_n = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_s = fig.add_subplot(gs[:, 1])

    bars(ax_n, pathways("N_metabolism"), share.loc["N_metabolism", "filterB"], "nitrogen metabolism")
    bars(ax_c, pathways("C_metabolism"), share.loc["C_metabolism", "filterB"], "carbon metabolism")
    ax_c.set_xlabel("liverwort-enriched orthogroups")
    # key outside the data, in the free corner below the pathway names of (b)
    leg = fs.shape_legend(ax_c, [("D", False, "expected from\nthe context")], loc="lower left",
                          anchor=(0.02, 0.02), color=fs.BLUE)
    leg.set_bbox_to_anchor((0.02, 0.02), transform=fig.transFigure)
    leg.set_frame_on(True)
    leg.get_frame().set_edgecolor(fs.GROUND)
    leg.get_frame().set_linewidth(0.6)
    # (a) and (c) share one height in the figure, with room above the letters
    top_y = lambda ax: (0.955 - ax.get_position().y0) / ax.get_position().height
    fs.panel_label(ax_n, "(a)", x=-1.15, y=top_y(ax_n))
    fs.panel_label(ax_c, "(b)", x=-1.15, y=1.12)

    # (c) share of Filter B also enriched in the moss, the lycophyte, or either
    ctxs = [("primary", "primary_metabolism"), ("carbon", "C_metabolism"),
            ("nitrogen", "N_metabolism")]
    series = [("Physcomitrium", "Ppatens", fs.GROUND),
              ("Selaginella", "Smoellendorffii", "#D9D5C7"),
              ("at least one of the two", "Ppatens_or_Smoellendorffii", fs.INK)]
    # baseline: share of the orthogroups of the context in which no angiosperm is enriched
    # (the condition Filter B imposes) that are enriched in the same species
    base = pd.read_csv(CTRL / "M_filterB_sharing_baseline.tsv", sep="\t")
    base = base[base.set == "no angiosperm enriched"].set_index("context")
    base_col = {"Ppatens": "Ppatens", "Smoellendorffii": "Smoellendorffii",
                "Ppatens_or_Smoellendorffii": "either"}
    h = 0.24
    for i, (lab, ctx) in enumerate(ctxs[::-1]):
        n = share.loc[ctx, "filterB"]
        for j, (name, col, c) in enumerate(series):
            v = 100 * share.loc[ctx, col] / n
            y = i + (1 - j) * h
            ax_s.barh(y, v, height=h * 0.9, color=c)
            b = 100 * base.loc[ctx, base_col[col]]
            ax_s.plot(b, y, linestyle="none", marker="D", markersize=4, markerfacecolor="white",
                      markeredgecolor=fs.BLUE, markeredgewidth=1.0, zorder=3)
            ax_s.text(max(v, b) + 4.5, y, f"{v:.0f}", va="center", fontsize=fs.LABEL_PT - 1, color=fs.INK)
        ax_s.text(-3, i, lab, ha="right", va="center", fontsize=fs.LABEL_PT, color=fs.INK)
    for j, (name, _, c) in enumerate(series):
        ky = len(ctxs) - 0.35 + 0.22 * (len(series) - 1 - j)
        ax_s.barh(ky, 6, left=0, height=0.13, color=c)
        ax_s.text(9, ky, name, fontsize=fs.LABEL_PT, va="center",
                  color=fs.INK, fontstyle="normal" if name.startswith("at least") else "italic")
    leg = fs.shape_legend(ax_s, [("D", False, "no angiosperm\nenriched")], loc="lower left",
                          anchor=(0.0, 1.0), color=fs.BLUE)
    ax_s.set_xlim(0, 112)
    ax_s.set_ylim(-0.5, len(ctxs) + 0.35)
    ax_s.set_xticks([0, 50, 100])
    ax_s.set_yticks([])
    ax_s.set_xlabel("liverwort-enriched orthogroups\nalso enriched (%)", fontsize=fs.LABEL_PT + 1)
    fs.bare(ax_s, keep=("bottom",))
    fs.panel_label(ax_s, "(c)", x=-0.2, y=top_y(ax_s))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_6_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_6_preliminary.pdf")


if __name__ == "__main__":
    main()
