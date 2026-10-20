#!/usr/bin/env python3
"""
Fig. 2 — the ranking and its consistency between annotation systems.

Story: RCA normalisation inverts the gene-count ordering and separates non-seed
from seed plants, and the separation is the same under KEGG and under GO. The
two liverworts carry the finding and are the only coloured marks; the rest of
the species stay ground.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
OUTDIR = ROOT / "plots" / "figures"
LIVERWORTS = ["Lcruciata", "Mpolymorpha"]
# label rows at the RCA end of panel (a), in data units of the ranking score
LABEL_Y = {"Lcruciata": 0.40, "Mpolymorpha": 0.32, "Creinhardtii": 0.24, "Smoellendorffii": 0.16,
           "Afiliculoides": 0.08, "Scucullata": 0.00, "Ppatens": -0.08, "Pabies": -0.30,
           "Athaliana": -0.37, "Slycopersicum": -0.47, "Pamericana": -0.54}
SHORT = {"Mpolymorpha": "Marchantia", "Lcruciata": "Lunularia",
         "Smoellendorffii": "Selaginella", "Ppatens": "Physcomitrium",
         "Afiliculoides": "Azolla", "Scucullata": "Salvinia",
         "Creinhardtii": "Chlamydomonas", "Athaliana": "Arabidopsis",
         "Slycopersicum": "Solanum", "Pamericana": "Persea", "Pabies": "Picea"}
KEGG_DIRS = {"primary": "results_primary_metabolism_v3",
             "carbon": "results_C_metabolism_v3",
             "nitrogen": "results_N_metabolism_v3",
             "regulatory": "results_regulatory_v3"}
# RCA score in the Gene Ontology contexts, written by 31_go_rankings.py
GO = pd.read_csv(ROOT / "data/results_go/MCi_rca_GO.tsv", sep="\t", index_col=0)


def colour(sp):
    return fs.PURPLE if sp in LIVERWORTS else fs.GROUND


def main():
    whole = pd.read_csv(ROOT / "data/results_v3/species_complexity.tsv",
                        sep="\t", index_col=0)
    ci = pd.read_csv(ROOT / "data/results_bootstrap/bootstrap_ci_whole_genome.tsv",
                     sep="\t", index_col=0)
    kegg = pd.DataFrame({k: pd.read_csv(ROOT / "data" / d / "species_complexity.tsv",
                                        sep="\t", index_col=0)["MCi_rca"]
                         for k, d in KEGG_DIRS.items()})
    order = whole.MCi_rca.sort_values(ascending=False).index.tolist()

    fig = plt.figure(figsize=(7.4, 6.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], width_ratios=[1, 1],
                          hspace=0.42, wspace=0.34,
                          left=0.13, right=0.99, top=0.93, bottom=0.16)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    # ── (a) raw to RCA ───────────────────────────────────────────────────
    for sp in whole.index:
        col = colour(sp)
        ax_a.plot([0, 1], [whole.MCi_raw[sp], whole.MCi_rca[sp]], color=col,
                  lw=1.6 if sp in LIVERWORTS else 0.9,
                  alpha=1.0 if sp in LIVERWORTS else 0.8, zorder=3 if sp in LIVERWORTS else 2)
        ax_a.scatter([0, 1], [whole.MCi_raw[sp], whole.MCi_rca[sp]], s=16,
                     color=col, zorder=4, edgecolors="none")
        # every species named at the RCA end, spread over evenly spaced rows
        # with a thin leader, since several end points nearly coincide
        ax_a.plot([1.03, 1.12], [whole.MCi_rca[sp], LABEL_Y[sp]], color=fs.INK, lw=0.6,
                  alpha=0.6, zorder=1)
        ax_a.text(1.14, LABEL_Y[sp], SHORT[sp], fontsize=fs.LABEL_PT - 1, style="italic",
                  color=col if sp in LIVERWORTS else fs.INK, va="center", ha="left")
    ax_a.set_xlim(-0.12, 1.75)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["presence\nscore", "RCA\nscore"])
    ax_a.set_ylabel("score, whole genome")
    fs.bare(ax_a, keep=("left", "bottom"))
    fs.panel_label(ax_a, "(a)", x=-0.20)

    # ── (b) bootstrap intervals ──────────────────────────────────────────
    ci = ci.reindex(order)
    ys = np.arange(len(order))[::-1]
    for y, sp in zip(ys, order):
        col = colour(sp)
        ax_b.plot([ci.MCi_RCA_lower95[sp], ci.MCi_RCA_upper95[sp]], [y, y],
                  color=col, lw=2.0 if sp in LIVERWORTS else 1.2, alpha=0.9,
                  solid_capstyle="round")
        ax_b.scatter(ci.MCi_RCA_observed[sp], y, s=14, color=col, zorder=4,
                     edgecolors="none")
        ax_b.text(0.62, y, SHORT[sp], fontsize=fs.LABEL_PT, va="center", style="italic",
                  color=col if sp in LIVERWORTS else fs.INK)
    ax_b.axvline(0, color=fs.INK, lw=0.7)
    ax_b.set_xlim(-0.62, 1.15)
    ax_b.set_xticks([-0.5, 0, 0.5])
    ax_b.set_xlabel("RCA score with bootstrap interval")
    fs.bare(ax_b, keep=("bottom",))
    fs.panel_label(ax_b, "(b)", x=-0.06)

    # ── (c) contexts under both annotations ──────────────────────────────
    contexts = ["nitrogen", "carbon", "primary"]   # bottom to top, so primary is on top as in Figs. 3c, 4a and 5a
    SP = 1.7          # vertical offset between contexts, in score units
    xs = np.arange(len(order))
    for i, ctx in enumerate(contexts):
        for j, (src, mark, lbl) in enumerate([(kegg, "o", "KEGG"), (GO, "s", "GO")]):
            vals = src[ctx].reindex(order)
            ax_c.scatter(xs + (j - 0.5) * 0.26, vals + i * SP,
                         s=[26 if sp in LIVERWORTS else 15 for sp in order],
                         color=[colour(sp) for sp in order],
                         marker=mark, edgecolors="none",
                         alpha=[1.0 if sp in LIVERWORTS else 0.8 for sp in order])
        ax_c.axhline(i * SP, color=fs.INK, lw=0.6, alpha=0.5)
        # white box so the zero line does not run through the label
        # context name above its zero line, at the right edge, clear of the tick labels
        from matplotlib.transforms import blended_transform_factory
        ax_c.text(0.995, i * SP + 0.04, ctx, fontsize=fs.LABEL_PT, color=fs.INK, ha="right",
                  va="bottom", zorder=4, transform=blended_transform_factory(ax_c.transAxes, ax_c.transData))
    fs.shape_legend(ax_c, [("o", True, "KEGG"), ("s", True, "GO")], ncol=2,
                    loc="lower right", anchor=(1.0, 1.0))
    ax_c.set_xticks(xs)
    ax_c.set_xticklabels([SHORT[s] for s in order], rotation=40, ha="right")
    for tick, sp in zip(ax_c.get_xticklabels(), order):
        tick.set_color(fs.PURPLE if sp in LIVERWORTS else fs.INK)
        tick.set_fontstyle("italic")
    ax_c.set_xlim(-2.1, len(order) - 0.4)
    ax_c.set_ylim(-0.85, 2 * SP + 0.65)
    # each context has its own zero line and scale
    ticks = [i * SP + t for i in range(len(contexts)) for t in (-0.5, 0, 0.5)]
    ax_c.set_yticks(ticks)
    ax_c.set_yticklabels(["\u22120.5", "0", "0.5"] * len(contexts), fontsize=fs.LABEL_PT - 1)
    ax_c.set_ylabel("RCA score within each\nfunctional context", labelpad=2)
    fs.bare(ax_c, keep=("left", "bottom"))
    fs.panel_label_aligned(ax_c, "(c)", ax_a, -0.20)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_2_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_2_preliminary.pdf")


if __name__ == "__main__":
    main()
