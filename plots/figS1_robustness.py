#!/usr/bin/env python3
"""
Fig. S1 — robustness of the ranking.

Story: the whole-genome result survives every control, and the regulatory
context fails all of them. Those two rows carry the finding; the intermediate
contexts stay ground.
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
CTRL = ROOT / "data" / "results_controls"
OUTDIR = ROOT / "plots" / "figures"
LIVERWORTS = ["Lcruciata", "Mpolymorpha"]
SHORT = {"Mpolymorpha": "Marchantia", "Lcruciata": "Lunularia",
         "Smoellendorffii": "Selaginella", "Ppatens": "Physcomitrium",
         "Afiliculoides": "Azolla", "Scucullata": "Salvinia",
         "Creinhardtii": "Chlamydomonas", "Athaliana": "Arabidopsis",
         "Slycopersicum": "Solanum", "Pamericana": "Persea", "Pabies": "Picea"}
CTX = ["whole_genome", "primary_metabolism", "C_metabolism", "N_metabolism",
       "regulatory_networks"]
CTX_LABEL = {"whole_genome": "whole genome", "primary_metabolism": "primary",
             "C_metabolism": "carbon", "N_metabolism": "nitrogen",
             "regulatory_networks": "regulatory"}
# symmetric design: every proteome subsampled below the smallest one
SYM_DEPTHS = [6000, 8000, 10000, 12000]
# depths equal to the gene totals of the proteomes, up to no subsampling
DEPTHS = ["12702", "14160", "17455", "24542", "34211", "none"]
DEPTH_LABEL = [f"{int(x):,}" if x != "none" else "no\nsubsampling" for x in DEPTHS]


def ctx_colour(ctx, value_ok):
    if ctx == "whole_genome":
        return fs.PURPLE
    if ctx == "regulatory_networks":
        return fs.RED
    return fs.GROUND


def main():
    depth = pd.read_csv(CTRL / "E_depth_series_liverworts.tsv", sep="\t")
    jack = pd.read_csv(CTRL / "C_jackknife.tsv", sep="\t")
    null = pd.read_csv(CTRL / "D_null_model.tsv", sep="\t")

    fig = plt.figure(figsize=(7.4, 6.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], wspace=0.38,
                          hspace=0.52, left=0.17, right=0.98,
                          top=0.90, bottom=0.13)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # ── (a) rarefaction depth series ─────────────────────────────────────
    d = depth[depth.set == "all_species"]
    sym = pd.read_csv(CTRL / "J_symmetric_rarefaction_summary.tsv", sep="\t")
    n_sym = len(SYM_DEPTHS)
    for i, ctx in enumerate(CTX[::-1]):
        col = ctx_colour(ctx, True)
        for j, dep in enumerate([str(x) for x in SYM_DEPTHS] + DEPTHS):
            if j < n_sym:
                row = sym[(sym.context == ctx) & (sym.depth.astype(str) == dep)]
                key = "liverworts_top2"
            else:
                row = d[(d.context == ctx) & (d.depth.astype(str) == dep)]
                key = "both_liverworts_top2"
            if row.empty:
                continue
            v = float(row[key].iloc[0])
            ax_a.scatter(j, i, s=12 + 190 * v, color=col,
                         alpha=1.0 if col != fs.GROUND else 0.85,
                         edgecolors="none")
            ax_a.text(j, i - 0.36, f"{v:.2f}", ha="center", va="top",
                      fontsize=fs.LABEL_PT, color=fs.INK)
        ax_a.text(-0.85, i, CTX_LABEL[ctx], ha="right", va="center",
                  fontsize=fs.LABEL_PT,
                  color=col if col != fs.GROUND else fs.INK)
    n_all = n_sym + len(DEPTHS)
    ax_a.axvline(n_sym - 0.5, color=fs.GROUND, lw=0.85, ls=":")
    ax_a.set_xticks(range(n_all))
    ax_a.set_xticklabels([f"{x:,}" for x in SYM_DEPTHS] + DEPTH_LABEL, fontsize=fs.LABEL_PT - 0.5)
    ax_a.set_xlabel("genes retained per proteome (rarefaction depth)")
    ax_a.set_xlim(-0.9, n_all - 0.4)
    ax_a.set_ylim(-0.75, len(CTX) - 0.4)
    ax_a.set_yticks([])
    fs.bare(ax_a, keep=("bottom",))
    fs.panel_label(ax_a, "(a)", x=-0.12)
    ax_a.set_title("proportion of replicates with both liverworts first and second",
                   fontsize=fs.AXIS_PT, color=fs.INK, loc="left", pad=8)

    # ── (b) jackknife ────────────────────────────────────────────────────
    j_lw = jack[jack.species.isin(LIVERWORTS)]
    rng = np.random.default_rng(1)
    for i, ctx in enumerate(CTX[::-1]):
        col = ctx_colour(ctx, True)
        worst = (j_lw[(j_lw.context == ctx) & (~j_lw.dropped.isin(["none"] + LIVERWORTS))]
                 .groupby("dropped")["rank"].max())
        jitter = rng.uniform(-0.16, 0.16, len(worst))
        ax_b.scatter(worst.values, i + jitter, s=20, color=col,
                     alpha=0.9 if col != fs.GROUND else 0.8, edgecolors="none")
        ax_b.text(0.35, i, CTX_LABEL[ctx], ha="right", va="center",
                  fontsize=fs.LABEL_PT,
                  color=col if col != fs.GROUND else fs.INK)
    # limits follow the data, so no removal is clipped (carbon reaches rank 5 without the moss)
    ax_b.set_xticks([1, 2, 3, 4, 5])
    ax_b.set_xlim(0.3, max(5, int(j_lw[j_lw.species.isin(LIVERWORTS)]["rank"].max())) + 0.6)
    ax_b.set_ylim(-0.7, len(CTX) - 0.4)
    ax_b.set_yticks([])
    ax_b.set_xlabel("worst liverwort rank")
    fs.bare(ax_b, keep=("bottom",))
    fs.panel_label(ax_b, "(b)", x=-0.30)
    ax_b.set_title("one dot per species removed", fontsize=fs.AXIS_PT,
                   color=fs.INK, loc="left", pad=8)

    # ── (c) null model ───────────────────────────────────────────────────
    n = null[null.species != "Pabies"]
    for i, ctx in enumerate(CTX[::-1]):
        sub = n[n.context == ctx]
        for _, row in sub.iterrows():
            is_lw = row.species in LIVERWORTS
            ax_c.scatter(row.z_vs_null, i, s=26 if is_lw else 14,
                         color=fs.PURPLE if is_lw else fs.GROUND,
                         alpha=1.0 if is_lw else 0.8, edgecolors="none",
                         zorder=3 if is_lw else 2)
        ax_c.text(-22, i, CTX_LABEL[ctx], ha="right", va="center",
                  fontsize=fs.LABEL_PT,
                  color=fs.RED if ctx == "regulatory_networks" else
                  fs.PURPLE if ctx == "whole_genome" else fs.INK)
    ax_c.axvline(0, color=fs.INK, lw=0.7)
    ax_c.axvline(1.96, color=fs.INK, lw=0.7, ls=(0, (2, 2)))
    ax_c.text(2.3, -0.62, "z = 1.96", fontsize=fs.LABEL_PT, color=fs.INK)
    ax_c.set_xscale("symlog", linthresh=5)
    ax_c.set_xlim(-24, 40)
    ax_c.set_xticks([-20, -5, 0, 5, 20])
    ax_c.set_xticklabels(["−20", "−5", "0", "5", "20"])
    ax_c.xaxis.set_minor_locator(plt.NullLocator())
    ax_c.set_xlabel("z score against a\ndegree-preserving null model")
    ax_c.set_ylim(-0.7, len(CTX) - 0.4)
    ax_c.set_yticks([])
    fs.bare(ax_c, keep=("bottom",))
    fs.panel_label(ax_c, "(c)", x=-0.30)
    ax_c.set_title("liverworts in colour", fontsize=fs.AXIS_PT,
                   color=fs.INK, loc="left", pad=8)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_S1_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_S1_preliminary.pdf")


if __name__ == "__main__":
    main()
