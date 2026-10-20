#!/usr/bin/env python3
"""
Fig. S2 — the index against simpler alternatives.

Story: once the correspondence axis is set aside as equivalent by construction
(Figure S4), no remaining statistic reproduces the ranking, and the ordering
depends on the RCA threshold through the orientation of the eigenvector. The
whole-genome series and the threshold at which the convention flips carry the
finding.
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
CTX = ["whole_genome", "primary_metabolism", "C_metabolism", "N_metabolism",
       "regulatory_networks"]
CTX_LABEL = {"whole_genome": "whole genome", "primary_metabolism": "primary",
             "C_metabolism": "carbon", "N_metabolism": "nitrogen",
             "regulatory_networks": "regulatory"}
BASE_LABEL = {"CA_axis1": "correspondence axis 1", "genes": "gene number",
              "diversity": "RCA diversity", "fitness": "Fitness-Complexity",
              "MCi_raw": "presence score", "low_ubiquity_share": "share of rare orthogroups"}


def main():
    base = pd.read_csv(CTRL / "G_baselines.tsv", sep="\t")
    sweep = pd.read_csv(CTRL / "G_threshold_sweep.tsv", sep="\t")
    vals = pd.read_csv(CTRL / "G_baseline_values_whole_genome.tsv",
                       sep="\t", index_col=0)

    sign = pd.read_csv(CTRL / "J_sign_rule_summary.tsv", sep="\t")

    fig = plt.figure(figsize=(6.4, 3.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.40,
                          left=0.25, right=0.98, top=0.86, bottom=0.21)
    ax_b, ax_c = (fig.add_subplot(gs[0, i]) for i in range(2))

    # ── (b) agreement of every baseline with the index ───────────────────
    order = ["genes", "diversity", "fitness", "low_ubiquity_share", "MCi_raw"]
    for i, bl in enumerate(order[::-1]):
        col = fs.GROUND
        for j, ctx in enumerate(CTX):
            row = base[(base.baseline == bl) & (base.context == ctx)]
            if row.empty:
                continue
            v = float(row.spearman_vs_MCi.iloc[0])
            # whole genome and regulatory context in colour, as in Fig. S1
            c = fs.PURPLE if ctx == "whole_genome" else fs.RED if ctx == "regulatory_networks" else col
            ax_b.scatter(v, i + (j - 2) * 0.13, s=18, color=c,
                         alpha=1.0 if c != fs.GROUND else 0.8,
                         edgecolors="none", zorder=3 if c != fs.GROUND else 2)
        # names in the left margin, clear of the dots
        ax_b.text(-0.05, i, BASE_LABEL[bl], ha="right", va="center",
                  fontsize=fs.LABEL_PT, color=col if col != fs.GROUND else fs.INK)
    ax_b.set_xlim(-0.02, 1.06)  # from zero, so the regulatory Fitness-Complexity value (0.06) is shown
    ax_b.set_xticks([0, 0.5, 1.0])
    ax_b.set_ylim(-0.6, len(order) - 0.15)
    ax_b.set_yticks([])
    ax_b.set_xlabel("|Spearman ρ| against the RCA score\n(one dot per context)")
    fs.bare(ax_b, keep=("bottom",))
    fs.panel_label(ax_b, "(a)", x=-0.42)

    # ── (c) RCA threshold sweep ──────────────────────────────────────────
    lw = sweep[sweep.species.isin(LIVERWORTS)]
    for ctx in CTX:
        col = fs.PURPLE if ctx == "whole_genome" else fs.GROUND
        worst = (lw[lw.context == ctx].groupby("threshold")["rank"].max()
                   .sort_index())
        ax_c.plot(worst.index, worst.values, color=col,
                  lw=1.6 if col != fs.GROUND else 1.0,
                  alpha=1.0 if col != fs.GROUND else 0.8,
                  zorder=3 if col != fs.GROUND else 2)
        if ctx == "whole_genome":
            ax_c.text(2.05, worst.values[-1] - 0.35, CTX_LABEL[ctx],
                      ha="right", va="bottom", fontsize=fs.LABEL_PT, color=col)
        if ctx == "primary_metabolism":
            ax_c.text(1.18, 8.0, "functional contexts",
                      ha="center", va="center", fontsize=fs.LABEL_PT, color=fs.INK)
    flip = sign[(sign.context == "whole_genome") & (~sign.sign_flip)].threshold
    if len(flip):
        ax_c.axvline(float(flip.max()), color=fs.INK, lw=1.3, ls=(0, (2, 2)),
                     zorder=1)
        ax_c.text(float(flip.max()) + 0.06, 9.9, "regulatory",
                  ha="left", va="center", fontsize=fs.LABEL_PT, color=fs.INK)
    ax_c.invert_yaxis()
    ax_c.set_xticks([0.5, 1.0, 1.5, 2.0])
    ax_c.set_yticks([1, 4, 7, 10])
    ax_c.set_xlim(0.4, 2.12)
    ax_c.set_xlabel("RCA threshold")
    ax_c.set_ylabel("worst liverwort rank")
    fs.bare(ax_c, keep=("left", "bottom"))
    fs.panel_label(ax_c, "(b)", x=-0.34)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_S2_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_S2_preliminary.pdf")


if __name__ == "__main__":
    main()
