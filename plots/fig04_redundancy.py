#!/usr/bin/env python3
"""
Fig. 4 — low redundancy: transcription factor copies and retained duplications.

Story: every ratio is below one. Liverworts carry a fraction of the angiosperm
transcription factor copies and a lower duplication density in every functional
context, so the breadth measured by the index is not built by duplication.
The ratio scale itself carries the finding, so the bars are coloured and the
reference line at a ratio of one is the only chrome kept.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
EXP = ROOT / "data" / "results_expansion_contraction"
TF = ROOT / "data" / "results_tf_redundancy"
OUTDIR = ROOT / "plots" / "figures"
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
ANGIOSPERMS = ["Athaliana", "Slycopersicum", "Pamericana"]
# Duplications retained in each species (both daughter clades carry its genes),
# per occupied orthogroup, from 27_regulatory_contrast.py. The dup_per_OG of
# 09c summed every duplication of the orthogroup whatever the species.
RATIOS = ROOT / "data" / "results_regulatory_contrast" / "liverwort_angiosperm_ratios.tsv"


def retained_ratios(annotation):
    r = pd.read_csv(RATIOS, sep="\t", index_col=0)["L/A retained_dup_per_OG"]
    return {ctx: r[f"{annotation} {ctx}"] for ctx in ("primary", "carbon", "nitrogen", "regulatory")}


def single_ratio(category):
    return pd.read_csv(RATIOS, sep="\t", index_col=0)["L/A retained_dup_per_OG"][category]


def main():
    tf = pd.read_csv(TF / "liverwort_vs_angiosperm_summary.tsv", sep="\t",
                     index_col=0).sort_values("ratio_L_A")
    par = pd.read_csv(TF / "tf_family_paralogs.tsv", sep="\t")
    # mean copies per species in each group, pooled over the twelve families
    pooled = ((par[par.species.isin(LIVERWORTS)].total_copies.sum() / len(LIVERWORTS)) /
              (par[par.species.isin(ANGIOSPERMS)].total_copies.sum() / len(ANGIOSPERMS)))
    kegg = retained_ratios("KEGG")
    go = retained_ratios("GO")

    fig = plt.figure(figsize=(7.2, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.35], wspace=0.45,
                          left=0.14, right=0.98, top=0.84, bottom=0.18)
    # duplications first, since the text cites them before the transcription factors
    ax_a = fig.add_subplot(gs[0, 1])
    ax_b = fig.add_subplot(gs[0, 0])

    # ── (b) transcription factor copies ──────────────────────────────────
    ys = range(len(tf))
    ax_a.barh(list(ys), tf.ratio_L_A.values, height=0.62, color=fs.PURPLE,
              alpha=0.92)
    for y, (fam, row) in zip(ys, tf.iterrows()):
        ax_a.text(row.ratio_L_A + 0.010, y, f"{row.ratio_L_A:.2f}", va="center",
                  fontsize=fs.LABEL_PT, color=fs.INK)
        ax_a.text(-0.012, y, fam.replace("_", "-"), va="center", ha="right",
                  fontsize=fs.LABEL_PT, color=fs.INK)
    ax_a.axvline(pooled, color=fs.INK, lw=1.3)
    ax_a.text(pooled + 0.006, len(tf) - 0.2,
              f"all twelve families pooled: {pooled:.2f}",
              fontsize=fs.LABEL_PT, color=fs.INK, va="center")
    ax_a.set_xlim(0, 0.42)  # room for the GRAS value and the pooled label at the right edge
    ax_a.set_ylim(-0.7, len(tf) + 0.1)
    ax_a.set_xticks([0, 0.1, 0.2, 0.3])
    ax_a.set_yticks([])
    ax_a.set_xlabel("transcription factor copies,\nliverworts / angiosperms")
    fs.bare(ax_a, keep=("bottom",))
    fs.panel_label(ax_a, "(b)", x=-0.14)

    # ── (a) duplication density ──────────────────────────────────────────
    # whole genome and transcription factor orthogroups have no annotation split
    contexts = ["whole genome", "primary", "carbon", "nitrogen", "regulatory",
                "transcription factors"]
    single = {"whole genome": (single_ratio("whole genome"), "D"),
              "transcription factors": (single_ratio("transcription factors"), "^")}
    for i, ctx in enumerate(contexts[::-1]):
        if ctx in single:
            v, mark = single[ctx]
            ax_b.scatter(v, i, s=30, marker=mark, color=fs.PURPLE, edgecolors="none")
        else:
            for j, (src, mark, lbl) in enumerate([(kegg, "o", "KEGG"),
                                                  (go, "s", "GO")]):
                ax_b.scatter(src[ctx], i + (j - 0.5) * 0.22, s=30, marker=mark,
                             color=fs.PURPLE, edgecolors="none")
        ax_b.text(-0.03, i, ctx.replace("transcription factors", "transcription\nfactors"), ha="right", va="center", fontsize=fs.LABEL_PT,
                  color=fs.INK)
    ax_b.axvline(1.0, color=fs.INK, lw=1.1, ls=(0, (4, 3)))
    fs.shape_legend(ax_b, [("o", True, "KEGG"), ("s", True, "GO"), ("D", True, "whole genome"),
                           ("^", True, "transcription factors")], ncol=2,
                    loc="lower left", anchor=(0.0, 1.0))
    ax_b.set_xlim(0.0, 1.12)
    ax_b.set_ylim(-0.75, len(contexts) - 0.1)
    ax_b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_b.set_yticks([])
    ax_b.set_xlabel("retained duplications per occupied\northogroup, liverworts / angiosperms")
    fs.bare(ax_b, keep=("bottom",))
    fs.panel_label(ax_b, "(a)", x=-0.28)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_4_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_4_preliminary.pdf",
          f"(pooled TF ratio {pooled:.3f}; KEGG ratios "
          + ", ".join(f"{k} {v:.2f}" for k, v in kegg.items()) + ")")


if __name__ == "__main__":
    main()
