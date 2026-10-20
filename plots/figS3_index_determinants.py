#!/usr/bin/env python3
"""
Fig. S3 — what the index measures.

Story: MCi-RCA tracks gene number, and once that is removed the excess belongs
to non-seed plants as a group rather than to the liverworts. The liverworts and
Selaginella carry the finding and are the only coloured marks.
"""

import importlib.util
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

STORY = {"Mpolymorpha": fs.PURPLE, "Lcruciata": fs.PURPLE,
         "Smoellendorffii": fs.TEAL}
SHORT = {"Mpolymorpha": "Marchantia", "Lcruciata": "Lunularia",
         "Smoellendorffii": "Selaginella", "Ppatens": "Physcomitrium",
         "Afiliculoides": "Azolla", "Scucullata": "Salvinia",
         "Creinhardtii": "Chlamydomonas", "Athaliana": "Arabidopsis",
         "Slycopersicum": "Solanum", "Pamericana": "Persea", "Pabies": "Picea"}
GROUPS = {"Smoellendorffii": "lycophyte", "Ppatens": "moss",
          "Afiliculoides": "ferns", "Scucullata": "ferns",
          "Mpolymorpha": "liverworts", "Lcruciata": "liverworts",
          "Pabies": "gymnosperm", "Creinhardtii": "alga",
          "Athaliana": "angiosperms", "Slycopersicum": "angiosperms",
          "Pamericana": "angiosperms"}


def load_counts():
    spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
    bm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bm)
    counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
    M = bm.filter_matrix(bm.build_rca_matrix(counts), 2)
    mci, _ = bm.compute_complexity(M)
    return counts.sum(axis=0), mci


def main():
    genes, mci = load_counts()
    fits = pd.read_csv(CTRL / "F_pgls_fits.tsv", sep="\t")
    res = pd.read_csv(CTRL / "F_pgls_residuals.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0],
                          hspace=0.70, wspace=0.28,
                          left=0.09, right=0.98, top=0.94, bottom=0.10)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # ── (a) index against gene number ────────────────────────────────────
    x = np.log(genes.values)
    y = mci.reindex(genes.index).values
    offsets = {"Mpolymorpha": (7, 4), "Lcruciata": (7, -10),
               "Creinhardtii": (0, -12), "Scucullata": (-7, -9),
               "Ppatens": (7, 2), "Afiliculoides": (7, -9),
               "Smoellendorffii": (7, 2), "Athaliana": (7, 2),
               "Slycopersicum": (7, 2), "Pamericana": (7, -9),
               "Pabies": (-7, 4)}
    for sp, xi, yi in zip(genes.index, x, y):
        col = STORY.get(sp, fs.GROUND)
        ax_a.scatter(xi, yi, s=46 if sp in STORY else 30, color=col,
                     alpha=0.85 if sp not in STORY else 1.0, zorder=3,
                     edgecolors="none")
        dx, dy = offsets.get(sp, (7, 2))
        ha = "center" if dx == 0 else ("right" if dx < 0 else "left")
        ax_a.annotate(SHORT[sp], (xi, yi), textcoords="offset points",
                      xytext=(dx, dy), fontsize=fs.LABEL_PT,
                      ha=ha, style="italic",
                      color=col if sp in STORY else fs.INK)
    row = fits[(fits.context == "whole_genome") & (fits.model == "log_genes")].iloc[0]
    xs = np.linspace(x.min() - 0.08, x.max() + 0.08, 50)
    intercept = y.mean() - row.beta_log_genes * x.mean()
    ax_a.plot(xs, intercept + row.beta_log_genes * xs, color=fs.INK,
              lw=1.25, zorder=2)
    ax_a.text(0.99, 0.97,
              (f"PGLS slope {row.beta_log_genes:.2f}, t = {row.t_log_genes:.1f}, "
               f"Pagel lambda = {row.lambda_ML:.0f}").replace("-", "\u2212"),
              transform=ax_a.transAxes, fontsize=fs.LABEL_PT, color=fs.INK,
              ha="right", va="top")
    ax_a.set_xlabel("genes assigned to orthogroups (natural log)")
    ax_a.set_ylabel("RCA score, whole genome")
    ax_a.set_ylim(min(y) - 0.12, max(y) + 0.1)
    fs.bare(ax_a, keep=("left", "bottom"))
    fs.panel_label(ax_a, "(a)", x=-0.07)

    # ── (b) variance explained ───────────────────────────────────────────
    # OLS R2 of the whole-genome score on log gene number (0.68), and on
    # log gene number + BUSCO completeness (0.81)
    parts = [("gene number", 0.68, fs.PURPLE),
             ("+ completeness", 0.13, fs.PURPLE_L),
             ("unexplained", 0.19, fs.GROUND)]
    left = 0.0
    for i, (label, val, col) in enumerate(parts):
        ax_b.barh(0, val, left=left, height=0.5, color=col,
                  alpha=0.85 if col == fs.GROUND else 1.0)
        ax_b.text(left + val / 2, -0.38, f"{val:.0%}", ha="center", va="top",
                  fontsize=fs.LABEL_PT, color=fs.INK)
        ax_b.text(left + val / 2, 0.30 + 0.30 * (i % 2), label, ha="center",
                  va="bottom", fontsize=fs.LABEL_PT, color=fs.INK)
        ax_b.plot([left + val / 2, left + val / 2],
                  [0.26, 0.28 + 0.30 * (i % 2)], color=fs.INK, lw=0.6)
        left += val
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(-0.9, 1.15)
    fs.bare(ax_b, keep=())
    fs.panel_label(ax_b, "(b)", x=-0.05, y=1.30)
    ax_b.set_title("variance of the RCA score explained", fontsize=fs.AXIS_PT,
                   color=fs.INK, loc="left", pad=16)

    # ── (c) residuals by group ───────────────────────────────────────────
    res = res[res.context == "whole_genome"].copy()
    res["group"] = res.species.map(GROUPS)
    means = res.groupby("group").residual.mean().sort_values(ascending=True)
    colors = [fs.PURPLE if g == "liverworts" else
              fs.TEAL if g == "lycophyte" else fs.GROUND for g in means.index]
    ax_c.barh(range(len(means)), means.values, color=colors, height=0.6,
              alpha=0.85)
    ax_c.axvline(0, color=fs.INK, lw=0.85)
    for i, (g, v) in enumerate(means.items()):
        if v >= 0:
            ax_c.text(v + 0.015, i, f"{v:+.2f}".replace("-", "\u2212"), va="center", ha="left",
                      fontsize=fs.LABEL_PT, color=fs.INK)
        else:
            ax_c.text(0.012, i, f"{v:+.2f}".replace("-", "\u2212"), va="center", ha="left",
                      fontsize=fs.LABEL_PT, color=fs.INK)
        ax_c.text(-0.38, i, g, va="center", ha="right", fontsize=fs.LABEL_PT,
                  color=fs.PURPLE if g == "liverworts" else fs.INK)
    ax_c.set_xlim(-0.95, 0.42)
    ax_c.set_ylim(-0.7, len(means) - 0.15)
    fs.bare(ax_c, keep=())
    fs.panel_label(ax_c, "(c)", x=-0.05, y=1.30)
    ax_c.set_title("residual after removing gene number",
                   fontsize=fs.AXIS_PT, color=fs.INK, loc="left", pad=16)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_S3_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_S3_preliminary.pdf")


if __name__ == "__main__":
    main()
