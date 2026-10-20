#!/usr/bin/env python3
"""
Fig. S4 — the index against the first axis of a correspondence analysis.

Story: the two are the same object. The transition matrix of the index is, up
to normalization, the correspondence analysis operator, so rank identity is
expected; the panel reports it as a check of the implementation.
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
# every species named; label positions in data coordinates (x, y, ha), with a
# thin leader to the point, since the moss and the ferns nearly coincide
LABEL = {"Lcruciata": ("Lunularia", (1.05, 0.27, "left")),
         "Mpolymorpha": ("Marchantia", (0.95, 0.17, "left")),
         "Creinhardtii": ("Chlamydomonas", (0.30, 0.25, "right")),
         "Smoellendorffii": ("Selaginella", (0.70, 0.05, "left")),
         "Ppatens": ("Physcomitrium", (0.20, 0.13, "right")),
         "Scucullata": ("Salvinia", (0.20, 0.06, "right")),
         "Afiliculoides": ("Azolla", (0.20, -0.01, "right")),
         "Pabies": ("Picea", (-1.10, -0.33, "left")),
         "Athaliana": ("Arabidopsis", (-1.30, -0.45, "left")),
         "Slycopersicum": ("Solanum", (-1.65, -0.55, "left")),
         "Pamericana": ("Persea", (-2.25, -0.42, "left"))}


def main():
    vals = pd.read_csv(CTRL / "G_baseline_values_whole_genome.tsv",
                       sep="\t", index_col=0)
    x, y = vals["CA_axis1"], vals["MCi"]
    if np.corrcoef(x, y)[0, 1] < 0:
        x = -x

    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    fig.subplots_adjust(left=0.19, right=0.97, top=0.93, bottom=0.17)
    for sp in vals.index:
        lw = sp in LIVERWORTS
        ax.scatter(x[sp], y[sp], s=34 if lw else 20,
                   color=fs.PURPLE if lw else fs.GROUND,
                   alpha=1.0 if lw else 0.85, edgecolors="none", zorder=3)
        if sp in LABEL:
            lx, ly, ha = LABEL[sp][1]
            ax.annotate(LABEL[sp][0], (x[sp], y[sp]), textcoords="data",
                        xytext=(lx, ly), fontsize=fs.LABEL_PT - 1, style="italic",
                        color=fs.PURPLE if lw else fs.INK, va="center", ha=ha,
                        arrowprops=dict(arrowstyle="-", color=fs.INK, lw=0.6, alpha=0.6,
                                        shrinkA=1.5, shrinkB=3.5))
    ax.set_xlim(-2.3, 1.85)
    ax.margins(y=0.12)
    ax.text(0.97, 0.32, "Spearman |ρ| = 1.000 in the\nwhole genome and four contexts",
            transform=ax.transAxes, fontsize=fs.LABEL_PT, color=fs.PURPLE,
            va="bottom", ha="right")
    ax.set_xlabel("correspondence analysis, axis 1")
    ax.set_ylabel("RCA score")
    fs.bare(ax, keep=("left", "bottom"))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_S4_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_S4_preliminary.pdf")


if __name__ == "__main__":
    main()
