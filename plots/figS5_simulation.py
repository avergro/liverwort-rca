#!/usr/bin/env python3
"""
Fig. S5 — simulation under the capability model.

Story: two simple statistics recover the planted capability endowment better
than the eigenvector index, and the eigenvector loses most of its accuracy when
capabilities are inherited along a tree. The three estimators that carry that
comparison are coloured; the rest stay ground.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
SIM = ROOT / "data" / "results_simulation"
OUTDIR = ROOT / "plots" / "figures"

STORY = {"MCi_rca": (fs.PURPLE, "RCA score"),
         "diversity": (fs.GREEN, "RCA diversity"),
         "fitness": (fs.TEAL, "Fitness-Complexity")}
GROUND_EST = {"MCi_raw": "presence score", "MCi_rarefied": "RCA score, rarefied",
              "CA_axis1": "correspondence axis 1"}
SCENARIOS = [("baseline", "independent"),
             ("phylogeny", "inherited"),
             ("incomplete", "incomplete"),
             ("neofunction", "neo-\nfunctionalization"),
             ("all_three", "all three")]


def main():
    g11 = pd.read_csv(SIM / "S1_grid_summary.tsv", sep="\t")
    g40 = pd.read_csv(SIM / "S1_grid_summary_S40.tsv", sep="\t")
    ext = pd.read_csv(SIM / "S2_extensions_summary.tsv", sep="\t")

    fig = plt.figure(figsize=(7.6, 3.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.28,
                          left=0.08, right=0.99, top=0.86, bottom=0.20)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # ── (a) recovery against duplication heterogeneity ───────────────────
    cell = 0.2
    for est, (col, label) in list(STORY.items()) + \
            [(k, (fs.GROUND, v)) for k, v in GROUND_EST.items()]:
        for df, ls, lw in ((g11, "-", 1.6), (g40, (0, (2.5, 1.5)), 1.1)):
            sub = (df[(df.cap_spread == cell) & (df.estimator == est)]
                     .sort_values("dup_spread"))
            if sub.empty:
                continue
            ax_a.plot(sub.dup_spread, sub.partial_signal_given_size,
                      color=col, ls=ls, lw=lw,
                      alpha=1.0 if est in STORY else 0.75, zorder=3 if est in STORY else 2)
        last = (df[(df.cap_spread == cell) & (df.estimator == est)]
                  .sort_values("dup_spread"))
        if est in STORY and not last.empty:
            sub11 = (g11[(g11.cap_spread == cell) & (g11.estimator == est)]
                       .sort_values("dup_spread"))
            nudge = {"diversity": 0.045, "fitness": -0.045}.get(est, 0.0)
            ax_a.text(sub11.dup_spread.iloc[-1] + 0.06,
                      sub11.partial_signal_given_size.iloc[-1] + nudge, label,
                      color=col, fontsize=fs.LABEL_PT, va="center")
    # names of the two ground series that stay visible (the correspondence axis
    # lies under the RCA score), at the end of the 11-species line
    for est, dy in (("MCi_raw", 0.0), ("MCi_rarefied", 0.0)):
        sub11 = (g11[(g11.cap_spread == cell) & (g11.estimator == est)]
                   .sort_values("dup_spread"))
        ax_a.text(sub11.dup_spread.iloc[-1] + 0.06,
                  sub11.partial_signal_given_size.iloc[-1] + dy, GROUND_EST[est],
                  color=fs.INK, fontsize=fs.LABEL_PT, va="center")
    # line-style key above the panel, outside the data
    from matplotlib.lines import Line2D
    ax_a.legend(handles=[Line2D([], [], color=fs.INK, lw=1.4, label="11 species"),
                         Line2D([], [], color=fs.INK, lw=1.4, ls=(0, (2.5, 1.5)), label="40 species")],
                loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False,
                fontsize=fs.LABEL_PT, labelcolor=fs.INK, handlelength=2.2, borderaxespad=0)
    ax_a.set_xlim(-0.08, 3.0)
    ax_a.set_ylim(-0.48, 1.05)  # the rarefied RCA score with 40 species falls to -0.43
    ax_a.set_xticks([0, 0.5, 1.0, 2.0])
    ax_a.set_yticks([-0.25, 0, 0.5, 1.0])
    ax_a.set_xlabel("heterogeneity of duplication across species")
    ax_a.set_ylabel("recovery of capability endowment")
    fs.bare(ax_a, keep=("left", "bottom"))
    fs.panel_label(ax_a, "(a)", x=-0.12)

    # ── (b) the four extensions ──────────────────────────────────────────
    width = 0.26
    order = list(STORY.items())
    for j, (est, (col, label)) in enumerate(order):
        vals = [ext[(ext.scenario == sc) & (ext.estimator == est)]
                .partial_signal_given_size.iloc[0] for sc, _ in SCENARIOS]
        xs = [i + (j - 1) * width for i in range(len(SCENARIOS))]
        ax_b.bar(xs, vals, width=width * 0.92, color=col,
                 alpha=1.0 if est == "MCi_rca" else 0.9)
        if est == "MCi_rca":
            for x, v in zip(xs, vals):
                ax_b.text(x, v + 0.03, f"{v:.2f}", ha="center",
                          fontsize=fs.LABEL_PT, color=col)
    for j, (est, (col, label)) in enumerate(order):
        ax_b.text(0.02 + j * 0.34, 1.04, label, transform=ax_b.transAxes,
                  color=col, fontsize=fs.LABEL_PT, va="bottom")
    ax_b.set_xticks(range(len(SCENARIOS)))
    ax_b.set_xticklabels([lbl for _, lbl in SCENARIOS], fontsize=fs.LABEL_PT)
    ax_b.set_xlabel("capability model, with each assumption relaxed")
    ax_b.set_ylim(0, 1.05)
    ax_b.set_yticks([0, 0.5, 1.0])
    ax_b.set_ylabel("recovery of capability endowment")
    fs.bare(ax_b, keep=("left", "bottom"))
    fs.panel_label(ax_b, "(b)", x=-0.10)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_S5_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_S5_preliminary.pdf")


if __name__ == "__main__":
    main()
