#!/usr/bin/env python3
"""
Fig. 1 — the framework, the dataset and the workflow.

Story: the mapping from the capability model to genome gene content. Panel (a)
carries it, so the species that has every required capability and the
orthogroups it reaches are the only coloured marks; panel (b) marks the two
liverworts because the rest of the figure refers to them, and panel (c) is
plain workflow.
"""

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyArrow, Rectangle

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "plots"))
import figure_style as fs  # noqa: E402

fs.apply()
OUTDIR = ROOT / "plots" / "figures"
LIVERWORTS = ["Mpolymorpha", "Lcruciata"]
SHORT = {"Mpolymorpha": "Marchantia", "Lcruciata": "Lunularia",
         "Smoellendorffii": "Selaginella", "Ppatens": "Physcomitrium",
         "Afiliculoides": "Azolla", "Scucullata": "Salvinia",
         "Creinhardtii": "Chlamydomonas", "Athaliana": "Arabidopsis",
         "Slycopersicum": "Solanum", "Pamericana": "Persea", "Pabies": "Picea"}
GENES = {"Creinhardtii": 12702, "Mpolymorpha": 14160, "Scucullata": 17083,
         "Lcruciata": 17455, "Ppatens": 24542, "Athaliana": 25554,
         "Afiliculoides": 18231, "Slycopersicum": 29730,
         "Smoellendorffii": 23017, "Pamericana": 34211, "Pabies": 38650}


def gradient_link(ax, p0, p1, c0, c1, lw=1.4, alpha=0.8, zorder=1):
    t = np.linspace(0, 1, 64)
    xs = p0[0] + (p1[0] - p0[0]) * t
    ys = p0[1] + (p1[1] - p0[1]) * t + 0.10 * np.sin(np.pi * t) * np.sign(p1[0] - p0[0])
    pts = np.array([xs, ys]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    tt = np.linspace(0, 1, len(segs))[:, None]
    cols = (1 - tt) * np.array(to_rgb(c0)) + tt * np.array(to_rgb(c1))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw,
                                     alpha=alpha, capstyle="round", zorder=zorder))


def panel_a(ax):
    caps = ["c1", "c2", "c3"]
    species = [("species A", fs.PURPLE, {0, 1, 2}),
               ("species B", fs.GROUND, {0, 1}),
               ("species C", fs.GROUND, {2})]
    ogs = [("OG 1", {0, 1}), ("OG 2", {2}), ("OG 3", {0, 1, 2})]

    ys = [2.4, 1.5, 0.6]
    for (name, col, held), y in zip(species, ys):
        ax.text(-0.05, y, name, ha="right", va="center", fontsize=fs.LABEL_PT,
                color=col if col != fs.GROUND else fs.INK)
        for j, cap in enumerate(caps):
            filled = j in held
            ax.scatter(0.22 + 0.19 * j, y, s=34,
                       color=col if filled else "none",
                       edgecolors=col, linewidths=1.0, zorder=3,
                       alpha=1.0 if col != fs.GROUND else 0.85)
    for j, cap in enumerate(caps):
        ax.text(0.22 + 0.19 * j, 2.85, cap, ha="center", fontsize=fs.LABEL_PT,
                color=fs.INK)
    ax.text(0.41, 3.15, "capabilities", ha="center", fontsize=fs.LABEL_PT,
            color=fs.INK)

    og_x = 1.45
    og_y = [2.4, 1.5, 0.6]
    for (name, req), y in zip(ogs, og_y):
        ax.text(og_x + 0.08, y, f"{name}  requires " + ", ".join(caps[i] for i in sorted(req)),
                ha="left", va="center", fontsize=fs.LABEL_PT, color=fs.INK)
        ax.scatter(og_x, y, s=30, color=fs.INK, zorder=3, edgecolors="none")

    for (name, col, held), sy in zip(species, ys):
        for (og_name, req), oy in zip(ogs, og_y):
            if req <= held:
                gradient_link(ax, (0.74, sy), (og_x - 0.02, oy), col, fs.INK,
                              lw=1.6 if col != fs.GROUND else 1.0,
                              alpha=0.85 if col != fs.GROUND else 0.55)
    ax.text(0.0, -0.15,
            "a species occupies an orthogroup only when it has every capability the orthogroup requires",
            fontsize=fs.LABEL_PT, color=fs.INK, ha="left")
    ax.set_xlim(-0.75, 3.2)
    ax.set_ylim(-0.45, 3.4)
    ax.axis("off")


def tree_positions():
    spec = importlib.util.spec_from_file_location("pg", ROOT / "14_pgls_mci.py")
    pg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pg)
    tree = pg.parse_newick(pg.NEWICK)
    tips, edges = [], []

    def walk(node, x_root):
        if node["children"] is None:
            x = x_root + node["length"]
            tips.append((node["name"], x))
            return x, len(tips) - 1
        x = x_root + node["length"]
        child_pos = []
        for ch in node["children"]:
            cx, cy = walk(ch, x)
            child_pos.append((cx, cy, ch))
        y = np.mean([c[1] for c in child_pos])
        for cx, cy, ch in child_pos:
            edges.append(((x, cy), (cx, cy)))          # horizontal
        edges.append(((x, min(c[1] for c in child_pos)),
                      (x, max(c[1] for c in child_pos))))  # vertical
        return x, y

    walk(tree, 0.0)
    return tips, edges


def panel_b(ax):
    tips, edges = tree_positions()
    height = max(x for _, x in tips)
    tree_w = 105.0

    def px(depth):
        """root at the left, tips at tree_w"""
        return depth / height * tree_w

    for (x0, y0), (x1, y1) in edges:
        ax.plot([px(x0), px(x1)], [y0, y1], color=fs.INK, lw=1.0,
                solid_capstyle="round", zorder=1)

    names = [t[0] for t in tips]
    gmax = max(GENES.values())
    for name, depth in tips:
        y = names.index(name)
        col = fs.PURPLE if name in LIVERWORTS else fs.GROUND
        ax.text(tree_w + 6, y, SHORT[name], fontsize=fs.LABEL_PT, va="center",
                style="italic", color=col if name in LIVERWORTS else fs.INK)
        w = 150 * GENES[name] / gmax
        ax.add_patch(Rectangle((215, y - 0.28), w, 0.56, color=col,
                               alpha=0.9, linewidth=0))
        ax.text(215 + w + 5, y, f"{GENES[name]:,}", fontsize=fs.LABEL_PT,
                va="center", color=fs.INK)
    ax.text(215, len(tips) - 0.25, "genes assigned to orthogroups",
            fontsize=fs.LABEL_PT, color=fs.INK)

    for age in (1000, 500, 0):
        x = px(height - age)
        ax.plot([x, x], [-1.05, -0.85], color=fs.INK, lw=0.85)
        ax.text(x, -1.24, f"{age}", ha="center", va="top",
                fontsize=fs.LABEL_PT, color=fs.INK)
    ax.plot([px(0), px(height)], [-0.95, -0.95], color=fs.INK, lw=0.85)
    ax.text(px(height / 2), -2.45, "million years ago", ha="center",
            fontsize=fs.LABEL_PT, color=fs.INK)
    ax.set_xlim(-6, 420)
    ax.set_ylim(-2.9, len(tips) + 0.2)
    ax.axis("off")


def panel_c(ax):
    """Upper row: the RCA ranking and its controls. Lower row: the two analyses that do not depend on RCA."""
    def box(x, y, w, text, colour):
        ax.add_patch(Rectangle((x, y), w, 0.76, facecolor="none", edgecolor=colour, linewidth=1.0))
        ax.text(x + w / 2, y + 0.38, text, ha="center", va="center", fontsize=fs.LABEL_PT, color=colour)

    def arrow(x0, y0, dx, dy):
        ax.add_patch(FancyArrow(x0, y0, dx, dy, width=0.012, head_width=0.07, head_length=0.09,
                                color=fs.INK, length_includes_head=True))

    steps = ["proteomes\n11 species", "OrthoFinder\n18,902 orthogroups",
             "RCA\nnormalization", "estimators\nRCA score,\nRCA diversity, CA,\nFitness-Complexity",
             "controls\nnull model, rarefaction,\nspecies removal,\nPGLS, simulation"]
    widths = (1.30, 1.55, 1.15, 1.70, 1.90)
    y_top, xs, x = 1.15, [], 0.0
    for i, text in enumerate(steps):
        w = widths[i]
        box(x, y_top, w, text, fs.INK if i < 3 else fs.PURPLE)
        xs.append((x, w))
        if i < len(steps) - 1:
            arrow(x + w + 0.04, y_top + 0.38, 0.18, 0)
        x += w + 0.30
    # analyses independent of RCA, below the step they start from
    (x0, w0), (x1, w1) = xs[0], xs[1]
    box(x0, 0.0, 1.30 + 0.30 + 1.55 - 1.85, "KOfam\nKEGG module\ncompleteness", fs.TEAL)
    box(x1 + 0.10, 0.0, 2.95, "decomposition of gene number\ncoverage and genes per\noccupied orthogroup", fs.TEAL)
    arrow(x0 + 0.5, y_top - 0.04, 0, -0.30)
    arrow(x1 + w1 / 2, y_top - 0.04, 0, -0.30)
    ax.set_xlim(-0.1, x)
    ax.set_ylim(-0.12, y_top + 0.88)
    ax.axis("off")


def main():
    fig = plt.figure(figsize=(7.4, 7.8))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.95, 1.35, 0.80], hspace=0.30,
                          left=0.09, right=0.99, top=0.95, bottom=0.04)
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[i, 0]) for i in range(3))
    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)
    for ax, lab in ((ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)")):
        fs.panel_label(ax, lab, x=0.0, y=1.02)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "Figure_1_preliminary.pdf", format="pdf")
    print("wrote", OUTDIR / "Figure_1_preliminary.pdf")


if __name__ == "__main__":
    main()
