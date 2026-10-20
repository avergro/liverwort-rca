"""Plain figure style for the public repository.

Same names as the style module used for the published figures, so every
plot script runs unchanged, but with matplotlib defaults: default fonts and
colour cycle, full axes, standard legends. The scripts reproduce the content
of every panel from the result tables; the published figures were then
refined for typography, colour and layout.
"""

import matplotlib as mpl

BLUE, RED, TEAL = "C0", "C3", "C2"
BLUE_L = "C9"
GROUND = "0.6"
INK = "black"
PAPER = "white"
PURPLE, PURPLE_L, GREEN, MAGENTA = BLUE, BLUE_L, RED, TEAL
LABEL_PT = 8
AXIS_PT = 10


def hue(deg, tier="dark"):
    return BLUE


def apply():
    mpl.rcParams.update({"font.size": LABEL_PT, "pdf.fonttype": 42})


def bare(ax, keep=("left",)):
    """Kept for compatibility; the plain style leaves the axes untouched."""


def panel_label(ax, text, x=-0.14, y=1.04):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=AXIS_PT, ha="left", va="bottom")


def panel_label_aligned(ax, text, ref_ax, ref_x, y=1.04):
    """Panel label placed so that its left edge lines up, in the figure, with the label of ref_ax at ref_x."""
    r, p = ref_ax.get_position(), ax.get_position()
    panel_label(ax, text, x=(r.x0 + ref_x * r.width - p.x0) / p.width, y=y)


def shape_legend(ax, items, ncol=1, loc="lower left", anchor=(0.0, 1.02), color=None):
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], linestyle="none", marker=m, color=color or INK,
                      markerfacecolor=(color or INK) if filled else "none", label=lab)
               for m, filled, lab in items]
    return ax.legend(handles=handles, loc=loc, bbox_to_anchor=anchor, ncol=ncol)
