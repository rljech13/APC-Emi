#!/usr/bin/env python
"""Supplementary figure — partner maps + soft proximity to nearest Emi2 substitution.

Distance metric: for each partner residue, median over models of the minimum
heavy-atom distance to the nearest Emi2 substitution site present in the AF3
construct (F549I and/or V590M). Other parental Emi2 sites lie outside 535–675.

Includes Apc10 / Apc1 as non-measurements (pose not trusted / not engaged).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
mpl.rcParams.update(
    {
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
    }
)
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "Figures"
DATA = json.loads((ROOT / "Interface_variability/partner_distance_to_emi2_muts.json").read_text())
# Apc1 distances = S3 all models (same CIF pool as heatmap row); still non-measurement.
if "deletion" in DATA.get("Apc1", {}):
    DATA["Apc1"]["deletion"] = tuple(DATA["Apc1"]["deletion"])
if "s3_window" in DATA.get("Apc1", {}):
    DATA["Apc1"]["s3_window"] = tuple(DATA["Apc1"]["s3_window"])

EMI2_SUBS_FULL = [int(x) for x in DATA.get("emi2_substitutions_full", [])]
EMI2_SUBS_IN = {int(x) for x in DATA.get("emi2_substitutions_in_construct", [549, 590])}
EMI2_SUB_LABELS = {
    76: "Q76H",
    92: "V92I",
    102: "A102T",
    171: "K171R",
    219: "E219D",
    391: "I391V",
    395: "D395N",
    425: "L425F",
    549: "F549I",
    590: "V590M",
}

BAR = "#E4E5E7"
MUT = "#B71C1C"
MUT_HALO = "#FFECB3"
CONTACT_LINE = "#C62828"
DIST_FACE = "#156B5E"
EMI2_IN = "#1565C0"
EMI2_OUT = "#90A4AE"
SIGMA = 4.5
PROX_FLOOR = 0.04


def proximity(d: float) -> float:
    return float(np.exp(-max(d, 0.0) / SIGMA))


def smooth(y: np.ndarray, win: int = 5) -> np.ndarray:
    if win < 2 or len(y) < win:
        return y
    k = np.ones(win) / win
    sm = np.convolve(y, k, mode="same")
    sm[y <= 0] = 0.0
    return sm


def panel(fig, rect, name, blob, letter, calls):
    """Proximity profile plus a sequence bar. Tall = close to F549I or V590M."""
    left, bottom, width, height = rect
    ax = fig.add_axes([left, bottom + height * 0.22, width, height * 0.72])
    ax_m = fig.add_axes([left, bottom, width, height * 0.16], sharex=ax)

    length = int(blob["length"])
    dist = {int(k): float(v) for k, v in blob.get("dist", {}).items()}
    xs = np.arange(1, length + 1, dtype=float)
    prox = np.zeros(length, dtype=float)
    for i, d in dist.items():
        p = proximity(d)
        if p >= PROX_FLOOR:
            prox[i - 1] = p
    prox_s = smooth(prox, win=5 if length > 100 else 3)
    h4 = proximity(4.0)

    ax.fill_between(xs, 0, prox_s, color="#2A9D8F", alpha=0.35, linewidth=0, zorder=2)
    ax.plot(xs, prox_s, color="#1B7A6E", lw=1.05, zorder=3)
    ax.axhline(h4, color=CONTACT_LINE, lw=0.8, ls=(0, (3, 1.8)), zorder=4)
    ax.text(8, h4 + 0.03, "4 Å", fontsize=8, color=CONTACT_LINE, va="bottom")

    for pos, lab, dx in calls:
        y = dist[pos]
        ax.axvline(pos, color=MUT, lw=0.7, ls=(0, (2, 1.4)), alpha=0.8, zorder=5)
        ax.plot([pos], [0.06], marker="o", ms=5, color=MUT, zorder=6, markeredgecolor="white", markeredgewidth=0.5)
        ax.annotate(
            f"{lab}, {y:.0f} Å",
            xy=(pos, 0.06),
            xytext=(pos + dx, 0.22),
            fontsize=8,
            color=MUT,
            ha="center",
            va="bottom",
            arrowprops=dict(arrowstyle="-", color=MUT, lw=0.6),
            annotation_clip=False,
            zorder=7,
        )

    ax.set_xlim(0.5, length + 0.5)
    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, h4, proximity(2.0), 1.0])
    ax.set_yticklabels(["far", "4 Å", "2 Å", "contact"], fontsize=8)
    ax.set_ylabel("Proximity to nearest\nEmi2 substitution", fontsize=9)
    ax.tick_params(axis="x", labelbottom=False, length=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(-0.07, 1.04, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="right", clip_on=False)
    ax.text(0.0, 1.04, name, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    for s in ax_m.spines.values():
        s.set_visible(False)
    ax_m.set_yticks([])
    ax_m.set_ylim(0, 1)
    ax_m.add_patch(Rectangle((0.5, 0.35), length, 0.45, facecolor=BAR, edgecolor="#8A9096", lw=0.6, zorder=1))
    if name == "Anapc2":
        ax_m.add_patch(Rectangle((0.5, 0.35), 219, 0.45, facecolor="#F2F2F2", edgecolor="none", zorder=2))
        ax_m.text(110, 0.57, "1–219 not in S2b", ha="center", va="center", fontsize=7.5, color="#777")
    for pos, _, _ in calls:
        ax_m.plot([pos], [0.57], marker="|", ms=9, mew=1.4, color=MUT, zorder=4)
    ax_m.set_xlabel(f"{name} residue", fontsize=9)
    ax_m.tick_params(labelsize=8, length=3)
    return ax


fig = plt.figure(figsize=(7.8, 6.2), facecolor="white")
panel(fig, [0.12, 0.54, 0.78, 0.40], "Cdc20", DATA["Cdc20"], "a", [(291, "V291I", 55)])
panel(fig, [0.12, 0.08, 0.78, 0.40], "Anapc2", DATA["Apc2"], "b", [(632, "V632A", -70), (654, "P654L", 70)])

out_png = FIG / "Figure_S_subs_vs_interface.png"
out_svg = FIG / "Figure_S_subs_vs_interface.svg"
fig.savefig(out_png, dpi=400, facecolor="white")
fig.savefig(out_svg, dpi=300, facecolor="white")
print("wrote", out_png)
print("wrote", out_svg)
