#!/usr/bin/env python
"""Figure 2 — story order: (a) 3D Emi2 → (b) pLDDT → (c) map + interface heatmap."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
import statistics as stats

import matplotlib as mpl
import numpy as np
from PIL import Image

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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import ConnectionPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "Figures"
IV = ROOT / "Interface_variability"
PDB = ROOT / "Fasta/Models/Emi2_Dval_interface400-675_AF3_model0.pdb"
MAT = json.loads((IV / "emi2_combined_distance_matrix.json").read_text())

# --- pLDDT from model ---
by_res: dict[int, list[float]] = defaultdict(list)
for line in PDB.read_text().splitlines():
    if line.startswith(("ATOM", "HETATM")):
        by_res[int(line[22:26])].append(float(line[60:66]))


def region_plddt(lo: int, hi: int):
    vals = [b for r, bs in by_res.items() if lo <= r <= hi for b in bs]
    vals_s = sorted(vals)
    n = len(vals_s)
    return (
        stats.mean(vals_s),
        vals_s[max(0, n // 10)],
        vals_s[min(n - 1, 9 * n // 10)],
    )


# rows tile 400–675; Skp1 motif indented under F-box
PLDDT_ROWS = [
    ("pre-F-box 400–451", 400, 451, False, "disordered", "—"),
    ("F-box 452–537", 452, 537, False, "fold", "no (Skp1)"),
    ("Skp1 479–519", 479, 519, True, "fold", "no (Skp1)"),
    ("538–547", 538, 547, False, "spacer", "—"),
    ("D-box 548–551", 548, 551, False, "degron", "yes (Cdc20)"),
    ("linker 552–599", 552, 599, False, "linker", "yes*"),
    ("ZBR 600–648", 600, 648, False, "fold", "yes (Apc2/11)"),
    ("649–657", 649, 657, False, "spacer", "—"),
    ("RL 658–675", 658, 675, False, "SLiM", "yes (Apc2)"),
]

# AlphaFold / AFDB official pLDDT legend colours
BAND_COLORS = {
    "very low": "#FF7D45",   # < 50
    "low": "#FFDB13",        # 50–70
    "confident": "#65CBF3",  # 70–90
    "very high": "#0053D6",  # > 90
}


def band(m: float) -> str:
    if m < 50:
        return "very low"
    if m < 70:
        return "low"
    if m < 90:
        return "confident"
    return "very high"


# --- figure geometry ---
FIG_W = 7.20
FIG_H = 9.40
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=400)

# Top: a | b(bars) | b(table) — keep b compact so it fits
gs_top = fig.add_gridspec(
    1,
    3,
    left=0.04,
    right=0.985,
    top=0.965,
    bottom=0.58,
    wspace=0.06,
    width_ratios=[1.35, 0.85, 0.42],
)
ax_a = fig.add_subplot(gs_top[0, 0])
ax_b = fig.add_subplot(gs_top[0, 1])
ax_btab = fig.add_subplot(gs_top[0, 2], sharey=ax_b)

# Bottom: c (map + heatmap)
gs_c = fig.add_gridspec(
    3,
    2,
    left=0.14,
    right=0.92,
    top=0.545,
    bottom=0.05,
    height_ratios=[0.38, 0.28, 1.55],
    width_ratios=[1, 0.028],
    hspace=0.30,
    wspace=0.04,
)
ax_map = fig.add_subplot(gs_c[0, 0])
ax_zoom = fig.add_subplot(gs_c[1, 0])
ax_heat = fig.add_subplot(gs_c[2, 0])
ax_cbar = fig.add_subplot(gs_c[2, 1])

# ========== (a) 3D model ==========
pil = Image.open(FIG / "CX_2b_v3_regions.png").convert("RGBA")
arr = np.asarray(pil)
# auto-trim near-white margins, then pad a little
rgb = arr[..., :3]
alpha = arr[..., 3] if arr.shape[2] == 4 else np.full(rgb.shape[:2], 255)
ink = (rgb.min(axis=2) < 245) & (alpha > 10)
rows = np.where(ink.any(axis=1))[0]
cols = np.where(ink.any(axis=0))[0]
if len(rows) and len(cols):
    r0, r1 = rows[0], rows[-1]
    c0, c1 = cols[0], cols[-1]
    pad = 18
    r0, r1 = max(0, r0 - pad), min(arr.shape[0] - 1, r1 + pad)
    c0, c1 = max(0, c0 - pad), min(arr.shape[1] - 1, c1 + pad)
    arr = arr[r0 : r1 + 1, c0 : c1 + 1]
ax_a.imshow(arr, aspect="equal", interpolation="lanczos")
# shift content slightly right inside the panel
ax_a.set_xlim(-0.02 * arr.shape[1], 1.02 * arr.shape[1])
ax_a.set_ylim(arr.shape[0] * 1.02, -0.02 * arr.shape[0])
ax_a.set_xticks([])
ax_a.set_yticks([])
for s in ax_a.spines.values():
    s.set_visible(False)
# compact region key under a
key_items = [
    ("#F8C9C4", "F-box domain"),
    ("#E8452C", "Skp1 motif"),
    ("#4EC3D8", "D-box"),
    ("#00A57F", "ZBR"),
    ("#F5A084", "RL"),
    ("#B0B4BA", "unannotated"),
]
for i, (col, lab) in enumerate(key_items):
    x = 0.08 + (i % 3) * 0.31
    y = -0.08 - (i // 3) * 0.055
    ax_a.add_patch(
        Rectangle(
            (x, y),
            0.045,
            0.035,
            transform=ax_a.transAxes,
            facecolor=col,
            edgecolor="#333",
            lw=0.4,
            clip_on=False,
        )
    )
    ax_a.text(
        x + 0.055,
        y + 0.017,
        lab,
        transform=ax_a.transAxes,
        fontsize=6.2,
        va="center",
        clip_on=False,
    )

fig.text(0.02, 0.972, "a", fontsize=14, fontweight="bold", va="top")

# ========== (b) pLDDT bars (0–100 only) + separate table ==========
labels = []
means = []
p10s = []
p90s = []
fills = []
chars = []
apcs = []
for lab, lo, hi, _ind, ch, ap in PLDDT_ROWS:
    m, p10, p90 = region_plddt(lo, hi)
    labels.append(lab)
    means.append(m)
    p10s.append(p10)
    p90s.append(p90)
    fills.append(BAND_COLORS[band(m)])
    chars.append(ch)
    apcs.append(ap)

y = np.arange(len(labels))[::-1]
for lo, hi, col in (
    (0, 50, BAND_COLORS["very low"]),
    (50, 70, BAND_COLORS["low"]),
    (70, 90, BAND_COLORS["confident"]),
    (90, 100, BAND_COLORS["very high"]),
):
    ax_b.axvspan(lo, hi, facecolor=col, alpha=0.22, zorder=0, linewidth=0)
ax_b.barh(
    y,
    means,
    color=fills,
    edgecolor="#333333",
    lw=0.35,
    height=0.68,
    xerr=[np.array(means) - np.array(p10s), np.array(p90s) - np.array(means)],
    error_kw=dict(ecolor="#444", lw=0.55, capsize=1.2),
    zorder=2,
)
for thr in (50, 70, 90):
    ax_b.axvline(thr, color="#666666", lw=0.4, ls=(0, (2, 2)), zorder=1)
ax_b.set_yticks(y)
ax_b.set_yticklabels(labels, fontsize=5.6)
ax_b.set_xlabel("mean pLDDT", fontsize=6.8, labelpad=1)
ax_b.tick_params(axis="x", labelsize=5.8, pad=1)
ax_b.set_xlim(0, 100)
ax_b.set_xticks([0, 50, 100])
ax_b.set_ylim(-0.55, len(labels) - 0.4)
for s in ("top", "right"):
    ax_b.spines[s].set_visible(False)

# table axis: compact, visually separated
ax_btab.set_xlim(0, 1)
ax_btab.set_ylim(ax_b.get_ylim())
ax_btab.set_xticks([])
ax_btab.tick_params(axis="y", left=False, labelleft=False)
ax_btab.set_facecolor("#F4F6F8")
for s in ax_btab.spines.values():
    s.set_visible(False)
ax_btab.spines["left"].set_visible(True)
ax_btab.spines["left"].set_color("#9AA3AD")
ax_btab.spines["left"].set_linewidth(1.0)
ax_btab.axvline(0.30, color="#D0D5DB", lw=0.5)
ax_btab.axvline(0.62, color="#D0D5DB", lw=0.5)
ax_btab.text(0.14, len(labels) - 0.18, "mean", fontsize=5.0, ha="center", color="#555", clip_on=False)
ax_btab.text(0.46, len(labels) - 0.18, "type", fontsize=5.0, ha="center", color="#555", clip_on=False)
ax_btab.text(0.82, len(labels) - 0.18, "APC/C", fontsize=5.0, ha="center", color="#555", clip_on=False)
# shorten APC labels for narrow table
apc_short = []
for ap in apcs:
    apc_short.append(
        {
            "—": "—",
            "no (Skp1)": "no",
            "yes (Cdc20)": "Cdc20",
            "yes*": "yes*",
            "yes (Apc2/11)": "Apc2/11",
            "yes (Apc2)": "Apc2",
        }.get(ap, ap)
    )
for yi, m, ch, ap in zip(y, means, chars, apc_short):
    ax_btab.text(0.14, yi, f"{m:.0f}", fontsize=5.5, ha="center", va="center", color="#222")
    ax_btab.text(0.46, yi, ch, fontsize=5.0, ha="center", va="center", color="#444")
    ax_btab.text(0.82, yi, ap, fontsize=4.9, ha="center", va="center", color="#222")

fig.text(0.50, 0.972, "b", fontsize=14, fontweight="bold", va="top")

# ========== (c) architecture + heatmap ==========
LEN = 675
FBOX_DOMAIN = (452, 537)
FBOX_TINT = "#F8C9C4"
FEATURES = [
    (479, 519, "FBOX", "#E8452C"),
    (548, 551, "DB", "#4EC3D8"),
    (600, 648, "ZBR", "#00A57F"),
    (658, 675, "RL", "#F5A084"),
]
MUTS = [76, 92, 102, 171, 219, 391, 395, 425, 549, 590]
ZLO, ZHI = 535, 675
COV_LO, COV_HI = 535, 675

ROWS = [
    ("Cdc20", "Cdc20", False),
    ("Apc10 (non-meas.)", "Apc10", True),
    ("Apc2", "Apc2", False),
    ("Apc11", "Apc11", False),
    ("Apc1 (non-meas.)", "Apc1", True),
]
cmap = LinearSegmentedColormap.from_list(
    "contact", ["#7F0000", "#C62828", "#EF6C00", "#FDD835", "#E8EAF0", "#F7F8FA"]
)
cmap.set_bad("#FFFFFF")


def strip(ax, lo, hi, *, ticks, mut_labels, title_muts=False):
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, 1)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_yticks([])
    ax.add_patch(Rectangle((lo, 0.30), hi - lo, 0.34, facecolor="#E4E5E7", edgecolor="none"))
    fa, fb = FBOX_DOMAIN
    if not (fb < lo or fa > hi):
        fa2, fb2 = max(fa, lo), min(fb, hi)
        ax.add_patch(
            Rectangle(
                (fa2, 0.30),
                fb2 - fa2,
                0.34,
                facecolor=FBOX_TINT,
                edgecolor="#E8452C",
                lw=0.7,
                zorder=2,
            )
        )
    for a, b, lab, col in FEATURES:
        if b < lo or a > hi:
            continue
        a2, b2 = max(a, lo), min(b, hi)
        ax.add_patch(
            Rectangle(
                (a2, 0.30),
                max(b2 - a2, LEN * 0.004),
                0.34,
                facecolor=col,
                edgecolor="none",
                zorder=3,
            )
        )
        cx = sum(FBOX_DOMAIN) / 2 if lab == "FBOX" else (a2 + b2) / 2
        ax.text(cx, 0.72, lab, ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=col)
    for m in MUTS:
        if lo <= m <= hi:
            ax.plot([m], [0.47], marker="o", ms=3.0, color="#111111", zorder=5)
            if mut_labels:
                ax.annotate(
                    str(m),
                    xy=(m, 0.30),
                    xytext=(m, 0.02),
                    fontsize=6.6,
                    ha="center",
                    va="top",
                    color="#111111",
                    arrowprops=dict(arrowstyle="-", lw=0.55, color="#111111", shrinkA=0, shrinkB=0),
                )
    ax.set_xticks(ticks)
    ax.tick_params(labelsize=6.8, length=2.2, pad=1.2)


strip(ax_map, 1, LEN, ticks=[1, 100, 200, 300, 400, 500, 600, 675], mut_labels=False)
ax_map.add_patch(
    Rectangle(
        (COV_LO, 0.21),
        COV_HI - COV_LO,
        0.47,
        facecolor="none",
        edgecolor="#1A237E",
        lw=1.0,
        ls=(0, (3, 2)),
        zorder=6,
    )
)
strip(ax_zoom, ZLO, ZHI, ticks=[], mut_labels=True)

# Zoom callout (full-width zoom). White pad under 300/675 wipes dashes there.
for x_full, x_zoom in ((COV_LO, ZLO), (COV_HI, ZHI)):
    fig.add_artist(
        ConnectionPatch(
            xyA=(x_full, 0.21),
            coordsA=ax_map.transData,
            xyB=(x_zoom, 0.95),
            coordsB=ax_zoom.transData,
            color="#1A237E",
            lw=0.9,
            ls=(0, (2.5, 1.8)),
            zorder=2,
        )
    )
for lab in ax_map.get_xticklabels():
    if lab.get_text() in {"300", "675"}:
        lab.set_bbox(dict(facecolor="white", edgecolor="none", pad=1.2, alpha=1.0))
        lab.set_zorder(50)

xs = list(range(ZLO, ZHI + 1))
data = np.full((len(ROWS), len(xs)), np.nan)
for r, (_, key, _) in enumerate(ROWS):
    series = MAT[key]
    for c, p2 in enumerate(xs):
        v = series.get(str(p2))
        if v is not None:
            data[r, c] = v

im = ax_heat.imshow(
    data,
    aspect="auto",
    cmap=cmap,
    vmin=2.0,
    vmax=20.0,
    extent=(ZLO - 0.5, ZHI + 0.5, len(ROWS) - 0.5, -0.5),
    interpolation="nearest",
)
ax_heat.set_yticks(range(len(ROWS)))
ax_heat.set_yticklabels([lab for lab, _, _ in ROWS], fontsize=7.4)
for i, (_, _, non) in enumerate(ROWS):
    if non:
        ax_heat.get_yticklabels()[i].set_color("#78838E")
ax_heat.set_xlim(ZLO - 0.5, ZHI + 0.5)
ax_heat.set_xticks([540, 560, 580, 600, 620, 640, 660, 675])
ax_heat.set_xlabel("Emi2 residue number", fontsize=8, labelpad=2)
ax_heat.tick_params(labelsize=7)

for i, (_, _, non) in enumerate(ROWS):
    if not non:
        continue
    ax_heat.add_patch(
        Rectangle(
            (ZLO - 0.5, i - 0.5),
            ZHI - ZLO + 1,
            1.0,
            facecolor="none",
            edgecolor="#B0B7C0",
            lw=0.4,
            hatch="///",
            zorder=2,
            alpha=0.55,
        )
    )

# allele markers — bright blue over warm heatmap
for x in (549, 590):
    ax_heat.axvline(x, color="#1565C0", lw=1.15, ls=(0, (2.8, 1.6)), zorder=5)

cb = fig.colorbar(im, cax=ax_cbar)
cb.set_label("min. distance (\u00c5)", fontsize=7.5)
cb.ax.tick_params(labelsize=6.5)
cb.ax.invert_yaxis()

fig.text(0.02, 0.555, "c", fontsize=14, fontweight="bold", va="top")

for ext in ("png", "svg"):
    out = FIG / f"Figure_2_model_and_map.{ext}"
    fig.savefig(out, dpi=400, bbox_inches="tight", pad_inches=0.03, facecolor="white")
    print("wrote", out)
