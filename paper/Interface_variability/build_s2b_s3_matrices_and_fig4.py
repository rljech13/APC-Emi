#!/usr/bin/env python
"""Build S2b/S3 distance matrices and regenerate Figure 4 (Emi2-native, multi-screen)."""
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import statistics as stats

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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
S2 = ROOT / "AF3_results_S2b"
S3 = ROOT / "AF3_results_S3"
IV = ROOT / "Interface_variability"
OUT = ROOT / "Figures"
OFF = 534


def parse_cif(path: Path):
    text = path.read_text(errors="replace").splitlines()
    i = 0
    atoms = []
    while i < len(text):
        if text[i].startswith("loop_"):
            i += 1
            cols = []
            while i < len(text) and text[i].startswith("_"):
                cols.append(text[i].strip())
                i += 1
            if any(c.startswith("_atom_site.") for c in cols):
                idx = {c.split(".")[-1]: n for n, c in enumerate(cols)}
                while (
                    i < len(text)
                    and text[i]
                    and not text[i].startswith("#")
                    and not text[i].startswith("loop_")
                    and not text[i].startswith("_")
                ):
                    line = text[i].strip()
                    i += 1
                    if not line:
                        continue
                    parts = []
                    cur = ""
                    inq = False
                    for ch in line:
                        if ch == '"':
                            inq = not inq
                            continue
                        if ch.isspace() and not inq:
                            if cur:
                                parts.append(cur)
                                cur = ""
                        else:
                            cur += ch
                    if cur:
                        parts.append(cur)
                    try:
                        atom = parts[idx["label_atom_id"]]
                        if atom.startswith("H"):
                            continue
                        atoms.append(
                            (
                                parts[idx["label_asym_id"]],
                                int(parts[idx["label_seq_id"]]),
                                atom,
                                float(parts[idx["Cartn_x"]]),
                                float(parts[idx["Cartn_y"]]),
                                float(parts[idx["Cartn_z"]]),
                            )
                        )
                    except Exception:
                        pass
                continue
        i += 1
    return atoms


def min_dist_maps(atoms, partner_chain: str):
    by_a: dict[int, list] = defaultdict(list)
    partner = []
    for ch, r, _at, x, y, z in atoms:
        if ch == "A":
            by_a[r].append((x, y, z))
        elif ch == partner_chain:
            partner.append((x, y, z))
    out = {}
    for r, pts in by_a.items():
        best = 1e9
        for x, y, z in pts:
            for x2, y2, z2 in partner:
                d = (x - x2) ** 2 + (y - y2) ** 2 + (z - z2) ** 2
                if d < best:
                    best = d
        out[r] = math.sqrt(best) if best < 1e9 else None
    return out


def build_matrices():
    s2tab = list(csv.DictReader((S2 / "AF3_S2b_model_table.csv").open()))
    eng = [(r["job"], int(r["model"])) for r in s2tab if float(r["pair_AC"] or 0) >= 0.5]
    print("S2b engaged", len(eng))

    acc_apc2: dict[int, list] = defaultdict(list)
    acc_apc11: dict[int, list] = defaultdict(list)
    for job, mid in eng:
        atoms = parse_cif(next((S2 / job).glob(f"*_model_{mid}.cif")))
        m2 = min_dist_maps(atoms, "B")
        m11 = min_dist_maps(atoms, "C")
        for r, d in m2.items():
            if d is not None:
                acc_apc2[r + OFF].append(d)
        for r, d in m11.items():
            if d is not None:
                acc_apc11[r + OFF].append(d)

    mat_s2 = {
        "source": "AF3_S2b engaged models (pair_iptm_emi2_apc11 >= 0.50)",
        "n_engaged": len(eng),
        "Apc2": {str(k): stats.median(v) for k, v in sorted(acc_apc2.items())},
        "Apc11": {str(k): stats.median(v) for k, v in sorted(acc_apc11.items())},
    }
    (IV / "emi2_s2b_distance_matrix.json").write_text(json.dumps(mat_s2, indent=2) + "\n")
    print(
        "S2b 590 Apc2",
        mat_s2["Apc2"].get("590"),
        "Apc11",
        mat_s2["Apc11"].get("590"),
    )

    s3tab = list(csv.DictReader((S3 / "AF3_S3_model_table.csv").open()))
    acc_apc1: dict[int, list] = defaultdict(list)
    for r in s3tab:
        job, mid = r["job"], int(r["model"])
        atoms = parse_cif(next((S3 / job).glob(f"*_model_{mid}.cif")))
        m1 = min_dist_maps(atoms, "B")
        for rr, d in m1.items():
            if d is not None:
                acc_apc1[rr + OFF].append(d)
    mat_s3 = {
        "source": "AF3_S3 ALL models (pair_iptm never >=0.5; NON-MEASUREMENT)",
        "n_models": len(s3tab),
        "Apc1": {str(k): stats.median(v) for k, v in sorted(acc_apc1.items())},
        "Apc1_note": "NON-MEASUREMENT: Emi2–Apc1 pair_iptm 0.07–0.21 in all 60 models.",
    }
    (IV / "emi2_s3_distance_matrix.json").write_text(json.dumps(mat_s3, indent=2) + "\n")
    print("S3 590 Apc1", mat_s3["Apc1"].get("590"))

    mat_s1 = json.loads((IV / "emi2_s1_distance_matrix.json").read_text())
    combined = {
        "Cdc20": mat_s1["Cdc20"],
        "Apc10": mat_s1["Apc10"],
        "Apc2": mat_s2["Apc2"],
        "Apc11": mat_s2["Apc11"],
        "Apc1": mat_s3["Apc1"],
        "notes": {
            "Cdc20": "S1 engaged (pair_iptm_emi2_cdc20>=0.5), n=45",
            "Apc10": "S1 NON-MEASUREMENT (co-receptor never occupied)",
            "Apc2": f"S2b engaged (pair_iptm_emi2_apc11>=0.5), n={len(eng)}",
            "Apc11": f"S2b engaged (pair_iptm_emi2_apc11>=0.5), n={len(eng)}",
            "Apc1": "S3 NON-MEASUREMENT (pair_iptm 0.07–0.21)",
        },
    }
    (IV / "emi2_combined_distance_matrix.json").write_text(json.dumps(combined, indent=2) + "\n")
    return combined, len(eng)


def make_fig4(combined, n_s2b_eng: int):
    rows = [
        ("Cdc20", "Cdc20", False),
        ("Apc10 (non-meas.)", "Apc10", True),
        ("Apc2", "Apc2", False),
        ("Apc11", "Apc11", False),
        ("Apc1 (non-meas.)", "Apc1", True),
    ]
    lo, hi = 2.0, 20.0
    length = 675
    fbox_domain = (452, 537)
    fbox_tint = "#F8C9C4"
    features = [
        (479, 519, "FBOX", "#E8452C"),
        (548, 551, "DB", "#4EC3D8"),
        (600, 648, "ZBR", "#00A57F"),
        (658, 675, "RL", "#F5A084"),
    ]
    muts = [76, 92, 102, 171, 219, 391, 395, 425, 549, 590]
    zlo, zhi = 535, length
    cov_lo, cov_hi = 535, 675
    cmap = LinearSegmentedColormap.from_list(
        "contact", ["#7F0000", "#C62828", "#EF6C00", "#FDD835", "#E8EAF0", "#F7F8FA"]
    )
    cmap.set_bad("#FFFFFF")

    fig_w = 7.09
    v_bot, v_band = 0.060, 4.55
    fig_h = v_bot + v_band + 0.330

    def band_y(offset_in: float) -> float:
        return (v_bot + offset_in) / fig_h

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=400)
    gs = fig.add_gridspec(
        4,
        2,
        height_ratios=[0.38, 0.28, 1.55, 0.05],
        width_ratios=[1, 0.022],
        hspace=0.38,
        wspace=0.03,
        left=0.24,
        right=0.955,
        top=band_y(v_band),
        bottom=band_y(0.0),
    )

    def strip(ax, lo_r, hi_r, *, ticks, mut_labels):
        ax.set_xlim(lo_r, hi_r)
        ax.set_ylim(0, 1)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_yticks([])
        ax.add_patch(
            Rectangle((lo_r, 0.30), hi_r - lo_r, 0.34, facecolor="#E4E5E7", edgecolor="none")
        )
        fa, fb = fbox_domain
        if not (fb < lo_r or fa > hi_r):
            fa2, fb2 = max(fa, lo_r), min(fb, hi_r)
            ax.add_patch(
                Rectangle(
                    (fa2, 0.30),
                    fb2 - fa2,
                    0.34,
                    facecolor=fbox_tint,
                    edgecolor="#E8452C",
                    lw=0.7,
                    zorder=2,
                )
            )
        for a, b, lab, col in features:
            if b < lo_r or a > hi_r:
                continue
            a2, b2 = max(a, lo_r), min(b, hi_r)
            ax.add_patch(
                Rectangle(
                    (a2, 0.30),
                    max(b2 - a2, length * 0.004),
                    0.34,
                    facecolor=col,
                    edgecolor="none",
                    zorder=3,
                )
            )
            cx = sum(fbox_domain) / 2 if lab == "FBOX" else (a2 + b2) / 2
            ax.text(
                cx,
                0.72,
                lab,
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="bold",
                color=col,
            )
        for m in muts:
            if lo_r <= m <= hi_r:
                ax.plot([m], [0.47], marker="o", ms=3.4, color="#111111", zorder=5)
                if mut_labels:
                    ax.annotate(
                        str(m),
                        xy=(m, 0.30),
                        xytext=(m, 0.02),
                        fontsize=7.4,
                        ha="center",
                        va="top",
                        color="#111111",
                        arrowprops=dict(
                            arrowstyle="-", lw=0.6, color="#111111", shrinkA=0, shrinkB=0
                        ),
                    )
        ax.set_xticks(ticks)
        ax.tick_params(labelsize=7.5, length=2.5, pad=1.5)

    axa = fig.add_subplot(gs[0, 0])
    strip(axa, 1, length, ticks=[1, 100, 200, 300, 400, 500, 600, 675], mut_labels=False)
    axa.add_patch(
        Rectangle(
            (cov_lo, 0.21),
            cov_hi - cov_lo,
            0.47,
            facecolor="none",
            edgecolor="#1A237E",
            lw=1.0,
            ls=(0, (3, 2)),
            zorder=6,
        )
    )
    axa.text(
        (cov_lo + cov_hi) / 2,
        0.015,
        f"{cov_lo}\u2013{cov_hi} = AF3 S1/S2b/S3 window, expanded in (b)",
        fontsize=6.6,
        color="#1A237E",
        ha="center",
        va="bottom",
        clip_on=False,
    )
    _t = axa.text(1, 1.30, "Emi2 (", fontsize=8.5, color="#37474F", va="bottom")
    for _s, _it in (("D. valentini", True), ("), residue", False)):
        _t = axa.annotate(
            _s,
            xycoords=_t,
            xy=(1, 0),
            va="bottom",
            fontsize=8.5,
            color="#37474F",
            fontstyle="italic" if _it else "normal",
        )

    axb = fig.add_subplot(gs[1, 0])
    strip(axb, zlo, zhi, ticks=[], mut_labels=True)

    xs = list(range(zlo, zhi + 1))
    data = np.full((len(rows), len(xs)), np.nan)
    for r, (_, key, _) in enumerate(rows):
        series = combined[key]
        for c, p2 in enumerate(xs):
            v = series.get(str(p2))
            if v is not None:
                data[r, c] = v

    ax = fig.add_subplot(gs[2, 0])
    im = ax.imshow(
        data,
        aspect="auto",
        cmap=cmap,
        vmin=lo,
        vmax=hi,
        extent=(zlo - 0.5, zhi + 0.5, len(rows) - 0.5, -0.5),
        interpolation="nearest",
    )
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for lab, _, _ in rows], fontsize=8.2)
    for i, (_, _, non) in enumerate(rows):
        if non:
            ax.get_yticklabels()[i].set_color("#78838E")
    ax.set_xlim(zlo - 0.5, zhi + 0.5)
    ax.set_xticks([540, 560, 580, 600, 620, 640, 660, 675])
    ax.set_xlabel("Emi2 residue number", fontsize=9, labelpad=2)
    ax.tick_params(labelsize=8)

    for i, (_, _, non) in enumerate(rows):
        if not non:
            continue
        ax.add_patch(
            Rectangle(
                (zlo - 0.5, i - 0.5),
                zhi - zlo + 1,
                1.0,
                facecolor="none",
                edgecolor="#B0B7C0",
                lw=0.4,
                hatch="///",
                zorder=2,
                alpha=0.55,
            )
        )

    ax.axvline(549, color="#1A237E", lw=0.8, ls=(0, (2.5, 1.8)), zorder=5)
    ax.axvline(590, color="#B71C1C", lw=0.8, ls=(0, (2.5, 1.8)), zorder=5)

    d549_cdc = float(combined["Cdc20"].get("549", float("nan")))
    d590_cdc = float(combined["Cdc20"].get("590", float("nan")))
    d590_a2 = float(combined["Apc2"].get("590", float("nan")))
    d590_a11 = float(combined["Apc11"].get("590", float("nan")))
    ax.text(
        0.02,
        -0.12,
        f"549: med {d549_cdc:.1f} \u00c5 to Cdc20 (side chain never \u22644 \u00c5);  "
        f"590: med {d590_cdc:.1f}/{d590_a2:.1f}/{d590_a11:.1f} \u00c5 to Cdc20/Apc2/Apc11 "
        f"(never \u22644 \u00c5; closest Apc11 6.5 \u00c5). S2b n={n_s2b_eng}.",
        transform=ax.transAxes,
        fontsize=6.4,
        color="#37474F",
        ha="left",
        va="top",
    )

    cax = fig.add_subplot(gs[2, 1])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("min. distance to partner (\u00c5)", fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    cb.ax.invert_yaxis()

    fig.text(0.012, band_y(v_band + 0.210), "a", fontsize=13, fontweight="bold", va="top")
    fig.text(0.012, band_y(3.55), "b", fontsize=13, fontweight="bold", va="top")

    for ext in ("svg", "png"):
        p = OUT / f"Figure_4_interface_map.{ext}"
        fig.savefig(p, dpi=400, bbox_inches="tight", pad_inches=0.02, facecolor="white")
        print("wrote", p)


if __name__ == "__main__":
    combined, n_eng = build_matrices()
    make_fig4(combined, n_eng)
