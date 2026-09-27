#!/usr/bin/env python3
"""Render Emi2 domain-colored structure figures for manuscript."""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

OUT_DIR = Path(__file__).resolve().parent
MODELS = Path(__file__).resolve().parents[1] / "Fasta" / "Models"
CLUSPRO_DNAI = Path("/Users/dmitrij/Downloads/cluspro.1372497/model.000.00.pdb")
CLUSPRO_DVAL = Path("/Users/dmitrij/Downloads/cluspro.1372499/model.004.00.pdb")

DOMAINS = [
    (400, 407, "N-proximal", "#BDBDBD"),
    (408, 458, "D-box (APC/C docking)", "#5C6BC0"),
    (459, 587, "Linker / middle", "#FFB74D"),
    (588, 636, "ZBR", "#E53935"),
    (637, 675, "RL-tail", "#26C6DA"),
]

HIGHLIGHT = {
    425: ("#00E676", 120, "425"),
    426: ("#69F0AE", 60, "426"),
}


def domain_color(resseq: int) -> str:
    for start, end, _, color in DOMAINS:
        if start <= resseq <= end:
            return color
    return "#9E9E9E"


def parse_pdb(path: Path) -> dict:
    atoms: dict[tuple[str, int, str], tuple[float, float, float]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            atom = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resseq = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            atoms[(chain, resseq, atom)] = (x, y, z)
            atoms.setdefault(("meta", resseq, "resname"), (0, 0, 0))
    return atoms


def load_chain(path: Path, chain_id: str) -> dict[int, dict[str, np.ndarray]]:
    residues: dict[int, dict[str, np.ndarray]] = defaultdict(dict)
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain_id:
                continue
            resseq = int(line[22:26])
            atom = line[12:16].strip()
            coord = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float
            )
            residues[resseq][atom] = coord
    return dict(residues)


def ca_trace(residues: dict[int, dict[str, np.ndarray]]) -> tuple[np.ndarray, list[int]]:
    order = sorted(residues)
    coords = np.array([residues[r]["CA"] for r in order if "CA" in residues[r]])
    return coords, order


def set_equal_3d(ax, points: np.ndarray, pad: float = 0.15) -> None:
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    center = (mins + maxs) / 2
    radius = max((maxs - mins).max() / 2, 1.0)
    radius *= 1 + pad
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)


def draw_backbone(ax, residues: dict[int, dict[str, np.ndarray]], lw: float = 2.8) -> None:
    order = sorted(residues)
    for i in range(len(order) - 1):
        r1, r2 = order[i], order[i + 1]
        if "CA" not in residues[r1] or "CA" not in residues[r2]:
            continue
        p1, p2 = residues[r1]["CA"], residues[r2]["CA"]
        color = domain_color(r1)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color=color, lw=lw, solid_capstyle="round")


def draw_residue_spheres(ax, residues: dict[int, dict[str, np.ndarray]], resseq: int, color: str, size: float, label: str | None = None) -> None:
    if resseq not in residues:
        return
    atoms = residues[resseq]
    if "CA" in atoms:
        ca = atoms["CA"]
        ax.scatter(*ca, s=size, c=color, edgecolors="black", linewidths=0.6, zorder=5)
        if label:
            ax.text(ca[0], ca[1], ca[2], label, fontsize=8, fontweight="bold", color="black")


def draw_sidechain(ax, residues: dict[int, dict[str, np.ndarray]], resseq: int, color: str) -> None:
    if resseq not in residues:
        return
    atoms = residues[resseq]
    if "CA" not in atoms:
        return
    ca = atoms["CA"]
    heavy = [atoms[a] for a in atoms if a not in {"H", "HA", "HB2", "HB3"} and not a.startswith("H")]
    if not heavy:
        return
    sc = np.mean(heavy, axis=0)
    ax.plot([ca[0], sc[0]], [ca[1], sc[1]], [ca[2], sc[2]], color=color, lw=3.5)
    ax.scatter(*sc, s=90, c=color, edgecolors="black", linewidths=0.5, zorder=6)


def domain_legend(ax) -> None:
    handles = [
        patches.Patch(facecolor=color, edgecolor="black", linewidth=0.4, label=label)
        for _, _, label, color in DOMAINS
    ]
    handles.append(patches.Patch(facecolor="#00E676", edgecolor="black", linewidth=0.4, label="aa425 (SNP)"))
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.95)


def plot_domain_map(out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 2.2), dpi=200)
    y = 0.35
    h = 0.35
    for start, end, label, color in DOMAINS:
        width = end - start + 1
        ax.add_patch(patches.Rectangle((start, y), width, h, facecolor=color, edgecolor="black", linewidth=0.6))
        if width > 25:
            ax.text((start + end) / 2, y + h / 2, label, ha="center", va="center", fontsize=8, fontweight="bold")
    ax.scatter([425], [y + h / 2], s=120, c="#00E676", edgecolors="black", zorder=5)
    ax.text(425, y + h + 0.12, "425\nL/F", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_xlim(398, 677)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Emi2 residue number (aa)", fontsize=10)
    ax.set_title("Emi2 functional domains (fragment 400–675)", fontsize=12, fontweight="bold")
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_standalone(path: Path, title: str, aa425: str, out_path: Path) -> None:
    residues = load_chain(path, "A")
    coords, _ = ca_trace(residues)

    fig = plt.figure(figsize=(10, 8), dpi=200)
    views = [(25, -55, "Overview"), (90, 0, "Side view"), (0, 0, "Front view")]
    for idx, (elev, azim, subtitle) in enumerate(views, 1):
        ax = fig.add_subplot(2, 2, idx, projection="3d")
        draw_backbone(ax, residues, lw=3.0)
        draw_sidechain(ax, residues, 425, "#00E676")
        draw_residue_spheres(ax, residues, 425, "#00E676", 140, aa425[0] + "425")
        set_equal_3d(ax, coords)
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(subtitle, fontsize=10, pad=0)

    ax_leg = fig.add_subplot(2, 2, 4)
    ax_leg.axis("off")
    domain_legend(ax_leg)
    ax_leg.text(
        0.02,
        0.55,
        f"{title}\nAlphaFold3 model_0 (aa 400–675)\n425 = {aa425}",
        transform=ax_leg.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_docking(path: Path, title: str, aa425: str, out_path: Path) -> None:
    emi2 = load_chain(path, "A")
    anapc2 = load_chain(path, "N")
    emi_coords, _ = ca_trace(emi2)
    apc_coords, _ = ca_trace(anapc2)

    fig = plt.figure(figsize=(12, 6), dpi=200)

    # Panel A: full complex
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    for resseq, atoms in sorted(anapc2.items()):
        if "CA" not in atoms:
            continue
        ax1.scatter(*atoms["CA"], s=8, c="#EF5350", alpha=0.85)
    draw_backbone(ax1, emi2, lw=2.5)
    draw_sidechain(ax1, emi2, 425, "#00E676")
    draw_residue_spheres(ax1, emi2, 425, "#00E676", 120, aa425[0] + "425")
    all_pts = np.vstack([emi_coords, apc_coords])
    set_equal_3d(ax1, all_pts, pad=0.05)
    ax1.view_init(elev=20, azim=-70)
    ax1.set_axis_off()
    ax1.set_title("Emi2–APC/C complex (ClusPro)", fontsize=10)

    # Panel B: D-box zoom
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    interface_emi = {r: atoms for r, atoms in emi2.items() if 408 <= r <= 458}
    nearby_apc = []
    for r, atoms in anapc2.items():
        if "CA" not in atoms:
            continue
        for er, e_atoms in interface_emi.items():
            if "CA" not in e_atoms:
                continue
            if np.linalg.norm(atoms["CA"] - e_atoms["CA"]) < 8.0:
                nearby_apc.append(r)
                break
    apc_subset = {r: anapc2[r] for r in set(nearby_apc)}
    for resseq, atoms in sorted(apc_subset.items()):
        if "CA" not in atoms:
            continue
        ax2.scatter(*atoms["CA"], s=35, c="#EF5350", alpha=0.95, edgecolors="darkred", linewidths=0.3)
    draw_backbone(ax2, interface_emi, lw=4.0)
    draw_sidechain(ax2, emi2, 425, "#00E676")
    draw_residue_spheres(ax2, emi2, 425, "#00E676", 160, aa425[0] + "425")
    if 426 in emi2:
        draw_residue_spheres(ax2, emi2, 426, "#69F0AE", 80, "426")
    zoom_pts = []
    for r in interface_emi:
        zoom_pts.append(interface_emi[r]["CA"])
    for r in apc_subset:
        zoom_pts.append(apc_subset[r]["CA"])
    zoom_pts = np.array(zoom_pts)
    set_equal_3d(ax2, zoom_pts, pad=0.25)
    ax2.view_init(elev=15, azim=-110)
    ax2.set_axis_off()
    ax2.set_title("D-box interface (408–458) × Anapc2", fontsize=10)

    fig.suptitle(f"{title}  |  Anapc2 = red  |  Emi2 domains colored  |  425 = {aa425}", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_combined_panel(out_path: Path) -> None:
    fig = plt.figure(figsize=(14, 12), dpi=200)

    standalone = [
        (MODELS / "Emi2_Dval_interface400-675_AF3_model0.pdb", "D. valentini", "Leu425"),
        (MODELS / "Emi2_Dnai_interface400-675_AF3_model0.pdb", "D. raddei nairensis", "Phe425"),
    ]
    for idx, (pdb, sp, aa425) in enumerate(standalone, 1):
        residues = load_chain(pdb, "A")
        coords, _ = ca_trace(residues)
        ax3d = fig.add_subplot(2, 2, idx, projection="3d")
        draw_backbone(ax3d, residues, lw=3.0)
        draw_sidechain(ax3d, residues, 425, "#00E676")
        draw_residue_spheres(ax3d, residues, 425, "#00E676", 130, aa425[0] + "425")
        set_equal_3d(ax3d, coords)
        ax3d.view_init(elev=25, azim=-55)
        ax3d.set_axis_off()
        ax3d.set_title(f"A{idx}. {sp} — AF3 Emi2", fontsize=10, fontweight="bold")

    docking = [
        (CLUSPRO_DVAL, "Dval Emi2 × Dnai APC/C", "Leu425"),
        (CLUSPRO_DNAI, "Dnai Emi2 × Dnai APC/C", "Phe425"),
    ]
    for i, (pdb, title, aa425) in enumerate(docking):
        ax3d = fig.add_subplot(2, 2, 3 + i, projection="3d")
        emi2 = load_chain(pdb, "A")
        anapc2 = load_chain(pdb, "N")
        for r, atoms in sorted(anapc2.items()):
            if "CA" in atoms:
                ax3d.scatter(*atoms["CA"], s=6, c="#EF5350", alpha=0.8)
        iface = {r: atoms for r, atoms in emi2.items() if 408 <= r <= 458}
        draw_backbone(ax3d, emi2, lw=2.2)
        draw_backbone(ax3d, iface, lw=4.5)
        draw_sidechain(ax3d, emi2, 425, "#00E676")
        draw_residue_spheres(ax3d, emi2, 425, "#00E676", 120, aa425[0] + "425")
        pts = np.vstack([ca_trace(emi2)[0], ca_trace(anapc2)[0]])
        set_equal_3d(ax3d, pts, pad=0.05)
        ax3d.view_init(elev=18, azim=-85)
        ax3d.set_axis_off()
        ax3d.set_title(f"B{i+1}. {title}", fontsize=10, fontweight="bold")

    handles = [patches.Patch(facecolor=c, edgecolor="black", linewidth=0.4, label=l) for _, _, l, c in DOMAINS]
    handles += [
        patches.Patch(facecolor="#EF5350", label="Anapc2 (APC/C)"),
        patches.Patch(facecolor="#00E676", label="aa425 SNP"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=True, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle("Emi2 domain architecture and APC/C docking", fontsize=14, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0.05, 1, 0.96])
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_domain_map(OUT_DIR / "Fig_Emi2_domain_map.png")
    plot_standalone(
        MODELS / "Emi2_Dval_interface400-675_AF3_model0.pdb",
        "Emi2 — D. valentini (paternal)",
        "Leu425",
        OUT_DIR / "Fig_Emi2_Dval_AF3_domains.png",
    )
    plot_standalone(
        MODELS / "Emi2_Dnai_interface400-675_AF3_model0.pdb",
        "Emi2 — D. raddei nairensis (maternal)",
        "Phe425",
        OUT_DIR / "Fig_Emi2_Dnai_AF3_domains.png",
    )
    plot_docking(
        CLUSPRO_DVAL,
        "Heterospecific: Dval Emi2 × Dnai APC/C",
        "Leu425",
        OUT_DIR / "Fig_Emi2_Dval_docking_domains.png",
    )
    plot_docking(
        CLUSPRO_DNAI,
        "Homospecific: Dnai Emi2 × Dnai APC/C",
        "Phe425",
        OUT_DIR / "Fig_Emi2_Dnai_docking_domains.png",
    )
    plot_combined_panel(OUT_DIR / "Fig_Emi2_combined_panel.png")
    print("Saved figures to", OUT_DIR)
    for p in sorted(OUT_DIR.glob("Fig_Emi2*.png")):
        print(" ", p.name)


if __name__ == "__main__":
    main()
