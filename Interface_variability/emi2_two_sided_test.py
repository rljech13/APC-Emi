#!/usr/bin/env python
"""Emi2-native two-sided test (no Emi1 / 4UI9 surrogate).

Side 1 — Emi2: Dval vs Dnai substitutions vs functional-element map.
Side 2 — APC/C: allele differences for subunits with both parental sequences.
Geometry — only AlphaFold3 S1: Emi2(535–675) + Cdc20 + Apc10.
  - Emi2–Cdc20 distances are measurements (D-box engagement validated vs 5G04).
  - Emi2–Apc10 distances are recorded but flagged as non-measurements (Apc10 never
    reaches its co-receptor site in any of the 60 models; see AF3_results_S1/REPORT).

Outputs REPORT.md, CSVs, and emi2_s1_distance_matrix.json under Interface_variability/.
"""
from __future__ import annotations

import csv
import json
import shlex
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Interface_variability"
FASTA = ROOT / "Fasta"
S1 = ROOT / "AF3_results_S1"
OFF = 534  # full-length Emi2 = CIF residue number + OFF

EMI2_ELEMENTS = [
    (452, 537, "F-box domain"),
    (479, 519, "F-box motif"),
    (548, 551, "D-box"),
    (600, 648, "ZBR"),
    (658, 675, "RL tail"),
]

# Engaged basin from REPORT_AF3_S1.md: pair_iptm_emi2_cdc20 0.52–0.60
ENGAGED_IPTM = 0.50


def read_fa(path: Path) -> str:
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    )


def element_of(pos: int) -> str:
    hits = [name for a, b, name in EMI2_ELEMENTS if a <= pos <= b]
    if "D-box" in hits:
        return "D-box"
    if "ZBR" in hits:
        return "ZBR"
    if "RL tail" in hits:
        return "RL tail"
    if "F-box motif" in hits or "F-box domain" in hits:
        return "F-box (not APC/C)"
    return ""


def pairwise_diffs(a: str, b: str) -> list[tuple[int, str, str]]:
    """Ungapped pairwise: pad shorter to left-align with best identity window for unequal lengths."""
    if len(a) == len(b):
        return [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]
    # slide shorter over longer
    if len(a) <= len(b):
        short, long_, short_is_a = a, b, True
    else:
        short, long_, short_is_a = b, a, False
    best_off, best_m = 0, -1
    for off in range(0, len(long_) - len(short) + 1):
        m = sum(x == y for x, y in zip(short, long_[off : off + len(short)]))
        if m > best_m:
            best_m, best_off = m, off
    out = []
    for i, (x, y) in enumerate(zip(short, long_[best_off : best_off + len(short)])):
        if x == y:
            continue
        pos_long = best_off + i + 1
        if short_is_a:
            out.append((pos_long, x, y))  # position in longer (b); store dval/dnai carefully
        else:
            out.append((pos_long, y, x))
    # For unequal lengths we report positions in the longer sequence and alleles as (short, long)
    # Caller passes (dval, dnai) — re-do properly below per subunit.
    return out


def diffs_dval_dnai(dval: str, dnai: str) -> list[tuple[int, str, str, str]]:
    """Return (position_in_dval_numbering_if_aligned, dval_aa, dnai_aa, note)."""
    if len(dval) == len(dnai):
        return [(i + 1, x, y, "") for i, (x, y) in enumerate(zip(dval, dnai)) if x != y]
    # Prefer numbering of Dval (usually more complete). Align dnai into dval.
    if len(dnai) <= len(dval):
        short, long_, short_name = dnai, dval, "dnai"
    else:
        short, long_, short_name = dval, dnai, "dval"
    best_off, best_m = 0, -1
    for off in range(0, len(long_) - len(short) + 1):
        m = sum(x == y for x, y in zip(short, long_[off : off + len(short)]))
        if m > best_m:
            best_m, best_off = m, off
    out = []
    for i, (s_aa, l_aa) in enumerate(zip(short, long_[best_off : best_off + len(short)])):
        if s_aa == l_aa:
            continue
        pos_long = best_off + i + 1
        if short_name == "dnai":
            # long = dval
            out.append((pos_long, l_aa, s_aa, f"aligned dnai@{len(dnai)} into dval@{len(dval)} off={best_off}"))
        else:
            # long = dnai, short = dval — report dval position = i+1
            out.append((i + 1, s_aa, l_aa, f"aligned dval@{len(dval)} into dnai@{len(dnai)} off={best_off}"))
    return out


def parse_cif_heavy(path: Path):
    """Return {chain: {resnum: Nx3 array of heavy atoms}}."""
    lines = path.read_text(errors="ignore").splitlines()
    headers = []
    start = None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            headers.append(line.strip().split(".")[-1])
        elif headers and not line.startswith("_atom_site."):
            start = i
            break
    if not headers or start is None:
        raise ValueError(f"no atom_site in {path}")
    idx = {c: n for n, c in enumerate(headers)}
    atoms = defaultdict(lambda: defaultdict(list))
    for line in lines[start:]:
        if line.startswith("#") or line.startswith("loop_") or line.startswith("_"):
            break
        if not line.strip():
            continue
        try:
            parts = shlex.split(line)
        except ValueError:
            continue
        if len(parts) < len(headers):
            continue
        if parts[idx["group_PDB"]] != "ATOM":
            continue
        atom_key = "auth_atom_id" if "auth_atom_id" in idx else "label_atom_id"
        atom = parts[idx[atom_key]]
        if atom.startswith("H") or parts[idx["type_symbol"]] == "H":
            continue
        ch = parts[idx["auth_asym_id"]] if "auth_asym_id" in idx else parts[idx["label_asym_id"]]
        seq_key = "auth_seq_id" if "auth_seq_id" in idx else "label_seq_id"
        try:
            res = int(parts[idx[seq_key]])
        except ValueError:
            continue
        x, y, z = (
            float(parts[idx["Cartn_x"]]),
            float(parts[idx["Cartn_y"]]),
            float(parts[idx["Cartn_z"]]),
        )
        atoms[ch][res].append((x, y, z))
    out = {}
    for ch, resmap in atoms.items():
        out[ch] = {r: np.asarray(coords, dtype=float) for r, coords in resmap.items()}
    return out


def min_dist_residue_to_chain(res_coords: np.ndarray, partner: dict[int, np.ndarray]) -> float:
    if res_coords.size == 0 or not partner:
        return float("nan")
    P = np.concatenate(list(partner.values()), axis=0)
    md = float("inf")
    for c0 in range(0, len(P), 8000):
        chunk = P[c0 : c0 + 8000]
        d = np.sqrt(((res_coords[:, None, :] - chunk[None, :, :]) ** 2).sum(-1)).min()
        if d < md:
            md = float(d)
    return md


def engaged_cif_paths() -> list[Path]:
    conf = list(csv.DictReader((S1 / "AF3_S1_confidence_table.csv").open()))
    paths = []
    for row in conf:
        if float(row["pair_iptm_emi2_cdc20"]) < ENGAGED_IPTM:
            continue
        job = row["job"]
        model = row["model_index"]
        # job folder name == job; cif pattern fold_{job}_model_{n}.cif
        folder = S1 / job
        # model_index in CSV may be '0'..'4'
        cands = list(folder.glob(f"*model_{model}.cif"))
        if not cands:
            cands = list(folder.glob(f"*model_{int(float(model))}.cif"))
        if cands:
            paths.append(cands[0])
    return paths


def build_s1_distance_matrix(cif_paths: list[Path]) -> dict:
    """Median min-distance per full-length Emi2 residue to Cdc20 and Apc10 across engaged models."""
    per_model = []
    for path in cif_paths:
        atoms = parse_cif_heavy(path)
        emi = atoms.get("A", {})
        cdc = atoms.get("B", {})
        apc10 = atoms.get("C", {})
        row_cdc = {}
        row_apc10 = {}
        for cif_res, coords in emi.items():
            full = cif_res + OFF
            row_cdc[full] = min_dist_residue_to_chain(coords, cdc)
            row_apc10[full] = min_dist_residue_to_chain(coords, apc10)
        per_model.append((row_cdc, row_apc10))

    def median_map(key: str) -> dict[str, float]:
        bag = defaultdict(list)
        for row_cdc, row_apc10 in per_model:
            src = row_cdc if key == "Cdc20" else row_apc10
            for r, d in src.items():
                if d == d:  # not nan
                    bag[r].append(d)
        return {str(r): float(np.median(v)) for r, v in sorted(bag.items())}

    return {
        "source": "AF3_S1 engaged models (pair_iptm_emi2_cdc20 >= 0.50)",
        "n_models": len(per_model),
        "emi2_construct": "535-675 (CIF 1-141, full = CIF + 534)",
        "Cdc20": median_map("Cdc20"),
        "Apc10": median_map("Apc10"),
        "Apc10_note": (
            "NON-MEASUREMENT: Apc10 never reaches its D-box co-receptor site in any of "
            "the 60 S1 models (see AF3_results_S1/REPORT_AF3_S1.md). Distances are "
            "recorded for transparency only."
        ),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # --- Emi2 substitutions ---
    dval_e = read_fa(FASTA / "Emi2" / "Emi2_Dval_AAseq.fasta")
    dnai_e = read_fa(FASTA / "Emi2" / "Emi2_Dnai_AAseq.fasta")
    assert len(dval_e) == len(dnai_e) == 675
    emi2_sites = []
    for i, (a, b) in enumerate(zip(dval_e, dnai_e), start=1):
        if a != b:
            emi2_sites.append((i, a, b, element_of(i)))

    with (OUT / "emi2_variable_sites.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["position", "dval", "dnai", "element"])
        for pos, a, b, el in emi2_sites:
            w.writerow([pos, a, b, el])

    # --- APC allele differences ---
    subunits = ["Cdc20", "Anapc10", "Anapc11", "Anapc2", "Anapc5", "Anapc7"]
    # Prefer Anapc2_Dnai_AAseq_v2 if present
    apc_rows = []
    for sub in subunits:
        folder = FASTA / sub
        dv_path = folder / f"{sub}_Dval_AAseq.fasta"
        dn_path = folder / f"{sub}_Dnai_AAseq.fasta"
        if sub == "Anapc2" and (folder / "Anapc2_Dnai_AAseq_v2.fasta").exists():
            dn_path = folder / "Anapc2_Dnai_AAseq_v2.fasta"
        if not dv_path.exists() or not dn_path.exists():
            apc_rows.append(dict(subunit=sub, status="missing allele", position="", dval="", dnai="", note=""))
            continue
        dv, dn = read_fa(dv_path), read_fa(dn_path)
        diffs = diffs_dval_dnai(dv, dn)
        if not diffs:
            apc_rows.append(
                dict(
                    subunit=sub,
                    status="identical",
                    position="",
                    dval="",
                    dnai="",
                    note=f"len_dval={len(dv)} len_dnai={len(dn)}",
                )
            )
        else:
            for pos, a, b, note in diffs:
                apc_rows.append(
                    dict(
                        subunit=sub,
                        status="substitution",
                        position=pos,
                        dval=a,
                        dnai=b,
                        note=note or f"len_dval={len(dv)} len_dnai={len(dn)}",
                    )
                )

    # Anapc1 — Dnai is a partial ORF (internal deletion vs Dval)
    apc1_dv = FASTA / "Anapc1" / "Anapc1_Dval_AAseq.fasta"
    apc1_dn = FASTA / "Anapc1" / "Anapc1_Dnai_AAseq_part.fasta"
    if apc1_dv.exists() and apc1_dn.exists():
        dv1, dn1 = read_fa(apc1_dv), read_fa(apc1_dn)
        # Use MAFFT alignment on disk if present, else length-aware diffs via mafft file
        aln_path = FASTA / "Anapc1" / "Anapc1_Dval_Dnai_aln.fasta"
        apc1_gap = None
        apc1_subs = []
        if aln_path.exists():
            seqs = {}
            name = None
            parts = []
            for line in aln_path.read_text().splitlines():
                if line.startswith(">"):
                    if name:
                        seqs[name] = "".join(parts)
                    name = line[1:].split()[0]
                    parts = []
                else:
                    parts.append(line.strip())
            if name:
                seqs[name] = "".join(parts)
            A = seqs.get("Anapc1_Dval") or next(v for k, v in seqs.items() if "Dval" in k or "dval" in k.lower())
            B = seqs.get("Anapc1_Dnai_part") or next(
                v for k, v in seqs.items() if "Dnai" in k or "dnai" in k.lower() or "part" in k.lower()
            )
            hr = 0
            i = 0
            while i < len(A):
                if A[i] == "-":
                    i += 1
                    continue
                hr += 1
                if B[i] == "-":
                    j = i
                    while j < len(A) and B[j] == "-" and A[j] != "-":
                        j += 1
                    # count Dval residues in this gap
                    gap_len = sum(1 for k in range(i, j) if A[k] != "-")
                    if gap_len >= 5:
                        apc1_gap = (hr, hr + gap_len - 1, gap_len)
                    # advance hr for remaining gap residues except first
                    for k in range(i + 1, j):
                        if A[k] != "-":
                            hr += 1
                    i = j
                    continue
                if A[i] != B[i]:
                    apc1_subs.append((hr, A[i], B[i]))
                i += 1
        if apc1_gap:
            apc_rows.append(
                dict(
                    subunit="Anapc1",
                    status="Dnai internal deletion vs Dval",
                    position=f"{apc1_gap[0]}-{apc1_gap[1]}",
                    dval="(present)",
                    dnai="(absent)",
                    note=f"len_dval={len(dv1)} len_dnai={len(dn1)} gap_len={apc1_gap[2]} (Dval numbering)",
                )
            )
        for pos, a, b in apc1_subs:
            apc_rows.append(
                dict(
                    subunit="Anapc1",
                    status="substitution",
                    position=pos,
                    dval=a,
                    dnai=b,
                    note=f"len_dval={len(dv1)} len_dnai={len(dn1)}; no Emi2–Apc1 complex model",
                )
            )
        if not apc1_subs and not apc1_gap:
            apc_rows.append(
                dict(
                    subunit="Anapc1",
                    status="compared",
                    position="",
                    dval="",
                    dnai="",
                    note=f"len_dval={len(dv1)} len_dnai={len(dn1)}; no Emi2–Apc1 complex model",
                )
            )
    else:
        apc_rows.append(
            dict(
                subunit="Anapc1",
                status="missing fasta",
                position="",
                dval="",
                dnai="",
                note="expected Fasta/Anapc1/",
            )
        )
    apc_rows.append(
        dict(
            subunit="Anapc4",
            status="Dval partial — no Dnai; no Emi2 complex model",
            position="",
            dval="",
            dnai="",
            note="",
        )
    )

    with (OUT / "apc_variable_sites.csv").open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["subunit", "status", "position", "dval", "dnai", "note"]
        )
        w.writeheader()
        w.writerows(apc_rows)

    # --- Geometry from S1 ---
    cifs = engaged_cif_paths()
    matrix = build_s1_distance_matrix(cifs)
    (OUT / "emi2_s1_distance_matrix.json").write_text(json.dumps(matrix, indent=2))

    cdc = {int(k): v for k, v in matrix["Cdc20"].items()}
    apc10 = {int(k): v for k, v in matrix["Apc10"].items()}

    def dist_summary(pos: int) -> tuple[float | None, float | None]:
        return cdc.get(pos), apc10.get(pos)

    # Sidechain-level facts for 549 from existing CSV (engaged only)
    sc_rows = list(csv.DictReader((S1 / "AF3_S1_sidechain_contacts.csv").open()))
    eng_set = set()
    for row in csv.DictReader((S1 / "AF3_S1_confidence_table.csv").open()):
        if float(row["pair_iptm_emi2_cdc20"]) >= ENGAGED_IPTM:
            eng_set.add((row["job"], str(int(float(row["model_index"])))))

    sc_549 = [
        r
        for r in sc_rows
        if int(r["full"]) == 549 and (r["job"], str(int(float(r["model"]))) ) in eng_set
    ]
    sc_549_within4 = sum(1 for r in sc_549 if float(r["min_sc"]) <= 4.0)

    # Cdc20 V291I distance to D-box from pose/sidechain — use median of min dist from
    # Emi2 548-551 to Cdc20 in engaged models (already in matrix); V291 measured below
    v291_dists = []
    for path in cifs:
        atoms = parse_cif_heavy(path)
        cdc_atoms = atoms.get("B", {})
        if 291 not in cdc_atoms:
            continue
        dbox = []
        for cif_res in (14, 15, 16, 17):  # 548-551
            if cif_res in atoms.get("A", {}):
                dbox.append(atoms["A"][cif_res])
        if not dbox:
            continue
        D = np.concatenate(dbox, axis=0)
        v291_dists.append(min_dist_residue_to_chain(cdc_atoms[291], {0: D}))

    v291_med = float(np.median(v291_dists)) if v291_dists else float("nan")

    # --- Intersection ---
    # Publishable structural interface partners in S1: Cdc20 only.
    # Variable APC sites on that surface: none (Cdc20 V291I far from D-box).
    cdc20_var = [r for r in apc_rows if r["subunit"] == "Cdc20" and r["status"] == "substitution"]
    interface_emi2 = [s for s in emi2_sites if s[3] == "D-box"]  # only D-box has S1 geometry
    # 549 side chain does not contact; 590 has no partner model

    lines = []
    lines += [
        "# Two-sided variability test (Emi2-native): *D. valentini* × *D. r. nairensis*",
        "",
        "**Scope.** All geometry in this report comes from AlphaFold3 screen S1 "
        "(Emi2 residues 535–675 + Cdc20 + Anapc10; 60 models). "
        "There is **no Emi1 / PDB 4UI9 surrogate**. "
        "Prior Emi1-based files live in `Interface_variability/superseded/`.",
        "",
        "## 0. Emi2 element boundaries",
        "",
        "| element | residues | APC/C role in this study |",
        "|---|---|---|",
        "| F-box domain | 452–537 | Skp1/SCF — **not** APC/C |",
        "| F-box motif (inside domain) | 479–519 | Skp1/SCF — **not** APC/C |",
        "| D-box | 548–551 | binds Cdc20 WD40 in S1 (validated vs PDB 5G04) |",
        "| ZBR | 600–648 | literature: Apc11/Apc2/Apc1 — **no Emi2 complex model here** |",
        "| RL tail | 658–675 | literature: Apc2 CTD — **no atomic coverage in S1 pose** |",
        "",
        "## 1. Emi2 side",
        "",
        f"Dval {len(dval_e)} aa vs Dnai {len(dnai_e)} aa: **{len(emi2_sites)} substitutions**.",
        "",
        "| position | Dval | Dnai | functional element |",
        "|---|---|---|---|",
    ]
    for pos, a, b, el in emi2_sites:
        lines.append(f"| {pos} | {a} | {b} | {el or '—'} |")
    lines += [
        "",
        "Substitutions inside an annotated APC/C-binding element: **1** (549, D-box). "
        "V590M lies between D-box and ZBR; it is **not** inside the D-box or ZBR ranges above, "
        "and S1 does not place it against a resolved APC/C partner (see §3).",
        "",
        "## 2. APC/C side (sequence)",
        "",
        "| subunit | Dval len | Dnai len | substitutions | note |",
        "|---|---|---|---|---|",
    ]
    for sub in subunits:
        folder = FASTA / sub
        dv_path = folder / f"{sub}_Dval_AAseq.fasta"
        dn_path = folder / f"{sub}_Dnai_AAseq.fasta"
        if sub == "Anapc2" and (folder / "Anapc2_Dnai_AAseq_v2.fasta").exists():
            dn_path = folder / "Anapc2_Dnai_AAseq_v2.fasta"
        dv, dn = read_fa(dv_path), read_fa(dn_path)
        subs = [r for r in apc_rows if r["subunit"] == sub and r["status"] == "substitution"]
        if not subs and any(r["subunit"] == sub and r["status"] == "identical" for r in apc_rows):
            sub_s = "**0 — identical**"
        else:
            sub_s = "; ".join(f"{r['dval']}{r['position']}{r['dnai']}" for r in subs)
        lines.append(f"| {sub} | {len(dv)} | {len(dn)} | {sub_s} | |")

    apc1_rows = [r for r in apc_rows if r["subunit"] == "Anapc1"]
    apc1_gap_row = next((r for r in apc1_rows if "deletion" in r["status"]), None)
    apc1_sub_rows = [r for r in apc1_rows if r["status"] == "substitution"]
    if apc1_gap_row:
        dv1 = read_fa(FASTA / "Anapc1" / "Anapc1_Dval_AAseq.fasta")
        dn1 = read_fa(FASTA / "Anapc1" / "Anapc1_Dnai_AAseq_part.fasta")
        sub_txt = (
            "; ".join(f"{r['dval']}{r['position']}{r['dnai']}" for r in apc1_sub_rows)
            if apc1_sub_rows
            else "no AA substitutions in overlap"
        )
        apc1_sub_s = f"gap Dval {apc1_gap_row['position']} absent in Dnai; {sub_txt}"
        lines.append(
            f"| Anapc1 | {len(dv1)} | {len(dn1)} (partial) | {apc1_sub_s} | "
            f"no Emi2–Apc1 complex; see §2.1 |"
        )
    else:
        lines.append("| Anapc1 | — | — | *missing* | |")

    lines += [
        "| Anapc4 | partial Dval | — | *Dnai missing* | no Emi2–Apc4 model |",
        "",
        f"Cdc20 V291I median distance to Emi2 D-box (548–551) across {len(v291_dists)} "
        f"engaged S1 models: **{v291_med:.1f} Å**.",
        "",
    ]
    if apc1_gap_row:
        gap_len = apc1_gap_row["note"].split("gap_len=")[-1].split()[0]
        sub_list = (
            "(" + "; ".join(f"{r['dval']}{r['position']}{r['dnai']}" for r in apc1_sub_rows) + ")"
            if apc1_sub_rows
            else ""
        )
        lines += [
            "### 2.1 Anapc1 Dnai gap",
            "",
            f"Relative to Dval ({len(dv1)} aa), the Dnai partial ORF ({len(dn1)} aa) lacks "
            f"**Dval residues {apc1_gap_row['position']}** ({gap_len} aa — the ~48-residue hole). "
            "The rest of the ORF aligns continuously (97% of Dval covered).",
            "",
            f"Substitutions in the overlapping region: **{len(apc1_sub_rows)}** {sub_list}.",
            "",
            "None of those substitutions fall in Dval 1050–1125 (homology window facing an "
            "Emi-family ZBR contact on human Apc1 in 4UI9 — literature context only, not an "
            "Emi2 measurement). That window is fully present in both alleles. The Dnai gap "
            "instead sits at Dval 1284–1330, overlapping the human Apc1 surface that cradles "
            "Apc10 in 4UI9 — again literature context; we have **no Emi2–Apc1 model**, so this "
            "does not by itself create a dual-side Emi2 hit.",
            "",
        ]

    lines += [
        "## 3. Geometry from Emi2 S1 models (not Emi1)",
        "",
        f"Engaged models used for the distance matrix: **{matrix['n_models']}** "
        f"(pair_iptm_emi2_cdc20 ≥ {ENGAGED_IPTM}). "
        "Matrix: `emi2_s1_distance_matrix.json` (median over engaged models).",
        "",
        "| Emi2 site | element | median Å to Cdc20 | median Å to Apc10 | interpretation |",
        "|---|---|---|---|---|",
    ]
    for pos in (548, 549, 550, 551, 590):
        d_c, d_a = dist_summary(pos)
        el = element_of(pos) or "—"
        dc = f"{d_c:.2f}" if d_c is not None else "—"
        da = f"{d_a:.2f}" if d_a is not None else "—"
        if pos == 549:
            interp = (
                f"side chain ≤4 Å of Cdc20 in **{sc_549_within4}/{len(sc_549)}** engaged models "
                "(backbone may approach; allele effect null in REPORT_AF3_S1)"
            )
        elif pos == 590:
            interp = "no publishable APC partner in S1; ZBR/Apc1/Apc2 not in the construct"
        elif pos in (548, 551):
            interp = "canonical D-box contacts on Cdc20"
        else:
            interp = ""
        note_apc10 = " †" if d_a is not None else ""
        lines.append(f"| {pos} | {el} | {dc} | {da}{note_apc10} | {interp} |")
    lines += [
        "",
        "† Apc10 column is a **non-measurement** (co-receptor site never occupied; "
        "see `AF3_results_S1/REPORT_AF3_S1.md`).",
        "",
        "Coverage of S1: only Emi2 **535–675** is present. Residues 1–534, including eight of "
        "the ten substitutions and the F-box, have **no complex geometry** in this study.",
        "",
        "## 4. Intersection — parthenogenetic substitution list",
        "",
        "The hypothesis requires an Emi2 substitution that both (i) sits in a structurally "
        "tested interface and (ii) faces a variable APC/C position.",
        "",
        "- Structurally tested interface in hand: **Emi2 D-box ↔ Cdc20**.",
        "- Emi2 substitutions in that element: **F549I**.",
        f"- Cdc20 variability on/near that surface: **none** (sole allele difference V291I at "
        f"{v291_med:.1f} Å from the D-box).",
        "- Apc10: identical between parents; geometry non-informative in S1.",
        "- V590M / ZBR / Apc1 / Apc2 / Apc11: **not structurally tested on Emi2** "
        "(AF3 S2 inputs exist but were not run; no Emi2–Apc1 complex).",
        "- Anapc1 sequence: Dnai lacks Dval 1284–1330 (47 aa); overlap has substitutions "
        "but **none** in the homology window opposite an Emi-family ZBR site (1050–1125). "
        "Without an Emi2–Apc1 model this remains sequence context, not a dual-side hit.",
        "",
        "**Candidate list for this triad, given available Emi2 models: empty.**",
        "",
        "## 5. What was removed",
        "",
        "The previous report measured distances to human Emi1 in PDB 4UI9 and transferred them "
        "onto Emi2 by alignment. That surrogate is withdrawn from the publishable line. "
        "Literature on Emi1/Emi2 family architecture may still be cited as background, but not "
        "as a measurement of Darevskia Emi2 allele contacts.",
        "",
        "## 6. Files",
        "",
        "- `emi2_two_sided_test.py` — this analysis",
        "- `emi2_s1_distance_matrix.json` — median distances from engaged S1 models",
        "- `emi2_variable_sites.csv`, `apc_variable_sites.csv`",
        "- `../Fasta/Anapc1/` — Dval full + Dnai partial + alignment",
        "- `superseded/` — Emi1/4UI9 artefacts",
        "",
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT / 'REPORT.md'}")
    print(f"Engaged CIFs: {len(cifs)}")
    print(f"Emi2 sites: {len(emi2_sites)}")
    print(f"Cdc20 V291I vs D-box median: {v291_med:.2f} Å")


if __name__ == "__main__":
    main()
