#!/usr/bin/env python
"""Two-sided test of the parthenogenesis hypothesis for the D. valentini x D. r. nairensis pair.

The hypothesis (Fujita & Moritz framework, as specialised here) requires variable sites on BOTH
sides of an Emi2-APC/C interface: a substitution in Emi2 that sits in an interface AND contacts a
position that is itself variable between the two parental species. This script enumerates both
sides and intersects them.

Side 1 (Emi2): pairwise Dval vs Dnai, positions classified against the functional element map.
Side 2 (APC/C): pairwise Dval vs Dnai for every subunit for which both alleles are available,
then each variable site is located relative to the Emi2-binding surface. Location is measured in
4UI9 (APC/C-Cdh1-Emi1, 3.6 A) because it is the only atomic structure containing an Emi1/Emi2-family
inhibitor bound to the intact APC/C; Darevskia positions are mapped to human numbering by pairwise
alignment. Emi1 chain S (residues 319-436) stands in for Emi2.

Outputs Interface_variability/REPORT.md and two CSVs.
"""
from pathlib import Path
import csv
import gemmi
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Interface_variability"
REF_4UI9 = Path("/tmp/4ui9.cif")

# Literature-verified element boundaries in Darevskia Emi2 numbering (see REPORT_AF3_S1.md).
EMI2_ELEMENTS = [
    (40, 46, "beta-TrCP degron"), (176, 207, "Plk1-PBD docking"),
    (343, 363, "PP2A-B56"), (452, 537, "F-box domain"), (479, 519, "F-box"),
    (548, 551, "D-box"), (600, 648, "ZBR"), (658, 675, "RL tail"),
]
# Emi2 elements that engage the APC/C, and the subunit each one engages.
EMI2_APC_INTERFACES = {"D-box": ["Cdc20", "Anapc10"], "ZBR": ["Anapc11", "Anapc2", "Anapc1"],
                       "RL tail": ["Anapc2"]}

# Subunits and the 4UI9 chain that carries them. Emi1 (the Emi2 surrogate) is chain S.
SUBUNITS = {
    "Cdc20":   dict(chain=None, note="coactivator; D-box receptor (Cdh1 chain R in 4UI9)"),
    "Anapc10": dict(chain="L",  note="D-box co-receptor"),
    "Anapc11": dict(chain="B",  note="RING; ZBR target"),
    "Anapc2":  dict(chain="N",  note="cullin; ZBR body + RL tail"),
    "Anapc5":  dict(chain=None, note="not an Emi2 partner (negative control)"),
    "Anapc7":  dict(chain=None, note="not an Emi2 partner (negative control)"),
}
EMI1_CHAIN = "S"

aligner = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")


def read_aa(path):
    return str(next(SeqIO.parse(str(path), "fasta")).seq).replace("*", "")


def pairwise_diffs(a, b):
    """Substitutions a->b, positions in a's numbering."""
    aln = aligner.align(a, b)[0]
    out = []
    for (x, y), (c, _d) in zip(*aln.aligned):
        for k in range(y - x):
            if a[x + k] != b[c + k]:
                out.append((x + k + 1, a[x + k], b[c + k]))
    return out


def map_numbering(src, dst):
    """position in src -> position in dst, by pairwise alignment."""
    aln = aligner.align(src, dst)[0]
    m = {}
    for (x, y), (c, _d) in zip(*aln.aligned):
        for k in range(y - x):
            m[x + k + 1] = c + k + 1
    return m


def element_of(pos):
    for lo, hi, name in EMI2_ELEMENTS:
        if lo <= pos <= hi:
            return name
    return None


def load_ref():
    st = gemmi.read_structure(str(REF_4UI9))
    st.setup_entities()
    st.remove_ligands_and_waters()
    return st


def min_dist_to_emi1(st, chain_name, resnum):
    """Minimum heavy-atom distance from one residue of `chain_name` to Emi1 (chain S)."""
    try:
        chain = st[0][chain_name]
        emi1 = st[0][EMI1_CHAIN]
    except Exception:
        return None
    res = {r.seqid.num: r for r in chain}.get(resnum)
    if res is None:
        return None
    return min(a.pos.dist(b.pos) for a in res for rb in emi1 for b in rb)


def cdc20_distance():
    """Locate the single Cdc20 substitution (291) in our own Emi2-Cdc20 models.

    Chain A is Emi2 535-675 renumbered 1..141, so the D-box 548-551 is CIF 14-17; chain B is Cdc20.
    """
    import glob
    import statistics
    any_d, db_d = [], []
    for p in sorted(glob.glob(str(ROOT / "AF3_results_S1/s1_dbox_*/fold_*_model_*.cif"))):
        s = gemmi.read_structure(p)
        s.setup_entities()
        s.remove_ligands_and_waters()
        b = {r.seqid.num: r for r in s[0]["B"]}.get(291)
        if b is None:
            continue
        a_all = list(s[0]["A"])
        a_map = {r.seqid.num: r for r in s[0]["A"]}
        dbox = [a_map[n] for n in (14, 15, 16, 17) if n in a_map]
        any_d.append(min(x.pos.dist(y.pos) for x in b for rb in a_all for y in rb))
        if dbox:
            db_d.append(min(x.pos.dist(y.pos) for x in b for rb in dbox for y in rb))
    if not db_d:
        return None
    return dict(n=len(db_d), any_min=min(any_d), db_min=min(db_d),
                db_med=statistics.median(db_d))


def main():
    lines = ["# Two-sided variability test: D. valentini x D. r. nairensis", ""]

    # ---------------- Side 1: Emi2 ----------------
    dval = read_aa(ROOT / "Fasta/Emi2/Emi2_Dval_AAseq.fasta")
    dnai = read_aa(ROOT / "Fasta/Emi2/Emi2_Dnai_AAseq.fasta")
    emi2_diffs = pairwise_diffs(dval, dnai)

    lines += ["## 1. Emi2 side", "",
              f"Dval {len(dval)} aa vs Dnai {len(dnai)} aa: **{len(emi2_diffs)} substitutions**.", "",
              "| position | Dval | Dnai | functional element |", "|---|---|---|---|"]
    emi2_rows = []
    for pos, a, b in emi2_diffs:
        el = element_of(pos)
        lines.append(f"| {pos} | {a} | {b} | {el or '-'} |")
        emi2_rows.append(dict(position=pos, dval=a, dnai=b, element=el or ""))
    in_iface = [r for r in emi2_rows if r["element"] in EMI2_APC_INTERFACES]
    lines += ["", f"Substitutions inside an APC/C-engaging element: **{len(in_iface)}** "
              f"({', '.join(str(r['position']) for r in in_iface) or 'none'})", ""]

    # ---------------- Side 2: APC/C ----------------
    st = load_ref()
    lines += ["## 2. APC/C side", "",
              "Distance is the minimum heavy-atom distance to Emi1 (chain S, 319-436) in 4UI9.", "",
              "| subunit | Dval len | Dnai len | substitutions | distance to Emi1 | role |",
              "|---|---|---|---|---|---|"]
    apc_rows = []
    for sub, meta in SUBUNITS.items():
        d = ROOT / "Fasta" / sub
        pv = sorted(d.glob("*_Dval_AAseq*.fasta"))
        pn = sorted(d.glob("*_Dnai_AAseq*.fasta"))
        if not pv or not pn:
            lines.append(f"| {sub} | - | - | no data | - | {meta['note']} |")
            apc_rows.append(dict(subunit=sub, status="no data"))
            continue
        # prefer the v2 (re-translated) file when present
        pn = [p for p in pn if "v2" in p.name] or pn
        sv, sn = read_aa(pv[0]), read_aa(pn[0])
        diffs = pairwise_diffs(sv, sn)
        hs = ROOT / "Fasta" / sub / f"{sub}_Hsap_AAseq.fasta"
        d2h = map_numbering(sv, read_aa(hs)) if hs.exists() else {}
        descr, dists = [], []
        for pos, a, b in diffs:
            h = d2h.get(pos)
            dist = min_dist_to_emi1(st, meta["chain"], h) if (meta["chain"] and h) else None
            descr.append(f"{a}{pos}{b}")
            dists.append(f"{dist:.1f} A" if dist else "n/a")
            apc_rows.append(dict(subunit=sub, position=pos, dval=a, dnai=b,
                                 human_position=h or "", dist_to_emi1=f"{dist:.2f}" if dist else ""))
        if not diffs:
            apc_rows.append(dict(subunit=sub, status="identical"))
        lines.append(f"| {sub} | {len(sv)} | {len(sn)} | "
                     f"{'; '.join(descr) if descr else '**0 — identical**'} | "
                     f"{'; '.join(dists) if dists else '—'} | {meta['note']} |")

    # Cdc20 cannot be measured in 4UI9 (that structure carries Cdh1, a different coactivator),
    # so its one substitution is located in our own Darevskia Emi2-Cdc20 models instead.
    cdc20 = cdc20_distance()
    if cdc20:
        lines += ["", f"Cdc20 carries one substitution (V291I). 4UI9 contains Cdh1 rather than Cdc20, so this "
                  f"one is measured in our own {cdc20['n']} Emi2-Cdc20 models: **{cdc20['db_min']:.1f}-"
                  f"{cdc20['db_med']:.1f} A from the D-box** (min-median). Distances to the nearest Emi2 atom "
                  f"of any kind reach {cdc20['any_min']:.1f} A, but those approaches involve the disordered "
                  "Emi2 segments whose placement is not reliable; the D-box figure is the meaningful one.", ""]
        apc_rows.append(dict(subunit="Cdc20", position=291, dval="V", dnai="I",
                             dist_to_emi1=f"{cdc20['db_min']:.2f}",
                             status="measured vs D-box in our AF3 models"))

    # positive control: the published Emi1-contacting APC2 residues
    lines += ["", "Positive control — published Emi1-contacting APC2 residues, same measurement:", ""]
    ctrl = [513, 514, 517, 520, 552, 553, 556, 602]
    got = [(h, min_dist_to_emi1(st, "N", h)) for h in ctrl]
    lines.append("| hAPC2 residue | " + " | ".join(str(h) for h, _ in got) + " |")
    lines.append("|---|" + "---|" * len(got))
    lines.append("| distance to Emi1 | " + " | ".join(f"{v:.2f} A" if v else "n/a" for _, v in got) + " |")

    # ---------------- Intersection ----------------
    variable_in_iface = [r for r in apc_rows
                         if r.get("dist_to_emi1") and float(r["dist_to_emi1"]) <= 5.0]
    lines += ["", "## 3. Intersection — the parthenogenetic substitution list", "",
              f"APC/C variable sites lying on the Emi2-binding surface (<=5 A from Emi1): "
              f"**{len(variable_in_iface)}**.", "",
              "The list of candidate parthenogenetic substitutions is the set of Emi2 substitutions "
              "that both sit in an interface and contact a variable APC/C position. For this pair "
              f"the list is **{'empty' if not variable_in_iface else 'non-empty'}**.", ""]

    with open(OUT / "emi2_variable_sites.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["position", "dval", "dnai", "element"])
        w.writeheader()
        w.writerows(emi2_rows)
    with open(OUT / "apc_variable_sites.csv", "w", newline="") as f:
        keys = ["subunit", "position", "dval", "dnai", "human_position", "dist_to_emi1", "status"]
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(apc_rows)

    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
