#!/usr/bin/env python
"""Task 4: is the D-box canonically engaged? Contacts + 5G04 validation + SASA."""
import json, re
from pathlib import Path
import numpy as np, pandas as pd, gemmi
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parent
REF = ROOT / "ref" / "5g04.cif"
AA3 = gemmi.ResidueInfo

def load(p):
    st = gemmi.read_structure(str(p)); st.setup_entities(); st.remove_alternative_conformations()
    return st

def seq_of(chain):
    return gemmi.one_letter_code([r.name for r in chain]).upper()

def heavy(res):
    return [a for a in res if a.element != gemmi.Element("H")]

def min_dist(r1, r2):
    return min(a.pos.dist(b.pos) for a in heavy(r1) for b in heavy(r2))

def contacts(st, ch1, ch2, cutoff=4.0):
    """residues of ch1 with any heavy atom within cutoff of ch2"""
    c1 = st[0][ch1]; c2 = st[0][ch2]
    ns = gemmi.NeighborSearch(st[0], st.cell, cutoff + 1).populate()
    out = {}
    for r1 in c1:
        best = 1e9; partners = set()
        for a in heavy(r1):
            for m in ns.find_atoms(a.pos, '\0', radius=cutoff):
                cra = m.to_cra(st[0])
                if cra.chain.name != ch2: continue
                if cra.atom.element == gemmi.Element("H"): continue
                d = a.pos.dist(cra.atom.pos)
                if d <= cutoff:
                    best = min(best, d); partners.add((cra.residue.seqid.num, cra.residue.name))
        if partners:
            out[(r1.seqid.num, r1.name)] = (round(best, 2), sorted(partners))
    return out

# ---------- ranking ----------
df = pd.read_csv(ROOT / "AF3_S1_confidence_table.csv")
top = (df.sort_values(["ranking_score", "pair_iptm_emi2_cdc20"], ascending=False)
         .groupby(["emi2_allele", "apc_source"]).head(1)
         .sort_values(["emi2_allele", "apc_source"]))
print("=== top-ranked model per combination ===")
print(top[["emi2_allele","apc_source","seed","model_index","ranking_score","iptm",
           "pair_iptm_emi2_cdc20","pair_iptm_emi2_apc10","pair_iptm_cdc20_apc10"]].to_string(index=False))

def path_of(row):
    j = f"s1_dbox_emi2{row.emi2_allele}_apc{row.apc_source}_seed{row.seed}"
    return ROOT / j / f"fold_{j}_model_{row.model_index}.cif"

# ---------- reference: 5G04 ----------
ref = load(REF)
refCdc20 = ref[0]["R"]; refHsl1 = ref[0]["S"]; refApc10 = ref[0]["L"]
print(f"\n5G04: Cdc20 chain R ({len(refCdc20)} res), Hsl1 D-box chain S "
      f"({[(r.seqid.num, r.name) for r in refHsl1]}), Apc10 chain L ({len(refApc10)} res)")

# canonical D-box receptor: Cdc20 residues within 4A of the Hsl1 peptide
ref_recv = contacts(ref, "R", "S", 4.0)
print(f"\n5G04 canonical D-box receptor = Cdc20 residues within 4 A of Hsl1 828-837 (n={len(ref_recv)}):")
print("  " + ", ".join(f"{n}{gemmi.find_tabulated_residue(nm).one_letter_code.upper()}"
                       for (n, nm) in sorted(ref_recv)))
ref_hsl1_recv = contacts(ref, "S", "R", 4.0)
print("  Hsl1 residues contacting Cdc20:", sorted(ref_hsl1_recv))
print("  Hsl1 residues contacting Apc10:", sorted(contacts(ref, "S", "L", 4.0)))
print("  Apc10(L) residues within 4 A of Cdc20(R):", len(contacts(ref, "L", "R", 4.0)))

# ---------- alignment our Cdc20 -> human Cdc20 ----------
aligner = Align.PairwiseAligner(scoring="blastp", mode="global",
                                open_gap_score=-11, extend_gap_score=-1)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

def map_residues(qchain, tchain):
    """returns dict qseqid -> tseqid from a global alignment"""
    qs, ts = seq_of(qchain), seq_of(tchain)
    aln = aligner.align(qs, ts)[0]
    qi = [r.seqid.num for r in qchain]; ti = [r.seqid.num for r in tchain]
    mp = {}
    for (qa, qb), (ta, tb) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(qb - qa):
            mp[qi[qa + k]] = ti[ta + k]
    ident = sum(1 for q, t in mp.items() if qs[qi.index(q)] == ts[ti.index(t)])
    return mp, aln.score, ident / max(1, len(mp))

def superpose(st, refst, mp, qch="B", tch="R", cycles=5):
    q = {r.seqid.num: r for r in st[0][qch]}
    t = {r.seqid.num: r for r in refst[0][tch]}
    pairs = [(q[a]["CA"][0].pos, t[b]["CA"][0].pos)
             for a, b in mp.items() if a in q and b in t and q[a].find_atom("CA", "*") and t[b].find_atom("CA", "*")]
    keep = list(range(len(pairs)))
    for _ in range(cycles):
        P = np.array([[pairs[i][0].x, pairs[i][0].y, pairs[i][0].z] for i in keep])
        Q = np.array([[pairs[i][1].x, pairs[i][1].y, pairs[i][1].z] for i in keep])
        pc, qc = P.mean(0), Q.mean(0)
        H = (P - pc).T @ (Q - qc)
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1, 1, d]) @ U.T
        allP = np.array([[p[0].x, p[0].y, p[0].z] for p in pairs])
        allQ = np.array([[p[1].x, p[1].y, p[1].z] for p in pairs])
        dev = np.linalg.norm((allP - pc) @ R.T + qc - allQ, axis=1)
        rms = np.sqrt((dev[keep] ** 2).mean())
        newkeep = [i for i in range(len(pairs)) if dev[i] < max(2.0, 2 * rms)]
        if newkeep == keep or len(newkeep) < 50: break
        keep = newkeep
    return R, pc, qc, rms, len(keep), len(pairs)

def apply_tf(pos, R, pc, qc):
    v = np.array([pos.x, pos.y, pos.z])
    w = (v - pc) @ R.T + qc
    return gemmi.Position(*w)

def sc_centroid(res, names):
    at = [a for a in res if a.name in names]
    if not at: return None
    return np.mean([[a.pos.x, a.pos.y, a.pos.z] for a in at], axis=0)

GUAN = {"CZ", "NH1", "NH2", "NE"}
LEUSC = {"CG", "CD1", "CD2"}

results = []
detail = {}
for row in top.itertuples():
    tag = f"emi2{row.emi2_allele}/apc{row.apc_source} seed{row.seed} model{row.model_index}"
    p = path_of(row)
    st = load(p)
    A, B, C = st[0]["A"], st[0]["B"], st[0]["C"]
    assert seq_of(A)[13:17] in ("RFAL", "RIAL"), seq_of(A)[13:17]
    print("\n" + "=" * 100)
    print(f"### {tag}   [{p.name}]   CIF 14-17 = {seq_of(A)[13:17]}")

    cAB = contacts(st, "A", "B"); cAC = contacts(st, "A", "C")
    cBC = contacts(st, "B", "C")
    def fmt(d):
        return ", ".join(f"{n}{gemmi.find_tabulated_residue(nm).one_letter_code.upper()}"
                         f"({n+534})" for (n, nm) in sorted(d))
    print(f"\n Emi2 residues within 4 A of Cdc20 (chain B): n={len(cAB)}")
    print("   " + (fmt(cAB) if cAB else "NONE"))
    print(f" Emi2 residues within 4 A of Apc10 (chain C): n={len(cAC)}")
    print("   " + (fmt(cAC) if cAC else "NONE"))
    print(f" Cdc20-Apc10 interface: n={len(cBC)} Cdc20 residues within 4 A of Apc10")

    # which Cdc20 residues does Emi2 touch, and do they include the canonical receptor?
    mp, score, ident = map_residues(B, refCdc20)          # our Cdc20 num -> human num
    inv = {v: k for k, v in mp.items()}
    print(f" [aln] our Cdc20 vs human Cdc20(5G04 R): {len(mp)} aligned, {ident*100:.0f}% identity")
    canon_ours = {inv[h] for (h, _) in ref_recv if h in inv}
    touched = set()
    for v in cAB.values():
        touched |= {n for n, _ in v[1]}
    print(f" canonical D-box-receptor residues of Cdc20 mapped into our numbering: n={len(canon_ours)}")
    ov = touched & canon_ours
    print(f" Cdc20 residues contacted by Emi2: n={len(touched)}; overlap with canonical receptor: "
          f"n={len(ov)} ({100*len(ov)/max(1,len(canon_ours)):.0f}% of receptor) -> {sorted(ov)}")
    # specifically who contacts R548(14) and L551(17)
    for cif, full in ((14, 548), (17, 551), (15, 549), (16, 550)):
        v = cAB.get((cif, A[cif-1].name))
        print(f"   CIF {cif} ({A[cif-1].name}{full}): "
              + (f"min d={v[0]} A to Cdc20 {[f'{n}{nm}' for n,nm in v[1]]}" if v else "no Cdc20 contact <4 A"))

    # ---- superpose our Cdc20 onto 5G04 Cdc20, then compare D-box placement
    R, pc, qc, rms, nk, npair = superpose(st, ref, mp)
    print(f" [superpose] our Cdc20 -> 5G04 Cdc20: CA RMSD {rms:.2f} A over {nk}/{npair} core residues")
    hsl1 = {r.seqid.num: r for r in refHsl1}
    ourA = {r.seqid.num: r for r in A}
    def dist_sc(ourres, refres, names):
        a = sc_centroid(ourres, names); b = sc_centroid(refres, names)
        if a is None or b is None: return None
        a = (a - pc) @ R.T + qc
        return float(np.linalg.norm(a - b))
    dR = dist_sc(ourA[14], hsl1[828], GUAN)
    dL = dist_sc(ourA[17], hsl1[831], LEUSC)
    dRca = dist_sc(ourA[14], hsl1[828], {"CA"})
    dLca = dist_sc(ourA[17], hsl1[831], {"CA"})
    print(f" >> D-box placement vs Hsl1 (after Cdc20 superposition):")
    print(f"      R548(CIF14) guanidinium centroid  displacement = {dR:.1f} A   (CA: {dRca:.1f} A)")
    print(f"      L551(CIF17) side-chain centroid   displacement = {dL:.1f} A   (CA: {dLca:.1f} A)")
    # backbone RMSD of the 4 D-box residues onto Hsl1 828-831
    bb = []
    for o, h in ((14, 828), (15, 829), (16, 830), (17, 831)):
        for an in ("N", "CA", "C", "O"):
            if ourA[o].find_atom(an, "*") and hsl1[h].find_atom(an, "*"):
                a = np.array([ourA[o][an][0].pos.x, ourA[o][an][0].pos.y, ourA[o][an][0].pos.z])
                a = (a - pc) @ R.T + qc
                b = hsl1[h][an][0].pos
                bb.append(np.linalg.norm(a - np.array([b.x, b.y, b.z])))
    dbox_rmsd = float(np.sqrt(np.mean(np.square(bb))))
    print(f"      D-box 548-551 backbone RMSD to Hsl1 828-831 (in-place, no local fit) = {dbox_rmsd:.1f} A")

    # where is our Apc10 relative to 5G04 Apc10?
    mpC, _, identC = map_residues(C, refApc10)
    ourC = {r.seqid.num: r for r in C}; refL = {r.seqid.num: r for r in refApc10}
    dd = []
    for a, b in mpC.items():
        if a in ourC and b in refL and ourC[a].find_atom("CA", "*") and refL[b].find_atom("CA", "*"):
            pa = ourC[a]["CA"][0].pos
            pa = (np.array([pa.x, pa.y, pa.z]) - pc) @ R.T + qc
            pb = refL[b]["CA"][0].pos
            dd.append(np.linalg.norm(pa - np.array([pb.x, pb.y, pb.z])))
    print(f" Apc10 placement: after Cdc20 superposition, our Apc10 CA are {np.mean(dd):.0f} A "
          f"(median {np.median(dd):.0f}) from their 5G04 counterparts ({identC*100:.0f}% id, n={len(dd)})")

    results.append(dict(emi2=row.emi2_allele, apc=row.apc_source, seed=row.seed, model=row.model_index,
                        ranking=row.ranking_score, iptm_AB=row.pair_iptm_emi2_cdc20,
                        n_contact_B=len(cAB), n_contact_C=len(cAC),
                        cdc20_recv_overlap=len(ov), cdc20_recv_n=len(canon_ours),
                        sup_rmsd=round(rms,2), dR548=round(dR,1), dL551=round(dL,1),
                        dbox_bb_rmsd=round(dbox_rmsd,1), apc10_disp=round(float(np.mean(dd)),0)))
    detail[tag] = dict(cAB=cAB, cAC=cAC, path=str(p))

print("\n\n=== TASK 4 SUMMARY ===")
rdf = pd.DataFrame(results)
print(rdf.to_string(index=False))
rdf.to_csv(ROOT / "AF3_S1_dbox_validation.csv", index=False)
json.dump({k: {"path": v["path"],
               "emi2_contacts_cdc20": {f"{n}{nm}": v["cAB"][(n,nm)] for (n,nm) in v["cAB"]},
               "emi2_contacts_apc10": {f"{n}{nm}": v["cAC"][(n,nm)] for (n,nm) in v["cAC"]}}
           for k, v in detail.items()}, open(ROOT / "AF3_S1_contacts.json", "w"), indent=1)
print("\nwrote AF3_S1_dbox_validation.csv, AF3_S1_contacts.json")
