#!/usr/bin/env python
"""Atom-level dissection of the Emi2 D-box residues against Cdc20, over every
D-box-engaged model rather than the four top-ranked ones.

For CIF residues 14 (R548), 15 (549 F/I), 17 (L551) and 56 (590 V/M) of chain A
we report, per model:
  * the minimum distance from any SIDE-CHAIN heavy atom to any Cdc20 heavy atom
    (side chain = everything except N, CA, C, O, OXT),
  * the minimum distance from any BACKBONE heavy atom to Cdc20, and the partner,
  * the full set of Cdc20 residues within CUT_CONTACT of the residue (any atom).

Engagement criterion is the one already used in pose_and_sasa.py: D-box backbone
RMSD to 5G04 Hsl1 828-831 < 2.0 A, read from AF3_S1_pose_table.csv.
"""
import re
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent

CUT_CONTACT = 4.0   # "in contact"
CUT_NEAR = 5.0      # "not even close" check
BACKBONE = {"N", "CA", "C", "O", "OXT"}
TARGETS = {14: 548, 15: 549, 17: 551, 56: 590}


def heavy_coords(res, which):
    """(names, Nx3 coords) of heavy atoms of a residue, filtered by group."""
    out_n, out_x = [], []
    for a in res:
        if a.element == gemmi.Element("H"):
            continue
        if which == "side" and a.name in BACKBONE:
            continue
        if which == "back" and a.name not in BACKBONE:
            continue
        out_n.append(a.name)
        out_x.append([a.pos.x, a.pos.y, a.pos.z])
    return out_n, (np.array(out_x) if out_x else np.zeros((0, 3)))


def chain_heavy(chain):
    """Flat arrays of heavy atoms of a chain: coords, residue numbers, names."""
    xyz, rnum, rname, aname = [], [], [], []
    for r in chain:
        for a in r:
            if a.element == gemmi.Element("H"):
                continue
            xyz.append([a.pos.x, a.pos.y, a.pos.z])
            rnum.append(r.seqid.num)
            rname.append(r.name)
            aname.append(a.name)
    return (np.array(xyz), np.array(rnum), np.array(rname, dtype=object),
            np.array(aname, dtype=object))


pose = pd.read_csv(ROOT / "AF3_S1_pose_table.csv")
pose["dbox_ok"] = pose.dbox_bb_rmsd < 2.0
engaged = pose[pose.dbox_ok]
print(f"engaged models: {len(engaged)}/{len(pose)}")

rows, pair_rows = [], []
for job in sorted(ROOT.glob("s1_dbox_*")):
    m = re.match(r"s1_dbox_emi2(dval|dnai)_apc(dval|dnai)_seed(\d)", job.name)
    emi2, apc, seed = m.group(1), m.group(2), int(m.group(3))
    for f in sorted(job.glob("*_model_*.cif")):
        idx = int(re.search(r"_(\d)\.cif$", f.name).group(1))
        sel = engaged[(engaged.job == job.name) & (engaged.model == idx)]
        if sel.empty:
            continue
        st = gemmi.read_structure(str(f))
        st.setup_entities()
        st.remove_alternative_conformations()
        A = {r.seqid.num: r for r in st[0]["A"]}
        Bx, Bnum, Bres, Batom = chain_heavy(st[0]["B"])

        for cif, full in TARGETS.items():
            res = A[cif]
            rec = dict(job=job.name, emi2=emi2, apc=apc, seed=seed, model=idx,
                       cif=cif, full=full, resname=res.name)
            for which, tag in (("side", "sc"), ("back", "bb"), ("all", "any")):
                names, X = heavy_coords(res, which)
                if len(X) == 0:                     # Gly has no side chain
                    rec[f"min_{tag}"] = np.nan
                    rec[f"min_{tag}_atom"] = ""
                    rec[f"min_{tag}_partner"] = ""
                    continue
                Dm = np.linalg.norm(X[:, None, :] - Bx[None, :, :], axis=2)
                i, j = np.unravel_index(np.argmin(Dm), Dm.shape)
                rec[f"min_{tag}"] = float(Dm[i, j])
                rec[f"min_{tag}_atom"] = names[i]
                rec[f"min_{tag}_partner"] = f"{Bres[j]}{Bnum[j]}:{Batom[j]}"
                if which == "all":
                    for cut, lbl in ((CUT_CONTACT, "4A"), (CUT_NEAR, "5A")):
                        hit = Dm.min(0) <= cut
                        part = sorted(set(zip(Bnum[hit], Bres[hit])))
                        rec[f"partners_{lbl}"] = ";".join(
                            f"{rn}{n}" for n, rn in part)
                        rec[f"n_partners_{lbl}"] = len(part)
                    # per-partner minimum distance, for the contact-detail table
                    for n in sorted(set(Bnum[Dm.min(0) <= CUT_NEAR])):
                        mask = Bnum == n
                        pair_rows.append(dict(
                            job=job.name, emi2=emi2, apc=apc, seed=seed,
                            model=idx, cif=cif, full=full, resname=res.name,
                            cdc20=f"{Bres[mask][0]}{n}",
                            dmin=float(Dm[:, mask].min())))
            rows.append(rec)
        print(".", end="", flush=True)
print()

df = pd.DataFrame(rows)
df.to_csv(ROOT / "AF3_S1_sidechain_contacts.csv", index=False)
pdf = pd.DataFrame(pair_rows)
pdf.to_csv(ROOT / "AF3_S1_sidechain_contacts_pairs.csv", index=False)

pd.set_option("display.width", 250, "display.max_columns", 60)
n = df.model.count() // len(TARGETS)
print(f"\n########## atom-level contacts, {n} D-box-engaged models ##########")
print(f"side chain = all heavy atoms except {sorted(BACKBONE)}")

for cif, full in TARGETS.items():
    d = df[df.cif == cif]
    print(f"\n===== CIF {cif} = Emi2 {full} "
          f"({'/'.join(sorted(d.resname.unique()))}) =====")
    for key, lbl in (("min_sc", "side chain -> Cdc20"),
                     ("min_bb", "backbone   -> Cdc20"),
                     ("min_any", "any atom   -> Cdc20")):
        s = d[key]
        print(f"  {lbl}: median {s.median():.2f} A, "
              f"range {s.min():.2f}-{s.max():.2f} A")
    print("  by allele (min side-chain distance, A):")
    for allele, g in d.groupby("emi2"):
        print(f"    Emi2^{allele} (n={len(g)}, {g.resname.iloc[0]}): "
              f"median {g.min_sc.median():.2f}, "
              f"range {g.min_sc.min():.2f}-{g.min_sc.max():.2f}; "
              f"backbone median {g.min_bb.median():.2f} "
              f"({g.min_bb.min():.2f}-{g.min_bb.max():.2f})")
    n_sc4 = int((d.min_sc <= CUT_CONTACT).sum())
    n_sc5 = int((d.min_sc <= CUT_NEAR).sum())
    print(f"  models with a side-chain atom within {CUT_CONTACT} A: {n_sc4}/{len(d)}")
    print(f"  models with a side-chain atom within {CUT_NEAR} A: {n_sc5}/{len(d)}")
    if n_sc5:
        print(d[d.min_sc <= CUT_NEAR][
            ["job", "model", "resname", "min_sc", "min_sc_atom",
             "min_sc_partner"]].to_string(index=False))
    print(f"  backbone closest partner (count over {len(d)} models):")
    print("    " + str(d.min_bb_partner.value_counts().to_dict()))
    print(f"  Cdc20 residues within {CUT_CONTACT} A (any atom of the residue), "
          "frequency:")
    from collections import Counter
    c = Counter(x for v in d.partners_4A for x in (v.split(";") if v else []))
    print("    " + ", ".join(f"{k} {v}/{len(d)}"
                             for k, v in sorted(c.items(), key=lambda t: -t[1])))
    print(f"  number of Cdc20 partners within {CUT_CONTACT} A: "
          f"median {d.n_partners_4A.median():.0f}, "
          f"range {d.n_partners_4A.min()}-{d.n_partners_4A.max()}")
    if len(pdf):
        q = pdf[(pdf.cif == cif) & (pdf.dmin <= CUT_CONTACT)]
        if len(q):
            print(f"  per-partner closest approach within {CUT_CONTACT} A:")
            print(q.groupby("cdc20").dmin.agg(["count", "median", "min", "max"])
                  .round(2).sort_values("count", ascending=False).to_string())
