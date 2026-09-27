#!/usr/bin/env python
"""Task 4 (SASA of 549) + Task 5 (pose reproducibility over all 60 models)."""
import json, re, warnings
from pathlib import Path
import numpy as np, pandas as pd, gemmi
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent
ref = gemmi.read_structure(str(ROOT / "ref" / "5g04.cif")); ref.setup_entities()
refCdc20 = ref[0]["R"]; refHsl1 = {r.seqid.num: r for r in ref[0]["S"]}
refApc10 = ref[0]["L"]

aligner = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1)
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

def seq_of(ch): return gemmi.one_letter_code([r.name for r in ch]).upper()

def map_residues(qchain, tchain):
    qs, ts = seq_of(qchain), seq_of(tchain)
    aln = aligner.align(qs, ts)[0]
    qi = [r.seqid.num for r in qchain]; ti = [r.seqid.num for r in tchain]
    return {qi[qa+k]: ti[ta+k] for (qa, qb), (ta, tb) in zip(*aln.aligned) for k in range(qb-qa)}

MP = None   # our Cdc20 numbering is identical across all models -> compute once

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P-pc).T @ (Q-qc); U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1,1,d]) @ U.T
    return R, pc, qc

def superpose_cdc20(st, mp, cycles=5):
    q = {r.seqid.num: r for r in st[0]["B"]}; t = {r.seqid.num: r for r in refCdc20}
    pr = [(q[a]["CA"][0].pos, t[b]["CA"][0].pos) for a, b in mp.items()
          if a in q and b in t and q[a].find_atom("CA","*") and t[b].find_atom("CA","*")]
    P = np.array([[p.x,p.y,p.z] for p,_ in pr]); Q = np.array([[p.x,p.y,p.z] for _,p in pr])
    keep = np.arange(len(P))
    for _ in range(cycles):
        R, pc, qc = kabsch(P[keep], Q[keep])
        dev = np.linalg.norm((P-pc) @ R.T + qc - Q, axis=1)
        rms = np.sqrt((dev[keep]**2).mean())
        nk = np.where(dev < max(2.0, 2*rms))[0]
        if len(nk) == len(keep) and (nk == keep).all(): break
        if len(nk) < 50: break
        keep = nk
    return R, pc, qc, rms

def tf(X, R, pc, qc): return (np.asarray(X) - pc) @ R.T + qc

def heavy(res): return [a for a in res if a.element != gemmi.Element("H")]

def contact_set(st, ch1, ch2, cutoff=4.0):
    ns = gemmi.NeighborSearch(st[0], st.cell, cutoff+1).populate()
    out = set()
    for r1 in st[0][ch1]:
        for a in heavy(r1):
            for m in ns.find_atoms(a.pos, '\0', radius=cutoff):
                cra = m.to_cra(st[0])
                if cra.chain.name == ch2 and cra.atom.element != gemmi.Element("H") \
                   and a.pos.dist(cra.atom.pos) <= cutoff:
                    out.add(r1.seqid.num); break
            else: continue
            break
    return out

GUAN = {"CZ","NH1","NH2","NE"}; LEUSC = {"CG","CD1","CD2"}
def cent(res, names):
    at=[a for a in res if a.name in names]
    return np.mean([[a.pos.x,a.pos.y,a.pos.z] for a in at],0) if at else None

rows=[]
for job in sorted(ROOT.glob("s1_dbox_*")):
    m=re.match(r"s1_dbox_emi2(dval|dnai)_apc(dval|dnai)_seed(\d)", job.name)
    emi2, apc, seed = m.group(1), m.group(2), int(m.group(3))
    for f in sorted(job.glob("*_model_*.cif")):
        idx=int(re.search(r"_(\d)\.cif$", f.name).group(1))
        st=gemmi.read_structure(str(f)); st.setup_entities(); st.remove_alternative_conformations()
        if MP is None: MP = map_residues(st[0]["B"], refCdc20)
        R,pc,qc,rms = superpose_cdc20(st, MP)
        A={r.seqid.num:r for r in st[0]["A"]}
        # D-box fidelity
        bb=[]
        for o,h in ((14,828),(15,829),(16,830),(17,831)):
            for an in ("N","CA","C","O"):
                if A[o].find_atom(an,"*") and refHsl1[h].find_atom(an,"*"):
                    p=A[o][an][0].pos; b=refHsl1[h][an][0].pos
                    bb.append(np.linalg.norm(tf([p.x,p.y,p.z],R,pc,qc)-np.array([b.x,b.y,b.z])))
        dbox_rmsd=float(np.sqrt(np.mean(np.square(bb))))
        cR=cent(A[14],GUAN); cL=cent(A[17],LEUSC)
        dR=float(np.linalg.norm(tf(cR,R,pc,qc)-cent(refHsl1[828],GUAN)))
        dL=float(np.linalg.norm(tf(cL,R,pc,qc)-cent(refHsl1[831],LEUSC)))
        # Emi2 COM in the Cdc20 frame
        allA=np.array([[a.pos.x,a.pos.y,a.pos.z] for r in st[0]["A"] for a in heavy(r)])
        com=tf(allA,R,pc,qc).mean(0)
        # Apc10 COM in the Cdc20 frame; and its displacement vs 5G04 Apc10 COM
        allC=np.array([[a.pos.x,a.pos.y,a.pos.z] for r in st[0]["C"] for a in heavy(r)])
        comC=tf(allC,R,pc,qc).mean(0)
        refC=np.array([[a.pos.x,a.pos.y,a.pos.z] for r in refApc10 for a in heavy(r)]).mean(0)
        setB=contact_set(st,"A","B"); setC=contact_set(st,"A","C")
        rows.append(dict(job=job.name,emi2=emi2,apc=apc,seed=seed,model=idx,
            sup_rmsd=rms, dbox_bb_rmsd=dbox_rmsd, dR548=dR, dL551=dL,
            emi2_com_x=com[0],emi2_com_y=com[1],emi2_com_z=com[2],
            apc10_com_disp=float(np.linalg.norm(comC-refC)),
            n_iface_B=len(setB), n_iface_C=len(setC),
            iface_B=";".join(map(str,sorted(setB))), iface_C=";".join(map(str,sorted(setC)))))
        print(".",end="",flush=True)
print()
pdf=pd.DataFrame(rows)
pdf.to_csv(ROOT/"AF3_S1_pose_table.csv",index=False)
conf=pd.read_csv(ROOT/"AF3_S1_confidence_table.csv").rename(columns={"emi2_allele":"emi2","apc_source":"apc","model_index":"model"})
pdf=pdf.merge(conf[["emi2","apc","seed","model","ranking_score","pair_iptm_emi2_cdc20","pair_iptm_emi2_apc10"]],
              on=["emi2","apc","seed","model"])
pd.set_option("display.width",250,"display.max_columns",60)

print("\n########## TASK 5: POSE REPRODUCIBILITY (all 60 models) ##########")
print("\n--- D-box fidelity vs 5G04 Hsl1, per job (mean over 5 models) ---")
g=pdf.groupby(["emi2","apc","seed"]).agg(
    iptm_AB=("pair_iptm_emi2_cdc20","mean"),
    dbox_rmsd=("dbox_bb_rmsd","mean"), dbox_rmsd_sd=("dbox_bb_rmsd","std"),
    dR548=("dR548","mean"), dL551=("dL551","mean"),
    n_iface_B=("n_iface_B","mean"), n_iface_C=("n_iface_C","mean"),
    apc10_disp=("apc10_com_disp","mean")).round(2)
print(g.to_string())

print("\n--- D-box engaged? (criterion: backbone RMSD to Hsl1 D-box < 2.0 A) ---")
pdf["dbox_ok"]=pdf.dbox_bb_rmsd<2.0
print(pdf.groupby(["emi2","apc","seed"])["dbox_ok"].agg(["sum","count"]).to_string())
print(f"\nOVERALL: {pdf.dbox_ok.sum()}/60 models place the D-box canonically")
print("\ndbox_bb_rmsd by whether the job is high- or low-ipTM:")
pdf["mode"]=np.where(pdf.pair_iptm_emi2_cdc20>0.4,"high-ipTM (>0.4)","low-ipTM (<0.4)")
print(pdf.groupby("mode")[["dbox_bb_rmsd","dR548","dL551","n_iface_B","n_iface_C",
                           "apc10_com_disp","pair_iptm_emi2_cdc20"]].agg(["mean","std","min","max"]).round(2).to_string())

print("\n--- Emi2 centre-of-mass spread in the Cdc20 frame (A) ---")
for (e,a),gg in pdf.groupby(["emi2","apc"]):
    X=gg[["emi2_com_x","emi2_com_y","emi2_com_z"]].values
    d=np.linalg.norm(X-X.mean(0),axis=1)
    print(f"emi2_{e}/apc_{a}: RMS deviation of Emi2 COM = {np.sqrt((d**2).mean()):.1f} A, max {d.max():.1f} A (n={len(gg)})")
X=pdf[["emi2_com_x","emi2_com_y","emi2_com_z"]].values
print(f"ALL 60: RMS COM deviation = {np.sqrt((np.linalg.norm(X-X.mean(0),axis=1)**2).mean()):.1f} A")

print("\n--- interface-set Jaccard overlap (Emi2 residues contacting Cdc20) ---")
def jac(a,b):
    A=set(a.split(';')) if a else set(); B=set(b.split(';')) if b else set()
    return len(A&B)/max(1,len(A|B))
for label,col in (("Emi2-Cdc20","iface_B"),("Emi2-Apc10","iface_C")):
    vals=pdf[col].values
    J=np.array([[jac(vals[i],vals[j]) for j in range(60)] for i in range(60)])
    iu=np.triu_indices(60,1)
    # within-job vs between-job
    jobs=pdf.job.values
    same=np.array([[jobs[i]==jobs[j] for j in range(60)] for i in range(60)])
    print(f"{label}: within-job J = {J[iu][same[iu]].mean():.2f}, "
          f"between-job J = {J[iu][~same[iu]].mean():.2f}")
    # core: residues present in >=90% of the high-ipTM models
    hi=pdf[pdf["mode"].str.startswith("high")]
    from collections import Counter
    c=Counter(x for v in hi[col] for x in (v.split(';') if v else []))
    core=sorted((int(k) for k,n in c.items() if n>=0.9*len(hi)))
    print(f"   residues in >=90% of high-ipTM models ({len(hi)}): {[f'{k}({k+534})' for k in core]}")

print("\n--- Apc10 placement (displacement of Apc10 COM from its 5G04 position, Cdc20 frame) ---")
print(pdf.groupby(["emi2","apc"])["apc10_com_disp"].agg(["mean","std","min","max"]).round(1).to_string())

# ================= SASA of residue 549 =================
print("\n\n########## TASK 4c: SASA of residue 549 (CIF 15) ##########")
# Tien et al. 2013 theoretical max ASA
MAXASA=dict(ALA=129,ARG=274,ASN=195,ASP=193,CYS=167,GLU=223,GLN=225,GLY=104,HIS=224,
            ILE=197,LEU=201,LYS=236,MET=224,PHE=240,PRO=159,SER=155,THR=172,TRP=285,TYR=263,VAL=174)
sr=ShrakeRupley(n_points=960)
p=MMCIFParser(QUIET=True)
sas=[]
top=(pd.read_csv(ROOT/"AF3_S1_confidence_table.csv")
     .sort_values(["ranking_score","pair_iptm_emi2_cdc20"],ascending=False)
     .groupby(["emi2_allele","apc_source"]).head(1))
for row in top.itertuples():
    j=f"s1_dbox_emi2{row.emi2_allele}_apc{row.apc_source}_seed{row.seed}"
    f=ROOT/j/f"fold_{j}_model_{row.model_index}.cif"
    s=p.get_structure("x",str(f))[0]
    for cid in [c.id for c in s]:
        if cid not in ("A","B","C"): s.detach_child(cid)
    sr.compute(s,level="R")
    complexed={r.id[1]:(r.sasa,r.get_resname()) for r in s["A"]}
    import copy
    s2=p.get_structure("y",str(f))[0]
    for cid in [c.id for c in s2]:
        if cid!="A": s2.detach_child(cid)
    sr.compute(s2,level="R")
    alone={r.id[1]:r.sasa for r in s2["A"]}
    for cif,full in ((14,548),(15,549),(16,550),(17,551),(56,590)):
        sa,nm=complexed[cif]; fa=alone[cif]
        sas.append(dict(emi2=row.emi2_allele,apc=row.apc_source,seed=row.seed,model=row.model_index,
                        cif=cif,full=full,resname=nm,
                        sasa_complex=round(sa,1),sasa_free=round(fa,1),d_sasa=round(fa-sa,1),
                        rel_sasa_complex=round(100*sa/MAXASA[nm],1),
                        rel_sasa_free=round(100*fa/MAXASA[nm],1),
                        buried_pct=round(100*(fa-sa)/max(fa,1e-6),1)))
sdf=pd.DataFrame(sas); sdf.to_csv(ROOT/"AF3_S1_sasa.csv",index=False)
print(sdf.to_string(index=False))
print("\n(rel_sasa in %, ref = Tien et al. 2013 theoretical max ASA; "
      "sasa_free = Emi2 chain A alone with the SAME conformation)")
