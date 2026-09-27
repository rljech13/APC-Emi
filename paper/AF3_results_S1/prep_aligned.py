#!/usr/bin/env python
"""Put the Dval model, the Dnai model and a trimmed 5G04 (Cdc20 R + Hsl1 S + Apc10 L)
into ONE common frame: the 5G04 Cdc20 frame. The AF3 Cdc20 construct has a ~90-residue
disordered N-terminus, so the fit is restricted to the ordered WD40 core that aligns to
5G04 Cdc20 residues 73-499."""
import numpy as np, gemmi
from pathlib import Path
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"aligned"; OUT.mkdir(exist_ok=True)
MODELS = {"dval": ROOT/"s1_dbox_emi2dval_apcdval_seed1/fold_s1_dbox_emi2dval_apcdval_seed1_model_0.cif",
          "dnai": ROOT/"s1_dbox_emi2dnai_apcdnai_seed1/fold_s1_dbox_emi2dnai_apcdnai_seed1_model_0.cif"}

def load(p):
    s=gemmi.read_structure(str(p)); s.setup_entities(); s.remove_alternative_conformations(); return s
def seq_of(ch): return gemmi.one_letter_code([r.name for r in ch]).upper()

aligner=Align.PairwiseAligner(mode="global",open_gap_score=-11,extend_gap_score=-1)
aligner.substitution_matrix=substitution_matrices.load("BLOSUM62")
def resmap(q,t):
    aln=aligner.align(seq_of(q),seq_of(t))[0]
    qi=[r.seqid.num for r in q]; ti=[r.seqid.num for r in t]
    return {qi[a+k]:ti[c+k] for (a,b),(c,d) in zip(*aln.aligned) for k in range(b-a)}

def fit(mob,ref,mp,cycles=8):
    q={r.seqid.num:r for r in mob}; t={r.seqid.num:r for r in ref}
    P=[];Q=[]
    for a,b in mp.items():
        if a in q and b in t and q[a].find_atom("CA","*") and t[b].find_atom("CA","*"):
            p1=q[a]["CA"][0].pos; p2=t[b]["CA"][0].pos
            P.append([p1.x,p1.y,p1.z]); Q.append([p2.x,p2.y,p2.z])
    P=np.array(P);Q=np.array(Q);keep=np.arange(len(P))
    for _ in range(cycles):
        pc,qc=P[keep].mean(0),Q[keep].mean(0)
        U,S,Vt=np.linalg.svd((P[keep]-pc).T@(Q[keep]-qc))
        R=Vt.T@np.diag([1,1,np.sign(np.linalg.det(Vt.T@U.T))])@U.T
        dev=np.linalg.norm((P-pc)@R.T+qc-Q,axis=1); rms=np.sqrt((dev[keep]**2).mean())
        nk=np.where(dev<max(1.5,2*rms))[0]
        if len(nk)<80 or (len(nk)==len(keep) and (nk==keep).all()): break
        keep=nk
    return R,pc,qc,rms,len(keep),len(P)

def transform(st,R,pc,qc):
    for m in st:
        for ch in m:
            for r in ch:
                for a in r:
                    v=np.array([a.pos.x,a.pos.y,a.pos.z]); a.pos=gemmi.Position(*((v-pc)@R.T+qc))

ref=load(ROOT/"ref"/"5g04.cif")
new=gemmi.Structure(); new.name="5g04_trim"; new.spacegroup_hm="P 1"
mo=gemmi.Model("1")
for cn in ("R","S","L"):
    ch=gemmi.Chain(cn)
    for r in ref[0][cn]: ch.add_residue(r)
    mo.add_chain(ch)
new.add_model(mo); new.setup_entities()
new.make_mmcif_document().write_file(str(OUT/"5g04_trim.cif"))

hsl1={r.seqid.num:r for r in ref[0]["S"]}
saved={}
for tag,p in MODELS.items():
    st=load(p)
    mp=resmap(st[0]["B"], ref[0]["R"])
    R,pc,qc,rms,nk,ntot=fit(st[0]["B"],ref[0]["R"],mp)
    print(f"{tag}: Cdc20 -> 5G04 Cdc20  CA RMSD {rms:.2f} A over {nk}/{ntot} core residues")
    transform(st,R,pc,qc)
    st.setup_entities()
    st.make_mmcif_document().write_file(str(OUT/f"{tag}_aligned.cif"))
    A={r.seqid.num:r for r in st[0]["A"]}
    saved[tag]=A
    for o,h in ((14,828),(17,831)):
        print(f"   CIF{o} CA vs Hsl1 {h} CA: {A[o]['CA'][0].pos.dist(hsl1[h]['CA'][0].pos):.2f} A")
d=saved["dval"][14]["CA"][0].pos.dist(saved["dnai"][14]["CA"][0].pos)
print(f"Dval vs Dnai R548 CA after common-frame alignment: {d:.2f} A")

# camera, computed in the common frame from the Dval model
st=load(OUT/"dval_aligned.cif"); A=st[0]["A"]
B=np.array([[a.pos.x,a.pos.y,a.pos.z] for r in st[0]["B"] for a in r if a.name=="CA"
            if r.seqid.num>90])
db=np.array([[a.pos.x,a.pos.y,a.pos.z] for r in A if 14<=r.seqid.num<=17 for a in r])
c,dd=B.mean(0),db.mean(0); z=dd-c; z/=np.linalg.norm(z)
v=np.array([A[16]["CA"][0].pos.x,A[16]["CA"][0].pos.y,A[16]["CA"][0].pos.z])-\
  np.array([A[13]["CA"][0].pos.x,A[13]["CA"][0].pos.y,A[13]["CA"][0].pos.z])
y=v-z*np.dot(v,z); y/=np.linalg.norm(y); x=np.cross(y,z)
R0=np.column_stack([x,y,z])
def cam(dist,tiltx=0.0,roll=0.0,centre=dd):
    R=R0.copy()
    if roll:
        t=np.radians(roll); R=R@np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])
    if tiltx:
        t=np.radians(tiltx); R=R@np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
    return "view matrix camera "+",".join(f"{q:.5f}" for q in np.column_stack([R,centre+R[:,2]*dist]).flatten())
(OUT/"cameras.txt").write_text("\n".join([
    "ZOOM_TOP  "+cam(80), "ZOOM_TILT "+cam(80,tiltx=-30), "ZOOM_SIDE "+cam(80,tiltx=-55,roll=100),
    "WIDE      "+cam(340), "WIDE_SIDE "+cam(340,tiltx=-75)]))
print((OUT/"cameras.txt").read_text())
