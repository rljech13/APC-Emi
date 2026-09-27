#!/usr/bin/env python
"""Emit ChimeraX 'view matrix camera' strings looking down the Cdc20 -> D-box axis."""
import numpy as np, gemmi, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
p = ROOT/"s1_dbox_emi2dval_apcdval_seed1/fold_s1_dbox_emi2dval_apcdval_seed1_model_0.cif"
st = gemmi.read_structure(str(p)); st.setup_entities()
B = np.array([[a.pos.x,a.pos.y,a.pos.z] for r in st[0]["B"] for a in r if a.name=="CA"])
A = st[0]["A"]
db = np.array([[a.pos.x,a.pos.y,a.pos.z] for r in A if 14<=r.seqid.num<=17 for a in r])
c, d = B.mean(0), db.mean(0)
z = d - c; z /= np.linalg.norm(z)
# up vector: along the D-box chain direction (14 CA -> 17 CA), orthogonalised
v = np.array([A[16]["CA"][0].pos.x,A[16]["CA"][0].pos.y,A[16]["CA"][0].pos.z]) - \
    np.array([A[13]["CA"][0].pos.x,A[13]["CA"][0].pos.y,A[13]["CA"][0].pos.z])
y = v - z*np.dot(v,z); y /= np.linalg.norm(y)
x = np.cross(y,z)
def cam(dist, extra_roll=0.0):
    R = np.column_stack([x,y,z])
    if extra_roll:
        th=np.radians(extra_roll); cth,sth=np.cos(th),np.sin(th)
        Rz=np.array([[cth,-sth,0],[sth,cth,0],[0,0,1]])
        R = R @ Rz
    pos = d + R[:,2]*dist
    M = np.column_stack([R, pos])
    return "view matrix camera " + ",".join(f"{v:.5f}" for v in M.flatten())
for dist,name in ((110,"near"),(190,"wide")):
    print(f"# {name}\n{cam(dist)}")
print("# side view (rotate 75 deg about the D-box axis y)")
th=np.radians(75); cth,sth=np.cos(th),np.sin(th)
Ry=np.array([[cth,0,sth],[0,1,0],[-sth,0,cth]])
R=np.column_stack([x,y,z])@Ry
pos=d+R[:,2]*190
print("view matrix camera " + ",".join(f"{v:.5f}" for v in np.column_stack([R,pos]).flatten()))
