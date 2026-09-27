#!/usr/bin/env python
"""Generate the D-box zoom + 5G04-validation ChimeraX scripts with explicit cameras."""
import numpy as np, gemmi
from pathlib import Path
ROOT = Path(__file__).resolve().parent
FIG = str(Path(__file__).resolve().parents[1] / "Figures")
DVAL = ROOT/"s1_dbox_emi2dval_apcdval_seed1/fold_s1_dbox_emi2dval_apcdval_seed1_model_0.cif"
DNAI = ROOT/"s1_dbox_emi2dnai_apcdnai_seed1/fold_s1_dbox_emi2dnai_apcdnai_seed1_model_0.cif"

st = gemmi.read_structure(str(DVAL)); st.setup_entities()
A = st[0]["A"]
B = np.array([[a.pos.x,a.pos.y,a.pos.z] for r in st[0]["B"] for a in r if a.name=="CA"])
db = np.array([[a.pos.x,a.pos.y,a.pos.z] for r in A if 14<=r.seqid.num<=17 for a in r])
c, d = B.mean(0), db.mean(0)
z = d-c; z/=np.linalg.norm(z)
v = np.array([A[16]["CA"][0].pos.x,A[16]["CA"][0].pos.y,A[16]["CA"][0].pos.z]) - \
    np.array([A[13]["CA"][0].pos.x,A[13]["CA"][0].pos.y,A[13]["CA"][0].pos.z])
y = v-z*np.dot(v,z); y/=np.linalg.norm(y); x=np.cross(y,z)
R0 = np.column_stack([x,y,z])

def cam(dist, tiltx=0.0, roll=0.0, centre=d):
    R = R0.copy()
    if roll:
        t=np.radians(roll); R = R@np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])
    if tiltx:
        t=np.radians(tiltx); R = R@np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
    return "view matrix camera " + ",".join(
        f"{q:.5f}" for q in np.column_stack([R, centre + R[:,2]*dist]).flatten())

TOP  = cam(80)
TILT = cam(80, tiltx=-30)
GROOVE = "131,132,133,135,156,158,163,165,172,422,424,425"

common = f"""set bgColor white
windowsize 1500 1300
hide atoms
hide cartoons
hide surfaces
clip off
set subdivision 8
lighting soft
lighting shadows false
graphics silhouettes true width 1.6
color name cval  rgb(13%,30%,70%)
color name cnai  rgb(88%,42%,8%)
color name cR548 rgb(10%,45%,85%)
color name c549  rgb(85%,10%,45%)
color name cL551 rgb(95%,70%,5%)
color name chsl1 rgb(30%,30%,33%)
color name cgroove rgb(80%,82%,86%)
color name cgroove2 rgb(58%,66%,78%)
"""

def panel(model, colr, tag):
    return f"""
surface {model}/B
color {model}/B cgroove
color {model}/B:{GROOVE} cgroove2
transparency {model}/B 0
show {model}/A:11-21 cartoons
cartoon style protein modeHelix tube sides 24
cartoon style coil thickness 0.5
color {model}/A {colr}
show {model}/A:14,15,16,17 atoms
style {model}/A:14,15,16,17 stick
size {model}/A:14,15,16,17 stickRadius 0.19
color {model}/A:14 cR548
color {model}/A:15 c549
color {model}/A:17 cL551
color {model}/A byhetero
cartoon style {model}/A arrows false modeHelix tube sides 24
label {model}/A:14 residues text "R548" height 2.3 color rgb(10%,45%,85%) offset 2.0,-2.6,4
label {model}/A:15 residues text "{{0.name[0]}}{{0.name[1:].lower()}}549" height 2.3 color rgb(85%,10%,45%) offset -6.5,1.6,4
label {model}/A:17 residues text "L551" height 2.3 color rgb(70%,50%,0%) offset 0.8,2.4,4
"""

s1 = common + f"open {DVAL}\nopen {DNAI}\nmatchmaker #2/B to #1/B\n" + common
s1 += panel("#1", "cval", "Dval") + panel("#2", "cnai", "Dnai")
s1 += f"""
hide #2 models
{TOP}
save "{FIG}/AF3_S1_dbox_zoom_Dval.png" width 2200 height 1900 supersample 4
{TILT}
save "{FIG}/AF3_S1_dbox_zoom_Dval_tilt.png" width 2200 height 1900 supersample 4
hide #1 models
show #2 models
{TOP}
save "{FIG}/AF3_S1_dbox_zoom_Dnai.png" width 2200 height 1900 supersample 4
{TILT}
save "{FIG}/AF3_S1_dbox_zoom_Dnai_tilt.png" width 2200 height 1900 supersample 4
exit
"""
(ROOT/"fig_zoom.cxc").write_text(s1)

s2 = common + f"open {DVAL}\nopen {DNAI}\nopen {ROOT}/ref/5g04.cif\n"
s2 += "matchmaker #2/B to #1/B\nmatchmaker #3/R to #1/B\n" + common
s2 += f"""
hide #3 models
surface #1/B
color #1/B cgroove
color #1/B:{GROOVE} cgroove2
transparency #1/B 60
show #1/A:11-21 cartoons
show #2/A:11-21 cartoons
cartoon style protein modeHelix tube sides 24
cartoon style coil thickness 0.5
color #1/A cval
color #2/A cnai
show #1,2/A:14,15,16,17 atoms
style #1,2/A:14,15,16,17 stick
size #1,2/A:14,15,16,17 stickRadius 0.19
show #3 models
hide #3 cartoons
hide #3 atoms
show #3/S:826-834 cartoons
show #3/S:828,829,830,831 atoms
style #3/S:828,829,830,831 stick
size #3/S:828,829,830,831 stickRadius 0.19
color #3/S chsl1
color #1,2,3 byhetero
label #1/A:14,17 residues text "{{0.name}}{{0.number}}" height 2.0 color rgb(13%,30%,70%) offset 1.2,1.2,3
label #3/S:828,831 residues text "Hsl1 {{0.name}}{{0.number}}" height 2.0 color rgb(30%,30%,33%) offset 1.2,-2.5,3
{TOP}
save "{FIG}/AF3_S1_dbox_vs_5G04_Hsl1.png" width 2200 height 1900 supersample 4
{TILT}
save "{FIG}/AF3_S1_dbox_vs_5G04_Hsl1_tilt.png" width 2200 height 1900 supersample 4
exit
"""
(ROOT/"fig_5g04.cxc").write_text(s2)
print("wrote fig_zoom.cxc and fig_5g04.cxc")
