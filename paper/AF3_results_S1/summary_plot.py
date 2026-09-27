#!/usr/bin/env python
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path(__file__).resolve().parent
FIG = Path(__file__).resolve().parents[1] / "Figures"
pose = pd.read_csv(ROOT/"AF3_S1_pose_table.csv")
conf = pd.read_csv(ROOT/"AF3_S1_confidence_table.csv").rename(
    columns={"emi2_allele":"emi2","apc_source":"apc","model_index":"model"})
df = pose.merge(conf, on=["emi2","apc","seed","model"])
df["engaged"] = df.dbox_bb_rmsd < 2.0
CV, CN = "#2150b3", "#e06b14"
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

# A: per-job ipTM, bimodality
jobs = df.groupby(["emi2","apc","seed"])
labels, xs, ys, cs = [], [], [], []
for i,(k,g) in enumerate(sorted(jobs, key=lambda kv:(kv[0][0],kv[0][1],kv[0][2]))):
    labels.append(f"{k[1]}\ns{k[2]}")
    xs += [i]*len(g); ys += list(g.pair_iptm_emi2_cdc20)
    cs += [CV if k[0]=="dval" else CN]*len(g)
ax[0].scatter(np.array(xs)+np.random.RandomState(0).uniform(-.13,.13,len(xs)), ys, c=cs, s=34,
              edgecolor="k", linewidth=.4, zorder=3)
ax[0].axhline(0.4, ls="--", c="grey", lw=1)
ax[0].set_xticks(range(12)); ax[0].set_xticklabels(labels, fontsize=7.5)
ax[0].set_ylabel("Emi2–Cdc20 chain-pair ipTM"); ax[0].set_ylim(0,0.7)
ax[0].set_title("A  Outcome is bimodal and set by the seed\n(blue = Emi2$^{Dval}$, orange = Emi2$^{Dnai}$;"
                " x-label = APC source)", fontsize=9.5)
ax[0].text(.02,.52,"D-box docked", transform=ax[0].transAxes, fontsize=8, color="grey")
ax[0].text(.02,.16,"D-box not docked", transform=ax[0].transAxes, fontsize=8, color="grey")
ax[0].grid(alpha=.25, zorder=0)

# B: ipTM vs D-box RMSD
for al,c,lab in (("dval",CV,"Emi2$^{Dval}$"),("dnai",CN,"Emi2$^{Dnai}$")):
    g=df[df.emi2==al]
    ax[1].scatter(g.dbox_bb_rmsd, g.pair_iptm_emi2_cdc20, c=c, s=34, edgecolor="k",
                  linewidth=.4, label=lab, alpha=.85)
ax[1].set_xscale("log"); ax[1].axvline(2.0, ls="--", c="grey", lw=1)
ax[1].set_xlabel("D-box backbone RMSD to 5G04 Hsl1 R828-L831 (Å, log)")
ax[1].set_ylabel("Emi2–Cdc20 chain-pair ipTM")
ax[1].set_title("B  ipTM reports D-box docking, nothing finer\n45/60 models: RMSD 0.46–0.99 Å", fontsize=9.5)
ax[1].legend(fontsize=8.5); ax[1].grid(alpha=.25)

# C: signal vs noise floor
e = df[df.engaged]
sig = e[e.emi2=="dval"].pair_iptm_emi2_cdc20.mean() - e[e.emi2=="dnai"].pair_iptm_emi2_cdc20.mean()
nul = [e[(e.emi2==al)&(e.apc=="dval")].pair_iptm_emi2_cdc20.mean() -
       e[(e.emi2==al)&(e.apc=="dnai")].pair_iptm_emi2_cdc20.mean() for al in ("dval","dnai")]
allsig = df[df.emi2=="dval"].pair_iptm_emi2_cdc20.mean() - df[df.emi2=="dnai"].pair_iptm_emi2_cdc20.mean()
allnul = [df[(df.emi2==al)&(df.apc=="dval")].pair_iptm_emi2_cdc20.mean() -
          df[(df.emi2==al)&(df.apc=="dnai")].pair_iptm_emi2_cdc20.mean() for al in ("dval","dnai")]
names = ["signal\n(all 60)","noise\n(all 60)","noise\n(all 60)",
         "signal\n(45 engaged)","noise\n(45 engaged)","noise\n(45 engaged)"]
vals  = [allsig]+allnul+[sig]+nul
cols  = ["#b0202a","#9aa4b0","#9aa4b0","#b0202a","#9aa4b0","#9aa4b0"]
ax[2].bar(range(6), vals, color=cols, edgecolor="k", linewidth=.6)
ax[2].axhline(0,c="k",lw=.8)
ax[2].set_xticks(range(6)); ax[2].set_xticklabels(names, fontsize=7.5)
ax[2].set_ylabel("Δ Emi2–Cdc20 ipTM")
ax[2].set_title("C  Allele effect vs the empirical noise floor\n(noise = identical-partner contrast)", fontsize=9.5)
for i,v in enumerate(vals):
    ax[2].text(i, v+(0.012 if v>=0 else -0.028), f"{v:+.3f}", ha="center", fontsize=8)
ax[2].grid(alpha=.25, axis="y")
plt.tight_layout()
out = FIG/"AF3_S1_summary_statistics.png"
plt.savefig(out, dpi=300); print("wrote", out)
