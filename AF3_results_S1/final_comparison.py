#!/usr/bin/env python
"""Re-run task 3 on the 45 D-box-engaged models only, with seed as the unit of replication,
plus interface buried-surface-area (BSA), a physically interpretable alternative to ipTM."""
import warnings, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent
pose = pd.read_csv(ROOT / "AF3_S1_pose_table.csv")
conf = pd.read_csv(ROOT / "AF3_S1_confidence_table.csv").rename(
    columns={"emi2_allele": "emi2", "apc_source": "apc", "model_index": "model"})
df = pose.merge(conf, on=["emi2", "apc", "seed", "model"])
df["dbox_ok"] = df.dbox_bb_rmsd < 2.0
pd.set_option("display.width", 250, "display.max_columns", 60)

print("########## A. Bimodality: the seed is the unit of replication ##########")
j = df.groupby(["emi2", "apc", "seed"]).agg(
    n_ok=("dbox_ok", "sum"), iptm_AB=("pair_iptm_emi2_cdc20", "mean"),
    dbox=("dbox_bb_rmsd", "mean")).reset_index()
j["outcome"] = np.where(j.n_ok == 5, "ENGAGED", np.where(j.n_ok == 0, "FAILED", "MIXED"))
print(j.round(2).to_string(index=False))
tab = pd.crosstab(j.emi2, j.outcome)
print("\nseed-level outcome by Emi2 allele:\n", tab)
a = int(tab.loc["dval"].get("ENGAGED", 0)); b = int(tab.loc["dval"].get("FAILED", 0))
c = int(tab.loc["dnai"].get("ENGAGED", 0)); d = int(tab.loc["dnai"].get("FAILED", 0))
print(f"Fisher exact (engaged vs failed seeds, dval {a}/{a+b} vs dnai {c}/{c+d}): "
      f"p = {stats.fisher_exact([[a, b], [c, d]])[1]:.3f}")
print("failed seeds:", j[j.outcome == 'FAILED'][['emi2', 'apc', 'seed']].to_dict('records'))

print("\n########## B. Comparison restricted to the 45 D-box-ENGAGED models ##########")
e = df[df.dbox_ok]
print(f"n engaged: dval={sum(e.emi2=='dval')}, dnai={sum(e.emi2=='dnai')}")
METRICS = ["pair_iptm_emi2_cdc20", "pair_iptm_emi2_apc10", "iptm", "ranking_score",
           "dbox_bb_rmsd", "dR548", "dL551", "n_iface_B"]
for m in METRICS:
    x = e[e.emi2 == "dval"][m]; y = e[e.emi2 == "dnai"][m]
    pooled = np.sqrt(((len(x)-1)*x.var(ddof=1)+(len(y)-1)*y.var(ddof=1))/(len(x)+len(y)-2))
    print(f"{m:24s} dval {x.mean():7.3f}+-{x.std(ddof=1):.3f}  dnai {y.mean():7.3f}+-{y.std(ddof=1):.3f}"
          f"  delta {x.mean()-y.mean():+.3f}  d={(x.mean()-y.mean())/pooled if pooled else np.nan:+5.2f}"
          f"  MWU p={stats.mannwhitneyu(x,y).pvalue:.3g}")

print("\n--- same, with SEED as the unit (n=6 seed-means per allele, engaged seeds only) ---")
sm = e.groupby(["emi2", "apc", "seed"])[METRICS].mean().reset_index()
for m in METRICS:
    x = sm[sm.emi2 == "dval"][m]; y = sm[sm.emi2 == "dnai"][m]
    print(f"{m:24s} dval {x.mean():7.3f}+-{x.std(ddof=1):.3f} (n={len(x)})  "
          f"dnai {y.mean():7.3f}+-{y.std(ddof=1):.3f} (n={len(y)})  delta {x.mean()-y.mean():+.3f}  "
          f"MWU p={stats.mannwhitneyu(x,y).pvalue:.3g}")

print("\n--- noise floor on the SAME engaged subset (partner-set contrast, must be null) ---")
for m in METRICS:
    ns = []
    for al in ("dval", "dnai"):
        p1 = e[(e.emi2 == al) & (e.apc == "dval")][m].mean()
        p2 = e[(e.emi2 == al) & (e.apc == "dnai")][m].mean()
        ns.append(p1 - p2)
    sig = e[e.emi2 == "dval"][m].mean() - e[e.emi2 == "dnai"][m].mean()
    print(f"{m:24s} signal {sig:+.3f}   null deltas {[round(v,3) for v in ns]}   "
          f"|signal|/max|null| = {abs(sig)/max(1e-9, max(abs(v) for v in ns)):.2f}")

print("\n########## C. Buried surface area of the Emi2-Cdc20 interface (engaged models) ##########")
sr = ShrakeRupley(n_points=500); P = MMCIFParser(QUIET=True)
rows = []
for r in e.itertuples():
    jn = f"s1_dbox_emi2{r.emi2}_apc{r.apc}_seed{r.seed}"
    f = ROOT / jn / f"fold_{jn}_model_{r.model}.cif"
    s = P.get_structure("x", str(f))[0]
    for cid in [c.id for c in s]:
        if cid not in ("A", "B"): s.detach_child(cid)
    sr.compute(s, level="C"); ab = sum(c.sasa for c in s)
    s2 = P.get_structure("y", str(f))[0]
    for cid in [c.id for c in s2]:
        if cid != "A": s2.detach_child(cid)
    sr.compute(s2, level="C"); sa = sum(c.sasa for c in s2)
    s3 = P.get_structure("z", str(f))[0]
    for cid in [c.id for c in s3]:
        if cid != "B": s3.detach_child(cid)
    sr.compute(s3, level="C"); sb = sum(c.sasa for c in s3)
    rows.append(dict(emi2=r.emi2, apc=r.apc, seed=r.seed, model=r.model, bsa=(sa + sb - ab) / 2))
    print(".", end="", flush=True)
print()
b = pd.DataFrame(rows); b.to_csv(ROOT / "AF3_S1_bsa.csv", index=False)
print(b.groupby(["emi2", "apc"])["bsa"].agg(["mean", "std", "min", "max", "count"]).round(1).to_string())
x = b[b.emi2 == "dval"].bsa; y = b[b.emi2 == "dnai"].bsa
print(f"\nEmi2-Cdc20 BSA: dval {x.mean():.0f}+-{x.std(ddof=1):.0f} A^2   dnai {y.mean():.0f}+-{y.std(ddof=1):.0f} A^2"
      f"   delta {x.mean()-y.mean():+.0f} A^2   MWU p={stats.mannwhitneyu(x,y).pvalue:.3g}")
nulls = [b[(b.emi2 == al) & (b.apc == 'dval')].bsa.mean() - b[(b.emi2 == al) & (b.apc == 'dnai')].bsa.mean()
         for al in ("dval", "dnai")]
print(f"null (partner-set) BSA deltas: {[round(v) for v in nulls]} A^2")
