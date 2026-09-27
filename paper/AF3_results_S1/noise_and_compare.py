#!/usr/bin/env python
"""Task 2 + 3: noise floor and the Dval-vs-Dnai comparison."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
df = pd.read_csv(ROOT / "AF3_S1_confidence_table.csv")
pd.set_option("display.width", 250, "display.max_columns", 60)

# ---- 0. verify APC sequences really are (near-)identical between sources
seqs = {}
for job in sorted(ROOT.glob("s1_dbox_*")):
    jr = json.loads(next(job.glob("*_job_request.json")).read_text())[0]
    s = [x["proteinChain"]["sequence"] for x in jr["sequences"] if "proteinChain" in x]
    seqs[job.name] = s
emi2 = {k: v[0] for k, v in seqs.items()}
cdc20 = {k: v[1] for k, v in seqs.items()}
apc10 = {k: v[2] for k, v in seqs.items()}

def diffs(a, b):
    assert len(a) == len(b), (len(a), len(b))
    return [(i + 1, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y]

ev = emi2["s1_dbox_emi2dval_apcdval_seed1"]; en = emi2["s1_dbox_emi2dnai_apcdnai_seed1"]
print("Emi2 fragment length:", len(ev))
print("Emi2 dval-vs-dnai diffs (CIF numbering / full-length):",
      [(i, x, y, i + 534) for i, x, y in diffs(ev, en)])
print("CIF 14-17 dval:", ev[13:17], " dnai:", en[13:17])
print("CIF 15 (SNP549) dval/dnai:", ev[14], en[14], "| CIF 56 (SNP590):", ev[55], en[55])
print("Cdc20 len", len(cdc20["s1_dbox_emi2dval_apcdval_seed1"]), "diffs:",
      diffs(cdc20["s1_dbox_emi2dval_apcdval_seed1"], cdc20["s1_dbox_emi2dval_apcdnai_seed1"]))
print("Apc10 len", len(apc10["s1_dbox_emi2dval_apcdval_seed1"]), "diffs:",
      diffs(apc10["s1_dbox_emi2dval_apcdval_seed1"], apc10["s1_dbox_emi2dval_apcdnai_seed1"]))
# sanity: all jobs with same tag share sequence
for tag in ("dval", "dnai"):
    assert len({v for k, v in emi2.items() if f"emi2{tag}" in k}) == 1
    assert len({v for k, v in cdc20.items() if f"apc{tag}" in k}) == 1

METRICS = ["pair_iptm_emi2_cdc20", "pair_iptm_emi2_apc10", "iptm", "ranking_score"]

print("\n########## TASK 2: NOISE FLOOR ##########")
print("\n--- per emi2 allele: distribution of Emi2-Cdc20 ipTM across the 2 partner sets x 3 seeds x 5 models (n=30) ---")
for al in ("dval", "dnai"):
    s = df[df.emi2_allele == al]["pair_iptm_emi2_cdc20"]
    print(f"emi2_{al}: mean={s.mean():.3f} SD={s.std(ddof=1):.3f} range=[{s.min():.2f}, {s.max():.2f}] n={len(s)}")

print("\n--- partner-set effect (should be ~0: identical contact residues) ---")
for al in ("dval", "dnai"):
    a = df[(df.emi2_allele == al) & (df.apc_source == "dval")]["pair_iptm_emi2_cdc20"]
    b = df[(df.emi2_allele == al) & (df.apc_source == "dnai")]["pair_iptm_emi2_cdc20"]
    t = stats.mannwhitneyu(a, b)
    print(f"emi2_{al}: apc_dval {a.mean():.3f}+-{a.std(ddof=1):.3f} vs apc_dnai {b.mean():.3f}+-{b.std(ddof=1):.3f}"
          f" | delta={a.mean()-b.mean():+.3f}  MWU p={t.pvalue:.3g}")

print("\n--- variance decomposition: seed-to-seed vs model-to-model within a seed ---")
for m in METRICS:
    rec = []
    for (al, ap), g in df.groupby(["emi2_allele", "apc_source"]):
        seedmeans = g.groupby("seed")[m].mean()
        within = g.groupby("seed")[m].std(ddof=1)
        rec.append((al, ap, seedmeans.std(ddof=1), within.mean(), seedmeans.min(), seedmeans.max()))
    r = pd.DataFrame(rec, columns=["emi2", "apc", "SD_between_seeds", "mean_SD_within_seed", "seedmean_min", "seedmean_max"])
    print(f"\n[{m}]")
    print(r.round(3).to_string(index=False))

print("\n--- per-seed means of Emi2-Cdc20 ipTM (all 12 jobs) ---")
piv = df.pivot_table(index=["emi2_allele", "seed"], columns="apc_source",
                     values="pair_iptm_emi2_cdc20", aggfunc=["mean", "std", "min", "max"])
print(piv.round(3))

print("\n--- raw per-job values ---")
for job, g in df.groupby("job"):
    v = g.sort_values("model_index")["pair_iptm_emi2_cdc20"].values
    print(f"{job:40s} {np.array2string(v, precision=2)}")

print("\n########## TASK 3: Emi2^Dval (n=30) vs Emi2^Dnai (n=30) ##########")
for m in METRICS:
    a = df[df.emi2_allele == "dval"][m]; b = df[df.emi2_allele == "dnai"][m]
    pooled = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2))
    d = (a.mean() - b.mean()) / pooled if pooled else np.nan
    mwu = stats.mannwhitneyu(a, b)
    tt = stats.ttest_ind(a, b, equal_var=False)
    print(f"\n[{m}] dval {a.mean():.3f}+-{a.std(ddof=1):.3f}  dnai {b.mean():.3f}+-{b.std(ddof=1):.3f}"
          f"\n   delta = {a.mean()-b.mean():+.3f}   Cohen d = {d:+.2f}   MWU p={mwu.pvalue:.3g}   Welch p={tt.pvalue:.3g}")
    # noise-floor comparison: same statistic computed on the null (partner-set) contrast
    nulls = []
    for al in ("dval", "dnai"):
        x = df[(df.emi2_allele == al) & (df.apc_source == "dval")][m]
        y = df[(df.emi2_allele == al) & (df.apc_source == "dnai")][m]
        nulls.append(x.mean() - y.mean())
    print(f"   null (partner-set) deltas for the SAME metric: {[round(v,3) for v in nulls]}"
          f"  -> |signal| {abs(a.mean()-b.mean()):.3f} vs max |null| {max(abs(v) for v in nulls):.3f}")

# paired-by-seed test (block on seed & partner set)
print("\n--- paired contrast: dval-minus-dnai within each (apc_source, seed, model_index) cell ---")
w = df.pivot_table(index=["apc_source", "seed", "model_index"], columns="emi2_allele",
                   values="pair_iptm_emi2_cdc20")
w["delta"] = w["dval"] - w["dnai"]
print(w["delta"].describe().round(3).to_string())
print("Wilcoxon:", stats.wilcoxon(w["dval"], w["dnai"]))
