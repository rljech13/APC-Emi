#!/usr/bin/env python
"""Parse all AF3 summary_confidences JSONs into a tidy CSV."""
import json, re, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent
CHAINS = ["A", "B", "C", "Zn1", "Zn2"]
LABEL = {"A": "emi2", "B": "cdc20", "C": "apc10", "Zn1": "zn1", "Zn2": "zn2"}

rows = []
for job in sorted(ROOT.glob("s1_dbox_*")):
    m = re.match(r"s1_dbox_emi2(dval|dnai)_apc(dval|dnai)_seed(\d)", job.name)
    if not m:
        sys.exit(f"unexpected job name {job.name}")
    emi2, apc, seed = m.group(1), m.group(2), int(m.group(3))
    for f in sorted(job.glob("*_summary_confidences_*.json")):
        idx = int(re.search(r"_(\d)\.json$", f.name).group(1))
        d = json.loads(f.read_text())
        row = dict(job=job.name, emi2_allele=emi2, apc_source=apc, seed=seed,
                   model_index=idx,
                   ranking_score=d["ranking_score"], iptm=d["iptm"], ptm=d["ptm"],
                   fraction_disordered=d["fraction_disordered"],
                   has_clash=d["has_clash"], num_recycles=d.get("num_recycles"))
        for i, c in enumerate(CHAINS):
            row[f"chain_iptm_{LABEL[c]}"] = d["chain_iptm"][i]
            row[f"chain_ptm_{LABEL[c]}"] = d["chain_ptm"][i]
        for i, ci in enumerate(CHAINS):
            for j, cj in enumerate(CHAINS):
                if i < j:
                    row[f"pair_iptm_{LABEL[ci]}_{LABEL[cj]}"] = d["chain_pair_iptm"][i][j]
        # PAE min matrix is not symmetric in AF3 output; keep both directions for protein chains
        for i, ci in enumerate(CHAINS[:3]):
            for j, cj in enumerate(CHAINS[:3]):
                if i != j:
                    row[f"pae_min_{LABEL[ci]}_to_{LABEL[cj]}"] = d["chain_pair_pae_min"][i][j]
        rows.append(row)

df = pd.DataFrame(rows).sort_values(["emi2_allele", "apc_source", "seed", "model_index"])
out = ROOT / "AF3_S1_confidence_table.csv"
df.to_csv(out, index=False)
print(f"wrote {out}  n={len(df)}")
print(df.groupby(["emi2_allele", "apc_source"]).size())
pd.set_option("display.width", 250, "display.max_columns", 100)
key = ["ranking_score", "iptm", "ptm", "pair_iptm_emi2_cdc20", "pair_iptm_emi2_apc10",
       "pair_iptm_cdc20_apc10", "fraction_disordered", "has_clash",
       "pae_min_emi2_to_cdc20", "pae_min_emi2_to_apc10", "pae_min_cdc20_to_apc10"]
print("\n=== overall (n=60) ===")
print(df[key].describe().T[["mean", "std", "min", "50%", "max"]].round(3))
print("\n=== by emi2 allele x apc source ===")
print(df.groupby(["emi2_allele", "apc_source"])[key].mean().round(3))
