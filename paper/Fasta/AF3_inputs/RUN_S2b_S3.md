# AF3 screens S2b (ZBR) and S3 (Apc1) — upload pack

## Priority
1. **S2b first** — Emi2(535–675) + Anapc2(matched catalytic) + Anapc11. Closes ZBR vs Apc2/Apc11.
2. **S3 second** — Emi2(535–675) + Anapc1(950–1250). Closes whether V590 / ZBR faces Apc1.

Do **not** upload full-length Anapc1 (1954 aa) for this test.

## Upload to AlphaFold Server
1. Open https://alphafoldserver.com/
2. Use **Upload JSON** (not manual paste).
3. Upload **`AF3_jobs_S2b_ZBR_withZn.json`** (12 jobs). Prefer this over noZn.
4. After S2b finishes, upload **`AF3_jobs_S3_Apc1_withZn.json`** (12 jobs).
5. Server limit is often ~20 jobs/batch — each file is 12, safe.

## Job design (matches S1)
- 4 allele combinations × 3 seeds × 5 diffusion samples = 60 models per screen
- Seeds 1,2,3 (same as published S1)
- `useStructureTemplate: true`
- dialect `alphafoldserver`, version 3
- S2b: Zn × 4 (Emi2 ZBR + Apc11 RING)
- S3: Zn × 2 (Emi2 ZBR)

## Constructs
### S2b (773 aa total)
| chain | construct |
|---|---|
| A | Emi2 535–675 (141) |
| B | Anapc2 Dval 220–767 / Dnai matched 1–548 (548) |
| C | Anapc11 1–84 (identical parents) |

### S3 (442 aa total)
| chain | construct |
|---|---|
| A | Emi2 535–675 (141) |
| B | Anapc1 Dval 950–1250 (301); Dnai = aligned equivalent (no gap 1284–1330) |

Apc1 window covers literature ZBR-facing homology **1050–1125** and stops before the Dnai deletion at 1284–1330.

## What to download
For every job folder, keep at least:
- `*_model_*.cif` (5)
- `*_summary_confidences_*.json` (5)
- `*_full_data_*.json` (5) if available
- `*_job_request.json`

Suggested drop folders:
- `AF3_results_S2b/`
- `AF3_results_S3/`

## How we will score (after you return results)
- Engagement filter: pair_ipTM(Emi2–Apc11) and/or Emi2–Apc2 (S2b); Emi2–Apc1 (S3)
- Distances from Emi2 **590** and ZBR **600–648** to partners
- Dual-side: only if a contact ≤4 Å faces a variable APC residue
- Seed = unit of replication (as in S1)

## Optional later
- Full Anapc2 S2 (non-matched, 767 aa) — only if S2b docks but looks truncated
- Cross-check noZn JSON if Zn jobs look over-constrained
