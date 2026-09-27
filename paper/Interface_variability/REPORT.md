# Two-sided variability test (Emi2-native): *D. valentini* × *D. r. nairensis*

**Scope.** All geometry in this report comes from AlphaFold3 screen S1 (Emi2 residues 535–675 + Cdc20 + Anapc10; 60 models). There is **no Emi1 / PDB 4UI9 surrogate**. Prior Emi1-based files live in `Interface_variability/superseded/`.

## 0. Emi2 element boundaries

| element | residues | APC/C role in this study |
|---|---|---|
| F-box domain | 452–537 | Skp1/SCF — **not** APC/C |
| F-box motif (inside domain) | 479–519 | Skp1/SCF — **not** APC/C |
| D-box | 548–551 | binds Cdc20 WD40 in S1 (validated vs PDB 5G04) |
| ZBR | 600–648 | literature: Apc11/Apc2/Apc1 — **no Emi2 complex model here** |
| RL tail | 658–675 | literature: Apc2 CTD — **no atomic coverage in S1 pose** |

## 1. Emi2 side

Dval 675 aa vs Dnai 675 aa: **10 substitutions**.

| position | Dval | Dnai | functional element |
|---|---|---|---|
| 76 | Q | H | — |
| 92 | V | I | — |
| 102 | A | T | — |
| 171 | K | R | — |
| 219 | E | D | — |
| 391 | I | V | — |
| 395 | D | N | — |
| 425 | L | F | — |
| 549 | F | I | D-box |
| 590 | V | M | — |

Substitutions inside an annotated APC/C-binding element: **1** (549, D-box). V590M lies between D-box and ZBR; it is **not** inside the D-box or ZBR ranges above, and S1 does not place it against a resolved APC/C partner (see §3).

## 2. APC/C side (sequence)

| subunit | Dval len | Dnai len | substitutions | note |
|---|---|---|---|---|
| Cdc20 | 454 | 454 | V291I | |
| Anapc10 | 185 | 185 | **0 — identical** | |
| Anapc11 | 84 | 84 | **0 — identical** | |
| Anapc2 | 767 | 568 | V632A; P654L | |
| Anapc5 | 743 | 722 | Y728H | |
| Anapc7 | 565 | 565 | **0 — identical** | |
| Anapc1 | 1954 | 1902 (partial) | gap Dval 1284-1330 absent in Dnai; I133V; M443L; G802A; I823V; F1844L; K1943N; S1944P; T1945Q; G1946A; V1948M; V1950L | no Emi2–Apc1 complex; see §2.1 |
| Anapc4 | partial Dval | — | *Dnai missing* | no Emi2–Apc4 model |

Cdc20 V291I median distance to Emi2 D-box (548–551) across 45 engaged S1 models: **31.8 Å**.

### 2.1 Anapc1 Dnai gap

Relative to Dval (1954 aa), the Dnai partial ORF (1902 aa) lacks **Dval residues 1284-1330** (47 aa — the ~48-residue hole). The rest of the ORF aligns continuously (97% of Dval covered).

Substitutions in the overlapping region: **11** (I133V; M443L; G802A; I823V; F1844L; K1943N; S1944P; T1945Q; G1946A; V1948M; V1950L).

None of those substitutions fall in Dval 1050–1125 (homology window facing an Emi-family ZBR contact on human Apc1 in 4UI9 — literature context only, not an Emi2 measurement). That window is fully present in both alleles. The Dnai gap instead sits at Dval 1284–1330, overlapping the human Apc1 surface that cradles Apc10 in 4UI9 — again literature context; we have **no Emi2–Apc1 model**, so this does not by itself create a dual-side Emi2 hit.

## 3. Geometry from Emi2 S1 models (not Emi1)

Engaged models used for the distance matrix: **45** (pair_iptm_emi2_cdc20 ≥ 0.5). Matrix: `emi2_s1_distance_matrix.json` (median over engaged models).

| Emi2 site | element | median Å to Cdc20 | median Å to Apc10 | interpretation |
|---|---|---|---|---|
| 548 | D-box | 2.24 | 11.03 † | canonical D-box contacts on Cdc20 |
| 549 | D-box | 3.72 | 17.39 † | side chain ≤4 Å of Cdc20 in **0/45** engaged models (backbone may approach; allele effect null in REPORT_AF3_S1) |
| 550 | D-box | 3.37 | 17.13 † |  |
| 551 | D-box | 2.62 | 17.79 † | canonical D-box contacts on Cdc20 |
| 590 | — | 11.94 | 7.24 † | no publishable APC partner in S1; ZBR/Apc1/Apc2 not in the construct |

† Apc10 column is a **non-measurement** (co-receptor site never occupied; see `AF3_results_S1/REPORT_AF3_S1.md`).

Coverage of S1: only Emi2 **535–675** is present. Residues 1–534, including eight of the ten substitutions and the F-box, have **no complex geometry** in this study.

## 4. Intersection — parthenogenetic substitution list

The hypothesis requires an Emi2 substitution that both (i) sits in a structurally tested interface and (ii) faces a variable APC/C position.

- Structurally tested interface in hand: **Emi2 D-box ↔ Cdc20**.
- Emi2 substitutions in that element: **F549I**.
- Cdc20 variability on/near that surface: **none** (sole allele difference V291I at 31.8 Å from the D-box).
- Apc10: identical between parents; geometry non-informative in S1.
- V590M / ZBR / Apc1 / Apc2 / Apc11: **not structurally tested on Emi2** (AF3 S2 inputs exist but were not run; no Emi2–Apc1 complex).
- Anapc1 sequence: Dnai lacks Dval 1284–1330 (47 aa); overlap has substitutions but **none** in the homology window opposite an Emi-family ZBR site (1050–1125). Without an Emi2–Apc1 model this remains sequence context, not a dual-side hit.

**Candidate list for this triad, given available Emi2 models: empty.**

## 5. What was removed

The previous report measured distances to human Emi1 in PDB 4UI9 and transferred them onto Emi2 by alignment. That surrogate is withdrawn from the publishable line. Literature on Emi1/Emi2 family architecture may still be cited as background, but not as a measurement of Darevskia Emi2 allele contacts.

## 6. Files

- `emi2_two_sided_test.py` — this analysis
- `emi2_s1_distance_matrix.json` — median distances from engaged S1 models
- `emi2_variable_sites.csv`, `apc_variable_sites.csv`
- `../Fasta/Anapc1/` — Dval full + Dnai partial + alignment
- `superseded/` — Emi1/4UI9 artefacts

