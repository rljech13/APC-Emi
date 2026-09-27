# AF3 screen S1 — Emi2(535–675) + Cdc20 + Anapc10, Dval vs Dnai

Analysis date 2026-07-27. Working dir `/Users/dmitrij/Desktop/Работа/APC:C-Emi/AF3_results_S1/`.
60 models = 4 combinations × 3 seeds × 5 diffusion samples. No manuscript file or AF3 input JSON was modified.

## 0. Construct verification (done against the CIFs, not assumed)
- CIF chains: A = Emi2 141 aa (renumbered 1–141), B = Cdc20 454 aa, C = Anapc10 185 aa, D/E = Zn. Confirmed.
- CIF 14–17 = **RFAL** in every Dval job and **RIAL** in every Dnai job. Mapping (full-length = CIF + 534) is correct.
- Within this 535–675 fragment the two alleles differ at exactly **two** positions, not ten:
  CIF 15 = F549I and CIF 56 = V590M. The other eight Emi2 SNPs lie outside the construct.
- APC partners: Cdc20 differs at **one** position (V291I, our numbering; ~30 Å from the D-box groove);
  Anapc10 sequences are **identical**. The partner-set contrast is therefore a genuine null.

## 1. Confidence table
`AF3_S1_confidence_table.csv` (60 rows). `has_clash` = 0 for all 60. `ptm` 0.45–0.49, `fraction_disordered` 0.28–0.33.

The informative columns are strongly **bimodal**, and the mode is set at the seed level:

| | Emi2–Cdc20 ipTM | Emi2–Apc10 ipTM | Cdc20–Apc10 ipTM | PAE_min Emi2→Apc10 |
|---|---|---|---|---|
| 9 "docked" jobs (45 models) | 0.52–0.60 | 0.06–0.13 | 0.09–0.32 | ~20 Å |
| 3 "failed" jobs (15 models) | 0.09–0.20 | 0.06–0.13 | 0.09–0.32 | ~20 Å |

Failed jobs: `emi2dnai_apcdval_seed3`, `emi2dval_apcdnai_seed2`, `emi2dval_apcdnai_seed3`.

**Apc10 never reaches its D-box co-receptor site** in any of the 60 models (ipTM ≈ 0.09 with both
partners, PAE_min 20–24 Å). It is not contact-free: 59/60 models put 1–39 Emi2 residues (median 14)
within 4 Å of Apc10, but the contact set is irreproducible and never uses the co-receptor surface —
see the caveat in section 4.

## 2. Noise floor (partner-set contrast; contact residues identical by construction)
- SD **within a seed** (5 diffusion samples): 0.005–0.026 ipTM. The 5 models per seed are near-duplicates.
- SD **between seeds**: 0.001–0.260 ipTM, entirely driven by whether the seed found the docked basin.
- Partner-set Δ on all 60: **+0.290** (Emi2^Dval) and **−0.127** (Emi2^Dnai) — both must be zero biologically.
- Partner-set Δ on the 45 engaged models: −0.007 and +0.023.

**Effective replication is 12 seeds, not 60 models.** Treating models as replicates inflates n five-fold.

## 3. Emi2^Dval vs Emi2^Dnai — null
All 60 models: Δ(Emi2–Cdc20 ipTM) = −0.069, Cohen d = −0.36, MWU p = 0.062.
This is **4× smaller than the 0.290 noise floor** and is entirely an artefact of 2/6 Dval seeds vs
1/6 Dnai seeds failing to dock (Fisher exact p = 1.00).

Restricted to the 45 D-box-engaged models, grouped into the 9 seeds that produced them
(4 Emi2^Dval, 5 Emi2^Dnai), every metric is null and every effect is smaller than the
corresponding identical-partner null. **The seed is the unit of replication throughout the
table** (see section 2), so entries are mean ± SD over seed means, n = 4 vs 5 seeds:

| metric | Dval | Dnai | Δ | max \|null\| | MWU p |
|---|---|---|---|---|---|
| Emi2–Cdc20 ipTM | 0.583 ± 0.004 | 0.577 ± 0.029 | +0.005 | 0.023 | 0.62 |
| ipTM | 0.464 ± 0.004 | 0.459 ± 0.021 | +0.005 | 0.020 | 0.81 |
| ranking_score | 0.624 ± 0.004 | 0.620 ± 0.018 | +0.003 | 0.014 | 0.71 |
| D-box RMSD vs Hsl1 (Å) | 0.65 ± 0.03 | 0.60 ± 0.14 | +0.04 | 0.08 | 0.19 |
| Emi2–Cdc20 BSA (Å²) | 779 ± 130 | 774 ± 69 | +4 | 93 | 0.90 |
| Emi2–Apc10 ipTM† | 0.084 ± 0.006 | 0.088 ± 0.011 | −0.004 | 0.008 | n/a |

Grouping by seed changes only the SDs and the p values; the means are identical either way,
because every engaged seed contributes all 5 of its diffusion samples. (`ranking_score` reads
0.624 rather than the 0.623 of earlier drafts only because its Dval mean is exactly 1247/2000 =
0.6235, a rounding tie now resolved upwards.) Per-model tests (n = 20 vs 25 models, for the five
measured rows in the order above) give p = 0.24, 0.69, 0.81, 0.03 and 0.94. Those
are pseudoreplication: the 5 diffusion samples of one seed are near-duplicates (within-seed SD
0.005–0.026 ipTM, section 2), so a per-model n inflates replication five-fold. The seed-level
column is the conservative one and is what Figure 3e shows.

\†Not a measurement. Apc10 does contact Emi2 in 59 of the 60 models, but never at its D-box
co-receptor site: its centre of mass sits 12.7–86.0 Å from its reference position in the
Cdc20 frame (median 64.9 Å, all 60 models), and the contacts it does make are irreproducible
(section 4). This ipTM therefore scores an arbitrary encounter, so no test is quoted for it.

**Conclusion: no detectable effect of F549I (or V590M) on Emi2–APC/C binding in this model set.**

## 4. Is the D-box canonically engaged? — YES
Reference: PDB **5G04** — cryo-EM of the intact human APC/C–Cdc20–Hsl1 at **3.9 Å**
(Zhang et al. 2016, Nature 533:260, doi:10.1038/nature17973); chain R = human Cdc20,
chain S = Hsl1 D-box peptide **R828-A-A-L831** (only 828–837 are modelled), chain L = Apc10.
Note the resolution when reading the RMSD values below: 0.5–0.6 Å is well inside the coordinate
error of the reference, so the defensible claim is the reproduced binding mode and contact set,
not the RMSD magnitude.

All nine canonical D-box receptor residues of Cdc20 are **identical** in Darevskia
(human → ours): I175→I131, L176→L132, D177→D133, P179→P135, Y207→Y163, W209→W165, I216→I172,
E465→E422, R468→R425. Emi2 contacts 8 of these 9 in all four top-ranked models.

After superposing our Cdc20 on 5G04 Cdc20 (CA RMSD 0.53–0.60 Å over ~300 core residues), with **no local fit**:

| | Dval/Dval | Dval/Dnai | Dnai/Dval | Dnai/Dnai |
|---|---|---|---|---|
| D-box 548–551 backbone RMSD to Hsl1 828–831 | 0.6 Å | 0.6 Å | 0.5 Å | 0.5 Å |
| L551 side-chain centroid vs Hsl1 L831 | 0.9 Å | 0.9 Å | 1.0 Å | 0.9 Å |
| R548 guanidinium centroid vs Hsl1 R828 | 2.2 Å | 2.8 Å | 2.9 Å | 2.9 Å |

The R548 guanidinium offset is a rotamer difference, not a placement failure: R548 makes the canonical
bidentate salt bridge to **E422 (= human E465)** at 2.34 Å plus a second to **D133 (= human D177)**,
in both alleles. In 5G04 itself Hsl1 R828 makes only the E465 contact at 3.59 Å, so AF3 is if anything
over-idealising an otherwise correct interaction.

**Caveat that must be reported: the D-box co-receptor never assembles.** Apc10 is in the input, and it
does touch Emi2 — 59/60 models place 1–39 Emi2 residues (median 14) within 4 Å of it. But it never
occupies its co-receptor position, and the contacts it makes are not the co-receptor contacts:

- **Position.** Displacement of the Apc10 centre of mass from its 5G04 position in the Cdc20 frame is
  **12.7–86.0 Å, median 64.9 Å** over all 60 models (`apc10_com_disp` in `AF3_S1_pose_table.csv`);
  the 12 per-job means span 47.2–76.1 Å. As a rigid body, with no local fit, Apc10 is 34–84 Å CA RMSD
  from its 5G04 placement in the 8 models checked in detail (the 4 closest plus the 4 top-ranked).
- **The three closest models are coincidence, not assembly.** Three models bring the COM within 20 Å
  (12.7, 17.8, 18.8 Å), all of them in `emi2dval_apcdnai` seeds 2–3, i.e. the failed low-ipTM jobs.
  In all three the D-box itself is 5–64 Å from its Hsl1 position (none passes the < 2 Å engagement
  criterion), the Apc10 CA RMSD is still 34 Å, and the Apc10 surface presented to Emi2 shares **zero**
  residues with the surface that grips Hsl1 833–837 in 5G04 (Apc10 87–89 and 146–149). The same zero
  overlap holds in the 4 top-ranked models. So no model assembles a native-like co-receptor.
- **What happens instead.** The Emi2 residues contacting Apc10 are drawn from 118 different positions
  across the whole 536–675 construct, 84 % of them in the 588–651 linker/ZBR region and only 3 % in the
  RL tail; 14/60 models happen to put some of 548–551 against Apc10. Contact-set Jaccard is 0.17 within
  a job vs 0.18 between jobs — chance level (Emi2–Cdc20, for contrast: 0.41 vs 0.34).

This is expected: Apc10 is held in place by Apc1 (which grips its core) and by Apc3/Cdc27 (which
receives its C-terminal IR tail) — Apc2 makes no contact with Apc10 at all. Both are absent from a three-chain
input. The consequence is that the C-terminal half of the degron (Hsl1 833–837 in 5G04) is unsupported
here, and every Emi2–Apc10 number in the dataset is noise, not a measurement.

### Where residue 549 points — solvent-exposed
Top-ranked models, Shrake–Rupley, Tien et al. 2013 reference ASA:

| residue | rel. SASA in complex | ΔSASA on binding | % buried |
|---|---|---|---|
| R548 | 24–41 % | 100–144 Å² | 51–69 % |
| **549 (F/I)** | **43–49 %** | **14–20 Å²** | **11–18 %** |
| A550 | 8–13 % | 78–81 Å² | 83–88 % |
| L551 | 1.0–1.3 % | 177–178 Å² | **99 %** |
| 590 (V/M) | 23–27 % | **0 Å²** | 0 % |

Residue 549 occupies the first "x" of the R-x-x-L motif and points **away from** the groove
(equivalent to Hsl1 A829). Residue 590 makes no intermolecular contact at all. Both polymorphic
positions in this construct are therefore structurally incapable of modulating D-box affinity —
which is the mechanistic reason the ipTM comparison is null.

### Atom level: the side chain of 549 is outside the interface — all 45 engaged models
The SASA table above is the four top-ranked models. Repeating the measurement atom by atom over
**all 45 D-box-engaged models** (`sidechain_contacts.py` → `AF3_S1_sidechain_contacts.csv` and
`AF3_S1_sidechain_contacts_pairs.csv`) shows where that small burial comes from. Side chain = every
heavy atom except N, CA, C, O, OXT; **"contact" = 4 Å heavy-atom distance throughout**, with 5 Å
used as a second, looser "not even close" check. Values are median (range) over the 45 models.

| Emi2 residue | min. side chain → Cdc20 | min. backbone → Cdc20 | Cdc20 residues within 4 Å |
|---|---|---|---|
| R548 | 2.24 Å (2.05–2.58) | 3.88 Å (3.54–4.27) | 3 (2–4): D133 45/45, E422 45/45, P135 41/45, A134 5/45 |
| **549 (F/I)** | **5.30 Å (4.21–6.11)** | **3.72 Å (3.46–3.95)** | **1 (1–1): D133 45/45 — nothing else, ever** |
| L551 | 3.56 Å (3.38–3.64) | 2.62 Å (2.40–2.72) | 6 (6–7): L132, D133, V156, Y163, W165, I172 all 45/45; L424 17/45 |
| 590 (V/M) | 12.70 Å (6.45–37.25) | 13.20 Å (4.94–36.32) | **0 (0–0)** |

- **No side-chain atom of residue 549 reaches Cdc20 in any model: 0/45 within 4 Å.** The 5.30 Å
  median above is bimodal and the two alleles do not overlap. Phe549 (20 Emi2^Dval models) never
  comes closer than **5.68 Å** (range 5.68–6.11, median 5.92), and its closest atom is always CB,
  the pivot — the ring itself is farther out; **0/20 within 5 Å**. Ile549 (25 Emi2^Dnai models)
  reaches **4.21–5.46 Å** (median 4.65) through one recurrent atom pair, CG2 to Ile172 CD1, which
  is the closest approach in 25/25 models; 20/25 fall inside 5 Å but **0/25 inside 4 Å**. The five
  that stay beyond 5 Å are the five models of `emi2dnai_apcdnai_seed3`, i.e. the small spread is
  seed-level, as everything else in this data set is; no seed or job brings the side chain into
  contact range.
- **The burial of 549 is entirely backbone.** In 45/45 models the closest Emi2 549 atom to Cdc20 is
  the backbone carbonyl O and its partner is always Asp133 CB (3.46–3.95 Å; Dval 3.46–3.77,
  Dnai 3.60–3.95). Asp133 is also the *only* Cdc20 residue within 4 Å of residue 549, in every one
  of the 45 models; no second partner appears anywhere in the set. The ΔSASA of 14–20 Å² in the
  table above is therefore backbone surface, not side-chain surface.
- **Consequence.** The backbone of position 549 is identical in the two alleles by definition, and
  it is the only part of the residue that touches Cdc20. The part that actually differs between
  *D. valentini* and *D. r. nairensis* — the side chain, Phe vs Ile — lies outside the interface in
  every engaged model. This is a structural reason for the null result that does not depend on the
  statistics of section 3.
- **Position 590 is not merely uninvolved, it is remote.** 0/45 models place any atom within 4 Å of
  Cdc20 (or within 5 Å, except one model at 4.94 Å); the closest approach of any atom is
  4.94–36.32 Å, median 11.9 Å. This is a direct distance confirmation of the ΔSASA = 0 Å² result.

The two flanking anchors behave exactly as the reference structure demands, which is the internal control for
the measurement: R548 grips through its guanidinium (side chain 2.05–2.58 Å; the closest partner is
E422 OE2 in 40/45 models and D133 OD1 in the other 5, with both residues within 4 Å in 45/45), and
L551 is fully engaged on both counts — backbone O to D133 at 2.40–2.72 Å plus a side chain buried
against Y163/W165/L424 at 3.38–3.64 Å, 6–7 Cdc20 partners in every model. Against those numbers the
single, backbone-only, single-partner contact of residue 549 is unambiguous.

## 5. Pose reproducibility — reproducible for the D-box only
- **D-box**: canonical in **45/60** models; within those, RMSD to Hsl1 = 0.62 ± 0.13 Å (range 0.46–0.99 Å)
  and L551 side chain 0.92 ± 0.05 Å. Extremely tight. All-or-nothing per seed (5/5 or 0/5).
- **Everything else**: not reproducible. Emi2 centre of mass, after superposing Cdc20, has an RMS spread
  of 18–25 Å within a single combination (max 43 Å). Interface-set Jaccard for Emi2–Cdc20 is 0.41 within a
  job vs 0.34 between jobs — barely above chance. For Emi2–Apc10 it is 0.17 vs 0.18, i.e. pure noise.
- Residues present in ≥90 % of the 45 engaged models: **CIF 14–20 only (Emi2 548–554)**. Zero residues
  qualify for the Emi2–Apc10 interface.

The linker, ZBR and RL tail are placed differently in essentially every model. Only the 7-residue D-box
segment is a real, reproducible structural observation.

## 6. Files written
Tables (this directory): `AF3_S1_confidence_table.csv`, `AF3_S1_pose_table.csv`,
`AF3_S1_dbox_validation.csv`, `AF3_S1_sasa.csv`, `AF3_S1_bsa.csv`, `AF3_S1_contacts.json`,
`AF3_S1_sidechain_contacts.csv`, `AF3_S1_sidechain_contacts_pairs.csv` (atom-level, 45 engaged
models; written by `sidechain_contacts.py`).
Aligned coordinates in the 5G04 frame: `aligned/`.
Figures in `/Users/dmitrij/Desktop/Работа/APC:C-Emi/Figures/`, all prefixed `AF3_S1_`.
