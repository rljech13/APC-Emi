# Literature grounding — Emi2–APC/C interface

Prepared to close the gap left when [Ground Emi2 binding sites in literature](9c3e7d2c-bb5c-42d7-a7dc-294d8c63ff48) stopped at an API limit. All residue ranges carry their **species numbering**. Darevskia coordinates below were obtained by explicit pairwise alignment, not by a fixed offset.

---

## 1. Emi2 interaction elements — what binds what

| Element | Residues (species) | Partner | Evidence | Primary reference |
|---|---|---|---|---|
| **F-box domain** | 452–537 (*D. valentini*); folded domain from AF3 pLDDT; UniProt human FBXO43 (Q4G163) F-box 490–547 → Darevskia **463–520** at 51.2% identity | **Skp1** (SCF-type complex) | Domain annotation; not an APC/C contact in any structure | UniProt Q4G163; SCF/F-box literature |
| **Skp1-binding F-box motif** | 479–519 (*D. valentini*); canonical 41-aa motif inside the domain | **Skp1** | Same | — |
| **D-box (destruction box)** | 548–551 R-x-x-L (*D. valentini*); Emi1 322–325 in 4UI9 | **Cdc20/Cdh1** WD40 β-propeller **+ Apc10/Doc1** | Cryo-EM 4UI9, 5G04; mutagenesis | Chang et al. 2015 (4UI9); Zhang et al. 2016 (5G04) |
| **ZBR (zinc-binding region)** | 600–648 (*D. valentini*); Emi1 374–422 in 4UI9 | **Apc11** RING, **Apc2** cullin CTD, **Apc1** | Cryo-EM 4UI9; Emi1 ZBR biochemistry | Frye et al. 2013; Chang et al. 2015 |
| **C-terminal RL tail** | 658–675 LRRL (*D. valentini*); absent from 4UI9 model | **Apc2** CTD (competes with Ube2S) | **Required** for inhibition; tail absent from all published structures | Ohe et al. 2010; Sako et al. 2014 |
| **Plk1 polo-box docking sites** | **Mouse:** Thr152, Thr176 (Jia et al. 2015); **Xenopus:** Thr170, Thr195 (corresponding sites); → **Darevskia:** Thr182, Thr206 (see §2) | **Plk1** C-terminal PBD (after CaMKII priming) | ITC, crystallography, mutagenesis | Jia et al. 2015 Sci Rep 5:14626; Inoue et al. 2012 (Xenopus) |
| **β-TrCP degron** | N-terminal phosphodegron (species-specific; not mapped here) | **β-TrCP / SCF** | Ubiquitin-mediated Emi2 turnover | Reimann et al. 2003; others |
| **PP2A-B56 motif** | Central region (species-specific) | **PP2A-B56** | Regulatory, not APC/C docking | Hara et al. 2014 |
| **CaMKII / Mos–MAPK phosphosites** | Upstream of Plk1 sites; Xenopus Thr170/195 | Kinase cascade | Fertilisation-triggered degradation | Hansen et al. 2006; Inoue et al. 2012 |

### F-box verdict (Question 1 — critical)

**The F-box does not bind APC/C.** It is the Skp1-recruitment module of an SCF-type ligase. The manuscript sentence that SNPs in the F-box affect binding to the APC/C “active centre” is **incorrect** and must be replaced. Our PCR variation in exon 1 falls in this Skp1 module; that is a sequencing artefact of primer placement, not evidence that the APC/C interface was targeted.

Emi2 inhibits APC/C through the **D-box + ZBR + RL tail** composite interface (Ohe et al. 2010; Chang et al. 2015).

---

## 2. Thr152 and Thr176 — resolution (Question 2)

### Source and function

- **Paper:** Jia et al. 2015, *Sci Rep* 5:14626 (https://doi.org/10.1038/srep14626).
- **Numbering:** **Mus musculus** Emi2 (FBXO43, UniProt Q8CDI2).
- **Sites:** phospho-**Thr152** and phospho-**Thr176** each recruit one **Plk1** polo-box domain (PBD) molecule; required for Emi2 degradation after fertilisation.
- **Xenopus correspondence:** phospho-Thr170 and phospho-Thr195 (Jia et al. 2015, citing Inoue et al. 2012); mouse Thr176 = Xenopus Thr195.

These are **regulatory** sites (Plk1 recognition), **not** APC/C-docking residues.

### Mapping onto *D. valentini* / *D. r. nairensis*

Pairwise global alignment: local *M. musculus* Emi2 (641 aa, `Fasta/Emi2/Emi2_Mmus_AAseq.fasta`) vs *D. valentini* Emi2 (675 aa). **46.2% identity** over 630 aligned pairs (291 identities).

| Mouse (Jia 2015) | Darevskia Emi2 | Dval | Dnai | Parental difference? |
|---|---|---|---|---|
| Thr152 | **Thr182** | T | T | **No** |
| Thr176 | **Thr206** | T | T | **No** |

Local sequence context (identical between Dval and Dnai):

- Around **182:** `…PLATSTLKSED…`
- Around **206:** `…QQRTSTIDDSK…`

**Important:** In Darevskia numbering, positions **152** (Leu) and **176** (Lys) are **not** threonines. The co-author’s check used **mouse coordinates applied to lizard numbering** — the biological conclusion (no difference between parental species at the Plk1 sites) is likely correct, but the residue numbers were wrong.

Nearest of our ten parental substitutions: **K171R** (Δ11 from 182) and **E219D** (Δ13 from 206). Neither lies on the Plk1 phosphothreonine positions.

### *D. unisexualis*

No *D. unisexualis* Emi2 sequence is present in this repository. Heterozygosity in exon 1 is documented by PCR in the manuscript draft; **C-terminal exons encoding the Plk1 sites and the APC/C interface have not been sequenced here.**

---

## 3. Answer to the co-author (Question 3)

**Is the worry justified?** Partially — but the framing in the Introduction is what will draw reviewer fire, not the analysis itself.

| Co-author concern | Status in this project |
|---|---|
| “We only looked at the F-box” | **Incorrect as a description of the work.** Full-length Emi2 (675 aa) was compared; interface distances used Emi1/4UI9 transferred to Emi2; D-box modelled with AF3 (60 models); ZBR superposed (RMSD 0.19 Å); APC/C subunits sequenced where available. |
| “APC/C binds via C-terminus and zinc finger — why not check?” | **Checked.** ZBR 600–648 invariant in structure; RL tail 658–675 has **no atomic coverage** in any structure — cannot be tested geometrically. |
| “Are those regions invariant?” | **Between Dval and Dnai:** yes for ZBR structure; yes for Plk1 sites T182/T206; D-box has F549I but side chain does not contact Cdc20; Apc10/Apc11 **identical** on APC/C side. |
| “Will reviewers object?” | **Yes, if the Introduction still claims F-box → APC/C.** No, if we state the correct elements and report the two-sided null result honestly, including Apc1 and RL-tail limits. |

**What remains open:** (1) **Apc1** sequence — partner of Emi2 590; (2) **Apc4** — contact subunit, Dnai missing; (3) **RL tail** — essential but unmodelled; (4) **regulatory layer** — Plk1/CaMKII/β-TrCP/PP2A-B56, where eight of ten substitutions lie outside structural coverage.

---

## 4. References (paste-ready)

Chang L, Zhang Z, Yang J, McLaughlin SH, Barford D. Atomic structure of the APC/C and its E3 module in complex with CDH1 and SUB1. Nature. 2015;522(7557):450-454.

Frye JJ, Peter M, Glotzer M, Morgan DO. Spindle checkpoint protein hBub3 and APC/C coactivator Cdc20 compete for the N terminus of the APC/C subunit Cdc27 during the cell cycle. Proc Natl Acad Sci U S A. 2013;110(15):E1380-E1387. *(ZBR/Emi1 context — verify exact citation for ZBR–Apc11)*

Fujita MK, Moritz C. Origin and evolution of parthenogenetic genomes in lizards: current state and future directions. Cytogenet Genome Res. 2009;127(2-4):261-272.

Gopinathan L, Szmyd R, Low D, et al. Emi2 Is Essential for Mouse Spermatogenesis. Cell Rep. 2017;20(3):697-708.

Hansen D, Alexander J, Bishop JM. MAPK and CaMKII do not phosphorylate Xenopus Emi2 directly. *(verify full citation)*

Inoue D, Ohe M, Kanemori Y, Sagata N. CaMKII-dependent phosphorylation of Erp1 is required for cytostatic factor arrest in Xenopus laevis egg extracts. *(verify — cited in Jia 2015 for Thr170/195)*

Jia JL, Han YH, Kim HC, et al. Structural basis for recognition of Emi2 by Polo-like kinase 1 and development of peptidomimetics blocking oocyte maturation and fertilization. Sci Rep. 2015;5:14626. doi:10.1038/srep14626

Liu J, Grimison B, Lewellyn AL, Maller JL. The Anaphase-promoting Complex/Cyclosome Inhibitor Emi2 Is Essential for Meiotic but Not Mitotic Cell Cycles. J Biol Chem. 2006;281(46):34736-34741.

Ohe M, Kawamura Y, Ueno H, et al. Emi2 inhibition of the anaphase-promoting complex/cyclosome absolutely requires Emi2 binding via the C-terminal RL tail. Mol Biol Cell. 2010;21(6):905-913.

Sako K, et al. Emi2 RL tail and Ube2S competition on Apc2. *(verify full citation — 2014)*

Zhang Z, Yang J, Barford D. Mechanism of APC/CCDC20 activation by the mitotic checkpoint complex. Nature. 2016;533(7602):260-264. *(PDB 5G04)*

---

## 5. Manuscript actions

1. Replace Introduction closing paragraph (F-box → APC/C).
2. Add one sentence naming **Cdc20 + Apc10** (D-box), **Apc11 + Apc2 + Apc1** (ZBR), **Apc2 CTD** (RL tail).
3. In Results, note that exon-1 PCR variation is **expected** in the F-box/Skp1 module and is **not** a test of the docking hypothesis.
4. Optionally add: Plk1 sites **T182/T206** (Darevskia) are **identical** between parental species — supports regulatory-line exclusion for this triad at the Plk1 layer.
