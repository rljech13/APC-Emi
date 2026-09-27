# Emi2 (FBXO43) regulatory-region analysis of the ten *D. valentini* / *D. r. nairensis* substitutions

Analysis directory: `SNP_regulatory_analysis/`. Nothing outside this directory was modified.

## 0. Sanity checks on the established results

| Claim checked | Result |
|---|---|
| Both sequences 675 aa, differing at exactly 10 positions | **Confirmed.** Q76H, V92I, A102T, K171R, E219D, I391V, D395N, L425F, F549I, V590M |
| β-TrCP degron DSGxxxS at ~40–46 | **Confirmed.** D40-S41-G42-Y43-N44-G45-S46 |
| D-box R548-x-x-L551 | **Confirmed.** R548-F549-A550-L551; F549I sits at D-box position +1 |
| ZBR essential Cys = Darevskia C607 (Xenopus C583) | **Confirmed** by alignment; C607 is fully buried (relSASA 0.0%) and coordinates Zn with C610/C612/C625/C630 in the AF3 model |
| RL tail ends Arg674-Leu675 | **Confirmed** |
| L425F maps to human Leu449 | **Confirmed** independently by MAFFT L-INS-i alignment of 56 orthologues |
| Human Emi2 = UniProt Q8NHZ8 | **Incorrect.** Q8NHZ8 is CDC26 (APC/C subunit 26). Human FBXO43/Emi2 is **Q4G163** (708 aa). Mouse is Q8CDI2; Xenopus Q8AXF4 is correct. The local `Emi2_Hsap_AAseq.fasta` does contain the right 708-aa Q4G163 sequence, so only the accession cited is wrong — worth fixing in the manuscript. |

Sequence-quality control: the study's own *Eublepharis macularius* Emi2 is **identical** to RefSeq XP_054841747.1 (0 differences), and the study's *Podarcis muralis* differs from RefSeq XP_028592231.2 at only 5 positions plus one indel. The in-house pipeline therefore reproduces reference sequences well, and the ten Darevskia differences are unlikely to be artefacts.

---

## (a) Regulatory-element map in Darevskia numbering

Built by MAFFT L-INS-i alignment of 56 vertebrate Emi2/FBXO43 orthologues (44 squamates including 14 lacertids, plus human, mouse, rat, chicken, alligator, green sea turtle, *Xenopus*, zebrafish). Every mapping below was cross-validated by direct motif matching in the raw sequences, and the two independent routes (via *Xenopus* and via mouse) agree exactly.

| Element | Darevskia | Mapped from | Source |
|---|---|---|---|
| **β-TrCP phosphodegron DSGxxxS**, Plk1-phosphorylated Ser | **40–46** (pS41, pS46) | Xen D32-S33…S38 (UniProt: "S33, phosphoserine; by PLK1") | Rauh 2005; Hansen 2006 |
| **Plk1 polo-box-domain docking module** | **176–207** | mouse Emi2 146–177 (Jia 2015, PDB 5DMV/5DMZ) | Jia 2015 |
|  · PBD phospho-Thr 1 | **T182** | mouse pT152 = Xen T170 | Jia 2015; Isoda 2011 |
|  · PBD acidic-patch contacts | **D187, P188, I189** | mouse D157/V158/V159 | Jia 2015 |
|  · PBD Tyr-pocket anchor | **F199** | mouse F169 (F169A almost abolishes Emi2 degradation) | Jia 2015 |
|  · PBD phospho-Thr 2 / **CaMKII site** | **T206** | mouse pT176 = Xen T195 (UniProt: "T195, phosphothreonine; by CaMK2") | Jia 2015; Rauh 2005; Hansen 2006 |
| **PP2A-B56 basic patch** | ~343–358 (K346,R347,K348,R350,K351,R354,R356,R357) | position −15 to −1 of the LxxIxE motif | Wang 2020 |
| **PP2A-B56 LxxIxE motif** | **L358-x-x-L361-x-E363** | Xen L334-S335-T336-L337-R338-E339; the Hertz 2016 Emi2 B56-binding peptide is "LSTLREQSSQS", K_D 41 µM | Hertz 2016; Wu 2007 |
| **p90rsk / Mos-MAPK sites** | **S359, T360** (inside the LxxIxE motif), **S366, S368** | Xen S335/T336, S342/S344 | Nishiyama 2007; Inoue 2007 |
| **F-box domain** (the folded domain) | **452–537** | human 490–547; Xen 424–499 | UniProt |
|  · **Skp1-binding F-box motif**, inside the domain | **479–519** | the canonical 41-residue F-box motif | F-box motif length; drawn nested in Figs 2, 4 |
| **D-box R-x-x-L** | **548–551** | Xen R529/L532 | Ohe 2010 |
| **Cdk1 site T-P** | **T574-P575** | Xen T551-P552 | Wu 2007; Isoda 2011 |
| **ZBR** (8 Cys: 607,610,612,625,630,635,638,648) | **603–651** | Xen 579–627 | Ohe 2010; Shoji 2014 |
| **post-ZBR region** (binds ANAPC2) | **652–661** | Xen 628–638 | Shoji 2014 |
| **RL tail** (obligatory APC/C docking) | **662–675** | Xen 629–651 | Ohe 2010 |

The coordinates in this table are as the cited sources give them, mapped onto *D. valentini*. The F-box entry is nested rather than a single range, and the figures draw it that way: the folded domain 452–537 with the canonical Skp1-binding motif 479–519 inside it. The enrichment test below counts the domain, because a *footprint* is a footprint of annotated domains. Two coordinates in this table still differ from the ones the figures use — the figures take ZBR 600–648 and a merged post-ZBR + RL tail 652–675 — and those are discussed under the enrichment test below and in `Interface_variability/REPORT.md` section 0.

Two notes. First, the Jia et al. 2015 sites resolve cleanly: that paper uses **mouse** Emi2 (Q8CDI2; PDB 5DMV DBREF confirms `FBX43_MOUSE 146-177`), and its pThr152/pThr176 are the *same two sites* as the classical *Xenopus* Thr170/Thr195, i.e. the CaMKII-primed Plk1 docking sites. Second, the Xenopus Cdk1 site Thr545 aligns to Darevskia His564 and is not conserved in lizards; the second Cdk1 site Thr551 is conserved as T574-P575.

---

## (b) Per-substitution table

`FI` = FoldIndex (Prilusky 2005) on the Darevskia sequence; `pLDDT` = AlphaFold-DB full-length models of human Q4G163 / Xenopus Q8AXF4 projected through the alignment (values < 50 indicate disorder).

| Sub | Nearest mapped element (distance) | Region | FI / pLDDT (hum/xen) | Δcharge | Δvolume (Å³) | Consensus created or destroyed |
|---|---|---|---|---|---|---|
| **Q76H** | β-TrCP degron, **30 aa** | disordered N-terminal linker, no known function | −0.20 / 38, 41 | +0.1 | +9 | none |
| **V92I** | β-TrCP degron, **46 aa** | disordered, no known function | +0.02 / 33, 34 | 0 | +27 | none |
| **A102T** | β-TrCP degron, **56 aa**; Plk1-PBD module 74 aa | disordered, no known function | +0.07 / 39, 33 | 0 | +28 | creates a Thr with acidic residue at n−3 (weak CK1-type context); **no S/T-P, no basophilic site, no degron** |
| **K171R** | Plk1-PBD docking module, **5 aa** (11 aa from pT182) | disordered, no known function | 0.00 / 40, 42 | 0 | +5 | **K-x-x-S → R-x-x-S at Ser174** (basophilic CaMKII/PKA/RSK consensus upgraded). Ser174 is not a documented phosphosite in any species, and 18/52 vertebrates including mouse, rat and chicken carry Lys here |
| **E219D** | Plk1-PBD docking module, **12 aa** | disordered, no known function | +0.06 / 37, 35 | 0 | −27 | none |
| **I391V** | p90rsk cluster, **23 aa**; LxxIxE 28 aa | disordered, no known function | −0.07 / 38, 44 | 0 | −27 | none |
| **D395N** | p90rsk cluster, **27 aa**; LxxIxE 32 aa | disordered, no known function | +0.01 / 30, 37 | **+1** (only charge change in the dataset) | +3 | none |
| **L425F** | F-box domain, **27 aa** | short isolated α-helix (416–429), solvent-exposed, no known function | +0.01 / **79, 77** | 0 | +23 | none |
| **F549I** | **inside the D-box** (R548-x-x-L551) | folded/marginal, at APC/C interface | +0.10 / 65, 63 | 0 | −23 | D-box preserved (RFAL → RIAL); already shown to be solvent-exposed and to have no ipTM effect |
| **V590M** | ZBR, **13 aa**; Cdk1 T574-P 15 aa | D-box–ZBR linker | +0.11 / **90, 91** | 0 | +23 | none |

**Global phospho-regulation checks (all negative):**

- All **15 proline-directed Cdk1/MAPK S/T-P sites are identical** between the two species (positions 19, 34, 84, 109, 169, 224, 234, 258, 273, 291, 304, 416, 574, 579, 650). Not one is created or destroyed.
- No Ser or Thr is lost; A102T adds one Thr (Dnai has 133 S+T vs 132 in Dval).
- No polo-box S-pS/pT-P consensus and no DSG-type degron is created or destroyed.
- Net charge differs by exactly one unit (D395N), i.e. +13 (Dval) vs +14 (Dnai).
- The **entire Mos-MAPK → p90rsk → PP2A-B56 module is invariant**: the LxxIxE motif (L358/L361/E363), the p90rsk sites (S359, T360, S366, S368) and the upstream B56 basic patch (343–358) are identical in the two species.
- The **entire Plk1-PBD docking module (176–207) is invariant**, including both docking threonines (T182, T206) and the Phe199 hydrophobic anchor whose mutation is the single most damaging lesion identified by Jia et al.
- The **β-TrCP degron (40–46) is invariant**.

**Enrichment test.** The mapped functional elements cover 215/675 residues (31.9% of the protein). One substitution (F549I) falls inside them; 3.19 are expected by chance. A 200 000-draw permutation test gives P(≥1 by chance) = 0.979 — the substitutions are, if anything, mildly *depleted* in functional elements, and certainly not enriched.

These figures use the F-box counted as the folded domain 452–537, together with ZBR 600–648 and post-ZBR and RL tail merged as 652–675 to stop the two overlapping. Counting the 41-residue Skp1-binding motif 479–519 in place of the domain instead gives a 170/675 footprint (25.2%) and 2.52 expected, P = 0.946. The observed count is 1 either way and the direction of the result does not depend on the choice.

**Resolved: the F-box is nested, and both ranges are now shown.** The project's two F-box definitions are not in conflict — they describe different objects, and the figures now draw one inside the other. The **folded domain is 452–537**. Two independent lines give that range: the UniProt F-box of human Emi2 (490–547) maps here to 452–537, and in the AF3 model of *D. valentini* the contiguous segment with Cα pLDDT ≥ 70 runs 454–535 (mean 92.4) — the confidently folded domain matches the wider range almost exactly at both ends. Inside it lies the **canonical 41-residue Skp1-binding F-box motif, 479–519**. Recomputed from the model's B-factor column as all-atom means (the convention used in Figure 2c), the domain scores 86.7 and the motif 93.4: the motif is the ordered core and the domain carries the frayed ends with it. The footprint above counts the domain; the figures label both and mark the motif as nested. No substitution falls in either range, so no conclusion in this report turns on the choice — but the two are not interchangeable and should not be quoted as if they were. Note also what the F-box does: it binds Skp1 as part of an SCF-type complex, not APC/C. The ZBR boundary is corroborated in the same way: the second ordered segment is 600–647 (Cα mean 77.3), and the last of the eight zinc-coordinating cysteines is Cys648.

---

## (c) Cross-vertebrate conservation at the ten positions

52 orthologues (44 squamates + 8 outgroups), Darevskia excluded from the tallies. "Percentile" ranks the column against all 670 scorable Emi2 columns — the median Emi2 column is 84% identical.

| Pos | Dval | Dnai | Residue distribution across 52 vertebrates | Column identity (percentile) | Class conserved? | Derived allele |
|---|---|---|---|---|---|---|
| 76 | Q | H | H:43, Y:3, D/G/S/R/E:1 each | 0.84 (50th) | **No** — 86% basic, Dval Gln unique | **Dval** |
| 92 | V | I | I:21, A:15, T:9, S:3, K:3, N:1 | 0.40 (4th) | yes | **Dval** |
| 102 | A | T | I:15, T:14, L:7, V:6, K:3, P:2, A:1, E/G/Q:1 | 0.29 (0th — least conserved decile) | no | **Dval** |
| 171 | K | R | R:33, K:18, D:1 | 0.63 (24th) | yes (98% basic) | **Dval** |
| 219 | E | D | **E:51, D:1** (only *Elgaria multicarinata*) | **0.98 (83rd)** | yes (100% acidic) | **Dnai** |
| 391 | I | V | I:39, V:4, L:3, S:2, A/T/F:1 | 0.76 (38th) | yes (92% aliphatic) | **Dnai** |
| 395 | D | N | N:29, S:14, Q:3, D:2, H:2, G/R:1 | 0.56 (17th) | no — 88% polar, Dval Asp rare | **Dval** |
| 425 | L | F | **L:27, I:18, V:6, M:1 — aliphatic in 100%, zero aromatics** | 0.52 (11th) | **No** — Dnai Phe is the only aromatic in the panel | **Dnai** |
| 549 | F | I | F:40, S:6, L:6 | 0.77 (38th) | **No** — Dnai Ile absent from panel | **Dnai** |
| 590 | V | M | A:27, V:13, E:5, Q:3, S/T/K/M:1 | 0.52 (11th) | yes | **Dnai** |

**Eight of the ten substitutions sit below the median conservation of Emi2**, and three (92, 102, 425/590) are in the bottom decile. Only **E219D** falls at a well-conserved position — and it is the chemically most conservative change possible (Glu→Asp, one methylene shorter, identical charge, no consensus altered, in a predicted-disordered segment).

**Polarity.** Against a 14-species lacertid outgroup (12 *Podarcis*, *Lacerta agilis*, *Zootoca vivipara*), the derived allele is carried by *D. valentini* at 5 positions and by *D. r. nairensis* at 5 — an exactly clock-like split with no lineage-specific excess.

**Divergence calibration.** Among the 66 pairwise comparisons of 12 congeneric *Podarcis* species, Emi2 differs by a median of 8 amino acids (range 3–14); 35% of pairs differ by ≥10. The 10 differences between *D. valentini* and *D. r. nairensis* are therefore **entirely typical background divergence for two congeneric lacertids**, not an unusual accumulation.

---

## (d) L425F assessment

**Annotated function: none.** Position 425 has no annotation in any species. Human Leu449 lies in the gap between the UniProt "Disordered" region (human 320–426 = Darevskia 292–403) and the F-box domain (human 490–547 = Darevskia 452–537). It is 27 residues N-terminal of the F-box and carries no modified-residue, binding-site or motif annotation in human, mouse or *Xenopus*. It is not part of the D-box, ZBR, RL tail, post-ZBR ANAPC2-binding region, the Plk1-PBD module, the degron, or the PP2A-B56 module.

**Folded or disordered: a short, isolated, solvent-exposed helix.** Full-length AlphaFold-DB models give pLDDT 79 (human) and 77 (*Xenopus*) — ordered, and markedly higher than the 30–44 seen at the other seven N-terminal substitutions. A CA(i)–CA(i+4) scan places it inside a discrete α-helix spanning Darevskia 416–429 (human 440–457), flanked by disorder on both sides. So it is a genuine secondary-structure element, but an isolated one.

**Burial: none.** In the study's AF3 Emi2(400–675) model the Leu425 side chain has 62.7 Å² SASA (31% relative); in the Dnai model Phe425 has 118.2 Å² (49% relative); in the full-length human AFDB model Leu449 is 62% exposed. Its only neighbours within 5 Å are its own helix partners (Ile421, Val422, Glu424, Phe426, Gln428, Asn429). A single additional contact to Leu451 appears in the Dval model but not the Dnai model — at pLDDT ~63 that difference is not interpretable. **The residue makes essentially no tertiary packing in any model.**

**Conservation: the most class-breaking substitution in the dataset.** Position 425 is aliphatic in 52/52 vertebrates (Leu 27, Ile 18, Val 6, Met 1) with no aromatic residue anywhere in the panel. *D. r. nairensis* Phe425 is unique. Against that, the column is only 52% identical overall (11th percentile), so the position tolerates substitution freely *within* the aliphatic class.

**Verdict on L425F.** It is the most chemically radical of the eight non-C-terminal substitutions and it breaks an otherwise invariant residue class, so it is the one worth mentioning. But it is solvent-exposed, unpacked, functionally unannotated, sits on an isolated helix 27 residues outside the nearest domain, and I421/V422/L425/F426 form a small amphipathic face that would only matter if the helix is a binding element — for which there is no evidence in any species. This supports, at most, a sentence flagging it as a hypothesis-generating observation. It does not support a mechanism.

---

## (e) Verdict and recommended framing

### No defensible alternative mechanism exists in these data.

Every documented regulatory element of Emi2 that could plausibly control its stability, activation or timing of degradation is **identical between the two species**:

- the β-TrCP phosphodegron and both Plk1-phosphorylated serines;
- the complete Plk1 polo-box docking module, including both CaMKII/Plk1 docking threonines and the Phe199 anchor;
- all 15 Cdk1/MAPK proline-directed sites across the whole protein;
- the p90rsk/Mos-MAPK sites, the PP2A-B56 LxxIxE motif, and the B56 basic patch;
- the F-box, the ZBR, the post-ZBR ANAPC2-binding region and the RL tail.

The eight non-C-terminal substitutions fall in regions with **no known function**, in predicted-disordered sequence (pLDDT 30–44 for seven of the eight), at positions that are mostly among the *least* conserved in the protein. They are not enriched in functional elements (1 observed vs 3.2 expected, P = 0.98). Their number and their 5:5 lineage split are exactly what two congeneric lacertids are expected to accumulate neutrally.

I looked hardest at the two candidates that could conceivably have carried a story, and both fail:

- **K171R** is the only substitution that changes a kinase consensus in a way one could write up — it upgrades Lys-x-x-Ser174 to the canonical basophilic Arg-x-x-Ser174, and it is the closest of the eight to a real element (5 residues from the crystallised Plk1-PBD peptide, 11 from pThr182). But Ser174 is not a phosphosite in any species; Lys and Arg are both acceptable at the −3 position; and 18 of 52 vertebrates — including mouse, rat and chicken — carry the "unfavourable" Lys. The *variant* is standing vertebrate polymorphism, not an innovation.
- **E219D** is the only substitution at a strongly conserved position (Glu in 51/52 vertebrates) but is the most conservative possible chemical change, alters no consensus, and lies in disorder.

Claiming a regulatory mechanism from any of this would require asserting function for an unannotated serine in a disordered loop, or biological meaning for a Glu→Asp change. Neither is defensible in *Scientific Reports* or anywhere else.

### Recommended framing: a well-controlled negative result with a positive methodological core

Do not reframe the paper around a second mechanism. Reframe it as **a systematic test of whether the reproductive biology of a hybrid parthenogen is reflected in the sequence of its APC/C-inhibitor, with a clear negative answer and a mechanistic explanation for why the answer is negative.** That is a legitimate, publishable, and unusually complete story, because the negative is *explained* rather than merely reported.

Suggested structure of the claim:

1. **The comparative dataset.** Emi2 and the Emi2-contacting APC/C subunits were sequenced/assembled for both parental species of a hybrid parthenogen. Ten Emi2 substitutions; all Emi2-contacting APC/C residues identical.
2. **The obligatory APC/C-binding machinery is invariant.** The ZBR and the RL tail — the elements without which Emi2 cannot bind or inhibit APC/C (Ohe 2010; Schmidt 2005) — are 100% identical. Add the new results: the F-box, the post-ZBR ANAPC2-binding region, and the *entire* known phosphoregulatory architecture are also identical. This is a strong purifying-selection statement backed by 52 orthologues.
3. **The one polymorphic residue at a functional site is structurally silent.** F549I lies at D-box position +1, but the validated AF3 model (0.5–0.6 Å backbone RMSD to the Hsl1 degron of 5G04) shows residue 549 is 43–49% solvent-exposed and buries only 14–20 Å²; the computed effect (+0.005 ipTM) is below the 0.023 empirical noise floor. Report the noise floor explicitly — a stated, measured null threshold is a methodological contribution in itself.
4. **The divergence is clock-like, not selected.** Substitutions are depleted rather than enriched in functional elements (1 vs 3.2 expected, P = 0.98); eight of ten sit below the median conservation of Emi2; the derived alleles split 5:5 between lineages; and 10 differences is the median for congeneric *Podarcis* pairs (median 8, range 3–14).
5. **The biological interpretation.** Because *D. unisexualis* is heterozygous and carries both parental Emi2 alleles, the relevant question was whether the two alleles are functionally divergent at the APC/C interface. They are not. Emi2 sequence divergence is therefore **not** a contributor to the parthenogenetic phenotype, and the CSF-arrest machinery appears to be under strong purifying selection even across a hybridisation event that produces a viable unisexual lineage. This directs future work to expression level, allele-specific expression, timing, or upstream Ca²⁺/CaMKII signalling rather than to Emi2 protein sequence.

The honest one-line abstract claim: *the obligatory APC/C-inhibitory elements of Emi2 are invariant between the parental species of a hybrid parthenogen, the single polymorphic residue in a functional motif is solvent-exposed and structurally inert, and the remaining divergence is indistinguishable from neutral background — so allelic divergence in Emi2 is not a candidate mechanism for parthenogenesis in* Darevskia.

If you want to keep one forward-looking sentence, L425F is the only substitution that earns it: unique among 52 vertebrates at a position that is aliphatic in all of them, but solvent-exposed, unpacked and unannotated — worth flagging, not worth claiming.

---

## Files

| File | Contents |
|---|---|
| `emi2_orthologues.fasta` | 56 curated Emi2/FBXO43 orthologues, one per species |
| `emi2_orthologues_linsi.aln.fasta` | MAFFT L-INS-i alignment |
| `build_dataset.py` | orthologue retrieval and curation |
| `map_landmarks.py` | landmark projection + conservation tallies |
| `context_analysis.py` | distances to elements + kinase-consensus scan |
| `disorder.py` | FoldIndex + pLDDT tracks |
| `l425_structure.py` | SASA / packing / helix analysis at 425 |
| `stats.py` | enrichment permutation test, conservation percentiles, QC |
| `outputs/01–05*.txt` | saved output of each script |
| `AF_Q4G163_v6.pdb`, `AF_Q8AXF4_v6.pdb` | AlphaFold-DB full-length reference models |
| `squamata_fbxo43_raw.fasta`, `outgroups.fasta`, `Emi2_reviewed_all.fasta` | raw downloads |

## References consulted

Rauh NR et al. (2005) *Nature* 437:1048. Hansen DV, Tung JJ, Jackson PK (2006) *PNAS* 103:608. Nishiyama T, Ohsumi K, Kishimoto T (2007) *Nature* 446:1096. Inoue D et al. (2007) *Nature* 446:1100. Wu Q et al. (2007) *PNAS* 104:16564. Wu Q et al. (2007) *Curr Biol* 17:213. Ohe M et al. (2010) *Mol Biol Cell* 21:905. Isoda M et al. (2011) *Dev Cell* 21:506. Shoji S et al. (2014) *FEBS Open Bio* 4:689 (PDB 2RT9). Jia J-L et al. (2015) *Sci Rep* 5:14626 (PDB 5DMV, 5DMZ). Hertz EPT et al. (2016) *Mol Cell* 63:686. Wang X et al. (2020) *eLife* 9:e55966. Prilusky J et al. (2005) *Bioinformatics* 21:3435 (FoldIndex). Tien MZ et al. (2013) *PLoS ONE* 8:e80635 (max ASA).
