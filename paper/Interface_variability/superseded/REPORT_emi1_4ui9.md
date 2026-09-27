# Two-sided variability test: D. valentini x D. r. nairensis

## 0. Emi2 element boundaries used throughout

These are the authoritative coordinates for this project, in *D. valentini* Emi2 numbering (675 aa). Figures 2, 2c, 3 and 4, `two_sided_test.py` and the enrichment test in `SNP_regulatory_analysis/stats.py` all use them, with the one exception noted at the end of this section:

| element | residues | binds |
|---|---|---|
| F-box domain (the folded domain) | 452-537 | Skp1 (SCF), **not** APC/C |
| · Skp1-binding F-box motif, **inside** the domain | 479-519 | Skp1 (SCF), **not** APC/C |
| D-box | 548-551 | Cdc20 WD40 + Apc10 |
| ZBR | 600-648 | Apc11 RING, Apc2, Apc1 |
| RL tail | 658-675 | Apc2 CTD (no atomic coverage available) |

The F-box is a nested pair, not a single range, and the figures draw it that way: 452-537 is the folded domain, delimited from the model itself by its contiguous C-alpha pLDDT >= 70 segment 454-535 (mean 92.4), and 479-519 is the canonical 41-residue Skp1-binding motif lying inside it. The UniProt F-box annotation of human Emi2 (FBXO43, Q4G163, F-box 490-547) is a third, intermediate range: aligned onto *D. valentini* at 51.2% identity it maps to 463-520, so it corroborates the nesting without coinciding with either boundary. An earlier version of this file, and of the figure captions, wrongly gave 452-537 as the image of that annotation.

Panels 2a and 4a draw the domain as a pale box with the motif solid inside it. Panel 2b carries no text labels at all: since it is coloured by region in the same palette as 2a, the nesting reads as a pale domain with a saturated motif inside it. Panel 2c lists both, the motif as an indented sub-row of the domain (all-atom mean pLDDT 86.7 for the domain, 93.4 for the motif). Where a single number is needed - the functional-element footprint - the domain is the unit used.

**Open, deliberately deferred.** The separate figure `Figure_2c_v3` shows the same model coloured by functional element in its left panel, but colours only the motif 479-519 and leaves the folded domain grey. Two panels of the same series therefore draw the F-box under different definitions, and the left panel is otherwise largely a differently-posed duplicate of panel 2b; its only unique content is the Phe549 and Leu425 spheres. Left as is by decision, to be resolved before submission.

An earlier set (ZBR 603-651, RL tail 662-675) is still present in the superseded scripts `fig2b_v2.py` and `measure_v2.py`. For the substitution set analysed here the differences are inconsequential: none of the ten sites falls in the symmetric difference of the old and new ranges (600-602, 649-651, 658-661), nor between the two F-box ranges (452-478, 520-537), so every `element` assignment below is identical under any of these choices. The one place the F-box choice does propagate is the enrichment test in `SNP_regulatory_analysis`, whose footprint counts the domain: the union of annotated elements is 215 of 675 residues with 3.19 substitutions expected by chance (P = 0.979 for the 1 observed). Counting the motif instead would give 170 of 675 and 2.52 expected (P = 0.946). The observed count stays at 1 and the direction of the result - no enrichment of substitutions in annotated elements - is unchanged.

`two_sided_test.py` still lists the F-box as 479-519 alone. Its output is unaffected, because no substitution lies in 452-478 or 520-537, but the entry should be brought into line with the nesting when that script is next touched.

## 1. Emi2 side

Dval 675 aa vs Dnai 675 aa: **10 substitutions**.

| position | Dval | Dnai | functional element |
|---|---|---|---|
| 76 | Q | H | - |
| 92 | V | I | - |
| 102 | A | T | - |
| 171 | K | R | - |
| 219 | E | D | - |
| 391 | I | V | - |
| 395 | D | N | - |
| 425 | L | F | - |
| 549 | F | I | D-box |
| 590 | V | M | - |

Substitutions inside an APC/C-engaging element: **1** (549)

## 2. APC/C side

Distance is the minimum heavy-atom distance to Emi1 (chain S, 319-436) in 4UI9.

| subunit | Dval len | Dnai len | substitutions | distance to Emi1 | role |
|---|---|---|---|---|---|
| Cdc20 | 454 | 454 | V291I | n/a | coactivator; D-box receptor (Cdh1 chain R in 4UI9) |
| Anapc10 | 185 | 185 | **0 — identical** | — | D-box co-receptor |
| Anapc11 | 84 | 84 | **0 — identical** | — | RING; ZBR target |
| Anapc2 | 767 | 568 | V632A; P654L | 19.5 A; 36.6 A | cullin; ZBR body + RL tail |
| Anapc5 | 743 | 722 | Y728H | n/a | not an Emi2 partner (negative control) |
| Anapc7 | 565 | 565 | **0 — identical** | — | not an Emi2 partner (negative control) |
| Anapc1 | — | — | *not sequenced* | 3.07 A | ZBR partner; **untested** |
| Anapc4 | 394 (partial) | — | *Dnai missing* | 4.03 A | contacts Emi1 432-436 (see section 4); **untested** |

Cdc20 carries one substitution (V291I). 4UI9 contains Cdh1 rather than Cdc20, so this one is measured in our own 60 Emi2-Cdc20 models: **18.6-31.8 A from the D-box** (min-median). Distances to the nearest Emi2 atom of any kind reach 8.1 A, but those approaches involve the disordered Emi2 segments whose placement is not reliable; the D-box figure is the meaningful one.


Positive control — published Emi1-contacting APC2 residues, same measurement:

| hAPC2 residue | 513 | 514 | 517 | 520 | 552 | 553 | 556 | 602 |
|---|---|---|---|---|---|---|---|---|
| distance to Emi1 | 2.85 A | 4.38 A | 2.65 A | 2.59 A | 2.98 A | 3.62 A | 3.20 A | 4.29 A |

## 3. Intersection — the parthenogenetic substitution list

APC/C variable sites lying on the Emi2-binding surface (<=5 A from Emi1): **0**.

The list of candidate parthenogenetic substitutions is the set of Emi2 substitutions that both sit in an interface and contact a variable APC/C position. For this pair the list is **empty**.

## 4. Full interface census — which subunits the inhibitor actually touches

Per-residue minimum heavy-atom distance from Emi1 (chain S) to every APC/C subunit in 4UI9, all 15 subunits scanned without assuming a candidate list. Chain identity was taken from `_atom_site.label_entity_id` (the `_struct_asym` table in this entry is misaligned and assigns zinc to 484-residue chains). Rendered as `Figures/Figure_4_interface_map.svg`; matrix archived as `emi1_4ui9_distance_matrix.json`.

| subunit | min. distance | residues <5 A | Emi1 segment engaged |
|---|---|---|---|
| Cdh1 / Cdc20 | 1.57 A | 8 | 321-328 (D-box) |
| Anapc2 | 2.59 A | 15 | 379-426 (ZBR body) |
| Anapc10 | 2.66 A | 9 | 320-331 (D-box + following turn) |
| Anapc11 | 2.67 A | 16 | 357-374 |
| Anapc1 | 3.07 A | 9 | 356, 361, 381, 398-403 |
| **Anapc4** | **4.03 A** | **4** | **432-436** |
| Anapc5 | 25.49 A | 0 | — |
| Anapc3/Cdc27 | 28.40 A | 0 | — |
| Anapc16 | 31.22 A | 0 | — |
| Anapc13 | 32.94 A | 0 | — |
| Anapc7 | 34.06 A | 0 | — |
| Anapc6/Cdc16 | 39.40 A | 0 | — |
| Anapc15 | 41.17 A | 0 | — |
| Anapc8/Cdc23 | 43.06 A | 0 | — |
| Cdc26 | 53.00 A | 0 | — |

Two results here revise earlier statements in this project.

**Apc4 is a contact, not a negative control.** It reaches 4.03 A of Emi1 432-436, the last modelled residues of the fragment, immediately upstream of where the RL tail continues. It was previously treated as an unrelated subunit. The contact is carried by hAPC4 F727, R728, K729, E745, D747 and I748 (4.0-7.1 A), with a peripheral patch at D33/L49/A50 at 6-8 A that does not qualify as contact.

Unlike Apc1, Apc4 is immediately testable: our *D. valentini* Apc4 is partial (394 aa, spanning human 411-804) but that span **covers all six contact residues**. Closing this gap needs only *D. r. nairensis* Apc4 over the same region, and the informative window is as narrow as human 720-760.

**Apc3/Cdc27 is 28.4 A away and contacts nothing.** This settles, on structural grounds, the proposal that the Emi2 RL tail docks on Apc3 by analogy with the Apc10 IR tail: over the whole resolved inhibitor there is no approach to Apc3.

Coverage limits that bound all of the above: Emi1 is modelled only over 319-436, with an internal break at **332-355** (linker unresolved), and the **C-terminal LRRL tail is absent from the model entirely**. Chains T and U are 21- and 24-residue poly-alanine stubs annotated only as "peptide" and cannot be assigned to Emi1. The RL-tail interface therefore has no atomic coverage in any available structure, ours or published.

## 5. Transfer onto Emi2 numbering, and a second interface substitution

Sections 2-4 are in Emi1 numbering. Transferring them onto Emi2 requires a residue-level alignment, not a fixed offset: the offset is **226 at the D-box but 229 at the ZBR**, because Emi2 carries a short linker insertion. A global BLOSUM62 alignment of *D. valentini* Emi2 against human Emi1 (Q9UKT4) over the C-terminal modules gives 195 aligned pairs and three unambiguous anchors - the D-box (Emi1 R322-T-P-L325 to Emi2 R548-F-A-L551), the ZBR cysteine ladder, and the tail itself (Emi1 KKNLRRL 441-447 to Emi2 KRNLKRL 669-675).

The modelled fragment maps to **Emi2 545-664: 120 of 675 residues**. Everything N-terminal of 545, the whole F-box included, has no structural coverage, as does 558-584 (the unresolved break) and 665-675 - the RL tail. Rendered as `Figures/Figure_4_interface_map.svg`.

This transfer corrects an earlier statement in this file. Emi2 590 was described as untestable on the assumption that it fell in the unresolved break; the alignment instead places it at **Emi1 361, which is modelled and lies 4.37 A from Apc1**. So the intersection in section 3 is not simply empty:

| Emi2 site | maps to | contacts | status of the partner |
|---|---|---|---|
| 549 F->I | Emi1 323 | Cdh1/Cdc20 3.78 A, Apc10 4.56 A | both **invariant**; and the 549 side chain does not reach Cdc20 in any engaged model |
| 590 V->M | Emi1 361 | **Apc1 4.37 A**, Apc11 7.25 A | Apc11 invariant; **Apc1 not sequenced** |

The 590 assignment is robust to alignment error. Across Emi1 356-370 every position lies within 3.1-9.7 A of Apc1 or Apc11 (minima 3.84 and 3.06 A), so a shift of several residues still places Emi2 590 on the same composite surface. What it cannot do is decide which of the two subunits carries the contact.

The conclusion of section 3 therefore stands but narrows to a single testable question: **does Emi2 590 face a variable position on Apc1?** Apc11, the alternative partner, is identical between the two species. This makes Apc1 the one sequence whose absence still bears on the hypothesis, and Apc4 (section 4) the one that is cheap to close.

