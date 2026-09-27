# Draft — Methods and Results (gene inference + structural modelling)

Status: 2026-08-07. Integrated co-author gene-inference narrative with our AF3 structural
part, in approximately the same plain first-person-plural style. Placeholders left where
IDs / metrics / citations are still missing: `(GenBank ID)`, `(SRA ID)`, `(assembly ID)`,
`(cite)`, `(metrics?)`, `(Fig. SX)`.

Main-text structural figure: **Figure 2**. D-box allele screen → Supplementary Figure S1.
Partner-side maps → Supplementary Figure S2. Numbers for the modelling sections trace to
`Figures/CAPTIONS.md`, `AF3_results_S1/REPORT_AF3_S1.md`,
`AF3_results_S2b/REPORT_AF3_S2b_S3.md`, `Interface_variability/REPORT.md`,
`emi2_combined_distance_matrix.json`, and
`Interface_variability/partner_distance_to_emi2_muts.json`.

**Naming.** We keep the co-author usage *D. r. nairensis* / *D. raddei* as written in the
gene-inference passages; both refer to the bisexual *raddei* parent lineage in this triad.
Introduction remains parked at the end and is not expanded here.

---

## Methods

### Gene inference and sequence assembly

Genomic sequences of Anapc1, Anapc2, Anapc5, Anapc7, Anapc10, Anapc11, Cdc20 and Emi2 of
*Darevskia valentini* were obtained from the genome assembly GCA_034642135.1. The Anapc1
sequence of *D. r. nairensis* was taken from the genome assembly `(assembly ID)`. Sequences
of the same genes from *D. mixta* were obtained from the assembly `(assembly ID)`.

In order to annotate exons, we used protein sequences of the corresponding genes from
*Podarcis muralis* as queries in tBLASTn searches against the *Darevskia* assemblies. For
Anapc1, exon positions were further checked by BLASTn. Illumina whole-genome sequencing
(WGS) reads of *D. mixta* were used to recover Anapc11 `(SRA ID)`; reads were mapped with
BWA `(cite)`. RNAseq data `(SRA ID)` were used for exon annotation with STAR `(cite)`.
Variants were called with bcftools `(cite)`. Multiple sequence alignments were performed
with MAFFT and Clustal Omega `(cite)`. Sequence inspection and manual curation were done
in Benchling.

### Domain annotation

`(placeholder — domain annotation section to be filled by co-authors)`

### Molecular modelling of Emi2 and Emi2–APC/C contacts

Protein sequences of Emi2, Cdc20, Anapc10, Anapc11, Anapc2 and Anapc1 used for modelling
were taken from the *D. valentini* and *D. r. nairensis* ORFs assembled as above. Complex
screens used the four reciprocal allele combinations of Emi2 with the partner set (parental
and both hybrid pairings), unless a subunit was identical between species.

Before the AlphaFold3 complex screens, we explored global docking of hybrid pairings
(APC/C of one parent with Emi2 of the other). The unstructured terminal region was removed
from the Emi2 model using Maestro (Schrödinger) `[ref: Schrödinger Release 2024-1: Maestro]`.
Global protein–protein docking of the prepared Emi2 model from *D. valentini* with the APC/C
complex from *D. r. nairensis*, and the reciprocal pairing, was performed with ClusPro
`[ref: Ashizawa et al. 2026, Proteins 94:183–191; DOI 10.1002/prot.70066]`. In each run,
70,000 mutual orientations of the receptor and ligand were analysed. A cluster for
structural analysis was selected based on the number of models it contained and the value
of the evaluation function. Docking was performed in two independent runs. Reproducibility
was assessed from the similarity of the interface geometry. ClusPro was used only as an
exploratory survey of hybrid interfaces. All distances, engagement calls and dual-side
tests reported below come from the AlphaFold3 screens, not from ClusPro poses.

We further modelled the *D. valentini* Emi2 segment 400–675 with AlphaFold3
`[ref: Abramson et al. 2024]` to obtain a monomer view of the APC/C-facing half of the
protein (Fig. 2a,b). Domain boundaries used throughout were: F-box domain 452–537, Skp1
motif 479–519, D-box 548–551, ZBR 600–648, RL 658–675. Mean per-atom pLDDT and 10th–90th
percentiles (Fig. 2b) were computed per annotated segment from the top-ranked model.

Three multimer screens were then run on the AlphaFold3 Server with default settings, twelve
seeds per screen (4 allele combinations × 3 seeds) and five diffusion samples per seed
(60 models each):

1. **S1** — Emi2(535–675) + full-length Cdc20 + Anapc10 (+ Zn where required).
2. **S2b** — Emi2(535–675) + Anapc2 catalytic fragment matched between alleles
   (*D. valentini* 220–767 / *D. r. nairensis* 1–548) + Anapc11 (+ Zn).
3. **S3** — Emi2(535–675) + Anapc1 residues 950–1250 (window covering the literature
   ZBR-facing homology 1050–1125 and ending before the *D. r. nairensis* deletion at
   1284–1330).

A model was called engaged for a given chain pair if the AlphaFold3 pair ipTM for that pair
was ≥ 0.50. Seed, not diffusion sample, is the unit of replication for allele contrasts.
D-box pose validation (S1) used superposition of modelled Cdc20 onto Cdc20 of PDB 5G04 and
backbone RMSD of Emi2 D-box atoms against the Hsl1 peptide in that structure
`[ref: Zhang et al. 2016]` (Supplementary Fig. S1).

For each engaged model, the minimum heavy-atom distance from every Emi2 residue in 535–675
to each partner chain was recorded. Fig. 2c shows the median of those minima across engaged
models of the relevant screen. Rows without engagement (Anapc10 in S1; Anapc1 in S3) are
reported as non-measurements. An Emi2 site was treated as an interface contact at the 4 Å
heavy-atom cut-off. For allele positions we separately checked whether any atom or only the
side chain reached that cut-off. A dual-side candidate required (i) an Emi2 parental
substitution at a tested contact and (ii) a parental difference on the facing APC/C surface
within the same modelled interface.

To place APC/C-side parental substitutions relative to the Emi2 substitutions present in the
modelled construct, we further computed for each residue of Cdc20, Anapc10, Anapc2 and
Anapc11 the minimum heavy-atom distance to the nearest Emi2 substitution site in the AF3
construct (positions 549 and 590). Distances were taken from engaged AlphaFold3 seeds and
summarised as the median across those seeds (Supplementary Fig. S2). Cdc20 and Anapc10
profiles used engaged S1 seeds; Anapc2 and Anapc11 used engaged S2b seeds. Anapc10
distances are retained for display but flagged as non-measurements, consistent with the
hatched Anapc10 row in Fig. 2c.

Figure 2 was assembled in Python (matplotlib). Molecular graphics for Fig. 2a were rendered
in ChimeraX.

---

## Results

### Inference of APC/C subunit and Emi2 coding sequences

Gene annotation for Anapc1, Anapc2, Anapc5, Anapc7, Anapc10, Anapc11, Cdc20 and Emi2 was
absent from the *Darevskia* assemblies used here. We therefore annotated exons by tBLASTn
using *P. muralis* proteins as queries, and refined exon boundaries with RNAseq (exons
supported by ≥100 reads). Coding sequences (CDSs) were inferred for all genes except
Anapc1. Comparison of *D. valentini* and *D. raddei* (*D. r. nairensis*) revealed the
mismatch patterns described below.

Anapc1 required special handling: the locus was split across contigs, leaving a 47 aa gap
`(metrics?)`. For *D. mixta*, CDSs were recovered by alignment to the parental CDSs in the
same manner as for Anapc1, except Anapc11, which was recovered from Illumina WGS. In each
case we retained the longest open reading frame. Notably, *D. unisexualis* yielded few Emi2
RNAseq reads `(metrics?)`, so a full hybrid Emi2 CDS was not inferred from that dataset
alone `(Fig. SX)`.

### Parental substitutions in Emi2

Full-length Emi2 of *D. valentini* and *D. raddei* (*D. r. nairensis*) differs at ten amino
acid positions `(Fig. SX / Fig. 2c)`. One of these, F549I, falls inside the destruction box
(D-box). The D-box of Emi2 is part of the APC/C-binding surface and has been implicated in
contacts that involve Cdc20 and the co-receptor Anapc10
`(cite: Suzuki et al. 2010, Development 137:3281; and/or related APC/C–D-box literature —
https://journals.biologists.com/dev/article/137/19/3281/44041/;
https://www.sciencedirect.com/science/article/pii/S2211546314000667#b0250)`. We therefore
treat F549I as a candidate interface substitution, while noting that Anapc10 geometry could
not be measured in our three-chain models (see below).

We further compared Emi2 of *D. mixta* with both parental alleles `(metrics?)`. Position 590
differs in all three species. The remaining parental differences and their distribution along
the protein are summarised in Fig. 2c.

### Survey of APC/C-side substitutions

In order to test whether parental Emi2 differences could face parental APC/C differences at
the same contact (a structural dual-side incompatibility), we surveyed substitutions in the
APC/C subunits that bind Emi2. Anapc10 and Anapc11 were identical between *D. valentini* and
*D. r. nairensis*. Anapc2, Anapc5, Cdc20 and Anapc1 carried differences `(list / Fig. SY)`.
None of these variable APC/C positions contact the variable Emi2 positions in the alignments
and models available to us; the structural evidence for that claim is given in the next
sections.

### Structural model of the APC/C-facing half of Emi2

To place the parental amino-acid differences relative to the APC/C-binding elements of Emi2,
we first asked what the C-terminal half of the protein looks like as a monomer and where that
monomer can be trusted (Fig. 2a,b).

AlphaFold3 modelling of *D. valentini* Emi2 residues 400–675 recovers the expected
architecture of an Emi-family inhibitor (Fig. 2a): a folded F-box domain (452–537) that
contains the Skp1-binding F-box motif (479–519); a short D-box (548–551); a zinc-binding
region (ZBR, 600–648); and the C-terminal RL tail (658–675). The F-box is a Skp1/SCF module
and is not part of the APC/C interface `[ref: Ohe et al. 2010; Chang et al. 2015]`.

Per-region mean pLDDT of the same monomer (Fig. 2b) shows that confidence is high only where
Emi2 has an autonomous fold: F-box domain 86.7, Skp1 motif 93.4, ZBR 72.0. The D-box (38.7),
the D-box–ZBR linker and the intervening spacers fall in the very-low band (<50). That
pattern is expected for short linear motifs that order on binding. Notably, APC/C contacts
cannot be read from the monomer alone: the degron and ZBR faces have to be resolved in
complex.

Of the ten parental substitutions, eight lie N-terminal of residue 535. Of the two that fall
in the APC/C-facing window, F549I sits inside the D-box (R548-x-x-L551) and V590M lies
between the D-box and the ZBR. We therefore restricted atomic complex modelling to Emi2
535–675 (dashed box in Fig. 2c), the segment that contains every annotated APC/C-binding
element and both interface-proximal substitutions. Residues 1–534, including the F-box and
eight parental substitutions, have no complex geometry in this study.

### Where Emi2 contacts APC/C in *Darevskia* models

We mapped median minimum heavy-atom distances from each Emi2 residue in 535–675 to five
APC/C subunits, using three AlphaFold3 screens (Fig. 2c; Methods). Engagement was defined
as pair_ipTM ≥ 0.50. In S1 (Emi2 + Cdc20 + Anapc10), 45 of 60 models were engaged for
Emi2–Cdc20; Anapc10 never reached its D-box co-receptor site in this three-chain construct
and is therefore a non-measurement (hatched row in Fig. 2c). In S2b (Emi2 + Anapc2 + Anapc11),
50 of 60 models were engaged for Emi2–Anapc11; Anapc2 and Anapc11 are measurements. In S3
(Emi2 + Anapc1 950–1250), 0 of 60 models engaged (pair_ipTM 0.07–0.21); Anapc1 is a
non-measurement, not evidence of absence of contact in the holocomplex.

In engaged S1 models the Emi2 D-box occupies the canonical Cdc20 degron site and reproduces
the Hsl1–Cdc20 contact geometry of PDB 5G04 `[ref: Zhang et al. 2016]` (backbone RMSD
0.5–0.6 Å; Supplementary Fig. S1). In engaged S2b models the ZBR packs against Anapc2 and
Anapc11 (median minima 2.41 Å and 2.94 Å). The RL-proximal C-terminus also approaches
Anapc2.

Allele-swap combinations of Emi2 with the partner set showed no meaningful structural
difference among engaged models (Supplementary Fig. S1).

### The two interface-proximal Emi2 substitutions are not dual-side hits

Blue lines in Fig. 2c mark parental positions 549 and 590.

**F549I (D-box).** Literature places the Emi2 D-box at a composite APC/C surface that can
involve Cdc20 and Anapc10 as co-receptor `(cite as above)`. In our engaged S1 models the
measurable contact is with Cdc20: any-atom median distance of residue 549 to Cdc20 is 3.7 Å,
but that approach is backbone. In 0 of 45 engaged S1 models does the side chain of residue
549 come within 4 Å of Cdc20. The flanking anchors R548 and L551 are the true contacts
(~2.2–2.6 Å). Allele statistics on the S1 screen show no detectable effect of Phe versus Ile
on Emi2–Cdc20 engagement (Supplementary Fig. S1). Anapc10 is identical between parents, and
Emi2–Anapc10 distances in S1 are non-informative because Anapc10 is not at its co-receptor
pose. We therefore cannot structurally confirm or reject an Anapc10-facing role for F549I
from the present construct; what we can say is that the Phe/Ile side chain does not contact
Cdc20, and that there is no parental Anapc10 difference to face it.

**V590M (pre-ZBR).** In engaged S2b models the closest approach of residue 590 is ~6.5–6.9 Å
to Anapc11 and ~13.5 Å to Anapc2; the side chain never reaches ≤4 Å (0/50). Neighbours on
the pre-ZBR stretch do contact Anapc11, but position 590 itself does not. Anapc11 is
identical between parents.

The intersection required by a structural dual-side mechanism for this triad — an Emi2
substitution inside a structurally tested contact that faces a variable APC/C position — is
therefore empty given the models in hand (Fig. 2c). Closing Anapc10 and Anapc1 geometry will
require larger constructs; until then those rows remain non-measurements.

### APC-side substitutions lie off the modelled Emi2 face

The same engaged models viewed from the partner side (Supplementary Fig. S2) place every
parental APC/C substitution far from the two in-construct Emi2 substitutions (F549I and
V590M). Median minimum heavy-atom distances to the nearest of those Emi2 sites are ~36 Å
for Cdc20 V291I, ~42 Å for Anapc2 V632A and ~59 Å for Anapc2 P654L — well outside the ≤4 Å
contact face. Anapc11 and Anapc10 are identical between parents; Anapc10 remains a
non-measurement in S1. These distances supply the structural support for the earlier claim
that none of the variable APC/C positions contact the variable Emi2 positions in the faces
we can currently measure.

---

## Discussion — stub (to expand after S4 / Anapc10 redo)

For the APC/C faces we can currently measure on *Darevskia* Emi2 — D-box↔Cdc20 and
ZBR↔Anapc2/Anapc11 — parental substitutions do not create a dual-side incompatibility.
F549I changes a side chain that never contacts Cdc20; V590M never reaches Anapc2/Anapc11
within 4 Å; Anapc10/Anapc11 are identical between parents; Cdc20 and Anapc2 differences lie
off the modelled faces. Anapc10 co-receptor geometry and Anapc1 remain unmeasured on Emi2.
Eight of ten Emi2 substitutions lie outside the complex window entirely and remain
candidates only for a regulatory, not docking, mechanism.

---

## Parked — Introduction (do not edit in this pass)

The Introduction text from co-authors is left unchanged for now. Open factual issues to
resolve in a later pass (not part of the present Results/Methods write-up):

1. Closing claim that F-box mutations bind the APC/C active centre — factually wrong
   (F-box → Skp1).
2. Inhibition-mechanism sentence should name receptors (D-box → Cdc20/Anapc10;
   ZBR → Anapc11/Anapc2/Anapc1; RL → Anapc2 CTD) without implying all are tested here.
3. Typos already listed previously (Moritz, nairensis, cleavage, etc.).

---

## Open tensions for co-author resolution

1. **F549I / Anapc10 wording.** Co-author draft: F549I “essential for interaction with
   Anapc10”. Literature: D-box is a composite Cdc20 (+ Anapc10 co-receptor) contact; the
   linked Development paper (Suzuki et al. 2010) maps D-box function in mouse Emi2 but does
   not by itself establish an Anapc10-only essential contact at this Phe. Our AF3 S1 screen
   measures Cdc20, not Anapc10. Current text keeps the literature framing and qualifies with
   our non-measurement — please choose final citation and soften or retain “essential”
   after checking the ScienceDirect target (`S2211546314000667`).
2. ***D. raddei* vs *D. r. nairensis*.** Both retained; standardise journal-wide if needed.
3. **Placeholders** still need real IDs / metrics / figure letters from the gene-inference
   side: `(assembly ID)`, `(SRA ID)`, `(GenBank ID)`, `(metrics?)`, `(Fig. SX/SY)`, `(cite)`.
