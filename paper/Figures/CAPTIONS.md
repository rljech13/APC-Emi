# Figure captions

Journal-ready captions. Numbers trace to `Interface_variability/REPORT.md`,
`AF3_results_S1/REPORT_AF3_S1.md`, `AF3_results_S2b/REPORT_AF3_S2b_S3.md`,
`emi2_combined_distance_matrix.json`, and the B-factor column of
`Fasta/Models/Emi2_Dval_interface400-675_AF3_model0.pdb`.

**Main-text figure:** Figure 2 only.  
**Supplementary:** former Figure 3 (D-box–Cdc20 screen).  
**Withdrawn from main:** standalone Figure 2c (ZBR superposition) and standalone
Figure 4 (heatmap now panel 2c).

---

## Figure 2 — `Figure_2_model_and_map.{svg,png}`  (script: `make_figure2_story.py`)

**Emi2 structure, where the monomer is trustworthy, and where parental substitutions
sit relative to modeled APC/C contacts.**

**(a)** AlphaFold3 monomer of Emi2 residues **400–675** (*D. valentini*), cartoon
coloured by functional region: pale pink, folded F-box domain 452–537; solid red,
Skp1-binding F-box motif 479–519 nested inside it; cyan, D-box 548–551; green, ZBR
600–648; salmon, RL tail 658–675; grey, unannotated. This panel answers: *what does
the APC/C-facing half of Emi2 look like?*

**(b)** Mean per-atom pLDDT of each region of the same monomer (whiskers: 10th–90th
percentile). Bar colours follow the AlphaFold/AFDB bands: <50 orange, 50–70 yellow,
70–90 light blue, >90 dark blue (keyed in this caption, not redrawn under the panel).
Columns at right: mean, structural character, and whether the region contacts APC/C in
the **complex** screens of this study (S1/S2b), not in the monomer. The model is
confident only where Emi2 has an autonomous fold — F-box domain **86.7**, Skp1 motif
**93.4**, ZBR **72.0**. D-box (**38.7**), linker (**39.5**) and spacers sit in the
very-low band: expected for SLiMs that order on binding. This panel answers: *where
may we trust the monomer, and why APC contacts require multimers?*

**(c)** Linear map of full-length Emi2 (675 aa) with the ten substitutions to
*D. r. nairensis* (black dots). A dashed navy box marks the complex window
**535–675**; callout lines link that box to the expanded strip below, then the
interface heatmap: median min. heavy-atom distance from each Emi2 residue in that window
to Cdc20 (S1 engaged), Apc10 (S1 non-measurement), Apc2 and Apc11 (S2b engaged), and
Apc1 (S3 non-measurement). Blue vertical lines mark substitutions **549** and **590**.
Heatmap cells are any-atom minima (backbone or side chain). For 549 the median to Cdc20
is **3.7 Å**, but the side chain never reaches ≤4 Å in engaged S1 models; neighbours
548/551 are closer (~2.2–2.6 Å), which is why the D-box band looks dark red. For 590
medians are ~12 / 13 / 6.9 Å to Cdc20 / Apc2 / Apc11 (never ≤4 Å). This panel answers:
*where are the parental differences, and do they sit on modeled APC/C contacts?*

*Not shown in main text.* Former Figure 2c (ZBR allele superposition, Cα RMSD 0.19 Å
over 600–648) is omitted as non-essential. Former Figure 4 is absorbed into (c).

---

## Supplementary Figure S1 (former Figure 3) — `Figure_S1.png`

Panel **e** (seed plot, allele table, noise-floor) is removed. The file to insert is `Figure_S1.png` (panels a–d only). The full sheet remains `Figure_3_v1.png`.


**The Emi2 D-box docks canonically on Cdc20, and the F549I substitution has no
detectable effect on that interaction.** All structural panels come from an AlphaFold3
screen of Emi2(535–675) + Cdc20 + Anapc10: 60 models = 4 allele combinations × 3 seeds
× 5 diffusion samples. Within this construct the two Emi2 alleles differ at exactly two
positions, F549I and V590M; the other eight substitutions lie outside it. Anapc10 is
identical between the two species and Cdc20 differs at one position (V291I), some 30 Å
from the D-box groove, so the partner-set contrast is a genuine null.

**(a)** Overview of the top-ranked Emi2–Cdc20 model: the Emi2 D-box 548–551 (red) lies
in the canonical D-box receptor site of the Cdc20 WD40 propeller (dark blue on a pale
propeller), with the ZBR (teal) and the RL tail (orange) projecting away from it.
Apc10 is omitted from the panel because it is never at its D-box co-receptor site. It
is not contact-free — 59 of the 60 models place 1–39 Emi2 residues (median 14) within
4 Å of Apc10 — but its centre of mass sits 12.7–86.0 Å (median 64.9 Å) from its
reference position in the Cdc20 frame, and the surface it presents to Emi2 shares zero
residues with the surface that grips the Hsl1 degron in the reference structure. Apc10
is held in place by Apc1 and Apc3/Cdc27, both absent from a three-chain input, so the
C-terminal half of the degron is unsupported here and every Emi2–Apc10 number in the
data set is noise rather than a measurement.

**(b)** Validation against the reference: the modelled pose reproduces the Cdc20–Hsl1
interaction of PDB 5G04, with the same contact set and a D-box backbone RMSD of
0.5–0.6 Å. The reference is cryo-EM of the intact human APC/C–Cdc20–Hsl1 at 3.9 Å
(Zhang et al. 2016, *Nature* 533:260); chain S is the Hsl1 D-box peptide R828-A-A-L831.
An RMSD of 0.5–0.6 Å is well inside the coordinate error of a 3.9 Å reference, so the
defensible claim is the reproduced binding mode and contact set, not the RMSD
magnitude. Our Cdc20 was superposed on 5G04 Cdc20 (Cα RMSD 0.53–0.60 Å over about 300
core residues) with no local fit. All nine canonical D-box receptor residues of Cdc20
are identical in *Darevskia*, and Emi2 contacts 8 of the 9 in all four top-ranked
models.

**(c, d)** The D-box of the two alleles at atomic detail: RFAL with Phe549 in
*D. valentini* (c) and RIAL with Ile549 in *D. r. nairensis* (d), with R548 and L551 as
the flanking anchors. In neither allele does the side chain of residue 549 enter the
interface: no side-chain atom reaches Cdc20 within 4 Å in any of the 45 D-box-engaged
models. Phe549 never comes closer than 5.68 Å (range 5.68–6.11 Å, median 5.92 Å;
0 of 20 models within 5 Å) and Ile549 reaches 4.21–5.46 Å (median 4.65 Å; 20 of 25
models within 5 Å, 0 of 25 within 4 Å) — hence the panel subtitles, which quote the
strongest cut-off the data support for each allele. The burial of position 549 is
entirely backbone: in 45 of 45 models the closest Emi2 atom is the backbone carbonyl
oxygen and its partner is always Asp133, the only Cdc20 residue within 4 Å of residue
549 anywhere in the set. The corresponding ΔSASA on binding is 14–20 Å² (11–18 % of the
residue buried), which is backbone rather than side-chain surface; residue 590 buries
0 Å². Because the backbone of position 549 is identical in the two alleles by
definition, the part that actually differs — Phe versus Ile — lies outside the interface
in every engaged model. This is a structural reason for the null result that does not
depend on the statistics of (e). For comparison, R548 grips through its guanidinium
(side chain 2.05–2.58 Å from Cdc20) and L551 is fully engaged with 6–7 Cdc20 partners
in every model (99 % buried, ΔSASA 177–178 Å²).

**(e)** The statistics, in three parts. *Left*, every AlphaFold3 seed: 12 seeds, of
which 3 never found the docked basin (open symbols). The seed-to-seed spread of the
Emi2–Cdc20 ipTM is 0.46 (between-seed SD up to 0.26), whereas the 5 diffusion samples
of one seed are near-duplicates (within-seed SD 0.005–0.026 ipTM). *Middle*, the allele
comparison over the 45 D-box-engaged models grouped into the 9 seeds that produced them
(4 *D. valentini*, 5 *D. r. nairensis*): no metric shows an allele effect. The last row
is the largest raw allele effect anywhere in the data set, Δ = −0.069 on the
Emi2–Cdc20 ipTM of all 60 models (p = 0.06), which is entirely an artefact of 2 of 6
*D. valentini* seeds versus 1 of 6 *D. r. nairensis* seeds failing to dock (Fisher
exact p = 1.00). *Right*, each effect divided by that metric's own empirical null, the
contrast between identical partner sets, which must be zero biologically; the shaded
band is ±1 null. Every effect falls inside its own noise floor, and the −0.069 raw
effect is about four times smaller than the 0.290 partner-set null on all 60 models.
The footnote retained in the panel states the test convention: mean ± SD over seed
means, Mann–Whitney *p* on seeds; per-model *n* = 20 vs 25 gives p = 0.24 / 0.69 /
0.81 / 0.03 / 0.94, but a seed's 5 diffusion samples are near-duplicates (SD
0.005–0.026 ipTM), so the seed is the unit of replication.

*Conclusion.* No detectable effect of F549I, or of V590M, on Emi2–APC/C binding in this
model set. Only the 7-residue D-box segment is a reproducible structural observation:
the linker, the ZBR and the RL tail are placed differently in essentially every model
(Emi2 centre-of-mass RMS spread 18–25 Å within one combination, maximum 43 Å).

---

## Figure 4 — absorbed into Figure 2c

Standalone `Figure_4_interface_map.{svg,png}` is kept on disk for provenance but is
**no longer a main-text figure**. Its content is panel **2c**. Older Emi1-transfer /
S1-only heatmaps remain under `Figures/superseded/`.

---

## Supplementary Figure S2 — `Figure_S_subs_vs_interface.{png,svg}`
(script: `make_supp_subs_vs_interface.py`;
data: `Interface_variability/partner_distance_to_emi2_muts.json`)

Two proximity profiles. Tall means close to the nearer of F549I and V590M.
The dashed line is the 4 Å contact cut-off. Parental substitutions are marked
in red and sit well below that line.

**(a) Cdc20**, engaged S1 models. V291I is 36 Å away.
**(b) Anapc2**, engaged S2b models. Residues 1–219 were not in the construct.
V632A and P654L are 42 Å and 59 Å away.

Apc10 and Apc1 are not drawn: both were excluded from the distance analysis.
Apc11 is identical between the parents.