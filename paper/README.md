# APC/C–Emi2 (Darevskia)

Code, sequences and AlphaFold3 models behind the manuscript on Emi2–APC/C compatibility in rock lizards.

The Google Doc export is in `article/manuscript_export.txt`.

## Layout

| Path | What it is |
|---|---|
| `Fasta/` | Inferred amino-acid sequences, AF3 job inputs, monomer and complex models |
| `AF3_results_S1/` | Screen S1: Emi2 535–675 + Cdc20 + Anapc10 (12 jobs × 5 models) |
| `AF3_results_S2b/` | Screen S2b: Emi2 + Anapc2 catalytic fragment + Anapc11 |
| `AF3_results_S3/` | Screen S3: Emi2 + Anapc1 950–1250 |
| `Interface_variability/` | Dual-side substitution test and distance matrices |
| `Figures/` | Figure scripts and rendered panels |
| `SNP_regulatory_analysis/` | Emi2 substitution context (not the main interface result) |
| `Manuscript/` | Earlier draft notes |

Scripts resolve the project root from `__file__`. They were checked with the local environment at `~/Desktop/emi2_viz_venv` (`biopython`, `pandas`, `numpy`, `scipy`, `gemmi`, `matplotlib`, `Pillow`).

## What is not in this tree

Left on the Desktop copy because they are large and not required to re-run the distance and confidence scripts:

- AlphaFold3 MSAs and template hits
- `*_full_data_*.json` (per-atom confidence dumps, about 1.1 GB)
- Raw unzip folders `AF3_results_incoming` and `AF3_results_incoming2` (duplicates of S2b/S3)
- Genome scaffold `DM1.ragtag.scaffold.fasta` (1.3 GB)
- Superseded figure renders

There is no script for the sequence inference in the methods (tBLASTn, STAR, bcftools, BWA). Those steps were done outside this tree; the resulting FASTA files are in `Fasta/`.

The older CombFold / ColabFold exploration (human APC/C, PDB 2RT9 and 4UI9) stays in the parent folder of this repository. It is not the AlphaFold3 analysis in the manuscript.

## Check

Re-run on the copied models (25 Sep 2026):

- S1 confidence table: Emi2–Cdc20 `pair_ipTM` median 0.58. The two-sided test counts 45 engaged models (`pair_ipTM` ≥ 0.50) and places Cdc20 V291I 31.8 Å from the D-box. The manuscript text says about 36 Å for that site; use the script number.
- S2b: 50 engaged models. Residue 590 is 6.9 Å from Anapc11 and 13.5 Å from Anapc2, as in the manuscript.
- Figure scripts rewrote `Figure_2_model_and_map`, `Figure_4_interface_map` and `Figure_S_subs_vs_interface`.
- SNP scripts run. `context_analysis.py` prints that F549I both creates and destroys an RxxL D-box (`RFAL` → `RIAL`). Both tetrapeptides match RxxL, so that line is a pattern-match artefact, not a lost motif.

`py_compile` passed for all 24 scripts. ChimeraX session scripts (`*.cxc`) were not executed.
