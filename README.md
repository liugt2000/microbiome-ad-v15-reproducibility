# Microbiome–AD v16 Reproducibility Package

Reproducibility package for:

> *Reproducible gut microbial–metabolite modules and cross-cohort functional signatures across the Alzheimer’s disease continuum* (v16, 2026).

**Current inference uses `results/corrected_discovery/` and `results/sensitivity/`.** The v15 scripts and outputs remain as a documented legacy audit trail; they contain a Spearman scaling error and an incorrectly swapped mediation-path summary. Exact regeneration of legacy files is not evidence that their formulas were correct.

No Zenodo DOI has been issued for this version. A DOI badge will be added only after a real GitHub Release and Zenodo Publish step.

## What the corrected study shows

- A 33-genus, SCFA-enriched module is associated with disease stage (r=−0.687) and MMSE (r=0.673). The stage association remains after correctly encoded age, sex and education adjustment (partial r=−0.485).
- M5–taurine is weak and does not pass FDR<0.05 (ρ=0.252, FDR=0.054). The Blautia–taurine–MMSE indirect effect is not statistically significant (0.0143; bootstrap 95% CI −0.0013 to 0.0346) and is near zero after covariate adjustment.
- Across five public cohorts (n=755), MaAsLin2 identifies 24 genera at FDR<0.05, but neither ANCOM-BC2 nor LinDA identifies any of 30 prespecified core genera at FDR<0.05. Single-genus signals are therefore method dependent.
- Aggregate SCFA capacity is lower in AD (pooled d=−0.524), with cohort heterogeneity. Fermentative H₂-producing guilds show a weak pooled decrease (d=−0.186, I²=0%), while the broad discovery-cohort association is null.
- LOCO-CV mean AUC=0.656 indicates modest transportability, not clinical diagnostic utility. Bidirectional MR remains secondary and instrument limited.

The data do **not** show that taurine supplementation or H₂-producing probiotics will slow AD.

## Current corrected workflow

```bash
pip install -r requirements.txt
python scripts/09_recompute_corrected_associations.py
python scripts/10_discovery_covariate_sensitivity.py
# Run the following two R scripts from the repository root:
Rscript scripts/11_compositional_sensitivity.R
Rscript scripts/12_summarize_compositional.R
```

Required R packages include ANCOMBC 2.14.0, MicrobiomeStat 1.4 and metafor 5.0.1. Inputs for all five cohorts are included in `data/compositional_prepared/`.

The independent 84-sample raw-read DADA2 script is `scripts/13_dada2_fresh.R`. It requires user-supplied `FASTQ_DIR` and `SILVA_TRAIN` environment variables; restricted raw reads are not committed to GitHub.

## Legacy workflow

Scripts 01–07 reproduce the earlier frozen v15 computational outputs. They are retained to explain prior irreproducibility and implementation errors, not to support current v16 claims. See the prominent warning and v16 run order in `docs/REPRODUCIBILITY_MANUAL.md`.

## Key directories

- `results/corrected_discovery/` — corrected Spearman matrices, network and taurine associations.
- `results/sensitivity/covariates/` — education encoding, adjusted associations and corrected mediation.
- `results/sensitivity/compositional/` — ANCOM-BC2/LinDA cohort outputs and REML+Hartung–Knapp summaries.
- `docs/figure_audit/` — 17-figure source hashes, closure checks and legacy inventory.
- `docs/V16_END_TO_END_AUDIT_CN.md` — Chinese end-to-end audit and evidence boundaries.
- `results/manuscript/Manuscript_Microbiome_v16_Microbiome.docx` — current manuscript.
- `results/manuscript/Supplementary_Material_Microbiome_v16.docx` — nine-figure supplement.

## Data boundaries

Processed de-identified discovery tables, public-cohort prepared inputs and derived MR tables are included. Public accessions and large GWAS sources are documented in `data/AVAILABILITY.md`. Raw discovery sequencing reads and raw mass-spectrometry spectra are not publicly deposited here; their archiving remains subject to ethics and consent requirements. Public 755-sample FASTQ files were not fully re-downloaded and rerun in this correction cycle.

## Citation and license

Cite the manuscript and this repository. A concept DOI will be added only after Zenodo issues it. Code is released under MIT; see `LICENSE`.

## Contact

Corresponding-author details will be completed before journal submission.
