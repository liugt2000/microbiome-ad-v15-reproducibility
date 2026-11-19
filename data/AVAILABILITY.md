# Data availability / 数据来源与获取

## 1. Public validation cohorts (16S rRNA, NCBI SRA)

| Cohort | Accession | n (AD/HC) | Region |
|---|---|---|---|
| Beijing | PRJNA946400 | 400 | V4-V5 |
| Hangzhou | PRJNA496408 | 61 | V3-V4 |
| Lishui | PRJNA633959 | 171 | V3-V4 |
| Shanghai | PRJNA489760 | 60 | V3-V4 |
| Turkey | PRJNA734525 | 63 | V3-V4 |

All cohorts were reprocessed with a unified DADA2 pipeline (SILVA v138.1; truncLen 270/220 for V3-V4 and 230/220 for V4-V5) and aggregated to genus level. Per-cohort effect sizes and downstream summaries used in v15 are archived in `data/meta_analysis` and `data/source_data`.

## 2. GWAS summary statistics

- **MiBioGen** (Kurilshikov et al. 2021, n=18,340): https://mibiogen.gcc.rug.nl — genus-level summary statistics. `data/mibiogen/MBG.allHits.p1e4.txt` (all hits at P<1e-4) is included for instrument re-selection. The multi-GB per-taxon zipped summary files (used only for an optional allele spot-check in the v14 audit) can be downloaded from the same site.
- **IGAP AD** (Kunkle et al. 2019, n=63,926): `Kunkle_etal_Stage1_results.txt` (543 MB) available from https://www.niagads.org/datasets/ng00036 (or the IGAP/EBI GWAS catalog mirror). The SNP effects required by `scripts/06_verify_mr_ivw.py` are already re-extracted and included as `data/mr_analysis/kunkle_reextracted.csv`; downloading the full file is therefore NOT necessary to reproduce v15.

## 3. Discovery cohort (this study, n=84)

Included in `data/discovery/` as **processed, de-identified** tables:
- `otu_taxon_Genus.percent.xls` — genus-level relative abundances (347 genera × 84 samples)
- `otu_taxon_Genus.full.xls` — taxonomic lineage table
- `metab_abund.txt` / `metab_desc.txt` — fecal metabolite intensities (14,678 features) and MSI-level descriptors
- `mediation_analysis_data.csv` — per-sample clinical metadata (diagnostic group, MMSE, age, sex, education years), sample-ID keyed, no direct identifiers

Ethics: discovery cohort approved by the local institutional review board (JAMHC-IRB-202201); all participants or legal guardians provided written informed consent. Raw sequencing reads and raw metabolomics spectra are NOT included in this repository; they will be deposited in public repositories upon manuscript acceptance, and are available from the corresponding author under the ethics approval framework.

## 4. What is NOT in this repository

- Raw sequencing reads / raw mass-spectrometry files (see above).
- Multi-GB MiBioGen per-taxon zip volumes (optional audit spot-check only).
- The 543-MB Kunkle GWAS text (pre-extracted SNP effects included instead).
- v14 audit scripts (delivered separately with the v14 audit report).
