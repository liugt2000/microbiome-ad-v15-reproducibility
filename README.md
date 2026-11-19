# Microbiome–AD v15 Reproducibility Package

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Reproducibility package for the manuscript:

> *Integrative Multi-Omics and Multi-Cohort Meta-Analysis Reveal Depleted Butyrate- and Hydrogen-Metabolizing Microbial Consortia and Reduced Fecal Taurine Across the Alzheimer's Disease Continuum* (v15, 2026).

**Goal:** any third party can regenerate every key number in the manuscript, without contacting the authors, by running the scripts in this repository. All discovery-cohort inputs, per-cohort validation effect sizes, and MR instrument tables are included; very large public GWAS files are documented in `data/AVAILABILITY.md`.

## What the study shows

- A protective, SCFA-enriched Firmicutes co-abundance module (33 genera) declines monotonically across the AD continuum (eigengene vs stage Spearman r = −0.687; MMSE r = +0.673).
- The acetogenic H₂-consuming guild (Blautia + Eubacterium hallii group) is depleted (r = −0.575); fermentative H₂ producers are depleted in the validation meta-analysis (pooled d = −0.186, I² = 0%).
- Fecal taurine is reduced in AD (p = 0.006), correlates with MMSE (r = +0.360), and mediates 19.5% of the Blautia–MMSE association (bootstrap 95% CI excludes zero; exploratory).
- Validation across five public cohorts (n = 755): 24 FDR-significant genera; total SCFA-producing capacity depleted (d = −0.524); LOCO-CV mean AUC = 0.656.
- Bidirectional MR: Bifidobacterium IVW OR = 1.306 (p = 0.0027); allele-harmonized reverse MR shows AD liability → lower Bacteroides (FDR = 0.014) and no reverse causation for Bifidobacterium (p = 0.626).

## Quick start

```bash
pip install -r requirements.txt
python scripts/01_recompute_discovery.py
python scripts/02_recompute_stage2.py
python scripts/03_reverse_harmonize.py
python scripts/04_verify_meta.py
python scripts/05_verify_meta2.py
python scripts/06_verify_mr_ivw.py
python scripts/07_verify_instruments.py
```

Outputs are written to `results/discovery` and `results/verification` and can be diffed against the frozen copies in this release. All scripts are path-independent and use fixed random seeds (`20261119`).

## Documentation

- `docs/REPRODUCIBILITY_MANUAL.md` — full manual: environment, frozen analysis parameters, run order, and a **manuscript-number → output-file → script** traceability table.
- `data/AVAILABILITY.md` — public data accessions (5 SRA cohorts, MiBioGen, IGAP) and notes on processed vs raw data.
- `results/manuscript/Manuscript_Microbiome_v15.docx` — the manuscript this package supports.

## Citation

If you use this package, please cite the manuscript and this archive (concept DOI above; release tag `v1.0`). See `CITATION.cff`.

## License

MIT (see `LICENSE`). Processed discovery-cohort tables are de-identified genus/metabolite-level measurements; no raw reads or directly identifying clinical information are included.

## Contact

Corresponding author: to be completed upon publication (see manuscript title page).
