# Reproducibility Manual / 复现手册（v16 当前版 + v15 legacy）

Current manuscript: *Reproducible gut microbial–metabolite modules and cross-cohort functional signatures across the Alzheimer’s disease continuum* (v16).

> **重要版本警示**：本文件第 1–7 节记录的是 v15 legacy 流程，其中第 4–5 节的 167 条网络边、19.5% 中介和部分旧数字已被 v16 更正，不能再用于当前稿件推断。v16 当前结果、执行顺序和纠错依据见文末第 8 节。保留 legacy 是为了让第三方准确追溯“为什么同一旧脚本能重复、但结论仍可能错误”。
>
> **Version warning**: Sections 1–7 document the legacy v15 workflow. The 167-edge network, 19.5% indirect effect and other superseded values in Sections 4–5 must not be used for v16 inference. Use Section 8 for the current corrected workflow.

本手册目标：让第三方能够区分并验证计算重复性、公式正确性、原始数据复现和独立生物学验证。所有 v16 关键数字均由仓库内的 corrected/sensitivity 结果支持；原始 reads 和原始质谱谱图的边界单独声明。

---

## 1. Environment / 运行环境

| Item | Version used for v15 |
|---|---|
| Python | 3.13.18 |
| numpy | 2.5.5 |
| pandas | 3.1.1 |
| scipy | 1.18.1 |
| scikit-learn | 1.9.1 |
| matplotlib | 3.11.2 |
| python-docx | (manuscript builder only) |
| R (upstream pipeline only) | 4.6.1 with metafor, MaAsLin2, vegan, glmnet, DADA2 |

Install: `pip install -r requirements.txt`. Python 3.10+ is expected to work; exact versions above are the ones used to freeze the results. No GPU required; full re-run takes ~10–20 min on a desktop (the 117×1,873 Spearman matrix is the slowest step).

**Random seeds / 随机种子**: discovery recomputation uses `numpy.random.default_rng(20261119)` (both `01_` and `02_`); the legacy v14-audit permutation used seed 20261103 (not part of v15 results). With these seeds all outputs are deterministic.

## 2. Directory map / 目录结构

```
data/
  discovery/      5 original discovery-cohort input files (genus % table, lineage table,
                  metabolite abundance + descriptor, mediation/clinical metadata; n=84)
  meta_analysis/  per-cohort effect sizes, feature importance, subgroup tables (22 files)
  mr_analysis/    instrument SNP tables (30 genera), AD clumped SNPs, reverse-MR SNP matches,
                  re-extracted Kunkle effects, Bifidobacterium MiBioGen SNP effects
  mibiogen/       MBG.allHits.p1e4.txt (MiBioGen all-hits, instrument re-selection)
  source_data/    v14 Source Data tables used for cross-checks
results/
  discovery/      frozen outputs of scripts 01–02 (CSVs + rho/p matrices + summary)
  verification/   frozen outputs of scripts 03–07 (repro_*.csv, reverse_mr_harmonized.csv, JSON reports)
  manuscript/     Manuscript_Microbiome_v15.docx (built by scripts/08_build_manuscript)
scripts/          numbered in execution order
docs/             this manual
```

## 3. Run order / 执行顺序

```bash
python scripts/01_recompute_discovery.py   # filters -> modules -> module-trait stats
python scripts/02_recompute_stage2.py      # rho/p matrices -> networks -> permutation -> guilds -> taurine -> mediation
python scripts/03_reverse_harmonize.py     # allele-harmonized reverse MR (105 genera)
python scripts/04_verify_meta.py           # genus-level meta-analysis re-pooled from per-cohort effects
python scripts/05_verify_meta2.py          # corrected BH + REML cross-check
python scripts/06_verify_mr_ivw.py         # forward IVW MR for 30 genera from instrument tables
python scripts/07_verify_instruments.py    # re-select instruments from MiBioGen all-hits (P<1e-5, 500-kb clump)
# optional: python scripts/08_build_manuscript/build_v15.py  # rebuild the .docx
```

Every script is path-independent (all paths are computed relative to the repository root) and writes into `results/`. To verify a fresh clone, delete `results/discovery` and `results/verification` outputs first, re-run, and diff against the frozen copies released with the Zenodo archive.

## 4. Frozen analysis parameters ("hidden knobs") / 冻结分析口径

These are the parameters that, if changed, change downstream results. They are fixed in the scripts and must not be altered when reproducing v15.

| Parameter | Frozen value | Where |
|---|---|---|
| Discovery samples | n=84 (HC 10 / SCD 22 / MCI 23 / AD 29); disease stage HC=0, SCD=1, MCI=2, AD=3 | `mediation_analysis_data.csv` (Group column, whitespace stripped; control group coded `CG`) |
| Genus filter | detected in ≥20% of samples (≥17/84) AND mean relative abundance ≥0.01% → 117 of 347 genera | script 01 |
| Metabolite filter | MSI level `B(i)` AND `compound_source != nonnatural` → 1,873 metabolites (1,966 − 93) | script 01 |
| Co-abundance clustering | Spearman distance (1−ρ)/2... (correlation of rank-transformed genera), average linkage, `fcluster(maxclust, k=6)`; eigengene = PC1 of Z-scored member abundances, sign aligned to mean | scripts 01–02 |
| Module robustness | re-cut at k=4…10 | script 02 |
| SCFA enrichment | Fisher exact, 11 of 16 curated SCFA genera in M5 (OR=7.9, p=3.3e-4) | script 01 |
| Genus–metabolite correlation | Spearman across 84 samples; Benjamini–Hochberg FDR over 117×1,873 pairs | script 02 |
| High-confidence network | FDR<0.05 AND \|ρ\|≥0.6 → 167 edges; sensitivity at 0.50/0.55/0.65/0.70 | script 02 |
| Hub permutation | 2,000 permutations of genus profiles; statistic = mean degree of top-3 hubs (observed 29.7, null 0.0, p<0.0005) | script 02 |
| Mediation | rank-based OLS: a=Blautia→taurine, b=taurine→MMSE adjusted for Blautia; indirect=a×b; 5,000 bootstrap resamples, percentile 95% CI | script 02 |
| MR instruments (forward) | P<1e-5, LD clumping within 500-kb windows; 370 SNPs across 30 genera (3 genera with <3 instruments excluded) | script 07 |
| Forward MR | IVW primary; weighted median, MR-Egger, Cochran Q, leave-one-SNP-out sensitivity | script 06 |
| Reverse MR | 45 LD-independent AD SNPs (P<1e-5, 500-kb clump); MiBioGen betas flipped when MiBioGen effect allele equals Kunkle non-effect allele; IVW; BH FDR across 105 genera | script 03 |
| Meta-analysis pooling | random effects; REML τ² (metafor-equivalent) with DL cross-check; BH FDR | scripts 04–05 |
| Alpha diversity / PERMANOVA / PICRUSt2 / LOCO-CV | computed upstream (R + PICRUSt2 on usegalaxy.eu, NSTI≤2.0); per-cohort inputs archived in `data/meta_analysis` | Methods 2.5–2.7 |

## 5. Manuscript number → output file → script map / 数字溯源表

| Manuscript item | Key numbers | Output file | Script |
|---|---|---|---|
| §3.1 M5 eigengene vs stage / MMSE | r=−0.687 (p=5.2e-13); r=+0.673 (p=2.3e-12); ANOVA F=23.1 (p=7.1e-11); adjusted r=−0.455 (p=1.4e-5) | `module_trait_correlations.csv`, `summary_stage2.md` | 01/02 |
| §3.1 M2 risk module | 64 genera; r=+0.358 (p=8.3e-4); MMSE r=−0.409 | `module_trait_correlations.csv`, `modules_membership.csv` | 01 |
| §3.1 SCFA enrichment / robustness | OR=7.9 (p=3.3e-4); k=4–10 r=−0.69…−0.71 | `module_robustness.csv`, summary | 01/02 |
| §3.1 guilds | SCFA r=−0.538; acetogens r=−0.575 (adj. −0.410); H2 producers r=−0.142 (n.s.); H2S r=+0.223 (adj. n.s.) | `guild_analysis.csv` | 02 |
| §3.1 Blautia & key genera | HC median 0.28% vs AD 0.05% (p<1e-4); stage r=−0.571; Akkermansia stage r=+0.323 (p=0.003) | `key_genera_ADvsCG.csv` | 02 |
| §3.2 meta-analysis | 24 genera FDR<0.05 (10 enriched/14 depleted); pooled coef/p/I² per genus | `repro_genus_meta_fixed.csv`, `verify_meta_report*.json` | 04/05 |
| §3.2 subgroups | V3-V4 n=355 / V4-V5 n=400 coefficients | `data/meta_analysis/subgroup_by_region.csv` | (upstream) |
| §3.5 LOCO-CV | mean AUC=0.656; top features Lactobacillus +0.276, Methanobrevibacter +0.231, Dialister −0.213, Oscillibacter +0.209, Romboutsia −0.186 | `data/meta_analysis/feature_importance.csv` | (upstream) |
| §3.6 forward MR | Bifidobacterium IVW OR=1.306 [1.097–1.555] p=0.0027; UCG010 OR=0.885 p=0.082; Dialister p=0.615 | `repro_mr_forward.csv` | 06 |
| §3.6 reverse MR (harmonized) | 105 genera; 13 nominal; Bacteroides β=−0.039 p=1.35e-4 FDR=0.014; Bifidobacterium p=0.626 | `reverse_mr_harmonized.csv` | 03 |
| §3.7 network | 167 edges (64 pos/103 neg); hubs R-7(47)/UCG-005(26)/Family XIII(16); permutation p<0.0005; M5 mean degree 1.0 vs 9.1 | `network_edges_strong.csv`, `network_edges_thr055.csv`, `hub_stage_assoc.csv` | 02 |
| §3.7 module–metabolite | 704 metabolites FDR<0.05; 5-methyluridine ρ=0.639; taurine ρ=0.255 (FDR=0.050) | `module_disease_metabolite_assoc.csv` | 02 |
| §3.7 taurine | AD vs HC p=0.006; medians 6.46/5.28/5.56/4.87; KW p=0.022; MMSE r=+0.360 (p=7.6e-4); 36 genera FDR<0.05 | `taurine_genus_assoc.csv`, `summary_stage2.md` | 02 |
| §3.7 mediation | a=0.214, b=0.146, c=0.161, indirect 19.5%, bootstrap 95% CI [0.0003, 0.0679] | `summary_stage2.md` | 02 |
| Table 1 / Abstract | all of the above | same files | same |

Numbers in §3.3–3.4 (guild/pathway/PICRUSt2 meta-analyses) and alpha/beta diversity are pooled from the per-cohort inputs archived in `data/meta_analysis` and `data/source_data`; upstream R/PICRUSt2 pipelines are described in Methods 2.3–2.7 and their inputs are public (see `data/AVAILABILITY.md`).

## 6. Known limitations of the package / 边界说明

- The discovery input tables are **processed** genus-level percentages and metabolite intensities, not raw reads; re-processing raw reads with other pipelines may yield slightly different genus tables (stated in Methods 2.2).
- `verify_mr_full.py`-style re-extraction from the 543-MB Kunkle GWAS text is optional; the re-extracted SNP effects needed by script 06 are included (`kunkle_reextracted.csv`).
- The Senegalimassilia raw-zip allele spot-check from the v14 audit is omitted (needs multi-GB MiBioGen volumes); it does not affect `reverse_mr_harmonized.csv`.
- v14-audit scripts (why specific v14 numbers were wrong) are delivered separately with the v14 audit report, not in this repository.

## 7. Verification evidence / 闭环核验

Before release, all seven scripts were re-run from inside this repository and every output CSV/NPY was diffed against the frozen results used to write v15. **Result: all 22 compared files matched exactly (maxabsdiff = 0, text identical), including the 117×1,873 rho/p matrices.** Machine-readable report: `docs/verification_diff.json`.

Additional closed-loop checks within the run: forward IVW MR for 30 genera — max |ΔOR| = 0 and max |Δp| = 0 vs archive; instrument re-selection from MiBioGen all-hits — 373 SNPs re-extracted with 0 mismatches (370 across the 30 genera with ≥3 instruments; Faecalibacterium, Romboutsia, Ruminococcus1 excluded with 1 SNP each); harmonized reverse MR — Bacteroides p = 1.35e-4 (FDR = 0.014), Bifidobacterium p = 0.626.

---

## 8. v16 corrected workflow / v16 当前更正流程

### 8.1 为什么 legacy 零差异仍不够

v15 的 Spearman 实现将 ddof=0 标准化后的交叉积除以 n−1，使全部相关系数放大 n/(n−1)=84/83。该错误可被原脚本稳定重复，所以“22 个输出零差异”不能证明公式正确。v16 改用秩变量的标准 Pearson 相关，交叉积除以 n；旧矩阵按 84/83 还原后，与 corrected 矩阵最大绝对差为 1.11×10⁻¹⁶。

另一个独立问题是教育字段为文本，旧脚本数值强制转换得到 NaN 后填 0，因此名义上的教育调整没有正确实现。v16 将学历映射为 6/9/12/16 年，并用分类哑变量、HAMD/HAMA、非 HC 和逐诊断组排除模型复核。

### 8.2 当前执行顺序

从仓库根目录运行：

```bash
python scripts/09_recompute_corrected_associations.py
python scripts/10_discovery_covariate_sensitivity.py
Rscript scripts/11_compositional_sensitivity.R
Rscript scripts/12_summarize_compositional.R
```

R 组成型分析需要 ANCOMBC 2.14.0、MicrobiomeStat 1.4 和 metafor 5.0.1。脚本读取 `data/compositional_prepared/*_prepared.rds`，输出到 `results/sensitivity/compositional/`。首轮逐队列模型和异构输出汇总分开保存，避免汇总失败后重复耗时建模。

原始 84 例 FASTQ 重跑使用 `scripts/13_dada2_fresh.R`。原始 reads 和 SILVA 数据库不随处理后人源数据包重复发布，运行前设置：

```text
FASTQ_DIR   = 84 例 R1/R2 raw FASTQ 所在目录
SILVA_TRAIN = silva_nr99_v138.1_train_set.fa.gz
DADA2_OUT   = 独立输出目录
```

脚本固定 set.seed(42)、280/220 截断、pool=FALSE、250–500 bp、ASV 总丰度≥5、consensus 去嵌合和 minBoot=80，并保存 FASTQ/参考库 MD5、过滤与合并轨道、ASV 序列、属计数及 sessionInfo。原始重跑的最终状态必须以输出目录和回执为准，不能仅凭脚本存在宣称成功。

### 8.3 v16 当前数字溯源

| 项目 | v16 当前值 | 当前结果文件 |
|---|---|---|
| M5–病程，年龄+性别+教育年限 | partial r=−0.485451，p=2.863×10⁻⁶ | `results/sensitivity/covariates/partial_association_sensitivity.csv` |
| 非 HC 同模型 | r=−0.431968，p=0.0001215 | `partial_association_nonHC_sensitivity.csv` |
| 未调整 Blautia–牛磺酸–MMSE | a=0.213648，b=0.066893，c′=0.146334，indirect=0.014292，8.8975%，95% CI −0.001266–0.034646，p=0.0808 | `mediation_covariate_sensitivity.csv` |
| M5–牛磺酸 | ρ=0.252101，p=0.020699，FDR=0.054452 | `results/corrected_discovery/corrected_association_summary.json` |
| Blautia–牛磺酸 | ρ=0.213648，p=0.051013，FDR=0.119371 | `taurine_genus_assoc_corrected.csv` |
| 强网络 | 137 边；52 正/85 负；20 属/73 代谢物 | `network_edges_strong_corrected.csv` |
| ANCOM-BC2/LinDA | 30 个核心属在两种 REML+Hartung–Knapp 合并中均无 FDR<0.05 | `results/sensitivity/compositional/core_genera_reml_hk.csv` |
| 图件闭环 | 8 主图+9 补图；源文件和 SHA-256；旧 51 图单独盘点 | `docs/figure_audit/` |
| Word 投稿包闭环 | 8 主图和 9 补图嵌入；旧数字禁用扫描通过 | `docs/manuscript_revision/submission_closure.json` |

### 8.4 当前证据边界

- 保留：M5 与病程/MMSE 的稳健关联、聚合 SCFA 功能信号、跨队列中等迁移性。
- 降级：牛磺酸仅弱名义相关；中介区间包含 0；单属跨队列结果具有模型依赖性。
- 不主张：补充牛磺酸或产氢益生菌可减缓 AD；PICRUSt2 等于实测功能；计算复现等于独立生物学验证。
- 尚未完成：公共 755 例 FASTQ 的本轮全量重跑、发现队列原始代谢组谱图重处理、原始数据公开入库及真实 Zenodo DOI。

完整中文审计见 `docs/V16_END_TO_END_AUDIT_CN.md`。
