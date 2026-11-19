# 07_v15_Recomputed_Results 说明

本文件夹包含 Microbiome v15 稿件的**全部发现队列重算结果与复现脚本**。重算目的：v14 审计发现发现队列侧 4 项实质问题（反向 MR 未做等位基因协调、网络 hub 标签错位、permutation p 值不实、牛磺酸 §4.6 无来源），v15 全部数字以本文件夹重算结果为准。

## 数据过滤口径
- 属：347 → 117（检出 ≥17/84 样本 且 平均相对丰度 ≥0.01%），84 样本（HC 10 / SCD 22 / MCI 23 / AD 29）
- 代谢物：14678 → 1873（B(i) 级 1966 个高置信内源代谢物 − 93 个非天然/外源）
- 原始数据：`..\..\..\粪便微生物多样性分析原始数据和深度分析数据\`（otu_taxon_Genus.percent.xls；深度分析文件\workflow_results\1.Preprocess\mix\metab_abund.txt + metab_desc.txt）

## 脚本（按顺序运行）
1. `recompute_discovery.py` — 过滤、模块聚类（k=6 maxclust）、模块-性状相关、模块成员
2. `recompute_stage2.py` — rho 行标准化修正后的属-代谢物相关、强阈值网络（|ρ|≥0.6 ∧ FDR<0.05 = 167 边；0.55 = 509 边存档）、置换检验（2000 次）、hub 分析、guild 分析、牛磺酸相关、Blautia→taurine→MMSE 中介（OLS + bootstrap 5000）、模块稳健性（k=4~10）、SCFA 富集 Fisher
3. `v15_content.py` + `build_v15.py` — v15 稿件内容与 docx 生成器（输出 `..\01_Manuscript\Manuscript_Microbiome_v15.docx`）

## 结果文件
| 文件 | 内容 |
|---|---|
| genus_117_filtered.csv / metab_1873_filtered.csv | 过滤后属/代谢物表 |
| modules_membership.csv / module_trait_correlations.csv | 模块成员与模块-性状相关（M5 r=−0.687 p=5.2e-13；MMSE +0.673 p=2.3e-12；M2 +0.358） |
| module_robustness.csv | k=4~10 模块稳健性（M5 核心 r=−0.69~−0.71） |
| module_disease_metabolite_assoc.csv | M5-代谢物相关（704 个 FDR<0.05；taurine ρ=0.255 FDR=0.050） |
| rho_matrix.npy / p_matrix.npy | 117×1873 Spearman 矩阵（修正版） |
| network_edges_fdr.csv / _nominal.csv / _strong.csv / _thr055.csv | 网络边表（22,658 / 167 / 509 边） |
| hub_stage_assoc.csv | hub 属-疾病分期相关 |
| guild_analysis.csv | SCFA/acetogens/H2/H2S guild 全套统计 |
| key_genera_ADvsCG.csv | 关键属组间比较 |
| taurine_genus_assoc.csv | 牛磺酸-属相关（36 属 FDR<0.05） |
| summary_stage2.md | 第二阶段结果汇总 |

## v15 相对 v14 的关键勘误
1. 亚组样本量 569/186 → 355/400（V3-V4 四队列 / V4-V5 北京队列）
2. Roseburia 荟萃 −0.306/p=0.073 → −0.303/p=0.112/I²=74.5%；Faecalibacterium I² 72.3% → 93.2%
3. PICRUSt2 n=661 → 631（Turkey 队列控件为正常对照而非 HC 标签）
4. MR 工具变量 P<1e-4 → P<1e-5、LD 500kb 窗口 clump；反向 MR 加入等位基因协调（105 属、13 nominal、Bacteroides FDR=0.014、Bifidobacterium p=0.626）
5. §3.5 top features 修正为 Lactobacillus +0.276 / Methanobrevibacter +0.231 / Dialister −0.213 / Oscillibacter +0.209 / Romboutsia −0.186
6. 网络 hub 改写为 Christensenellaceae_R-7_group(47)/UCG-005(26)/Family_XIII_AD3011_group(16) 等，全部位于风险模块 M2；置换 p<0.0005
7. 牛磺酸段落全部以重算结果改写（AD vs HC p=0.006、MMSE r=+0.360 p=7.6e-4、中介 19.5% bootstrap CI [0.0003, 0.0679]）
8. taurine 引用编号 [178-187] → [47-53]；H₂ 抗氧化 → Ohsawa 2007 [54]；GF 淀粉样 → Harach 2017 [55]
