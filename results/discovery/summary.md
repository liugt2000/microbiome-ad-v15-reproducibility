# v15 发现队列重算结果汇总（n=84，117 属 x 1873 代谢物）

## 过滤口径
- 属：检出率>=20%（>=17/84 样本）且平均相对丰度>=0.01%，347->117 属
- 代谢物：MSI level B(i) 且非 nonnatural（药物/外源），14678->1873
- 样本：全部 84 例（CG 10 / SCD 22 / MCI 23 / AD 29）

## 共丰度模块（average linkage on 1-Spearman, maxclust k=6）

 module  n_genera  r_stage  p_stage  r_MMSE  p_MMSE  anova_F  anova_p  firmicutes_frac  mean_CG  mean_SCD  mean_MCI  mean_AD
      1         8   0.1296   0.2401 -0.0954  0.3878   1.5596   0.2057           0.5000  -0.0289   -0.0049   -0.0248   0.0333
      2        64   0.3578   0.0008 -0.4086  0.0001   5.6295   0.0015           0.6562  -0.0189   -0.0316   -0.0399   0.0622
      4         3  -0.0023   0.9838  0.0435  0.6947   1.1284   0.3427           0.6667  -0.0126   -0.0058    0.0348  -0.0189
      5        33  -0.6873   0.0000  0.6730  0.0000  23.1159   0.0000           0.8788   0.1152    0.0716   -0.0096  -0.0864
      6         8   0.1635   0.1372 -0.0937  0.3967   1.1742   0.3249           1.0000  -0.0301   -0.0079   -0.0166   0.0296

- 疾病关联模块 M5：33 属，r_stage=-0.687 (p=5.2e-13)，r_MMSE=0.673 (p=2.4e-12)，ANOVA p=7.1e-11
- 调整年龄/性别/教育后 r_stage=-0.455 (p=1.37e-05)

## 网络
- 主口径：FDR<0.05，共 22658 条边（正 11040 / 负 11618），涉及 117 属 / 1698 代谢物

- Top hub：Christensenellaceae_R-7_group(685), norank_f__Ruminococcaceae(607), UCG-005(598), norank_f__Oscillospiraceae(581), Family_XIII_AD3011_group(566), norank_f__UCG-010(508), Ruminococcus_gnavus_group(500), Blautia(494)

- 置换检验：top3 平均度 p=0.000；SCFA 边占比 p=0.5385；Fisher(hub SCFA 富集) p=0.5964

## guild

         guild  n_members                                                                                                                                                                                                                                                        members  r_stage  p_stage  r_MMSE  p_MMSE   KW_p  ADvsCG_MWU_p  r_stage_adjusted  p_stage_adjusted
SCFA_producers         16 Faecalibacterium;Roseburia;Butyricicoccus;Coprococcus;Eubacterium_hallii_group;Anaerostipes;Agathobacter;Subdoligranulum;Lachnospira;Lachnospiraceae_NK4A136_group;Anaerotruncus;Odoribacter;Dorea;Ruminococcus_torques_group;Lachnoclostridium;Flavonifractor  -0.5375   0.0000  0.5424  0.0000 0.0000        0.0001           -0.3272            0.0024
  H2_producers         13                                       Ruminococcus;Ruminococcus_gnavus_group;Ruminococcus_torques_group;Blautia;Dorea;Clostridium_sensu_stricto_1;Eubacterium_hallii_group;Bifidobacterium;Bacteroides;Prevotella;Enterobacter;Klebsiella;Escherichia-Shigella  -0.1424   0.1964  0.0694  0.5306 0.0611        0.0796           -0.1288            0.2428
 H2S_producers          3                                                                                                                                                                                                                          Desulfovibrio;Bilophila;Fusobacterium   0.2226   0.0418 -0.2316  0.0340 0.0095        0.1400            0.0219            0.8432
     acetogens          2                                                                                                                                                                                                                               Blautia;Eubacterium_hallii_group  -0.5752   0.0000  0.5446  0.0000 0.0000        0.0000           -0.4103            0.0001

## 牛磺酸
- 分组中位数 CG=6.463 SCD=5.280 MCI=5.564 AD=4.870；KW p=0.022；AD vs CG p=0.006；r_stage=-0.324 (p=0.003)；r_MMSE=0.360 (p=0.001)

- FDR<0.05 的牛磺酸相关属 36 个；top：unclassified_o__Oscillospirales(-0.50); Alistipes(-0.49); norank_f__Oscillospiraceae(-0.48); Christensenella(-0.46); NK4A214_group(-0.45)
