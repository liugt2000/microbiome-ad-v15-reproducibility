# v15 重算·第二阶段结果

## 模块稳健性（不同k下 Faecalibacterium 所在模块）

 k           anchor  module_size  n_scfa_in  r_stage   p
 4 Faecalibacterium           41         11  -0.7105 0.0
 5 Faecalibacterium           41         11  -0.7105 0.0
 6 Faecalibacterium           33         11  -0.6873 0.0
 7 Faecalibacterium           33         11  -0.6873 0.0
 8 Faecalibacterium           33         11  -0.6873 0.0
10 Faecalibacterium           33         11  -0.6873 0.0

## SCFA 在 M5 的富集
Fisher p=0.00033 (OR=7.9), 11/16 SCFA 属位于 M5

## 模块M5-代谢物 top（修正rho）

正相关 top10:
   metab_id                                                                     name    rho   p  FDR
metab_16558                                                          5-methyluridine 0.6387 0.0  0.0
metab_22818                                               2,3-dinor Prostaglandin E1 0.6376 0.0  0.0
 metab_1020                                                           Nicotinic acid 0.6283 0.0  0.0
metab_26932                                                                  Uridine 0.6217 0.0  0.0
metab_11702 Ethyl 4-[3-(4-benzylpiperazin-1-yl)propylamino]quinazoline-2-carboxylate 0.6207 0.0  0.0
metab_20363                                  '5,8,11-trihydroxyoctadec-9-enoic acid' 0.6184 0.0  0.0
metab_13098                                                         5-methylcytosine 0.6160 0.0  0.0
metab_25821                                                  '1,7-dimethyluric acid' 0.6020 0.0  0.0
metab_27191                                                           Picolinic acid 0.5961 0.0  0.0
 metab_1847                 3-isobutyl-6-(1-hydroxy-2-methylpropyl)-2(1h)-pyrazinone 0.5856 0.0  0.0
负相关 top10:
   metab_id                                                                           name     rho      p    FDR
 metab_8571                                           N-cyclohexanecarbonyltetradecylamine -0.5903 0.0000 0.0000
metab_11884                                                                    Lys Leu Gln -0.4924 0.0000 0.0001
metab_21157                                                               Sulfocholic Acid -0.4561 0.0000 0.0003
  metab_782                                                     Methylimidazoleacetic acid -0.4503 0.0000 0.0004
metab_11892 1-[2-(4-fluorophenyl)pyrimidin-4-yl]-N-(3-methylbutyl)piperidine-3-carboxamide -0.4491 0.0000 0.0004
 metab_1040                                                                    Eritadenine -0.4397 0.0000 0.0005
metab_11376                                                                        Peg n12 -0.4247 0.0001 0.0007
 metab_5599                                                                     Pyrogallol -0.4228 0.0001 0.0008
metab_14490                                                                        Peg n14 -0.4217 0.0001 0.0008
 metab_5584                                                          '3,5-dimethoxyphenol' -0.4203 0.0001 0.0008

## 强阈值网络 (|rho|>=0.6 & FDR<0.05)
167 边, top hub: Christensenellaceae_R-7_group(47), UCG-005(26), Family_XIII_AD3011_group(16), norank_f__Oscillospiraceae(15), norank_f__Ruminococcaceae(13), unclassified_o__Oscillospirales(13), NK4A214_group(7), UCG-002(5)
置换 p=0.0000

## 牛磺酸
- 正相关 top: Ruminococcus_gnavus_group(0.37); Dorea(0.23); Granulicatella(0.23); Anaerostipes(0.23); Erysipelotrichaceae_UCG-003(0.22)
- 负相关 top: unclassified_o__Oscillospirales(-0.50); Alistipes(-0.49); norank_f__Oscillospiraceae(-0.48); Christensenella(-0.46); NK4A214_group(-0.45)
- Blautia-taurine r=0.214 p=0.051
- Blautia->taurine->MMSE: 间接效应 0.0313 (95%CI [0.0003,0.0679]), 占 19.5%
