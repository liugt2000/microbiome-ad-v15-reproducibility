# -*- coding: utf-8 -*-
"""反向MR等位基因协调重算（v15 主口径）
输入: data/mr_analysis/{reverse_mr_snp_matches.csv, ad_gwas_sig_snps_clumped.csv,
      ad_gwas_suggestive_snps_clumped.csv, bifido_mibio_full_snps.csv}
输出: results/verification/reverse_mr_harmonized.csv
注: 原审计脚本中的 Senegalimassilia 原始 zip 抽查（需 GB 级 MiBioGen 分卷）不影响本输出，此处略去。
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mr = os.path.join(REPO, 'data', 'mr_analysis')
OUTDIR = os.path.join(REPO, 'results', 'verification')
os.makedirs(OUTDIR, exist_ok=True)
mat = pd.read_csv(os.path.join(mr, 'reverse_mr_snp_matches.csv'))
bio = pd.read_csv(os.path.join(mr, 'bifido_mibio_full_snps.csv'), index_col=0)
# ===== 用 Bifidobacterium 文件提供的等位基因做协调，重算全部104属反向MR =====
mbg_eff = bio['eff'].str.upper()
mbg_ref = bio['ref'].str.upper()
ad_sig = pd.read_csv(mr + r'\ad_gwas_sig_snps_clumped.csv')
ad_sugg = pd.read_csv(mr + r'\ad_gwas_suggestive_snps_clumped.csv')
ad_all = pd.concat([ad_sig, ad_sugg]).drop_duplicates(subset='MarkerName').set_index('MarkerName')

merged = mat[mat['genus'].notna()].copy()
merged['beta_ad'] = merged['snp'].map(ad_all['Beta'])
merged['se_ad'] = merged['snp'].map(ad_all['SE'])
merged['eff_ad'] = merged['snp'].map(ad_all['Effect_allele'].str.upper())
merged['ref_ad'] = merged['snp'].map(ad_all['Non_Effect_allele'].str.upper())
merged['eff_mbg'] = merged['snp'].map(mbg_eff)
merged = merged.dropna(subset=['beta_ad', 'se_ad', 'beta_microbiome', 'se_microbiome', 'eff_mbg', 'eff_ad'])

# 协调: 若 MiBioGen eff == Kunkle non-effect, 翻转 MiBioGen beta
need_flip = (merged['eff_mbg'] == merged['ref_ad'])
merged['beta_harm'] = np.where(need_flip, -merged['beta_microbiome'], merged['beta_microbiome'])
print('需翻转的 SNP-genus 行:', int(need_flip.sum()), '/', len(merged))

def ivw(bx, by, byse):
    w = 1 / byse ** 2
    beta = np.sum(bx * by * w) / np.sum(bx ** 2 * w)
    se = np.sqrt(1 / np.sum(bx ** 2 * w))
    p = 2 * (1 - stats.norm.cdf(abs(beta / se)))
    return beta, p

res_orig, res_harm = [], []
for genus, g in merged.groupby('genus'):
    if len(g) < 3: continue
    b1, p1 = ivw(g['beta_ad'].values, g['beta_microbiome'].values, g['se_microbiome'].values)
    b2, p2 = ivw(g['beta_ad'].values, g['beta_harm'].values, g['se_microbiome'].values)
    res_orig.append({'genus': genus, 'p_orig': p1})
    res_harm.append({'genus': genus, 'p_harm': p2, 'beta_harm': b2})

o = pd.DataFrame(res_orig); h = pd.DataFrame(res_harm)
comp = o.merge(h, on='genus')
comp['p_adj_harm'] = stats.false_discovery_control(comp['p_harm'], method='bh')
print()
print('协调后 FDR<0.05 属数:', int((comp['p_adj_harm'] < 0.05).sum()))
print('协调后 nominal p<0.05 属数:', int((comp['p_harm'] < 0.05).sum()))
top = comp.nsmallest(10, 'p_harm')
print('协调后 top10:')
print(top[['genus', 'p_orig', 'p_harm', 'p_adj_harm']].to_string(index=False))
b = comp[comp['genus'] == 'Bifidobacterium']
print()
print('Bifidobacterium: 未协调 p=%.4f -> 协调后 p=%.4f (beta=%.4f)' % (b['p_orig'].iloc[0], b['p_harm'].iloc[0], b['beta_harm'].iloc[0]))
comp.to_csv(os.path.join(OUTDIR, 'reverse_mr_harmonized.csv'), index=False)
print('已保存 reverse_mr_harmonized.csv')
