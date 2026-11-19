# -*- coding: utf-8 -*-
"""用重提取的 Kunkle SNP 重跑 30 属 IVW，与存档比对"""
import pandas as pd
import numpy as np
from scipy import stats
import glob
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ws = os.path.join(REPO, 'results', 'verification')
os.makedirs(ws, exist_ok=True)
mr_dir = os.path.join(REPO, 'data', 'mr_analysis')

ad_df = pd.read_csv(os.path.join(mr_dir, 'kunkle_reextracted.csv'))

def mr_ivw(bx, byse, by):
    w = 1 / byse ** 2
    beta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1 / np.sum(w * bx ** 2))
    p = 2 * (1 - stats.norm.cdf(abs(beta / se)))
    return beta, se, p

results = []
for inst_file in sorted(glob.glob(os.path.join(mr_dir, 'instruments_genus_*.csv'))):
    name = os.path.basename(inst_file).replace('instruments_', '').replace('.csv', '')
    exp = pd.read_csv(inst_file)
    if len(exp) < 3:
        continue
    merged = pd.merge(exp, ad_df, left_on='rsID', right_on='MarkerName', how='inner')
    if len(merged) < 3:
        continue
    bx = merged['beta'].values.astype(float)
    by = merged['Beta'].values.astype(float)
    byse = merged['SE_y'].values.astype(float) if 'SE_y' in merged.columns else merged['SE'].values.astype(float)
    exp_eff = merged['eff.allele'].str.upper().values
    exp_ref = merged['ref.allele'].str.upper().values
    out_eff = merged['Effect_allele'].str.upper().values
    out_ref = merged['Non_Effect_allele'].str.upper().values
    flip = np.array([(exp_eff[i] == out_ref[i] and exp_ref[i] == out_eff[i]) for i in range(len(merged))])
    by = np.where(flip, -by, by)
    beta, se, p = mr_ivw(bx, byse, by)
    q = np.sum((1 / byse ** 2) * (by - beta * bx) ** 2)
    q_p = 1 - stats.chi2.cdf(q, len(bx) - 1)
    results.append({'genus': name, 'n_snps': len(merged), 'OR': np.exp(beta),
                    'lo': np.exp(beta - 1.96 * se), 'hi': np.exp(beta + 1.96 * se),
                    'p': p, 'Q_p': q_p})

re_df = pd.DataFrame(results)
arch = pd.read_csv(os.path.join(mr_dir, 'mr_results_final.csv'))
re_df['key'] = re_df['genus'].str.replace('genus_', '', regex=False)
arch['key'] = arch['genus'].str.replace(' ', '', regex=False).str.replace('(id.', '_id_', regex=False).str.replace(')', '', regex=False)
m = arch.merge(re_df, on='key', suffixes=('_arch', '_re'))
print('比对属数:', len(m), '/', len(arch))
print('OR 最大绝对差:', np.abs(m['IVW_OR'] - m['OR']).max())
print('p 最大绝对差:', np.abs(m['IVW_p'] - m['p']).max())
print('n_snps 不一致:', int((m['n_snps_arch'] != m['n_snps_re']).sum()), '行')
b = m[m['key'].str.contains('Bifidobacterium')].iloc[0]
print('Bifidobacterium 重跑: OR=%.4f [%.4f, %.4f] p=%.4g Q_p=%.4f n=%d' % (b['OR'], b['lo'], b['hi'], b['p'], b['Q_p'], b['n_snps_re']))
print('Bifidobacterium 存档: OR=%.4f [%.4f, %.4f] p=%.4g Q_p=%.4f n=%d' % (b['IVW_OR'], b['IVW_OR_lower'], b['IVW_OR_upper'], b['IVW_p'], b['Q_pval'], b['n_snps_arch']))
sig = m[m['p'] < 0.05]
print('IVW p<0.05 (重跑):', sig[['key', 'OR', 'p']].values.tolist())
re_df.to_csv(os.path.join(ws, 'repro_mr_forward.csv'), index=False)
print('已保存 repro_mr_forward.csv')
