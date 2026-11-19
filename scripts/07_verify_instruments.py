# -*- coding: utf-8 -*-
"""从 MiBioGen 原始文件重提全部 33 属工具变量，与存档比对"""
import pandas as pd
import numpy as np
import glob
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base = REPO
mbg = pd.read_csv(os.path.join(REPO, 'data', 'mibiogen', 'MBG.allHits.p1e4.txt'), sep='\t')

def ld_clump(snps_df, window_kb=500):
    selected = []
    for chr_id in snps_df['chr'].unique():
        chr_snps = snps_df[snps_df['chr'] == chr_id].sort_values('P.weightedSumZ')
        kept_positions = []
        for _, row in chr_snps.iterrows():
            bp = row['bp']
            if not kept_positions or all(abs(bp - kp) > window_kb * 1000 for kp in kept_positions):
                kept_positions.append(bp)
                selected.append(row)
    return pd.DataFrame(selected)

arch = {}
for f in glob.glob(os.path.join(REPO, 'data', 'mr_analysis', 'instruments_genus_*.csv')):
    taxon = os.path.basename(f).replace('instruments_', '').replace('.csv', '').replace('_', '.')
    arch[taxon] = pd.read_csv(f)

mismatch = 0
total_re = 0
for taxon, adf in arch.items():
    taxon_df = mbg[mbg['bac'] == taxon]
    sig = taxon_df[taxon_df['P.weightedSumZ'] < 5e-8]
    thr = '5e-8'
    if len(sig) == 0:
        sig = taxon_df[taxon_df['P.weightedSumZ'] < 1e-5]
        thr = '1e-5'
    clumped = ld_clump(sig.copy())
    total_re += len(clumped)
    only_re = sorted(set(clumped['rsID']) - set(adf['rsID']))
    only_ar = sorted(set(adf['rsID']) - set(clumped['rsID']))
    if only_re or only_ar:
        mismatch += 1
        print('MISMATCH %s: 重提%d vs 存档%d; 仅重提:%s; 仅存档:%s' % (taxon, len(clumped), len(adf), only_re, only_ar))

print('33属工具变量重提取完成: 不匹配 %d 个; 重提SNP总数 %d (存档373)' % (mismatch, total_re))
