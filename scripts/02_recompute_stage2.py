# -*- coding: utf-8 -*-
"""v15 重算第二阶段：修正模块-代谢物相关、强阈值网络、牛磺酸深挖、模块稳健性"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import warnings
warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(REPO, 'data', 'discovery')
OUT = os.path.join(REPO, 'results', 'discovery')
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261119)

def sp(x, y):
    r = stats.spearmanr(x, y)
    return r.statistic, r.pvalue

def bh(p):
    p = np.asarray(p, dtype=float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)

# ---- 载入 ----
G = pd.read_csv(os.path.join(OUT, 'genus_117_filtered.csv'), index_col=0)
M = pd.read_csv(os.path.join(OUT, 'metab_1873_filtered.csv'), index_col=0)
samples = list(G.columns)
md = pd.read_csv(os.path.join(RAW, 'metab_desc.txt'), sep='\t')
name_map = dict(zip(md['metab_id'], md['Metabolite']))
med = pd.read_csv(os.path.join(RAW, 'mediation_analysis_data.csv'))
med['Group'] = med['Group'].str.strip()
med['stage'] = med['Group'].map({'CG': 0, 'SCD': 1, 'MCI': 2, 'AD': 3})
med = med.set_index('Sample Analysis Name').loc[samples]
stage = med['stage'].values; mmse = med['MMSE'].values.astype(float); grp = med['Group'].values
age = med['Age'].values.astype(float); male = (med['Gender'] == 'Male').astype(float).values
edu = pd.to_numeric(med['Education'], errors='coerce').values

RHO = np.load(os.path.join(OUT, 'rho_matrix.npy'))
P = np.load(os.path.join(OUT, 'p_matrix.npy'))
Q = bh(P.ravel()).reshape(P.shape)

# ---- 模块（k=6 主体 + 稳健性） ----
Rg = np.clip(np.corrcoef(stats.rankdata(G.values, axis=1)), -1, 1)
Z = linkage(squareform(1 - Rg, checks=False), method='average')
lab6 = fcluster(Z, t=6, criterion='maxclust')
Gv = G.values.astype(float)
Gz = (Gv - Gv.mean(axis=1, keepdims=True)) / Gv.std(axis=1, ddof=0, keepdims=True)

def eig_of(idx):
    sub = Gz[idx, :]
    U, S, Vt = np.linalg.svd(sub, full_matrices=False)
    pc1 = Vt[0]
    if np.corrcoef(pc1, sub.mean(axis=0))[0, 1] < 0:
        pc1 = -pc1
    return pc1

idx5 = np.where(lab6 == 5)[0]
eig5 = eig_of(idx5)
m1_genera = list(G.index[idx5])
scfa_all = ['Faecalibacterium', 'Roseburia', 'Butyricicoccus', 'Coprococcus', 'Eubacterium_rectale_group',
            'Eubacterium_hallii_group', 'Anaerostipes', 'Agathobacter', 'Subdoligranulum', 'Lachnospira',
            'Lachnospiraceae_NK4A136_group', 'Anaerotruncus', 'Odoribacter', 'Dorea', 'Butyrivibrio',
            'Ruminococcus_torques_group', 'Lachnoclostridium', 'Flavonifractor']
scfa_set = set(x for x in scfa_all if x in G.index)

print('===== 稳健性：不同切割下的 Faecalibacterium/Blautia 所在模块 r_stage =====')
rob_rows = []
for k in [4, 5, 6, 7, 8, 10]:
    labk = fcluster(Z, t=k, criterion='maxclust')
    for anchor in ['Faecalibacterium']:
        mk = labk[list(G.index).index(anchor)]
        idxk = np.where(labk == mk)[0]
        if len(idxk) < 3:
            continue
        ek = eig_of(idxk)
        r_, p_ = sp(ek, stage)
        rob_rows.append(dict(k=k, anchor=anchor, module_size=len(idxk),
                             n_scfa_in=len(scfa_set & set(G.index[idxk])), r_stage=r_, p=p_))
print(pd.DataFrame(rob_rows).round(4).to_string())
pd.DataFrame(rob_rows).to_csv(os.path.join(OUT, 'module_robustness.csv'), index=False)

# SCFA 在 M5 中的富集（Fisher）
inm5 = set(m1_genera)
tab = [[len(scfa_set & inm5), len(scfa_set - inm5)],
       [len(inm5 - scfa_set), len(set(G.index) - scfa_set - inm5)]]
orf, pf = stats.fisher_exact(tab)
print('SCFA 属在保护模块 M5 内: %d 个; 模块外: %d 个; 非SCFA 模块内: %d; 非SCFA 模块外: %d; Fisher p=%.2g OR=%.1f' % (tab[0][0], tab[0][1], tab[1][0], tab[1][1], pf, orf))

# ---- 修正：模块-代谢物相关（行标准化） ----
print('===== 模块M5特征基因 vs 代谢物（修正rho） =====')
Mr = stats.rankdata(M.values, axis=1).astype(float)
Mrs = (Mr - Mr.mean(axis=1, keepdims=True)) / Mr.std(axis=1, ddof=0, keepdims=True)
er = stats.rankdata(eig5); er = (er - er.mean()) / er.std()
rho_mm = Mrs @ er / (84 - 1)
rho_mm = np.clip(rho_mm, -1, 1)
t_mm = rho_mm * np.sqrt(82 / np.maximum(1 - rho_mm ** 2, 1e-12))
p_mm = 2 * stats.t.sf(np.abs(t_mm), 82)
q_mm = bh(p_mm)
mm_tab = pd.DataFrame({'metab_id': M.index, 'name': [name_map.get(x, x) for x in M.index],
                       'rho': rho_mm, 'p': p_mm, 'FDR': q_mm}).sort_values('p')
mm_tab.to_csv(os.path.join(OUT, 'module_disease_metabolite_assoc.csv'), index=False)
print('FDR<0.05:', (q_mm < 0.05).sum(), '; rho 范围 [%.3f, %.3f]' % (rho_mm.min(), rho_mm.max()))
print('-- Top15 正相关 --'); print(mm_tab[mm_tab.rho > 0].head(15).round(4).to_string())
print('-- Top15 负相关 --'); print(mm_tab[mm_tab.rho < 0].head(15).round(4).to_string())
print('牛磺酸:', mm_tab[mm_tab.metab_id == 'metab_15244'].round(4).to_dict('records'))

# ---- 强阈值网络 ----
print('===== 强阈值网络 =====')
for thr in [0.5, 0.55, 0.6, 0.65, 0.7]:
    n = int(((Q < 0.05) & (np.abs(RHO) >= thr)).sum())
    print('|rho|>=%.2f & FDR<0.05 -> %d 边' % (thr, n))
THR = 0.6
mask = (Q < 0.05) & (np.abs(RHO) >= THR)
ii, jj = np.where(mask)
E = pd.DataFrame({'genus': G.index[ii], 'metab_id': M.index[jj],
                  'name': [name_map.get(x, x) for x in M.index[jj]],
                  'rho': RHO[ii, jj], 'p': P[ii, jj], 'FDR': Q[ii, jj]}).sort_values('rho')
E.to_csv(os.path.join(OUT, 'network_edges_strong.csv'), index=False)
deg = E.groupby('genus').size().sort_values(ascending=False)
print('强网络: %d 边 (正 %d/负 %d), %d 属, %d 代谢物' % (len(E), (E.rho > 0).sum(), (E.rho < 0).sum(), E.genus.nunique(), E.metab_id.nunique()))
print('Top10 hub:'); print(deg.head(10).to_string())
inm5_deg = deg.reindex([g for g in deg.index if g in inm5])
print('M5 模块内属的平均度: %.1f vs 模块外: %.1f' % (inm5_deg.mean(), deg[~deg.index.isin(inm5)].mean()))

# hub 的疾病关联 + 双阈值比较
print('\n-- 强网络(0.6) hub 的 r_stage --')
hub_stage = []
for gname_ in deg.head(10).index:
    r_, p_ = sp(G.loc[gname_].values, stage)
    hub_stage.append(dict(genus=gname_, degree=int(deg[gname_]), r_stage=r_, p_stage=p_))
print(pd.DataFrame(hub_stage).round(4).to_string())
pd.DataFrame(hub_stage).to_csv(os.path.join(OUT, 'hub_stage_assoc.csv'), index=False)

mask55 = (Q < 0.05) & (np.abs(RHO) >= 0.55)
ii5, jj5 = np.where(mask55)
E55 = pd.DataFrame({'genus': G.index[ii5], 'metab_id': M.index[jj5],
                    'name': [name_map.get(x, x) for x in M.index[jj5]],
                    'rho': RHO[ii5, jj5], 'p': P[ii5, jj5], 'FDR': Q[ii5, jj5]})
E55.to_csv(os.path.join(OUT, 'network_edges_thr055.csv'), index=False)
deg55 = E55.groupby('genus').size().sort_values(ascending=False)
print('\n-- 阈值0.55: %d 边 (%d 属), top10 hub --' % (len(E55), E55.genus.nunique()))
hub55 = []
for gname_ in deg55.head(10).index:
    r_, p_ = sp(G.loc[gname_].values, stage)
    inM5 = gname_ in inm5
    hub55.append(dict(genus=gname_, degree=int(deg55[gname_]), r_stage=r_, p_stage=p_, in_M5=inM5))
print(pd.DataFrame(hub55).round(4).to_string())

# 置换（强阈值）
Grank = stats.rankdata(G.values, axis=1).astype(float)
Grank = (Grank - Grank.mean(axis=1, keepdims=True)) / Grank.std(axis=1, ddof=0, keepdims=True)
Mrn = Mrs.T  # 84 x 1873? Mrs is 1873x84 -> transpose
top3_obs = deg.sort_values(ascending=False).head(3).mean()
nperm = 2000
null_top3 = []
for _ in range(nperm):
    Rp = Grank[:, rng.permutation(84)]
    RHO_p = Rp @ Mrn / (84 - 1)
    Tp = RHO_p * np.sqrt(82 / np.maximum(1 - RHO_p ** 2, 1e-12))
    Pp = 2 * stats.t.sf(np.abs(Tp), 82)
    Qp = bh(Pp.ravel()).reshape(RHO_p.shape)
    degp = pd.Series(((Qp < 0.05) & (np.abs(RHO_p) >= THR)).sum(axis=1), index=G.index)
    null_top3.append(degp.sort_values(ascending=False).head(3).mean())
p_perm = np.mean(np.array(null_top3) >= top3_obs)
print('置换: 观测 top3 平均度=%.1f, 零假设=%.1f, p=%.4f (n=%d)' % (top3_obs, np.mean(null_top3), p_perm, nperm))

# ---- 牛磺酸深挖 ----
print('===== 牛磺酸深挖 =====')
tau = M.loc['metab_15244'].values
tau_tab = pd.read_csv(os.path.join(OUT, 'taurine_genus_assoc.csv'))
pos = tau_tab[tau_tab.rho_taurine > 0].head(10)
neg = tau_tab[tau_tab.rho_taurine < 0].head(10)
print('-- 正相关 top10 --'); print(pos.round(4).to_string())
print('-- 负相关 top10 --'); print(neg.round(4).to_string())
bl = G.loc['Blautia'].values
r_bt, p_bt = sp(bl, tau)
r_bm, p_bm = sp(bl, mmse)
r_tm, p_tm = sp(tau, mmse)
print('Blautia-taurine r=%.3f p=%.3g; Blautia-MMSE r=%.3f p=%.3g; taurine-MMSE r=%.3f p=%.3g' % (r_bt, p_bt, r_bm, p_bm, r_tm, p_tm))

# 诚实中介：Blautia -> taurine -> MMSE（OLS + bootstrap 5000）
X = stats.rankdata(bl); Mm_ = stats.rankdata(tau); Y = mmse
Xc = np.column_stack([X - X.mean(), np.ones(84)])
Mc = np.column_stack([Mm_ - Mm_.mean(), np.ones(84)])
a = np.linalg.lstsq(Xc, Mm_, rcond=None)[0][0]
b = np.linalg.lstsq(Mc, Y - Xc @ np.linalg.lstsq(Xc, Y, rcond=None)[0] if False else Y, rcond=None)[0][0]  # placeholder
# 标准三路径
a_coef = np.linalg.lstsq(Xc, Mm_, rcond=None)[0][0]
XcY = np.column_stack([X - X.mean(), Mm_ - Mm_.mean(), np.ones(84)])
coefs = np.linalg.lstsq(XcY, Y, rcond=None)[0]
c_total = np.linalg.lstsq(Xc, Y, rcond=None)[0][0]
c_prime, b_coef = coefs[1], coefs[0]
indirect = a_coef * b_coef
prop = indirect / c_total if c_total != 0 else np.nan
boot = []
for _ in range(5000):
    bidx = rng.integers(0, 84, 84)
    Xb, Mb, Yb = X[bidx], Mm_[bidx], Y[bidx]
    Xb_ = np.column_stack([Xb - Xb.mean(), np.ones(84)])
    ab = np.linalg.lstsq(Xb_, Mb, rcond=None)[0][0]
    XbY = np.column_stack([Xb - Xb.mean(), Mb - Mb.mean(), np.ones(84)])
    cb = np.linalg.lstsq(XbY, Yb, rcond=None)[0]
    boot.append(ab * cb[0])
boot = np.array(boot)
ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
print('Blautia->taurine->MMSE 中介: a=%.3f b=%.3f c=%.3f c\'=%.3f 间接=%.4f 占总效应=%.1f%%' % (a_coef, b_coef, c_total, c_prime, indirect, 100 * prop))
print('bootstrap 95%%CI 间接效应: [%.4f, %.4f] (p=%s)' % (ci_lo, ci_hi, '不含0' if ci_lo * ci_hi > 0 else '含0'))

# M5 特征基因与 taurine
r_et, p_et = sp(eig5, tau)
print('M5特征基因-taurine r=%.3f p=%.3g' % (r_et, p_et))

# ---- 汇总更新 ----
print('===== 写 summary2 =====')
with open(os.path.join(OUT, 'summary_stage2.md'), 'w', encoding='utf-8') as f:
    f.write('# v15 重算·第二阶段结果\n\n')
    f.write('## 模块稳健性（不同k下 Faecalibacterium 所在模块）\n\n')
    f.write(pd.DataFrame(rob_rows).round(4).to_string(index=False) + '\n\n')
    f.write('## SCFA 在 M5 的富集\nFisher p=%.2g (OR=%.1f), %d/%d SCFA 属位于 M5\n\n' % (pf, orf, tab[0][0], len(scfa_set)))
    f.write('## 模块M5-代谢物 top（修正rho）\n\n')
    f.write('正相关 top10:\n' + mm_tab[mm_tab.rho > 0].head(10).round(4).to_string(index=False) + '\n')
    f.write('负相关 top10:\n' + mm_tab[mm_tab.rho < 0].head(10).round(4).to_string(index=False) + '\n\n')
    f.write('## 强阈值网络 (|rho|>=0.6 & FDR<0.05)\n%d 边, top hub: %s\n置换 p=%.4f\n\n' % (
        len(E), ', '.join(f'{g}({d})' for g, d in deg.head(8).items()), p_perm))
    f.write('## 牛磺酸\n- 正相关 top: %s\n- 负相关 top: %s\n- Blautia-taurine r=%.3f p=%.3g\n- Blautia->taurine->MMSE: 间接效应 %.4f (95%%CI [%.4f,%.4f]), 占 %.1f%%\n' % (
        '; '.join(f'{r_.genus}({r_.rho_taurine:.2f})' for r_ in pos.head(5).itertuples()),
        '; '.join(f'{r_.genus}({r_.rho_taurine:.2f})' for r_ in neg.head(5).itertuples()),
        r_bt, p_bt, indirect, ci_lo, ci_hi, 100 * prop))
print('完成')
