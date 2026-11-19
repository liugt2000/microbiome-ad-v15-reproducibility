# -*- coding: utf-8 -*-
"""v15 发现队列全量重算：模块、网络、置换检验、guild、牛磺酸
输入全部为 SHA256 已验证的原始数据；输出到 07_v15_Recomputed_Results/
"""
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
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261119)

def sp(x, y):
    r = stats.spearmanr(x, y)
    return r.statistic, r.pvalue

def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)

# ============ 1. 数据准备 ============
print('=' * 60)
print('1. 数据准备')
gf = pd.read_csv(os.path.join(RAW, 'otu_taxon_Genus.full.xls'), sep='\t')
lineage = gf[['phylum', 'class', 'order', 'family', 'genus']].copy()
lineage['genus'] = lineage['genus'].str.replace('g__', '', regex=False)
phyl_map = dict(zip(lineage['genus'], lineage['phylum'].str.replace('p__', '', regex=False)))

g = pd.read_csv(os.path.join(RAW, 'otu_taxon_Genus.percent.xls'), sep='\t', index_col=0)
g = g.loc[[x for x in g.index if isinstance(x, str)]]
g.index = [x.replace('g__', '') for x in g.index]
samples = list(g.columns)
prev = (g > 0).sum(axis=1); meanab = g.mean(axis=1)
keep = g.index[(prev >= 0.2 * 84) & (meanab >= 0.0001)]
G = g.loc[keep]                       # 117 x 84
print('过滤后属数:', G.shape[0], '; 样本:', G.shape[1])

md = pd.read_csv(os.path.join(RAW, 'metab_desc.txt'), sep='\t')
ma = pd.read_csv(os.path.join(RAW, 'metab_abund.txt'), sep='\t', index_col=0)
ma = ma[[c for c in ma.columns if not str(c).startswith('QC')]]
met_ids = md[(md['level'] == 'B(i)') & (md['compound_source'] != 'nonnatural')]['metab_id']
met_ids = [m for m in met_ids if m in ma.index]
M = ma.loc[met_ids]                   # 1873 x 84
M = M[samples]                        # 对齐样本顺序
name_map = dict(zip(md['metab_id'], md['Metabolite']))
kegg_map = dict(zip(md['metab_id'], md['KEGG Compound ID']))
print('过滤后代谢物:', M.shape[0])

med = pd.read_csv(os.path.join(RAW, 'mediation_analysis_data.csv'))
med['Group'] = med['Group'].str.strip()
stage_map = {'CG': 0, 'SCD': 1, 'MCI': 2, 'AD': 3}
med['stage'] = med['Group'].map(stage_map)
med = med.set_index('Sample Analysis Name').loc[samples]
stage = med['stage'].values
mmse = med['MMSE'].values.astype(float)
grp = med['Group'].values
age = med['Age'].values.astype(float)
male = (med['Gender'] == 'Male').astype(float).values
edu = pd.to_numeric(med['Education'], errors='coerce').values
print('分组:', pd.Series(grp).value_counts().to_dict())

# 保存过滤矩阵
G.to_csv(os.path.join(OUT, 'genus_117_filtered.csv'))
M.to_csv(os.path.join(OUT, 'metab_1873_filtered.csv'))

# ============ 2. 共丰度模块 ============
print('=' * 60)
print('2. 共丰度模块（117 属，average linkage on 1-r, maxclust k=6, n=84）')
Rg = np.corrcoef(stats.rankdata(G.values, axis=1))  # 117x117 spearman
Rg = np.clip(Rg, -1, 1)
Dm = squareform(1 - Rg, checks=False)
Z = linkage(Dm, method='average')
# 自动选切割高度：在 0.30-0.80 之间取"模块数最多且各模块>=5 属"的高度（平局取较小高度）
best = None
for h in np.arange(0.30, 0.801, 0.05):
    lab_h = fcluster(Z, t=h, criterion='distance')
    sizes = pd.Series(lab_h).value_counts()
    n_mod = (sizes >= 5).sum()
    if n_mod >= 2 and sizes.min() >= 5:
        cand = (n_mod, -h)
        if best is None or cand > best[0]:
            best = (cand, h, lab_h)
    print('  h=%.2f sizes=%s' % (h, sorted(sizes.tolist(), reverse=True)))
if best is None:
    lab = fcluster(Z, t=6, criterion='maxclust')
    cut_rule = 'maxclust k=6 (fallback)'
else:
    lab = best[2]; cut_rule = 'average linkage, cut height h=%.2f' % best[1]
sizes = pd.Series(lab).value_counts().sort_index()
print('采用切割规则:', cut_rule, '-> 模块大小:', sizes.tolist())

Gv = G.values.astype(float)
Gz = (Gv - Gv.mean(axis=1, keepdims=True)) / Gv.std(axis=1, ddof=0, keepdims=True)
mod_rows, eigengenes = [], {}
for m in sorted(set(lab)):
    idx = np.where(lab == m)[0]
    if len(idx) < 3:
        continue
    sub = Gz[idx, :]                        # k x 84
    U, S, Vt = np.linalg.svd(sub, full_matrices=False)
    pc1 = Vt[0]
    # 方向：与模块成员平均 z 同号
    if np.corrcoef(pc1, sub.mean(axis=0))[0, 1] < 0:
        pc1 = -pc1
    eig = pc1
    r_st, p_st = sp(eig, stage)
    r_mm, p_mm = sp(eig, mmse)
    F, p_an = stats.f_oneway(*[eig[grp == gp_] for gp_ in ['CG', 'SCD', 'MCI', 'AD']])
    means = {gp_: float(np.mean(eig[grp == gp_])) for gp_ in ['CG', 'SCD', 'MCI', 'AD']}
    phyla = [phyl_map.get(x, '?') for x in G.index[idx]]
    firm_frac = np.mean([p_ == 'Firmicutes' for p_ in phyla])
    mod_rows.append(dict(module=m, n_genera=len(idx), r_stage=r_st, p_stage=p_st,
                         r_MMSE=r_mm, p_MMSE=p_mm, anova_F=F, anova_p=p_an,
                         firmicutes_frac=firm_frac, **{f'mean_{gp_}': means[gp_] for gp_ in ['CG', 'SCD', 'MCI', 'AD']}))
    eigengenes[m] = eig
modtab = pd.DataFrame(mod_rows)
modtab.to_csv(os.path.join(OUT, 'module_trait_correlations.csv'), index=False)
print(modtab.round(4).to_string())

# 疾病模块 = |r_stage| 最大
m1 = modtab.loc[modtab['r_stage'].abs().idxmax(), 'module']
idx1 = np.where(lab == m1)[0]
m1_genera = list(G.index[idx1])
print('\n疾病关联模块 = M%s (%d 属, r_stage=%.3f, p=%.2g)' % (m1, len(idx1), modtab.loc[modtab.module == m1, 'r_stage'].iloc[0], modtab.loc[modtab.module == m1, 'p_stage'].iloc[0]))
print('构成:', m1_genera)
print('Firmicutes 占比: %.0f%%' % (100 * np.mean([phyl_map.get(x, '?') == 'Firmicutes' for x in m1_genera])))

# Tukey HSD（疾病模块）
try:
    tk = stats.tukey_hsd(*[eigengenes[m1][grp == gp_] for gp_ in ['CG', 'SCD', 'MCI', 'AD']])
    gnames = ['CG', 'SCD', 'MCI', 'AD']
    tukey_p = {f'{gnames[i]}-{gnames[j]}': float(tk.pvalue[i, j]) for i in range(4) for j in range(i + 1, 4)}
    print('Tukey:', {k_: round(v, 4) for k_, v in tukey_p.items()})
except Exception as e:
    tukey_p = {}
    print('Tukey 不可用:', e)

# 模块成员表（含谱系）
mem = pd.DataFrame({'genus': G.index, 'phylum': [phyl_map.get(x, '?') for x in G.index], 'module': lab})
mem.to_csv(os.path.join(OUT, 'modules_membership.csv'), index=False)

# 协变量调整（年龄/性别/教育）的偏 Spearman：秩残差
def partial_spearman(y, x, covs):
    ry = stats.rankdata(y); rx = stats.rankdata(x)
    C = np.column_stack([stats.rankdata(c) if c.dtype.kind == 'f' else c for c in covs] + [np.ones(len(y))])
    by = np.linalg.lstsq(C, ry, rcond=None)[0]
    bx = np.linalg.lstsq(C, rx, rcond=None)[0]
    ey, ex = ry - C @ by, rx - C @ bx
    return sp(ey, ex)

covs = [age, male, np.nan_to_num(edu.astype(float))]
r_adj, p_adj = partial_spearman(eigengenes[m1], stage.astype(float), covs)
print('调整年龄/性别/教育后 r_stage=%.3f (p=%.3g)' % (r_adj, p_adj))

# ============ 3. 疾病模块-代谢物关联 ============
print('=' * 60)
print('3. 疾病模块特征基因 vs 1873 代谢物')
Mr = stats.rankdata(M.values, axis=1)
er = stats.rankdata(eigengenes[m1])
er = (er - er.mean()) / er.std()
Mc = Mr - Mr.mean(axis=1, keepdims=True)
rho_mm = (Mc @ er) / (84 - 1)
t_mm = rho_mm * np.sqrt(82 / np.maximum(1 - rho_mm ** 2, 1e-12))
p_mm = 2 * stats.t.sf(np.abs(t_mm), 82)
q_mm = bh(p_mm)
mm_tab = pd.DataFrame({'metab_id': M.index, 'name': [name_map.get(x, x) for x in M.index],
                       'rho': rho_mm, 'p': p_mm, 'FDR': q_mm})
mm_tab.sort_values('p').to_csv(os.path.join(OUT, 'module_disease_metabolite_assoc.csv'), index=False)
print('FDR<0.05 关联代谢物数:', (q_mm < 0.05).sum())
print(mm_tab.sort_values('p').head(15).round(4).to_string())
tau_row = mm_tab[mm_tab['metab_id'] == 'metab_15244']
print('牛磺酸与疾病模块特征基因:', tau_row.round(4).to_dict('records'))

# ============ 4. 微生物-代谢物网络 ============
print('=' * 60)
print('4. 网络（117 属 x 1873 代谢物，Spearman，n=84）')
Gr = stats.rankdata(G.values, axis=1)
Gr = (Gr - Gr.mean(axis=1, keepdims=True)) / Gr.std(axis=1, ddof=0, keepdims=True)
Mrn = (Mr - Mr.mean(axis=1, keepdims=True)) / Mr.std(axis=1, ddof=0, keepdims=True)
RHO = Gr @ Mrn.T / (84 - 1)             # 117 x 1873
RHO = np.clip(RHO, -1, 1)
T = RHO * np.sqrt(82 / np.maximum(1 - RHO ** 2, 1e-12))
P = 2 * stats.t.sf(np.abs(T), 82)
Q = bh(P.ravel()).reshape(P.shape)
print('FDR<0.05 边数:', int((Q < 0.05).sum()))
print('nominal p<0.05 且 |rho|>=0.3 边数:', int(((P < 0.05) & (np.abs(RHO) >= 0.3)).sum()))
np.save(os.path.join(OUT, 'rho_matrix.npy'), RHO)
np.save(os.path.join(OUT, 'p_matrix.npy'), P)

for tag, mask in [('fdr', Q < 0.05), ('nominal', (P < 0.05) & (np.abs(RHO) >= 0.3))]:
    ii, jj = np.where(mask)
    edges = pd.DataFrame({'genus': G.index[ii], 'metab_id': M.index[jj],
                          'name': [name_map.get(x, x) for x in M.index[jj]],
                          'rho': RHO[ii, jj], 'p': P[ii, jj], 'FDR': Q[ii, jj]})
    edges = edges.sort_values('p')
    edges.to_csv(os.path.join(OUT, f'network_edges_{tag}.csv'), index=False)

edges_f = pd.read_csv(os.path.join(OUT, 'network_edges_fdr.csv'))
if len(edges_f) >= 30:
    E = edges_f; tag_use = 'FDR<0.05'
else:
    E = pd.read_csv(os.path.join(OUT, 'network_edges_nominal.csv')); tag_use = 'nominal p<0.05 & |rho|>=0.3'
print('主网络口径:', tag_use, '边数:', len(E))
deg = E.groupby('genus').size().sort_values(ascending=False)
pos = E[E['rho'] > 0].groupby('genus').size()
print('度数 top10:')
print(pd.DataFrame({'degree': deg, 'positive': pos}).head(10).to_string())
print('正边:', (E['rho'] > 0).sum(), '负边:', (E['rho'] < 0).sum(), '; 属数:', E['genus'].nunique(), '; 代谢物数:', E['metab_id'].nunique())

# ============ 5. SCFA/H2/H2S guild ============
print('=' * 60)
print('5. guild 分析')
scfa_all = ['Faecalibacterium', 'Roseburia', 'Butyricicoccus', 'Coprococcus', 'Eubacterium_rectale_group',
            'Eubacterium_hallii_group', 'Anaerostipes', 'Agathobacter', 'Subdoligranulum', 'Lachnospira',
            'Lachnospiraceae_NK4A136_group', 'Anaerotruncus', 'Odoribacter', 'Dorea', 'Butyrivibrio',
            'Ruminococcus_torques_group', 'Lachnoclostridium', 'Flavonifractor']
h2_prod = ['Ruminococcus', 'Ruminococcus_gnavus_group', 'Ruminococcus_torques_group', 'Blautia', 'Dorea',
           'Clostridium_sensu_stricto_1', 'Eubacterium_rectale_group', 'Eubacterium_hallii_group',
           'Bifidobacterium', 'Bacteroides', 'Prevotella', 'Enterobacter', 'Klebsiella', 'Escherichia-Shigella']
h2s_prod = ['Desulfovibrio', 'Bilophila', 'Fusobacterium']
acetogen = ['Blautia', 'Eubacterium_hallii_group', 'Acetobacterium']

def guild_score(lst):
    members = [x for x in lst if x in G.index]
    return G.loc[members].sum(axis=0).values, members

guild_rows = []
for gname, gl in [('SCFA_producers', scfa_all), ('H2_producers', h2_prod), ('H2S_producers', h2s_prod), ('acetogens', acetogen)]:
    s, members = guild_score(gl)
    r_st, p_st = sp(s, stage)
    r_mm, p_mm = sp(s, mmse)
    H, p_kw = stats.kruskal(*[s[grp == gp_] for gp_ in ['CG', 'SCD', 'MCI', 'AD']])
    U, p_mwu = stats.mannwhitneyu(s[grp == 'AD'], s[grp == 'CG'])
    r_adjg, p_adjg = partial_spearman(s, stage.astype(float), covs)
    guild_rows.append(dict(guild=gname, n_members=len(members), members=';'.join(members),
                           r_stage=r_st, p_stage=p_st, r_MMSE=r_mm, p_MMSE=p_mm,
                           KW_p=p_kw, ADvsCG_MWU_p=p_mwu, r_stage_adjusted=r_adjg, p_stage_adjusted=p_adjg))
guild_tab = pd.DataFrame(guild_rows)
guild_tab.to_csv(os.path.join(OUT, 'guild_analysis.csv'), index=False)
print(guild_tab.round(4).to_string())

# 关键单个属（AD vs CG）
print('\n关键属 AD vs CG（Mann-Whitney）:')
key_genera = ['Faecalibacterium', 'Roseburia', 'Coprococcus', 'Butyricicoccus', 'Blautia',
              'Bifidobacterium', 'Desulfovibrio', 'Bilophila', 'Akkermansia']
rows_kg = []
for gn in key_genera:
    if gn in G.index:
        v = G.loc[gn].values
        U, p = stats.mannwhitneyu(v[grp == 'AD'], v[grp == 'CG'])
        r_st, p_st = sp(v, stage)
        rows_kg.append(dict(genus=gn, median_CG=np.median(v[grp == 'CG']), median_AD=np.median(v[grp == 'AD']),
                            ADvsCG_p=p, r_stage=r_st, p_stage=p_st))
pd.DataFrame(rows_kg).to_csv(os.path.join(OUT, 'key_genera_ADvsCG.csv'), index=False)
print(pd.DataFrame(rows_kg).round(4).to_string())

# ============ 6. 牛磺酸分析 ============
print('=' * 60)
print('6. 牛磺酸（metab_15244）')
tau = M.loc['metab_15244'].values
H, p_kw = stats.kruskal(*[tau[grp == gp_] for gp_ in ['CG', 'SCD', 'MCI', 'AD']])
U, p_mwu = stats.mannwhitneyu(tau[grp == 'AD'], tau[grp == 'CG'])
r_st, p_st = sp(tau, stage)
r_mm, p_mm = sp(tau, mmse)
med_by = {gp_: float(np.median(tau[grp == gp_])) for gp_ in ['CG', 'SCD', 'MCI', 'AD']}
print('分组中位数:', {k_: round(v, 3) for k_, v in med_by.items()})
print('KW p=%.4f; AD vs CG MWU p=%.4f; r_stage=%.3f (p=%.3f); r_MMSE=%.3f (p=%.3f)' % (p_kw, p_mwu, r_st, p_st, r_mm, p_mm))

# 牛磺酸与 117 属的相关（BH）
tau_r = Grn_saved = None
er2 = stats.rankdata(tau); er2 = (er2 - er2.mean()) / er2.std()
rho_t = Mrn_saved = (Mc_gr := (stats.rankdata(G.values, axis=1) - stats.rankdata(G.values, axis=1).mean(axis=1, keepdims=True)))
rho_t = (rho_t / rho_t.std(axis=1, ddof=0, keepdims=True)) @ er2 / (84 - 1)
t_t = rho_t * np.sqrt(82 / np.maximum(1 - rho_t ** 2, 1e-12))
p_t = 2 * stats.t.sf(np.abs(t_t), 82)
q_t = bh(p_t)
tau_tab = pd.DataFrame({'genus': G.index, 'rho_taurine': rho_t, 'p': p_t, 'FDR': q_t}).sort_values('p')
tau_tab.to_csv(os.path.join(OUT, 'taurine_genus_assoc.csv'), index=False)
print('FDR<0.05 属数:', (q_t < 0.05).sum())
print(tau_tab.head(10).round(4).to_string())

# ============ 7. 置换检验 ============
print('=' * 60)
print('7. 置换检验（2000 次属表样本标签置换）')
top3_obs = deg.head(3).mean()
scfa_set = set(x for x in scfa_all if x in G.index)
scfa_edge_frac_obs = E['genus'].isin(scfa_set).mean()
nperm = 2000
null_top3, null_frac = [], []
Grank = stats.rankdata(G.values, axis=1).astype(float)
Grank = (Grank - Grank.mean(axis=1, keepdims=True)) / Grank.std(axis=1, ddof=0, keepdims=True)
for _ in range(nperm):
    perm = rng.permutation(84)
    Rp = Grank[:, perm]
    RHO_p = Rp @ Mrn.T / (84 - 1)
    Tp = RHO_p * np.sqrt(82 / np.maximum(1 - RHO_p ** 2, 1e-12))
    Pp = 2 * stats.t.sf(np.abs(Tp), 82)
    maskp = (Pp < 0.05) & (np.abs(RHO_p) >= 0.3)
    degp = pd.Series(maskp.sum(axis=1), index=G.index)
    null_top3.append(degp.head(3).mean() if False else degp.sort_values(ascending=False).head(3).mean())
    ii_p, jj_p = np.where(maskp)
    if len(ii_p) > 0:
        null_frac.append(np.mean([G.index[i] in scfa_set for i in ii_p]))
    else:
        null_frac.append(0.0)
null_top3 = np.array(null_top3); null_frac = np.array(null_frac)
p_perm_top3 = np.mean(null_top3 >= top3_obs)
p_perm_frac = np.mean(null_frac >= scfa_edge_frac_obs)
print('观测 top3 平均度=%.1f; 置换均值=%.1f; p=%.3f' % (top3_obs, null_top3.mean(), p_perm_top3))
print('观测 SCFA 边占比=%.3f; 置换均值=%.3f; p=%.4f' % (scfa_edge_frac_obs, null_frac.mean(), p_perm_frac))

# Fisher：hub 属中 SCFA 产生菌富集
top_hubs = deg.head(8).index.tolist()
in_scfa = [x in scfa_set for x in top_hubs]
all_in = [x in scfa_set for x in G.index]
table = [[sum(in_scfa), len(in_scfa) - sum(in_scfa)],
         [sum(all_in) - sum(in_scfa), (len(all_in) - sum(all_in)) - (len(in_scfa) - sum(in_scfa))]]
orv, p_fisher = stats.fisher_exact(table)
print('Top8 hub 中 SCFA 属: %d/%d; 全体 %d/%d; Fisher p=%.4f' % (sum(in_scfa), len(in_scfa), sum(all_in), len(all_in), p_fisher))

# ============ 8. 汇总 ============
print('=' * 60)
print('8. 全部关键数字汇总（写入 summary.md）')
m1row = modtab[modtab.module == m1].iloc[0]
lines = []
lines.append('# v15 发现队列重算结果汇总（n=84，117 属 x 1873 代谢物）\n')
lines.append('## 过滤口径\n- 属：检出率>=20%（>=17/84 样本）且平均相对丰度>=0.01%，347->117 属\n- 代谢物：MSI level B(i) 且非 nonnatural（药物/外源），14678->1873\n- 样本：全部 84 例（CG 10 / SCD 22 / MCI 23 / AD 29）\n')
lines.append('## 共丰度模块（average linkage on 1-Spearman, maxclust k=6）\n')
lines.append(modtab.round(4).to_string(index=False) + '\n')
lines.append('- 疾病关联模块 M%s：%d 属，r_stage=%.3f (p=%.2g)，r_MMSE=%.3f (p=%.2g)，ANOVA p=%.2g' % (m1, len(idx1), m1row.r_stage, m1row.p_stage, m1row.r_MMSE, m1row.p_MMSE, m1row.anova_p))
lines.append('- 调整年龄/性别/教育后 r_stage=%.3f (p=%.3g)\n' % (r_adj, p_adj))
lines.append('## 网络\n- 主口径：%s，共 %d 条边（正 %d / 负 %d），涉及 %d 属 / %d 代谢物\n' % (tag_use, len(E), (E['rho'] > 0).sum(), (E['rho'] < 0).sum(), E['genus'].nunique(), E['metab_id'].nunique()))
lines.append('- Top hub：' + ', '.join(f'{g_}({d_})' for g_, d_ in deg.head(8).items()) + '\n')
lines.append('- 置换检验：top3 平均度 p=%.3f；SCFA 边占比 p=%.4f；Fisher(hub SCFA 富集) p=%.4f\n' % (p_perm_top3, p_perm_frac, p_fisher))
lines.append('## guild\n')
lines.append(guild_tab.round(4).to_string(index=False) + '\n')
lines.append('## 牛磺酸\n- 分组中位数 CG=%.3f SCD=%.3f MCI=%.3f AD=%.3f；KW p=%.3f；AD vs CG p=%.3f；r_stage=%.3f (p=%.3f)；r_MMSE=%.3f (p=%.3f)\n' % (med_by['CG'], med_by['SCD'], med_by['MCI'], med_by['AD'], p_kw, p_mwu, r_st, p_st, r_mm, p_mm))
lines.append('- FDR<0.05 的牛磺酸相关属 %d 个；top：%s\n' % ((q_t < 0.05).sum(), '; '.join(f'{r_.genus}({r_.rho_taurine:.2f})' for r_ in tau_tab.head(5).itertuples())))
with open(os.path.join(OUT, 'summary.md'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print('完成。输出目录:', OUT)
