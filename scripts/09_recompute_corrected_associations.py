# -*- coding: utf-8 -*-
"""Correct rank-correlation scaling in legacy discovery outputs.
Legacy code divided ddof=0 standardized cross-products by n-1, multiplying every
correlation by n/(n-1). This script uses /n and preserves all legacy files.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

REPO=Path(__file__).resolve().parents[1]
RES=REPO/'results'/'discovery'; RAW=REPO/'data'/'discovery'
OUT=REPO/'results'/'corrected_discovery'; OUT.mkdir(parents=True,exist_ok=True)
G=pd.read_csv(RES/'genus_117_filtered.csv',index_col=0)
M=pd.read_csv(RES/'metab_1873_filtered.csv',index_col=0)
md=pd.read_csv(RAW/'metab_desc.txt',sep='\t'); names=dict(zip(md.metab_id,md.Metabolite))
med=pd.read_csv(RAW/'mediation_analysis_data.csv'); med['Group']=med.Group.str.strip(); med=med.set_index('Sample Analysis Name').loc[G.columns]
stage=med.Group.map({'CG':0,'SCD':1,'MCI':2,'AD':3}).to_numpy(); mmse=med.MMSE.to_numpy(float)
n=G.shape[1]
def bh(p):
    p=np.asarray(p,float); m=p.size; o=np.argsort(p); q=np.empty(m)
    q[o]=np.minimum.accumulate((p[o]*m/np.arange(1,m+1))[::-1])[::-1]
    return np.minimum(q,1)
def zp(a):
    a=stats.rankdata(a,axis=1).astype(float)
    return (a-a.mean(axis=1,keepdims=True))/a.std(axis=1,ddof=0,keepdims=True)
def corrp(r):
    t=r*np.sqrt((n-2)/np.maximum(1-r*r,1e-12)); return 2*stats.t.sf(np.abs(t),n-2)
Gr=zp(G.to_numpy()); Mr=zp(M.to_numpy())
RHO=np.clip(Gr@Mr.T/n,-1,1); P=corrp(RHO); Q=bh(P.ravel()).reshape(P.shape)
np.save(OUT/'rho_matrix_corrected.npy',RHO); np.save(OUT/'p_matrix_corrected.npy',P); np.save(OUT/'fdr_matrix_corrected.npy',Q)
# Exact legacy comparison after removing the known n/(n-1) multiplier.
oldR=np.load(RES/'rho_matrix.npy')
maxdiff=float(np.max(np.abs(oldR*(n-1)/n-RHO)))
# Modules and trait associations.
mem=pd.read_csv(RES/'modules_membership.csv'); Gv=G.to_numpy(float); Gz=(Gv-Gv.mean(1,keepdims=True))/Gv.std(1,ddof=0,keepdims=True)
eigs={}; trait=[]
for mod in sorted(mem.module.unique()):
    gs=mem.loc[mem.module==mod,'genus'].tolist(); sub=Gz[[list(G.index).index(g) for g in gs]]
    _,_,vt=np.linalg.svd(sub,full_matrices=False); e=vt[0]
    if np.corrcoef(e,sub.mean(0))[0,1]<0:e=-e
    eigs[mod]=e; rs,ps=stats.spearmanr(e,stage); rm,pm=stats.spearmanr(e,mmse)
    trait.append({'module':mod,'n_genera':len(gs),'r_stage':rs,'p_stage':ps,'r_MMSE':rm,'p_MMSE':pm})
pd.DataFrame(trait).to_csv(OUT/'module_trait_correlations_corrected.csv',index=False)
e5=eigs[5]; er=stats.rankdata(e5); er=(er-er.mean())/er.std(ddof=0)
# Both rank vectors must be Z-standardized for a Pearson correlation of ranks.
rmm=np.clip(Mr@er/n,-1,1); pmm=corrp(rmm); qmm=bh(pmm)
mm=pd.DataFrame({'metab_id':M.index,'name':[names.get(x,x) for x in M.index],'rho':rmm,'p':pmm,'FDR':qmm}).sort_values('p')
mm.to_csv(OUT/'module_disease_metabolite_assoc_corrected.csv',index=False)
# Correct taurine-genus associations.
tau=M.loc['metab_15244'].to_numpy(float); tr=stats.rankdata(tau); tr=(tr-tr.mean())/tr.std(ddof=0)
# Independent SciPy identity check prevents silent matrix/indexing errors.
tau_pos=list(M.index).index('metab_15244')
assert abs(rmm[tau_pos]-stats.spearmanr(e5,tau).statistic)<1e-12
rt=np.clip(Gr@tr/n,-1,1); pt=corrp(rt); qt=bh(pt)
pd.DataFrame({'genus':G.index,'rho_taurine':rt,'p':pt,'FDR':qt}).sort_values('p').to_csv(OUT/'taurine_genus_assoc_corrected.csv',index=False)
# Strong network and all thresholds.
thresholds=[]
for th in [.5,.55,.6,.65,.7]:
    mask=(Q<.05)&(np.abs(RHO)>=th)
    thresholds.append({'threshold':th,'n_edges':int(mask.sum()),'n_genera':int((mask.sum(1)>0).sum()),'n_metabolites':int((mask.sum(0)>0).sum())})
pd.DataFrame(thresholds).to_csv(OUT/'network_threshold_counts.csv',index=False)
mask=(Q<.05)&(np.abs(RHO)>=.6); ii,jj=np.where(mask)
edges=pd.DataFrame({'genus':G.index[ii],'metab_id':M.index[jj],'name':[names.get(x,x) for x in M.index[jj]],'rho':RHO[ii,jj],'p':P[ii,jj],'FDR':Q[ii,jj]}).sort_values('p')
edges.to_csv(OUT/'network_edges_strong_corrected.csv',index=False)
deg=edges.groupby('genus').size().sort_values(ascending=False)
# Corrected 2000-permutation hub test, same seed and permutation statistic.
rng=np.random.default_rng(20261119); top3obs=deg.head(3).mean(); null=[]
for _ in range(2000):
    rp=Gr[:,rng.permutation(n)]; rr=np.clip(rp@Mr.T/n,-1,1); pp=corrp(rr); qq=bh(pp.ravel()).reshape(rr.shape)
    dp=pd.Series(((qq<.05)&(np.abs(rr)>=.6)).sum(1),index=G.index)
    null.append(dp.sort_values(ascending=False).head(3).mean())
null=np.asarray(null); pperm=float(np.mean(null>=top3obs)); np.save(OUT/'hub_permutation_null.npy',null)
summary={'n_samples':n,'correction':'ddof=0 standardized cross-products divided by n, not n-1','legacy_scaling_factor':n/(n-1),'max_abs_rho_difference_after_rescaling_legacy':maxdiff,
 'module5_taurine':mm[mm.metab_id=='metab_15244'].iloc[0].to_dict(),
 'Blautia_taurine':{'rho':float(rt[list(G.index).index('Blautia')]),'p':float(pt[list(G.index).index('Blautia')]),'FDR':float(qt[list(G.index).index('Blautia')])},
 'strong_edges':len(edges),'positive_edges':int((edges.rho>0).sum()),'negative_edges':int((edges.rho<0).sum()),'genera':int(edges.genus.nunique()),'metabolites':int(edges.metab_id.nunique()),
 'hub_top3_observed_mean_degree':float(top3obs),'hub_permutation_p':pperm,'n_permutations':2000,'seed':20261119}
(OUT/'corrected_association_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print(pd.DataFrame(thresholds).to_string(index=False))
