# -*- coding: utf-8 -*-
"""
Level-2 复现验证：从分队列效应量重算荟萃分析，与 results/meta_analysis 及 v14 Source Data 比对
覆盖：属级荟萃(REML/DL)、α多样性荟萃、功能通路荟萃、氢代谢guild荟萃
"""
import numpy as np
import pandas as pd
from scipy import stats
import os, sys, json

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "data", "source_data")
RES = os.path.join(REPO, "data", "meta_analysis")
OUT = os.path.join(REPO, "results", "verification")
os.makedirs(OUT, exist_ok=True)

def dl_meta(yi, sei):
    """DerSimonian-Laird随机效应，返回(pooled, se, p, i2, Q, Qp)"""
    yi = np.asarray(yi, float); sei = np.asarray(sei, float)
    w = 1.0 / sei**2
    fe = np.sum(w*yi) / np.sum(w)
    Q = np.sum(w*(yi-fe)**2)
    df = len(yi)-1
    c = np.sum(w) - np.sum(w**2)/np.sum(w)
    tau2 = max(0.0, (Q-df)/c) if c > 0 else 0.0
    wr = 1.0/(sei**2 + tau2)
    mu = np.sum(wr*yi)/np.sum(wr)
    se_mu = np.sqrt(1.0/np.sum(wr))
    z = mu/se_mu
    p = 2*(1-stats.norm.cdf(abs(z)))
    i2 = max(0.0, (Q-df)/Q*100) if Q > 0 else 0.0
    Qp = 1-stats.chi2.cdf(Q, df) if df > 0 else np.nan
    return mu, se_mu, p, i2, Q, Qp

def reml_tau2(yi, vi):
    """metafor式REML tau2（Fisher scoring等价的方程求根）"""
    lo, hi = 0.0, 1e6
    def eq(t):
        w = 1.0/(vi+t)
        mu = np.sum(w*yi)/np.sum(w)
        return np.sum(w**2*((yi-mu)**2 - (vi+t)))
    f_lo, f_hi = eq(lo), eq(hi)
    if f_lo <= 0: return 0.0
    if f_hi > 0: return hi
    for _ in range(200):
        mid = (lo+hi)/2
        if eq(mid) > 0: lo = mid
        else: hi = mid
    return (lo+hi)/2

def reml_meta(yi, sei):
    yi = np.asarray(yi, float); sei = np.asarray(sei, float)
    vi = sei**2
    t2 = reml_tau2(yi, vi)
    w = 1.0/(vi+t2)
    mu = np.sum(w*yi)/np.sum(w)
    se_mu = np.sqrt(1.0/np.sum(w))
    z = mu/se_mu
    p = 2*(1-stats.norm.cdf(abs(z)))
    # metafor I2 基于DL的Q
    wf = 1.0/vi
    fe = np.sum(wf*yi)/np.sum(wf)
    Q = np.sum(wf*(yi-fe)**2); df = len(yi)-1
    i2 = max(0.0, (Q-df)/Q*100) if Q > 0 else 0.0
    return mu, se_mu, p, i2

def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); rk = np.empty(n); rk[o] = np.arange(1, n+1)
    q = p*n/rk
    q = np.minimum.accumulate(q[o])[np.argsort(o)]
    return np.minimum(q, 1.0)

report = {}

# ============ 1. 属级荟萃分析 ============
print("="*70); print("1. 属级荟萃分析: maaslin2_all_results -> genus_meta_analysis")
ma = pd.read_csv(os.path.join(SRC, "maaslin2_all_results.csv"))
gm = pd.read_csv(os.path.join(SRC, "genus_meta_analysis.csv"))
rows = []
for g, sub in ma.groupby("feature"):
    if len(sub) < 2: continue
    mu_dl, se_dl, p_dl, i2_dl, Q, Qp = dl_meta(sub["coef"], sub["stderr"])
    mu_rm, se_rm, p_rm, i2_rm = reml_meta(sub["coef"], sub["stderr"])
    rows.append(dict(genus=g, n=len(sub), dl_coef=mu_dl, dl_p=p_dl, reml_coef=mu_rm, reml_p=p_rm, i2_dl=i2_dl, i2_rm=i2_rm))
rep = pd.DataFrame(rows)
rep["dl_q"] = bh(rep["dl_p"]); rep["reml_q"] = bh(rep["reml_p"])
m = gm.merge(rep, on="genus", how="left")
m["absdiff_dl"] = (m["pooled_coef"]-m["dl_coef"]).abs()
m["absdiff_reml"] = (m["pooled_coef"]-m["reml_coef"]).abs()
n_sig_dl = int((rep["dl_q"] < 0.05).sum()); n_sig_reml = int((rep["reml_q"] < 0.05).sum())
n_sig_file = int((gm["p_adj"] < 0.05).sum())
print(f"文件中FDR<0.05属数: {n_sig_file}; 复算DL: {n_sig_dl}; 复算REML: {n_sig_reml}")
print(f"|coef差| 均值: DL={m['absdiff_dl'].mean():.2e}, REML={m['absdiff_reml'].mean():.2e}; 最大: DL={m['absdiff_dl'].max():.2e}, REML={m['absdiff_reml'].max():.2e}")
best = "REML" if m['absdiff_reml'].mean() < m['absdiff_dl'].mean() else "DL"
print(f"==> 存档数据更接近: {best}")
report["genus_meta"] = dict(n_genera=len(rep), n_sig_file=n_sig_file, n_sig_dl=n_sig_dl, n_sig_reml=n_sig_reml,
                            max_absdiff_dl=float(m['absdiff_dl'].max()), max_absdiff_reml=float(m['absdiff_reml'].max()),
                            best_match=best,
                            mean_absdiff_dl=float(m['absdiff_dl'].mean()), mean_absdiff_reml=float(m['absdiff_reml'].mean()))
# 显著属名单对比
sig_file = set(gm.loc[gm["p_adj"]<0.05, "genus"]); sig_dl = set(rep.loc[rep["dl_q"]<0.05, "genus"]); sig_rm = set(rep.loc[rep["reml_q"]<0.05, "genus"])
print(f"文件 vs DL 复算 显著属差异: 文件有而DL无 {sorted(sig_file-sig_dl)}; DL有而文件无 {sorted(sig_dl-sig_file)}")
print(f"文件 vs REML 复算 显著属差异: 文件有而REML无 {sorted(sig_file-sig_rm)}; REML有而文件无 {sorted(sig_rm-sig_file)}")
report["genus_meta"]["only_in_file_vs_dl"] = sorted(sig_file-sig_dl)
report["genus_meta"]["only_in_dl_vs_file"] = sorted(sig_dl-sig_file)
report["genus_meta"]["only_in_file_vs_reml"] = sorted(sig_file-sig_rm)
report["genus_meta"]["only_in_reml_vs_file"] = sorted(sig_rm-sig_file)
rep.to_csv(os.path.join(OUT, "repro_genus_meta.csv"), index=False)

# 关键属逐一核对（稿件点名的属）
key_genera = {
 "Lactobacillus": (0.529, 1.46e-7, 13.2), "Incertae.Sedis": (0.434, 1.80e-7, 0.0),
 "Collinsella": (0.508, 8.44e-6, 6.5), "Ruminococcus.gnavus.group": (0.566, 1.56e-4, 44.8),
 "Enterococcus": (0.403, 3.98e-4, 46.4), "Akkermansia": (0.549, 0.029, 71.6),
 "X.Eubacterium..eligens.group": (-0.355, 6.48e-7, 0.0), "Romboutsia": (-0.617, 3.75e-6, 21.3),
 "Dialister": (-0.796, 4.37e-5, 50.6), "Lachnospiraceae.NK4A136.group": (-0.329, 1.17e-4, 0.0),
 "Butyricicoccus": (-0.457, 0.007, 54.7), "Agathobacter": (-0.617, 0.005, 68.6),
 "Lachnospira": (-0.344, 0.008, 72.9), "Clostridium.sensu.stricto.1": (-0.404, 0.011, 38.9),
 "Coprococcus": (-0.357, 0.045, 72.2), "Bifidobacterium": (0.577, 0.017, None),
 "Faecalibacterium": (-0.557, 0.189, 72.3), "Roseburia": (-0.306, 0.073, None)}
print("\n稿件点名列 vs 存档文件 vs 复算(REML/DL):")
kg_rows = []
for g, (c_ms, p_ms, i2_ms) in key_genera.items():
    fr = gm[gm["genus"] == g]
    rr = rep[rep["genus"] == g]
    if len(fr) == 0 or len(rr) == 0:
        print(f"  {g:35s} !! 未找到"); continue
    fr = fr.iloc[0]; rr = rr.iloc[0]
    ok_c = abs(fr.pooled_coef - c_ms) < 0.006; ok_p = abs(fr.p_value - p_ms)/max(p_ms,1e-12) < 0.06
    ok_reml_c = abs(rr.reml_coef - c_ms) < 0.006
    kg_rows.append(dict(genus=g, ms=(c_ms,p_ms), file=(round(fr.pooled_coef,4), f"{fr.p_value:.3g}"),
                        reml=(round(rr.reml_coef,4), f"{rr.reml_p:.3g}"), dl=(round(rr.dl_coef,4), f"{rr.dl_p:.3g}"),
                        match_file=bool(ok_c and ok_p)))
    print(f"  {g:35s} 稿件={c_ms:+.3f}/p={p_ms:.3g}  文件={fr.pooled_coef:+.4f}/p={fr.p_value:.3g}  REML={rr.reml_coef:+.4f}/p={rr.reml_p:.3g}  DL={rr.dl_coef:+.4f}/p={rr.dl_p:.3g}  {'OK' if ok_c and ok_p else 'MISMATCH'}")
report["genus_meta"]["key_genera_all_match_file"] = all(r["match_file"] for r in kg_rows)

# ============ 2. α多样性荟萃 ============
print("\n" + "="*70); print("2. α多样性荟萃")
ab = pd.read_csv(os.path.join(SRC, "alpha_diversity_by_cohort.csv"))
am = pd.read_csv(os.path.join(SRC, "alpha_diversity_meta.csv"))
for idx in ["richness", "shannon", "simpson", "chao1"]:
    sub = ab[ab["index"] == idx]
    mu_dl, se_dl, p_dl, i2_dl, _, _ = dl_meta(sub["hedges_g"], sub["se_g"])
    mu_rm, se_rm, p_rm, i2_rm = reml_meta(sub["hedges_g"], sub["se_g"])
    f = am[am["index"] == idx].iloc[0]
    print(f"  {idx:9s} 文件: d={f.pooled_es:+.4f} p={f.p_value:.4f} I2={f.i2:.1f} | 复算DL: {mu_dl:+.4f} p={p_dl:.4f} I2={i2_dl:.1f} | 复算REML: {mu_rm:+.4f} p={p_rm:.4f} I2={i2_rm:.1f}")
    report.setdefault("alpha", {})[idx] = dict(file=(float(f.pooled_es), float(f.p_value)), dl=(mu_dl, p_dl), reml=(mu_rm, p_rm))

# ============ 3. 通路荟萃 ============
print("\n" + "="*70); print("3. 功能通路荟萃 (Cohen's d, SE= sqrt(1/n1+1/n2+d^2/(2(n1+n2))))")
pw = pd.read_csv(os.path.join(SRC, "pathway_differential_by_cohort.csv"))
pm = pd.read_csv(os.path.join(SRC, "pathway_meta_analysis.csv"))
for pathway in pm["pathway"]:
    sub = pw[pw["pathway"] == pathway].dropna(subset=["cohens_d"])
    sub = sub[sub["p_value"].notna()]
    if len(sub) < 2: continue
    d = sub["cohens_d"].values
    se = np.sqrt(1/sub["n_ad"].values + 1/sub["n_hc"].values + d**2/(2*(sub["n_ad"].values+sub["n_hc"].values)))
    mu, se_mu, p, i2, Q, Qp = dl_meta(d, se)
    f = pm[pm["pathway"] == pathway].iloc[0]
    print(f"  {pathway:22s} 文件: d={f.pooled_d:+.4f} p={f.p_value:.4f} I2={f.I2:.1f} | 复算DL: {mu:+.4f} p={p:.4f} I2={i2:.1f} | Δd={abs(mu-f.pooled_d):.1e}")
    report.setdefault("pathway", {})[pathway] = dict(file=(float(f.pooled_d), float(f.p_value)), repro_dl=(float(mu), float(p)))

# ============ 4. 氢代谢guild荟萃 ============
print("\n" + "="*70); print("4. 氢代谢guild荟萃")
hb = pd.read_csv(os.path.join(RES, "hydrogen_metabolism_by_cohort.csv"))
hm = pd.read_csv(os.path.join(RES, "hydrogen_metabolism_meta.csv"))
for pathway in hm["pathway"]:
    sub = hb[hb["pathway"] == pathway].dropna(subset=["cohens_d"])
    sub = sub[sub["p_value"].notna()]
    if len(sub) < 2: continue
    d = sub["cohens_d"].values
    se = np.sqrt(1/sub["n_ad"].values + 1/sub["n_hc"].values + d**2/(2*(sub["n_ad"].values+sub["n_hc"].values)))
    mu, se_mu, p, i2, Q, Qp = dl_meta(d, se)
    f = hm[hm["pathway"] == pathway].iloc[0]
    print(f"  {pathway:28s} 文件: d={f.pooled_cohens_d:+.4f} p={f.p_value:.4g} I2={f.i2:.1f} | 复算DL: {mu:+.4f} p={p:.4g} I2={i2:.1f} | Δd={abs(mu-f.pooled_cohens_d):.1e}")
    report.setdefault("hydrogen", {})[pathway] = dict(file=(float(f.pooled_cohens_d), float(f.p_value)), repro_dl=(float(mu), float(p)))

with open(os.path.join(OUT, "verify_meta_report.json"), "w", encoding="utf-8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2, default=str)
print("\n报告已保存 verify_meta_report.json")
