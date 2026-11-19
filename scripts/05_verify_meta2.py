# -*- coding: utf-8 -*-
"""修正版：BH修正方向 + 正确REML方程(Σw²(y-μ̂)²=c) + NaN行纳入 + 详细定位不一致属"""
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
import os, json

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "data", "source_data")
RES = os.path.join(REPO, "data", "meta_analysis")
OUT = os.path.join(REPO, "results", "verification")
os.makedirs(OUT, exist_ok=True)

def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p)
    ps = p[o]
    q = ps * n / np.arange(1, n+1)
    q = np.minimum.accumulate(q[::-1])[::-1]   # 反向累积最小值（正确的BH）
    out = np.empty(n); out[o] = np.minimum(q, 1.0)
    return out

def dl_meta(yi, sei):
    yi = np.asarray(yi, float); sei = np.asarray(sei, float)
    w = 1.0/sei**2
    fe = np.sum(w*yi)/np.sum(w)
    Q = np.sum(w*(yi-fe)**2); df = len(yi)-1
    c = np.sum(w) - np.sum(w**2)/np.sum(w)
    tau2 = max(0.0, (Q-df)/c) if c > 0 else 0.0
    wr = 1.0/(sei**2+tau2)
    mu = np.sum(wr*yi)/np.sum(wr)
    se_mu = np.sqrt(1.0/np.sum(wr))
    z = mu/se_mu; p = 2*(1-stats.norm.cdf(abs(z)))
    i2 = max(0.0, (Q-df)/Q*100) if Q > 0 else 0.0
    return mu, se_mu, p, i2

def reml_tau2(yi, vi):
    """正确REML: 解 Σw²(y-μ̂)² = c, w=1/(v+τ²)"""
    def f(t):
        w = 1.0/(vi+t)
        mu = np.sum(w*yi)/np.sum(w)
        lhs = np.sum(w**2*(yi-mu)**2)
        c = np.sum(w) - np.sum(w**2)/np.sum(w)
        return lhs - c
    if f(0.0) <= 0: return 0.0
    hi = 1.0
    while f(hi) > 0 and hi < 1e8: hi *= 4
    if f(hi) > 0: return hi
    return brentq(f, 0.0, hi, xtol=1e-14)

def reml_meta(yi, sei):
    yi = np.asarray(yi, float); sei = np.asarray(sei, float)
    t2 = reml_tau2(yi, sei**2)
    w = 1.0/(sei**2+t2)
    mu = np.sum(w*yi)/np.sum(w)
    se_mu = np.sqrt(1.0/np.sum(w))
    z = mu/se_mu; p = 2*(1-stats.norm.cdf(abs(z)))
    wf = 1.0/sei**2; fe = np.sum(wf*yi)/np.sum(wf)
    Q = np.sum(wf*(yi-fe)**2); df = len(yi)-1
    i2 = max(0.0, (Q-df)/Q*100) if Q > 0 else 0.0
    return mu, se_mu, p, i2

report = {}
# ===== 1. 属级荟萃 =====
ma = pd.read_csv(os.path.join(SRC, "maaslin2_all_results.csv"))
gm = pd.read_csv(os.path.join(SRC, "genus_meta_analysis.csv"))
rows = []
for g, sub in ma.groupby("feature"):
    if len(sub) < 2: continue
    mu_dl, se_dl, p_dl, i2_dl = dl_meta(sub["coef"], sub["stderr"])
    mu_rm, se_rm, p_rm, i2_rm = reml_meta(sub["coef"], sub["stderr"])
    rows.append(dict(genus=g, n=len(sub), dl_coef=mu_dl, dl_p=p_dl, dl_i2=i2_dl,
                     reml_coef=mu_rm, reml_p=p_rm, reml_i2=i2_rm))
rep = pd.DataFrame(rows)
rep["dl_q"] = bh(rep["dl_p"]); rep["reml_q"] = bh(rep["reml_p"])
m = gm.merge(rep, on="genus", how="left")
m["d_dl"] = (m["pooled_coef"]-m["dl_coef"]).abs()
m["d_rm"] = (m["pooled_coef"]-m["reml_coef"]).abs()
n_file = int((gm["p_adj"]<0.05).sum()); n_dl = int((rep["dl_q"]<0.05).sum()); n_rm = int((rep["reml_q"]<0.05).sum())
print(f"属级荟萃: 文件FDR显著={n_file}, 复算DL={n_dl}, 复算REML={n_rm}")
print(f"|coef-文件|: DL均值={m['d_dl'].mean():.2e}/最大={m['d_dl'].max():.2e}; REML均值={m['d_rm'].mean():.2e}/最大={m['d_rm'].max():.2e}")
print("DL差异最大的8个属:"); print(m.nlargest(8, "d_dl")[["genus","pooled_coef","dl_coef","reml_coef","p_value","dl_p","reml_p","i2","dl_i2","reml_i2"]].to_string(index=False))
print("REML差异最大的8个属:"); print(m.nlargest(8, "d_rm")[["genus","pooled_coef","dl_coef","reml_coef","p_value","dl_p","reml_p","i2","dl_i2","reml_i2"]].to_string(index=False))
sig_f = set(gm.loc[gm["p_adj"]<0.05,"genus"]); sig_d = set(rep.loc[rep["dl_q"]<0.05,"genus"]); sig_r = set(rep.loc[rep["reml_q"]<0.05,"genus"])
print(f"\n显著属集合: 文件\\DL={sorted(sig_f-sig_d)} DL\\文件={sorted(sig_d-sig_f)}")
print(f"            文件\\REML={sorted(sig_f-sig_r)} REML\\文件={sorted(sig_r-sig_f)}")
report["genus"] = dict(n_file=n_file, n_dl=n_dl, n_reml=n_rm, only_file_vs_dl=sorted(sig_f-sig_d), only_dl_vs_file=sorted(sig_d-sig_f),
                       only_file_vs_reml=sorted(sig_f-sig_r), only_reml_vs_file=sorted(sig_r-sig_f),
                       mean_dl=float(m['d_dl'].mean()), max_dl=float(m['d_dl'].max()),
                       mean_reml=float(m['d_rm'].mean()), max_reml=float(m['d_rm'].max()))
rep.to_csv(os.path.join(OUT, "repro_genus_meta_fixed.csv"), index=False)

# 方向拆分核对（稿件: 10富集/14耗竭）
gd = gm[gm["p_adj"]<0.05]
print(f"\n文件中24个显著属方向: enriched={int((gd['direction']=='AD_enriched').sum())}, depleted={int((gd['direction']=='AD_depleted').sum())}")

# ===== 2. α多样性(用正确REML) =====
print("\nα多样性(文件 vs DL vs REML):")
ab = pd.read_csv(os.path.join(SRC, "alpha_diversity_by_cohort.csv"))
am = pd.read_csv(os.path.join(SRC, "alpha_diversity_meta.csv"))
report["alpha"] = {}
for idx in ["richness","shannon","simpson","chao1"]:
    s = ab[ab["index"]==idx]
    mu_dl,_,p_dl,i2_dl = dl_meta(s["hedges_g"], s["se_g"])
    mu_rm,_,p_rm,i2_rm = reml_meta(s["hedges_g"], s["se_g"])
    f = am[am["index"]==idx].iloc[0]
    print(f"  {idx:9s} 文件 d={f.pooled_es:+.4f} p={f.p_value:.4f} I2={f.i2:.1f} | DL {mu_dl:+.4f} p={p_dl:.4f} I2={i2_dl:.1f} | REML {mu_rm:+.4f} p={p_rm:.4f} I2={i2_rm:.1f}")
    report["alpha"][idx] = dict(file=(float(f.pooled_es),float(f.p_value),float(f.i2)), dl=(mu_dl,p_dl,i2_dl), reml=(mu_rm,p_rm,i2_rm))

# ===== 3&4. 通路/氢代谢(纳入p=NaN但d=0的行) =====
def pathway_repro(by_cohort_file, meta_file, label):
    pw = pd.read_csv(by_cohort_file); pm = pd.read_csv(meta_file)
    print(f"\n{label} (文件 vs 复算DL, 含NaN-p行):")
    out = {}
    for pathway in pm["pathway"]:
        sub = pw[pw["pathway"]==pathway]
        sub = sub[sub["cohens_d"].notna()]
        d = sub["cohens_d"].values
        se = np.sqrt(1/sub["n_ad"].values + 1/sub["n_hc"].values + d**2/(2*(sub["n_ad"].values+sub["n_hc"].values)))
        mu, se_mu, p, i2 = dl_meta(d, se)
        f = pm[pm["pathway"]==pathway].iloc[0]
        ncf = getattr(f, "n_cohorts", None)
        print(f"  {pathway:28s} 文件 d={f.iloc[1] if False else ''}", end="")
        fd = float(f["pooled_d"] if "pooled_d" in f.index else f["pooled_cohens_d"])
        fp = float(f["p_value"]); fi2 = float(f["I2"] if "I2" in f.index else f["i2"])
        print(f" d={fd:+.4f} p={fp:.4g} I2={fi2:.1f} n={ncf} | 复算 d={mu:+.4f} p={p:.4g} I2={i2:.1f} n={len(sub)} | Δ={abs(mu-fd):.1e}")
        out[pathway] = dict(file=(fd,fp,fi2), repro=(float(mu),float(p),float(i2)), n_file=int(ncf), n_repro=int(len(sub)))
    return out

report["pathway"] = pathway_repro(os.path.join(SRC,"pathway_differential_by_cohort.csv"), os.path.join(SRC,"pathway_meta_analysis.csv"), "功能通路荟萃")
report["hydrogen"] = pathway_repro(os.path.join(RES,"hydrogen_metabolism_by_cohort.csv"), os.path.join(RES,"hydrogen_metabolism_meta.csv"), "氢代谢guild荟萃")

# ===== 5. Roseburia / Faecalibacterium 溯源计算 =====
print("\nRoseburia 精确复算:")
sub = ma[ma["feature"]=="Roseburia"]
print(sub[["feature","coef","stderr","cohort"]].to_string(index=False))
for name, fn in [("DL", dl_meta), ("REML", reml_meta)]:
    mu, se, p, i2 = fn(sub["coef"], sub["stderr"])
    print(f"  {name}: coef={mu:+.4f} p={p:.4f} I2={i2:.1f}")
print("  文件: coef=-0.3035 p=0.1123 I2=74.5 | 稿件: coef=-0.306 p=0.073")

with open(os.path.join(OUT, "verify_meta_report2.json"), "w", encoding="utf-8") as f:
    json.dump(report, f, ensure_ascii=False, indent=2, default=str)
print("\n已保存 verify_meta_report2.json")
