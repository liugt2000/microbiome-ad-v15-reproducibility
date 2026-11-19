# -*- coding: utf-8 -*-
"""Discovery-cohort covariate sensitivity with explicit legacy-bug reproduction.
The unadjusted rank-based mediation is preserved; adjusted models are new sensitivity analyses.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / 'results' / 'sensitivity' / 'covariates'
OUT.mkdir(parents=True, exist_ok=True)
RAW = REPO / 'data' / 'discovery'
RES = REPO / 'results' / 'discovery'
G = pd.read_csv(RES / 'genus_117_filtered.csv', index_col=0)
samples = list(G.columns)
med = pd.read_csv(RAW / 'mediation_analysis_data.csv')
med['Group'] = med['Group'].str.strip()
med['stage'] = med['Group'].map({'CG': 0, 'SCD': 1, 'MCI': 2, 'AD': 3})
med = med.set_index('Sample Analysis Name').loc[samples]
membership = pd.read_csv(RES / 'modules_membership.csv')
m5 = membership.loc[membership['module'] == 5, 'genus'].tolist()
sub = G.loc[m5].to_numpy(float)
z = (sub - sub.mean(axis=1, keepdims=True)) / sub.std(axis=1, ddof=0, keepdims=True)
_, _, vt = np.linalg.svd(z, full_matrices=False)
eig5 = vt[0]
if np.corrcoef(eig5, z.mean(axis=0))[0, 1] < 0:
    eig5 = -eig5

age = med['Age'].to_numpy(float)
male = (med['Gender'] == 'Male').to_numpy(float)
mmse = med['MMSE'].to_numpy(float)
stage = med['stage'].to_numpy(float)
hamd = med['HAMD'].to_numpy(float)
hama = med['HAMA'].to_numpy(float)
blautia = G.loc['Blautia'].to_numpy(float)
metab = pd.read_csv(RES / 'metab_1873_filtered.csv', index_col=0)
taurine = metab.loc['metab_15244', samples].to_numpy(float)
edu_map = {'Elementary School': 6.0, 'Middle School': 9.0, 'High School': 12.0, 'College': 16.0}
edu_years = med['Education'].map(edu_map)
if edu_years.isna().any():
    raise ValueError('Unmapped education labels: ' + repr(med.loc[edu_years.isna(), 'Education'].unique().tolist()))
edu_years = edu_years.to_numpy(float)
edu_cat = pd.get_dummies(med['Education'], prefix='Edu', drop_first=True, dtype=float).to_numpy()
legacy_edu_zero = np.nan_to_num(pd.to_numeric(med['Education'], errors='coerce').to_numpy(float))

covariates = {
    'unadjusted': np.empty((len(med), 0)),
    'legacy_age_sex_education_zero': np.column_stack([age, male, legacy_edu_zero]),
    'age_sex': np.column_stack([age, male]),
    'age_sex_education_years': np.column_stack([age, male, edu_years]),
    'age_sex_education_categorical': np.column_stack([age, male, edu_cat]),
    'age_sex_education_years_HAMD_HAMA': np.column_stack([age, male, edu_years, hamd, hama]),
}


def residual(y, C):
    D = np.column_stack([np.ones(len(y)), C]) if C.shape[1] else np.ones((len(y), 1))
    return y - D @ np.linalg.lstsq(D, y, rcond=None)[0]


def partial_assoc(y, x, C, legacy_rank_residual=False):
    ey = residual(stats.rankdata(y), C)
    ex = residual(stats.rankdata(x), C)
    if legacy_rank_residual:
        ans = stats.spearmanr(ey, ex)
        return float(ans.statistic), float(ans.pvalue)
    ans = stats.pearsonr(ey, ex)
    return float(ans.statistic), float(ans.pvalue)


assoc_rows = []
variables = {
    'M5_eigengene': eig5,
    'Blautia': blautia,
    'taurine': taurine,
}
pairs = [
    ('M5_eigengene', 'stage', eig5, stage),
    ('M5_eigengene', 'MMSE', eig5, mmse),
    ('Blautia', 'MMSE', blautia, mmse),
    ('Blautia', 'taurine', blautia, taurine),
    ('taurine', 'MMSE', taurine, mmse),
]
for name, C in covariates.items():
    for xname, yname, x, y in pairs:
        r0, p0 = partial_assoc(y, x, C, legacy_rank_residual=False)
        rl, pl = partial_assoc(y, x, C, legacy_rank_residual=True)
        assoc_rows.append({
            'model': name, 'x': xname, 'y': yname, 'n': len(x),
            'partial_r_pearson_residuals': r0, 'p_pearson_residuals': p0,
            'partial_r_legacy_ranked_residuals': rl, 'p_legacy_ranked_residuals': pl,
        })
assoc = pd.DataFrame(assoc_rows)
assoc.to_csv(OUT / 'partial_association_sensitivity.csv', index=False)

# Disease-range-restricted sensitivity, excluding HC/CG.
mask = med['Group'].to_numpy() != 'CG'
restricted = []
for name in ['unadjusted', 'age_sex', 'age_sex_education_years', 'age_sex_education_categorical']:
    C = covariates[name][mask]
    for xname, yname, x, y in pairs:
        r, p = partial_assoc(y[mask], x[mask], C)
        restricted.append({'model': name, 'x': xname, 'y': yname, 'n': int(mask.sum()), 'partial_r': r, 'p': p})
pd.DataFrame(restricted).to_csv(OUT / 'partial_association_nonHC_sensitivity.csv', index=False)

# Same rank-based OLS mediation definition as v15, with covariates added to every path model.
X = stats.rankdata(blautia)
M = stats.rankdata(taurine)
Y = mmse.astype(float)


def design(parts, C):
    cols = [p - p.mean() for p in parts]
    if C.shape[1]:
        cols.extend([C[:, j] - C[:, j].mean() for j in range(C.shape[1])])
    cols.append(np.ones(len(Y)))
    return np.column_stack(cols)


def mediation_point(Xv, Mv, Yv, Cv):
    n = len(Yv)
    intercept = np.ones(n)
    def dm(vectors):
        cols = [v - v.mean() for v in vectors]
        if Cv.shape[1]:
            cols.extend([Cv[:, j] - Cv[:, j].mean() for j in range(Cv.shape[1])])
        cols.append(intercept)
        return np.column_stack(cols)
    a = np.linalg.lstsq(dm([Xv]), Mv, rcond=None)[0][0]
    total = np.linalg.lstsq(dm([Xv]), Yv, rcond=None)[0][0]
    coefs = np.linalg.lstsq(dm([Xv, Mv]), Yv, rcond=None)[0]
    c_prime, b = coefs[0], coefs[1]
    return {'a': a, 'b': b, 'c_total': total, 'c_prime': c_prime,
            'indirect': a * b, 'proportion': (a * b / total if total != 0 else np.nan)}

model_names = ['unadjusted', 'age_sex', 'age_sex_education_years',
               'age_sex_education_categorical', 'age_sex_education_years_HAMD_HAMA']
B = 5000
# Shared bootstrap indices make between-model comparisons less Monte-Carlo noisy.
boot_indices = np.random.default_rng(20261119).integers(0, len(Y), size=(B, len(Y)))
summary = []
boot_store = {}
for name in model_names:
    C = covariates[name]
    point = mediation_point(X, M, Y, C)
    boots = np.empty(B)
    for i, idx in enumerate(boot_indices):
        boots[i] = mediation_point(X[idx], M[idx], Y[idx], C[idx])['indirect']
    lo, hi = np.percentile(boots, [2.5, 97.5])
    pboot = min(1.0, 2 * min(np.mean(boots <= 0), np.mean(boots >= 0)))
    summary.append({'model': name, 'n': len(Y), 'bootstrap': B, **point,
                    'ci_low': lo, 'ci_high': hi, 'p_bootstrap_two_sided': pboot,
                    'ci_excludes_zero': bool(lo * hi > 0)})
    boot_store[name] = boots
    np.save(OUT / f'mediation_bootstrap_{name}.npy', boots)
summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / 'mediation_covariate_sensitivity.csv', index=False)

# Leave-one-diagnostic-group-out mediation under the education-years model.
loo = []
C_full = covariates['age_sex_education_years']
for group in ['CG', 'SCD', 'MCI', 'AD']:
    keep = med['Group'].to_numpy() != group
    pt = mediation_point(X[keep], M[keep], Y[keep], C_full[keep])
    boots = np.empty(B)
    n = int(keep.sum())
    local_indices = np.random.default_rng(20261119).integers(0, n, size=(B, n))
    for i, idx in enumerate(local_indices):
        boots[i] = mediation_point(X[keep][idx], M[keep][idx], Y[keep][idx], C_full[keep][idx])['indirect']
    lo, hi = np.percentile(boots, [2.5, 97.5])
    loo.append({'excluded_group': group, 'n': n, **pt, 'ci_low': lo, 'ci_high': hi,
                'ci_excludes_zero': bool(lo * hi > 0)})
pd.DataFrame(loo).to_csv(OUT / 'mediation_leave_one_group_out.csv', index=False)

availability = pd.DataFrame([
    {'dataset': 'discovery', 'variable': v, 'n_available': int(med[v].notna().sum()), 'n_total': len(med)}
    for v in ['Age', 'Gender', 'Education', 'MMSE', 'HAMD', 'HAMA', 'Course of Disease']
] + [
    {'dataset': 'discovery', 'variable': v, 'n_available': 0, 'n_total': len(med)}
    for v in ['BMI', 'APOE', 'antibiotics', 'proton_pump_inhibitors', 'statins']
] + [
    {'dataset': 'public_755', 'variable': v, 'n_available': 0, 'n_total': 755}
    for v in ['Age', 'Gender', 'Education', 'BMI', 'MMSE', 'HAMD', 'HAMA', 'APOE', 'medications']
])
availability.to_csv(OUT / 'covariate_availability.csv', index=False)

notes = {
    'education_mapping_years': edu_map,
    'legacy_issue': 'String education labels were coerced to NaN and then replaced by zero; the historical adjusted estimate therefore adjusted only age and sex.',
    'association_primary_sensitivity': 'Partial Spearman calculated as Pearson correlation of rank residuals. Legacy ranked-residual output is also retained.',
    'mediation': 'Blautia and taurine are rank transformed before OLS; MMSE remains on its original scale. Covariates enter a, b and total-effect models. Bootstrap uses 5000 shared resamples per model.',
    'seed': 20261119,
    'm5_n_genera': len(m5),
    'education_counts': med['Education'].value_counts().to_dict(),
}
(OUT / 'analysis_notes.json').write_text(json.dumps(notes, ensure_ascii=False, indent=2), encoding='utf-8')
print(summary_df.round(6).to_string(index=False))
print('\nM5-stage associations:')
print(assoc[(assoc.x == 'M5_eigengene') & (assoc.y == 'stage')].round(6).to_string(index=False))
print('\nDONE covariate sensitivity')
