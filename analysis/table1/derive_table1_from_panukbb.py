#!/usr/bin/env python3
"""Derives the Table 1 sample sizes and EUR SNP heritabilities from the public Pan-UK Biobank phenotype manifest
and compares them field by field with the accepted Table 1 (results/main_tables/table1_accepted_values.tsv, text of the
released Table 1). Written for the release (P5); nothing is edited. Manifest (not redistributed):
  https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz
  (16 Nov 2026: Last-Modified 26 Jun 2024; SHA-256 of the .bgz 58343c50b06be5cee7a56a1d1a38c30065aa4321a8d8e6b3471d971e0a5b812c)
Usage: python3 analysis/table1/derive_table1_from_panukbb.py <phenotype_manifest.tsv(.bgz)>"""
import sys, os, re, pandas as pd
R = os.path.abspath(os.environ.get('CTAS_REPO', '.')); m = pd.read_csv(sys.argv[1], sep='\t', dtype=str, compression='gzip' if sys.argv[1].endswith(('.bgz', '.gz')) else None)
T = pd.read_csv(os.path.join(R, 'results/main_tables/table1_accepted_values.tsv'), sep='\t', dtype=str)
num = lambda x: float(x.replace(',', '')); sp = lambda x: tuple(map(float, re.match(r'(-?[\d.]+) \((-?[\d.]+)\)', x).groups()))
COL = {'cases_all_ancestries': 'n_cases_full_cohort_both_sexes', 'EUR_cases': 'n_cases_EUR', 'EUR_controls': 'n_controls_EUR', 'h2_observed_EUR': 'sldsc_25bin_h2_observed_EUR',
       'h2_observed_se_EUR': 'sldsc_25bin_h2_observed_se_EUR', 'h2_liability_EUR': 'sldsc_25bin_h2_liability_EUR', 'h2_liability_se_EUR': 'sldsc_25bin_h2_liability_se_EUR'}
out = []
for _, t in T.iterrows():
    tt, pc = ('icd10' if t.pan_ukbb_phenotype.startswith('ICD-10') else 'phecode'), t.pan_ukbb_phenotype.split()[-1]
    r = m[(m.trait_type == tt) & (m.phenocode == pc) & (m.pheno_sex == 'both_sexes')]; assert len(r) == 1, (t.trait, len(r)); r = r.iloc[0]
    o, ose = sp(t.h2_observed_EUR_se); l, lse = sp(t.h2_liability_EUR_se)
    for f, tv, dp in [('cases_all_ancestries', num(t.cases_all_ancestries), 0), ('EUR_cases', num(t.EUR_cases), 0), ('EUR_controls', num(t.EUR_controls), 0),
                      ('h2_observed_EUR', o, 4), ('h2_observed_se_EUR', ose, 4), ('h2_liability_EUR', l, 4), ('h2_liability_se_EUR', lse, 4)]:
        mv = float(r[COL[f]]); mr = round(mv, dp) if dp else mv
        out.append(dict(trait=t.trait, pan_ukbb_code=f'{tt}-{pc}-both_sexes', field=f, table1_value=tv, manifest_column=COL[f], manifest_value_full=mv, manifest_rounded=mr, match=abs(tv - mr) < 1e-12))
O = pd.DataFrame(out); os.makedirs(os.path.join(R, 'work'), exist_ok=True); O.to_csv(os.path.join(R, 'work/TABLE1_VS_PANUKBB_MANIFEST.tsv'), sep='\t', index=False)
print(f'{len(T)} traits x 7 fields = {len(O)} comparisons; matches {int(O.match.sum())}; mismatches {int((~O.match).sum())}'); sys.exit(0 if O.match.all() else 1)
