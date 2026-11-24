# CTAS V2 public release copy of outputs/finngen_corrected_assets_20261104/code/133_st11_families.py. Changes: paths only.
#!/usr/bin/env python3
"""ST9 colocalisation (Task 3), ST10 prior-study convergence (signed four-study matrix subset to the final sets;
coordinate-Foley primary), ST11 FinnGen external support (Task 2)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stlib import *
g = gene_table(); l = locus_table(); li = l.set_index('canonical_locus_id')
rd = lambda p, **k: pd.read_csv(p, sep='\t', **k)

# ---------------- ST11
gs = rd(T2 + 'gene_nested_set_summary.tsv'); gd = rd(T2 + 'gene_support_distribution.tsv'); ge = rd(T2 + 'gene_enrichment_tests.tsv'); led = rd(T2 + 'ctas_gene_finngen_ledger.tsv', dtype={'endpoints_p01': str}, keep_default_na=False)
miss = rd(T2 + 'gene_missingness.tsv'); lr = rd(T2 + 'locus_representative_finngen.tsv', keep_default_na=False); ln = rd(T2 + 'locus_nested_set_summary.tsv', keep_default_na=False)
core = rd(T2 + 'strict_core_21_locus_finngen.tsv', keep_default_na=False); v2 = rd(T2 + 'validation_checks.tsv'); N = CT + 'inputs/finngen_r12/derived/'; tier = rd(N + 'enrichment_with_CI.tsv'); tier = tier[[c for c in tier.columns if 'old' not in c]]; ts = rd(N + 'CORRECTED_TIER_SUMMARY.tsv', nrows=4)
s = gs.set_index('set')
head = [(k, int(s.at[k, 'supported_ge3']), int(s.at[k, 'n_evaluable']), s.at[k, 'pct_ge3_exact'], int(s.at[k, 'pct_ge3_rounded'])) for k in ['q<0.05', 'q<0.01', 'q<0.001', 'background_q>=0.05']]
headline = pd.DataFrame([dict(set=k, supported=a, evaluable=n, exact_pct=p, rounded_pct=r, statement=f'{a:,}/{n:,} = {p:.4f}% -> {r}%') for k, a, n, p, r in head])
EXPC = {'SUPPORT_q<0.05': '155/223', 'SUPPORT_q<0.01': '80/99', 'SUPPORT_q<0.001': '34/38', 'SUPPORT_background_q>=0.05': '1861/18701', 'ROUNDED_PERCENTAGES_42_52_51_VS_10': '70/81/89 vs 10'}
v2['expected'] = [EXPC.get(c_, e_) for c_, e_ in zip(v2.check_id, v2.expected)]
v2['check_id'] = v2.check_id.replace({'ROUNDED_PERCENTAGES_42_52_51_VS_10': 'ROUNDED_PERCENTAGES_70_81_89_VS_10'})
ok11 = [(a, n) for _, a, n, _, _ in head] == [(155, 223), (80, 99), (34, 38), (1861, 18701)]
FD = CT + 'work/finngen/families/'
fts = rd(FD + 'FAMILY_TIER_SUMMARY.tsv'); fte = rd(FD + 'FAMILY_ENRICHMENT_TESTS.tsv'); fbd = rd(FD + 'FAMILY_BANDS.tsv'); fdi = rd(FD + 'FAMILY_DISTRIBUTION.tsv'); fpf = rd(FD + 'FAMILY_PER_FAMILY_RATES.tsv')
fls = rd(FD + 'FAMILY_LOCUS_SUMMARY.tsv'); flo = rd(FD + 'FAMILY_LOCUS_SUPPORT.tsv', keep_default_na=False); fgl = rd(FD + 'FAMILY_GENE_LEDGER.tsv', keep_default_na=False)
fsel = fts.set_index('set')
fhead = [(k, int(fsel.at[k, 'supported_ge3']), int(fsel.at[k, 'n_evaluable']), fsel.at[k, 'pct_ge3']) for k in ['q<0.05', 'q<0.01', 'q<0.001', 'background_q>=0.05']]
fheadline = pd.DataFrame([dict(set=k, supported_3_of_6_families=a, evaluable=n, exact_pct=p, rounded_pct=round(p), CI95_low=fsel.at[k, 'pct_ge3_CI95_low'], CI95_high=fsel.at[k, 'pct_ge3_CI95_high'], statement=f'{a:,}/{n:,} = {p:.2f}%') for k, a, n, p in fhead])
okf = [(a, n) for _, a, n, _ in fhead] == [(109, 223), (56, 99), (23, 38), (1154, 18701)]
led2 = led.merge(fgl[['gene', 'n_families_p01', 'n_noncoronary_families', 'families', 'supported_ge3_families']], on='gene', how='left')
lr2 = lr.merge(flo[['canonical_locus_id', 'n_families_p01', 'families', 'supported_ge3_families']], on='canonical_locus_id', how='left')
core2 = core.merge(flo[['canonical_locus_id', 'n_families_p01', 'families']], on='canonical_locus_id', how='left')
b = Book('ST11')
b.readme('Supplementary Table 11. FinnGen external support for the final CTAS genes and loci', headline_pairs() + [
    ('Primary measure', 'FinnGen R12 endpoints were grouped into six clinical families so that correlated endpoints (e.g. CAD, MI and angina) were not counted as independent support: coronary (Angina, CAD = I9_CORATHER, MI), arrhythmia (AF), myocardial (Cardiomyopathy, HF), pressure (Hypertension), cerebrovascular (Stroke, TIA) and venous (VTE). A family is supported if any member endpoint has MAGMA gene-level P < 0.01. Primary: a gene is supported if at least 3 of the 6 families are supported.'),
    ('Gene-level headline (families)', '; '.join(f'{k}: {a:,}/{n:,} = {p:.2f}%' for k, a, n, p in fhead) + '.'),
    ('Secondary measures', 'Support in at least 1 and at least 2 families, mean number of supported families, and support in at least 2 non-coronary families (family_* sheets).'),
    ('Sensitivity', 'The earlier endpoint-based counts (at least 1, 2 and 3 of 10 endpoints) are retained in the sheets prefixed sens_endpoint_ and are labelled as a sensitivity analysis: ' + '; '.join(f'{k}: {a:,}/{n:,} = {p:.1f}%' for k, a, n, p, r in head) + ' (>=3 of 10 endpoints).'),
    ('Trend', 'Support increases across the nested sets (q<0.05, q<0.01, q<0.001); the sets are nested, so this is descriptive. In disjoint q-value bands support is 42.7% (0.01<=q<0.05), 54.1% (0.001<=q<0.01) and 60.5% (q<0.001); see family_bands.'),
    ('Statistics', 'One-sided Fisher exact tests for each nested set vs background (q>=0.05), Benjamini-Hochberg across the nine set x threshold tests, Woolf 95% CIs for odds ratios, Wilson 95% CIs for proportions; same framework for the endpoint sensitivity.'),
    ('Missing genes', 'Genes absent from the FinnGen gene outputs are not evaluable and are never scored as zero (gene_missingness); all 223 CTAS genes are evaluable and 156 background genes are not.'),
    ('Coordinate build', 'FinnGen R12 summary statistics (GRCh38) were converted to GRCh37 with the UCSC hg38ToHg19 chain through a validated Python implementation (not the liftOver binary) before MAGMA annotation (NCBI37.3, window 0 kb/0 kb, 1000 Genomes EUR LD, N = 500,348); 99.99% agreement with the 1000 Genomes BIM among uniquely matched rsIDs; 0.9% of variants failed conversion.'),
    ('Locus-representative analysis', 'Companion analysis (one highest-CTAS representative per locus) on the family measure (family_locus_* sheets) and the endpoint sensitivity (sens_endpoint_locus_*). In the Table 3 evidence integration a representative counts as FinnGen-supported with gene-level P < 0.01 in at least one endpoint.'),
    ('Interpretation', 'Independent-cohort external support for UK Biobank-selected CTAS genes and loci. CTAS was not rerun in FinnGen, so this does not re-estimate the CTAS statistic.'),
    ('Source', 'target_pipeline/finngen_corrected_rerun_01/ (corrected FinnGen rerun, 4 November 2026). An earlier FinnGen analysis used mismatched coordinate builds and is superseded.')])
b.sheet('family_headline', fheadline, 'PRIMARY HEADLINE: support in >=3 of 6 clinical families, nested CTAS gene sets vs background.')
b.sheet('family_nested_set_summary', fts, 'Family-based support at >=1, >=2 and >=3 families, mean families and non-coronary >=2 families, with Wilson CIs.')
b.sheet('family_enrichment_tests', fte, 'One-sided Fisher tests vs fixed background (BH over nine tests) with Woolf 95% CIs for odds ratios.')
b.sheet('family_bands', fbd, 'Disjoint q-value bands (support in >=3 families) and two-sided band comparisons.')
b.sheet('family_distribution', fdi, 'Distribution of supported-family counts 0-6 by set.')
b.sheet('family_per_family_rates', fpf, 'Percentage of genes supported in each clinical family by set.')
b.sheet('family_locus_summary', fls, 'LOCUS-REPRESENTATIVE COMPANION (families): nested 91/46/21 locus sets.')
b.sheet('family_locus_detail', flo, 'Per final locus representative: supported families and endpoint count.')
b.sheet('family_strict_core_21', core2, 'Strict-core 21-locus FinnGen evidence with supported families.')
b.sheet('gene_missingness', miss, 'Genes lacking FinnGen coverage (not evaluable).')
b.sheet('gene_ledger', led2, 'All 19,080 genes: family and endpoint support fields and all ten endpoint P values (blank = not available).')
b.sheet('sens_endpoint_headline', headline, 'SENSITIVITY (endpoint-based): support in >=3 of 10 endpoints, nested CTAS gene sets vs background.')
b.sheet('sens_endpoint_nested_summary', gs, 'SENSITIVITY: gene-level support at >=1, >=2 and >=3 of 10 endpoints.')
b.sheet('sens_endpoint_enrichment', ge, 'SENSITIVITY: one-sided Fisher tests vs background (BH over nine tests), disjoint-band checks.')
b.sheet('sens_endpoint_CIs', tier, 'SENSITIVITY: endpoint-based enrichment with Woolf 95% CIs.')
b.sheet('sens_endpoint_distribution', gd, 'SENSITIVITY: distribution of endpoint counts 0-10 by set.')
b.sheet('sens_endpoint_locus_summary', ln, 'SENSITIVITY: locus-representative companion on endpoints.')
b.sheet('sens_endpoint_locus_detail', lr2, 'SENSITIVITY: per representative endpoint counts (family counts added).')
b.sheet('sens_endpoint_strict_core', core, 'SENSITIVITY: strict-core 21 endpoint evidence.')
il = pd.DataFrame([dict(statement_id='I1_primary', statement='FinnGen external support (MAGMA P<0.01 in at least 3 of 6 clinical families) was found for ' + ', '.join(f"{a:,}/{n:,} ({p:.1f}%) CTAS genes at {k}" for k, a, n, p in fhead[:3]) + f' versus {fhead[3][1]:,}/{fhead[3][2]:,} ({fhead[3][3]:.2f}%) background genes (q>=0.05).', source_file='family_nested_set_summary'),
    dict(statement_id='I2_trend', statement='Support increases across the nested sets; the sets are nested and not independent groups. Disjoint bands: 42.7%, 54.1%, 60.5%.', source_file='family_bands'),
    dict(statement_id='I3_sensitivity', statement='Sensitivity analysis counting individual endpoints (at least 3 of 10): support was ' + ', '.join(f"{p:.0f}%" for k, a, n, p, r in head[:3]) + f' versus {head[3][3]:.0f}% background (sens_endpoint sheets).', source_file='sens_endpoint_headline'),
    dict(statement_id='I4_independent_cohort_boundary', statement='FinnGen is independent-cohort external support for UK Biobank-selected CTAS genes and loci. CTAS is not rerun in FinnGen, so this does not re-estimate the CTAS statistic, validate the correlation transformation or the breadth term independently, or identify causal genes; representatives remain statistical labels.', source_file='')])
b.sheet('interpretation', il, 'Downstream-safe statements.')
b.sheet('validation', v2, 'Analysis validation gates (corrected rerun, endpoint measure).')
for k, a, n, p in fhead: b.num('README / family_headline', f'FinnGen family support {k}', f'{a:,}/{n:,} ({round(p)}%)', p, f'{n:,} evaluable genes', '% genes', FD + 'FAMILY_TIER_SUMMARY.tsv')
b.save('ST11_CTAS_FinnGen_support.xlsx')
print('ST11 ok', ok11, okf)
