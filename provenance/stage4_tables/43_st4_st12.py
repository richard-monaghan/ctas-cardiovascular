#!/usr/bin/env python3
"""ST4 behaviour simulations (signed values carried over unchanged; terminology/method text updated only) and
ST12 trait missingness (regenerated from the unchanged 19,080 x 16 matrix and reconciled to the signed workbook)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stlib import *
ST = PAPER + 'Supplementary Tables/'

# ---------------- ST4
src = ST + 'ST4_CTAS_behaviour_simulations.xlsx'
b = Book('ST4'); sheets = {}
for s in ['calibration_typeI', 'fixed_total_signal', 'breadth_by_strength', 'parameters_provenance', 'validation_null_reproduction']:
    sheets[s] = pd.read_excel(src, sheet_name=s, header=0)
pp = sheets['parameters_provenance'].copy()
pp['item'] = pp['item'].replace({'sigma_md5': 'sigma_md5_cached_reference_null', 'M_null': 'M_reference_null'})
pp = pd.concat([pp, pd.DataFrame([dict(item='note', value='Signal and null are specified in whitened space (z_white ~ N(0, I_16)), so the simulated CTAS behaviour does not depend on which trait correlation matrix is whitened. The cached reference null (raw-Spearman working matrix, 10^6 draws) is used only for the null-reproduction check.')])])
cal = sheets['calibration_typeI']
b.readme('Supplementary Table 4. CTAS behaviour simulations', [
    ('Scope', 'Controlled simulations of the implemented CTAS statistic. CTAS depends on the 16-trait evidence only through the whitened vector z_white = Sigma^(-1/2) z (ZCA-cor whitening). Under the null z ~ N(0, Sigma), z_white ~ N(0, I_16) for any valid trait correlation matrix, so evidence is specified directly in whitened-axis space.'),
    ('fixed_total_signal', 'Total non-centrality lambda = sum(mu^2) held constant while spread over k active axes (mu = sqrt(lambda / k)); shows that CTAS rewards signal distributed across more whitened axes.'),
    ('breadth_by_strength', 'Grid of k active axes x per-axis strength mu, distinguishing weak-broad, strong-narrow and strong-broad evidence.'),
    ('calibration_typeI', f'Nominal versus observed type-I error on independent null genes (N = {int(cal.n_test.iloc[0]):,}); all {len(cal)} levels within Monte Carlo error.'),
    ('Selection probabilities', 'Columns sel_prob_a05/a01/a001 use per-gene nominal null-survival cut-offs (alpha = 0.05 / 0.01 / 0.001; CTAS 9.96 / 11.92 / 14.25); these are single-test nominal alphas, not the genome-wide BH-FDR gene cut-offs.'),
    ('Values', 'Simulation values, seeds, grids, replicate counts and checks are carried over unchanged from the signed simulation workbook; nothing was rerun. Only titles and method wording were updated.'),
    ('Source', 'CTAS V2 paper/Supplementary Tables/ST4_CTAS_behaviour_simulations.xlsx (signed values); parameters_provenance lists seeds and checksums.')])
b.sheet('calibration_typeI', cal, 'Nominal vs observed type-I error under the null.')
b.sheet('fixed_total_signal', sheets['fixed_total_signal'], 'Fixed total signal spread over k whitened axes.')
b.sheet('breadth_by_strength', sheets['breadth_by_strength'], 'Breadth (k axes) by per-axis strength (mu).')
b.sheet('parameters_provenance', pp, 'Seeds, grids, replicate counts and checksums.')
b.sheet('validation_null_reproduction', sheets['validation_null_reproduction'], 'Simulated N(0, I_16) null quantiles versus the cached reference null.')
b.num('README', 'calibration null genes', f'{int(cal.n_test.iloc[0]):,}', int(cal.n_test.iloc[0]), '', 'genes', src)
b.num('README', 'calibration levels within MC error', f'{int(cal.within_MC.astype(str).str.upper().eq("TRUE").sum())}/{len(cal)}', '', '', 'levels', src)
b.num('README', 'nominal CTAS cut-offs a=0.05/0.01/0.001', '9.96 / 11.92 / 14.25', '9.9595/11.9246/14.2534', '', 'CTAS', src)
same = all(pd.read_excel(src, sheet_name=s, header=0).equals(sheets[s]) for s in ['calibration_typeI', 'fixed_total_signal', 'breadth_by_strength', 'validation_null_reproduction'])
b.sheet('validation', pd.DataFrame([dict(check='data sheets identical to the signed workbook', result=same), dict(check='calibration within MC at all levels', result=bool(cal.within_MC.astype(str).str.upper().eq('TRUE').all()))]), 'Checks.')
b.save('ST4_CTAS_behaviour_simulations.xlsx')

# ---------------- ST12
raw = pd.read_csv(MATRIX, sep='\t', dtype=str, keep_default_na=False).set_index('GENE')
assert raw.shape == (19080, 16)
def classify(s):
    if s in ('', 'NA'): return 'missing'
    try: v = float(s)
    except ValueError: return 'invalid'
    if not np.isfinite(v): return 'nonfinite'
    return 'valid' if 0 < v <= 1 else 'invalid'
C = raw.map(classify); valid = C == 'valid'
pt = pd.DataFrame([dict(trait=t, total_genes=len(C), valid_genes=int(valid[t].sum()), missing_count=int((C[t] == 'missing').sum()), missing_pct=round(100 * (C[t] == 'missing').mean(), 2),
                        nonfinite_count=int((C[t] == 'nonfinite').sum()), invalid_p_count=int((C[t] == 'invalid').sum())) for t in raw.columns])
V = valid.astype(int).values; pc = pd.DataFrame(V.T @ V, index=raw.columns, columns=raw.columns)
pcp = (100 * pc / len(raw)).round(2)
avail = valid.sum(axis=1); pg = pd.DataFrame(dict(gene=raw.index, n_traits_available=avail.values))
ad = avail.value_counts().sort_index().rename_axis('n_traits_available').reset_index(name='n_genes')
summ = pd.DataFrame([('Total genes', len(raw)), ('Total traits', 16), ('Total cells (genes x traits)', raw.size), ('Genes complete across all 16 traits', int((avail == 16).sum())),
    ('Genes missing at least one trait', int((avail < 16).sum())), ('Total missing cells', int((C == 'missing').values.sum())), ('Total non-finite p-values', int((C == 'nonfinite').values.sum())),
    ('Total invalid p-values (outside 0,1]', int((C == 'invalid').values.sum())), ('Minimum pairwise-complete n', int(pc.values.min())), ('Minimum pairwise-complete %', round(100 * pc.values.min() / len(raw), 2))], columns=['metric', 'value'])
old = ST + 'ST12_V2_trait_missingness.xlsx'
o_pt = pd.read_excel(old, sheet_name='per_trait'); o_sum = pd.read_excel(old, sheet_name='summary'); o_ad = pd.read_excel(old, sheet_name='availability_distribution')
o_pc = pd.read_excel(old, sheet_name='pairwise_complete_n', index_col=0)
rec = [dict(check='per_trait identical to signed workbook', result=bool(o_pt.reset_index(drop=True).astype(str).equals(pt.astype(str)))),
       dict(check='summary identical to signed workbook', result=bool((o_sum.value.astype(float).values == summ.value.astype(float).values).all())),
       dict(check='availability distribution identical', result=bool((o_ad.values == ad.values).all())),
       dict(check='pairwise-complete n identical', result=bool((o_pc.loc[pc.index, pc.columns].values == pc.values).all())),
       dict(check='matrix 19,080 x 16; cells reconcile', result=bool(pt.eval('valid_genes + missing_count + nonfinite_count + invalid_p_count').eq(19080).all())),
       dict(check='per-gene availability sums to 19,080', result=int(ad.n_genes.sum()) == 19080)]
b = Book('ST12')
s = summ.set_index('metric').value
b.readme('Supplementary Table 12. Trait-matrix missingness QC (19,080 genes x 16 traits)', [
    ('Scope', 'Quality control of the 19,080 x 16 gene-level MAGMA P-value matrix used by CTAS.'),
    ('Key result', f"{int(s['Genes complete across all 16 traits']):,} genes complete across all 16 traits; {int(s['Genes missing at least one trait'])} missing at least one; {int(s['Total missing cells']):,} missing cells; no invalid or non-finite P values. Minimum pairwise-complete n = {int(s['Minimum pairwise-complete n']):,} ({s['Minimum pairwise-complete %']}%)."),
    ('Pattern', 'Per-trait missingness is concentrated in AF, CAD, HF and hypertension (about 3.9% each) versus 0.01-0.18% for the remaining 12 traits.'),
    ('Handling in CTAS', 'Missing trait evidence is set to z = 0 before whitening (no evidence); missingness is never imputed as significance.'),
    ('Definitions', "Valid P: finite and in (0, 1]. Missing: empty or 'NA'. Non-finite: NaN or +/-Inf. Invalid: finite but outside (0, 1]. Pairwise-complete: both traits valid for the same gene."),
    ('Regeneration', 'All values regenerated directly from the matrix and reconciled exactly to the signed missingness workbook (validation sheet); only titles and provenance wording changed.'),
    ('Source', 'analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv')])
b.sheet('per_trait', pt, 'Valid / missing / non-finite / invalid counts per trait.')
b.sheet('pairwise_complete_n', pc.reset_index().rename(columns={'index': 'trait'}), 'Genes with valid P values for each trait pair.')
b.sheet('pairwise_complete_pct', pcp.reset_index().rename(columns={'index': 'trait'}), 'Pairwise-complete percentage of 19,080 genes.')
b.sheet('per_gene_availability', pg, 'Number of traits with valid P values per gene.')
b.sheet('availability_distribution', ad, 'Genes by number of available traits.')
b.sheet('summary', summ, 'QC headline.')
b.sheet('validation', pd.DataFrame(rec), 'Regeneration and reconciliation checks.')
for k in ['Genes complete across all 16 traits', 'Genes missing at least one trait', 'Total missing cells', 'Minimum pairwise-complete n']:
    b.num('README / summary', k, f'{int(s[k]):,}', int(s[k]), '19,080 genes', 'genes' if 'cells' not in k else 'cells', 'UKBB_16trait_gene_matrix.tsv')
b.num('README', 'matrix dimensions', '19,080 x 16', '19080x16', '', '', 'UKBB_16trait_gene_matrix.tsv')
b.save('ST12_CTAS_trait_missingness.xlsx')
print('ST4 identical', same); print(pd.DataFrame(rec).to_string())
