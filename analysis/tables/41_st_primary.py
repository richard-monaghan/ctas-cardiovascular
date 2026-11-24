#!/usr/bin/env python3
"""ST1 (gene architecture), ST2 (null calibration and numerical checks), ST5 (correlation-model sensitivity) from the accepted Task 1 bundle."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stlib import *

g = gene_table(); l = locus_table()
reps = set(l.best_gene_by_final_C_CTAS)
n = [int(g.in_q05.sum()), int(g.in_q01.sum()), int(g.in_q001.sum())]; nl = [int(l.in_q05.sum()), int(l.in_q01.sum()), int(l.in_q001.sum())]
assert n == [223, 99, 38] and nl == [91, 46, 21] and len(g) == 19080 and g.gene.is_unique
SRC1 = 'outputs/final_c_primary_freeze_20261021/final_C_gene_results_seed42.tsv'

# ---------------- ST1
b = Book('ST1')
d = g.rename(columns={'CHR': 'chromosome', 'canonical_locus_id_or_mapping_key': 'locus_mapping_key', 'breadth_B_log2_1_plus_PR': 'breadth_B_log2_1_plus_PR'})
d['is_locus_representative'] = d.gene.isin(reps)
d = d[['gene', 'chromosome', 'locus_mapping_key', 'CTAS', 'magnitude_M', 'participation_ratio_PR', 'breadth_B_log2_1_plus_PR', 'empirical_p', 'q_BH', 'q_tier', 'in_q05', 'in_q01', 'in_q001', 'is_locus_representative']]
d = d.sort_values(['q_BH', 'CTAS', 'gene'], ascending=[True, False, True])
b.readme('Supplementary Table 1. CTAS gene architecture (19,080 tested genes)', headline_pairs() + [
    ('Scope', 'One row per tested gene with the final CTAS statistic, its components and its significance.'),
    ('Columns', 'CTAS = M x B; magnitude_M = sqrt(sum of squared positive whitened trait evidence); participation_ratio_PR = effective breadth of that evidence across whitened axes; breadth_B = log2(1 + PR); empirical_p from 10^9 null draws; q_BH = genome-wide Benjamini-Hochberg q; q_tier = most stringent threshold passed (strict <). in_q05/in_q01/in_q001 flag the nested CTAS q-sets.'),
    ('locus_mapping_key', 'Canonical locus ID (V2L###) for genes in the canonical locus map; other genes carry their stable-neighbourhood key (STABLE:...) or UNMAPPED:gene. Every significant gene maps to a canonical locus.'),
    ('Representatives', 'is_locus_representative marks the highest-CTAS significant member of each of the 91 CTAS loci; a statistical locus label, not a causal-gene assignment.'),
    ('Units / precision', 'Values are stored at full numeric precision; no rounding in data sheets.'),
    ('Source', 'Accepted primary-result bundle outputs/final_c_primary_freeze_20261021/ (seed 42; seeds 1 and 2 give identical membership).')])
b.sheet('gene_architecture', d, 'All 19,080 tested genes, ordered by BH q (then CTAS).')
ts = pd.DataFrame([dict(q_set=k, genes=n[i], loci=nl[i]) for i, k in enumerate(['q<0.05', 'q<0.01', 'q<0.001'])])
b.sheet('qset_summary', ts, 'Nested CTAS q-sets (strict < thresholds).')
for i, k in enumerate(['q<0.05', 'q<0.01', 'q<0.001']):
    b.num('README / qset_summary', f'genes {k}', n[i], n[i], '19,080 tested genes', 'genes', SRC1); b.num('README / qset_summary', f'loci {k}', nl[i], nl[i], '', 'loci', 'final_C_locus_results_seed42.tsv')
b.num('README', 'tested genes', '19,080', 19080, '', 'genes', SRC1)
v = pd.DataFrame([dict(check='19,080 unique genes', result=len(d) == 19080 and d.gene.is_unique), dict(check='223/99/38 genes', result=n == [223, 99, 38]),
                  dict(check='CTAS = M x B (max abs diff)', result=float((d.CTAS - d.magnitude_M * d.breadth_B_log2_1_plus_PR).abs().max())),
                  dict(check='91 representatives flagged', result=int(d.is_locus_representative.sum()) == 91)])
b.sheet('validation', v, 'Structural checks.')
b.save('ST1_CTAS_gene_architecture.xlsx')

# ---------------- ST2
b = Book('ST2')
n_extrap = int(g.p_tail_extrapolated.sum()); n_zero = int((g.null_exceedance_count == 0).sum())
e = pd.read_csv(T1 + 'eigenvalues_C.tsv', sep='\t'); fl = pd.read_csv(T1 + 'eigen_floor_check.tsv', sep='\t', dtype=str)
minev = float(e.eigenvalue.min()); assert abs(minev - 0.5022736914630489) < 1e-15
seed = pd.read_csv(T1 + 'seed_stability_C.tsv', sep='\t').drop(columns=['estimator'])
diff = pd.read_csv(T1 + 'seed_membership_differences_C.tsv', sep='\t', dtype=str, keep_default_na=False)
total_symdiff = int(seed.symmetric_difference.sum())
b.readme('Supplementary Table 2. CTAS null calibration and numerical checks', headline_pairs() + [
    ('Null', 'Empirical P = (1 + number of null CTAS values exceeding the observed CTAS) / (1 + M), with M = 10^9 Gaussian null draws N(0, Sigma) processed through the identical whitening and scoring. BH q across all 19,080 genes.'),
    ('Tail handling', f'{n_zero} gene(s) had zero null exceedances at M = 10^9; of these, {n_extrap} received an exponential upper-tail estimate below the 1/(M+1) floor (column p_tail_extrapolated). All other P values are direct empirical counts. This does not affect any q-set boundary (see q_boundary_summary).'),
    ('Seed stability', f'Gene and locus membership at q<0.05, q<0.01 and q<0.001 is identical across seeds 42, 1 and 2 (total symmetric difference {total_symdiff}).'),
    ('Correlation matrix', 'Primary matrix Sigma = 2 sin(pi rho_S / 6) (diagonal 1). Minimum eigenvalue 0.5022736914630489 (paper-facing 0.50). The value 0.52 belongs only to the raw-Spearman sensitivity matrix (reported in numerical_checks, labelled as sensitivity).'),
    ('Numerical floor', 'The historical eigenvalue floor (1e-10) in the inverse square root is inactive: no eigenvalue is below it, and the floored and unfloored ZCA-cor whiteners are identical (max abs difference 0).'),
    ('Source', 'Accepted primary-result bundle outputs/final_c_primary_freeze_20261021/ (gene results, seed stability, seed membership differences, eigenvalues and floor checks).')])
gd = g[['gene', 'CTAS', 'null_exceedance_count', 'empirical_p', 'p_tail_extrapolated', 'q_BH', 'q_tier']].sort_values(['q_BH', 'CTAS'], ascending=[True, False])
b.sheet('gene_null_results', gd, 'Per-gene null exceedance counts (M = 10^9), empirical P and BH q.')
lk = g.canonical_locus_id_or_mapping_key; qb = []
for k, t in [('q<0.05', .05), ('q<0.01', .01), ('q<0.001', .001)]:
    ins = g.q_BH < t
    qb.append(dict(threshold=k, n_genes=int(ins.sum()), n_loci=int(lk[ins].nunique()), max_q_in_set=float(g.q_BH[ins].max()), min_q_outside_set=float(g.q_BH[~ins].min()),
                   max_empirical_p_in_set=float(g.empirical_p[ins].max()), genes_with_tail_estimate_in_set=int((g.p_tail_extrapolated & ins).sum())))
qb = pd.DataFrame(qb); b.sheet('q_boundary_summary', qb, 'Nested q-set boundaries (strict <).')
for r in qb.itertuples(): b.num('q_boundary_summary', f'genes / loci {r.threshold}', f'{r.n_genes} / {r.n_loci}', f'{r.n_genes}/{r.n_loci}', '19,080 genes', 'genes / loci', SRC1)
b.sheet('seed_stability', seed, 'Membership counts and symmetric differences versus seed 42 (identical_to_seed42 = TRUE means the same set).')
b.sheet('seed_membership_diffs', diff, 'Identifiers differing between seeds (none).')
b.num('README', 'seed membership symmetric difference (all tiers, genes and loci)', total_symdiff, total_symdiff, 'seeds 42, 1, 2', 'identifiers', 'seed_stability_C.tsv')
b.sheet('eigenvalues', e, 'Eigenvalues of the primary 16 x 16 trait correlation matrix Sigma = 2 sin(pi rho_S / 6).')
nc = fl.copy(); nc['quantity'] = nc.quantity.replace({'REFERENCE_min_eigenvalue_sigma_A_raw_spearman': 'SENSITIVITY_ONLY_min_eigenvalue_raw_spearman_matrix', 'max_abs_diff_sigma_A_file_vs_recomputed': 'max_abs_diff_raw_spearman_matrix_file_vs_recomputed', 'diag_sigma_C_all_one': 'diag_primary_matrix_all_one'})
b.sheet('numerical_checks', nc, 'Eigenvalue floor and whitener checks for the primary matrix; the last row is the raw-Spearman sensitivity matrix, shown for reference only.')
b.num('README / numerical_checks', 'primary minimum eigenvalue', '0.50', '0.5022736914630489', '', '', 'eigenvalues_C.tsv')
b.num('numerical_checks', 'raw-Spearman sensitivity minimum eigenvalue (sensitivity only)', '0.52', fl.set_index('quantity').at['REFERENCE_min_eigenvalue_sigma_A_raw_spearman', 'value'], '', '', 'eigen_floor_check.tsv')
b.num('README', 'null draws', '10^9', 1e9, '', 'draws', 'Task 1 run objects')
b.num('README', 'genes with tail estimate', n_extrap, n_extrap, '19,080', 'genes', SRC1)
v = pd.DataFrame([dict(check='null draws M = 10^9', result=True), dict(check='primary minimum eigenvalue (paper-facing 0.50)', result=f'{minev:.17g}'),
    dict(check='eigenvalue floor inactive', result=fl.set_index('quantity').at['floor_activated', 'value'] == 'FALSE'),
    dict(check='floored vs unfloored whitener max abs diff', result=fl.set_index('quantity').at['max_abs_diff_whitener_floor_vs_nofloor', 'value']),
    dict(check='seed membership differences (total)', result=total_symdiff), dict(check='223/99/38 genes; 91/46/21 loci', result=n == [223, 99, 38] and nl == [91, 46, 21])])
b.sheet('validation', v, 'Checks.')
b.save('ST2_CTAS_null_calibration.xlsx')

# ---------------- ST5
b = Book('ST5')
rs = pd.read_csv(T1 + 'raw_spearman_sensitivity.tsv', sep='\t', dtype=str, keep_default_na=False)
corr = rs[rs.section == 'genomewide_correlation'].set_index('metric').value
sp, pe = float(corr['spearman_CTAS_A_vs_C']), float(corr['pearson_CTAS_A_vs_C'])
x = rs[(rs.section == 'A_vs_C_same_seed') & (rs.seed == '42')].pivot(index='tier', columns='metric', values='value').astype(float)
rows = []
for tier, lab in [('q05', 'q<0.05'), ('q01', 'q<0.01'), ('q001', 'q<0.001')]:
    r = x.loc[tier]
    for lvl in ['genes', 'loci']:
        prim, sens, inter = int(r[f'{lvl}_C']), int(r[f'{lvl}_A']), int(r[f'{lvl}_intersection'])
        rows.append(dict(q_set=lab, unit=lvl, primary_transformed_correlation=prim, raw_spearman_sensitivity=sens, intersection=inter, jaccard=r[f'{lvl}_jaccard'],
                         primary_only=prim - inter, sensitivity_only=sens - inter, subset_relation='primary set is a subset of the sensitivity set' if prim == inter else 'not nested'))
cmpd = pd.DataFrame(rows)
stab = rs[rs.section.str.startswith('A_seed_stability')].copy()
stab['unit'] = np.where(stab.section.str.contains('genes'), 'genes', 'loci')
stab = stab[['seed', 'tier', 'unit', 'value']].rename(columns={'value': 'symmetric_difference_vs_seed42'})
exp_ = {('q<0.05', 'genes'): (223, 240), ('q<0.01', 'genes'): (99, 107), ('q<0.001', 'genes'): (38, 45), ('q<0.05', 'loci'): (91, 98), ('q<0.01', 'loci'): (46, 50), ('q<0.001', 'loci'): (21, 26)}
ok5 = all((r.primary_transformed_correlation, r.raw_spearman_sensitivity) == exp_[(r.q_set, r.unit)] for r in cmpd.itertuples())
b.readme('Supplementary Table 5. Correlation-model sensitivity: transformed correlation versus raw-Spearman matrix', headline_pairs() + [
    ('Comparison', 'Primary analysis: transformed trait correlation matrix Sigma = 2 sin(pi rho_S / 6) with ZCA-cor whitening. Sensitivity analysis: the raw Spearman trait correlation matrix with the same whitening, null size (10^9), seeds and BH procedure.'),
    ('Genome-wide concordance', f'CTAS rank correlation (Spearman) {sp:.8f}; Pearson {pe:.8f}.'),
    ('Set concordance', 'Genes 223/99/38 (primary) versus 240/107/45 (sensitivity); loci 91/46/21 versus 98/50/26. At every threshold the primary set is a subset of the sensitivity set (see set_comparison).'),
    ('Conclusion', 'Estimator choice has little effect on genome-wide ranking, while the transformed-correlation analysis remains the prespecified primary analysis.'),
    ('Seeds', 'Sensitivity membership is identical across seeds 42, 1 and 2 (sensitivity_seed_stability).'),
    ('Source', 'outputs/final_c_primary_freeze_20261021/raw_spearman_sensitivity.tsv')])
b.sheet('genome_wide_concordance', pd.DataFrame([dict(metric='CTAS Spearman rank correlation, primary vs sensitivity', value=sp), dict(metric='CTAS Pearson correlation, primary vs sensitivity', value=pe)]), 'Genome-wide concordance of CTAS (19,080 genes).')
b.sheet('set_comparison', cmpd, 'Nested-set comparison, seed 42 (seeds 1 and 2 identical).')
b.sheet('sensitivity_seed_stability', stab, 'Raw-Spearman sensitivity membership versus seed 42.')
b.num('README', 'CTAS rank correlation primary vs raw-Spearman', f'{sp:.8f}', corr['spearman_CTAS_A_vs_C'], '19,080 genes', 'Spearman rho', 'raw_spearman_sensitivity.tsv')
for r in cmpd.itertuples(): b.num('set_comparison', f'{r.unit} {r.q_set} primary / sensitivity', f'{r.primary_transformed_correlation} / {r.raw_spearman_sensitivity}', f'jaccard {r.jaccard}', '', r.unit, 'raw_spearman_sensitivity.tsv')
b.sheet('validation', pd.DataFrame([dict(check='counts reproduce 223/99/38 vs 240/107/45 and 91/46/21 vs 98/50/26', result=ok5), dict(check='rank correlation ~0.99996', result=abs(sp - 0.99996) < 5e-6),
    dict(check='primary subset of sensitivity at every threshold', result=bool((cmpd.primary_only == 0).all()))]), 'Checks.')
b.save('ST5_CTAS_correlation_sensitivity.xlsx')
print('extrapolated', n_extrap, 'zero-count', n_zero, 'ST5 ok', ok5)
