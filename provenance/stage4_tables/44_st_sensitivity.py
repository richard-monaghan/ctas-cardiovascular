#!/usr/bin/env python3
"""ST6 effective-N, ST7 leave-one-trait-out and ST8 obstructive-cluster sensitivity: deterministic rebase of the stored
complete gene-level sensitivity results onto the final 223/99/38 CTAS genes and 91/46/21 CTAS loci. Nothing is rerun."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stlib import *
from scipy.stats import spearmanr, pearsonr
g = gene_table(); l = locus_table(); mem = final_members(g); gi = g.set_index('gene')
TH = [('q<0.05', .05, 'in_q05'), ('q<0.01', .01, 'in_q01'), ('q<0.001', .001, 'in_q001')]
reps = l.set_index('canonical_locus_id').best_gene_by_final_C_CTAS
def locus_present(qser, lid, t): return any(qser.get(x, 1.0) < t for x in mem[lid])

# ---------------- ST6 effective N
EF = CT + 'panel_v2_effective_n_sensitivity/v2_gene_comparison_all.csv'
e = pd.read_csv(EF).set_index('gene'); assert len(e) == 19080 and set(e.index) == set(g.gene)
tr = pd.read_excel(PAPER + 'Supplementary Tables/ST6_V2_effective_N_sensitivity.xlsx', sheet_name='Trait_effective_N', header=1)
conc = pd.DataFrame([dict(metric='CTAS Pearson, flat-N vs effective-N', value=pearsonr(e.CTAS_flat, e.CTAS_eff)[0]), dict(metric='CTAS Spearman, flat-N vs effective-N', value=spearmanr(e.CTAS_flat, e.CTAS_eff)[0]),
    dict(metric='maximum |delta CTAS|', value=float((e.CTAS_eff - e.CTAS_flat).abs().max())), dict(metric='median |delta CTAS|', value=float((e.CTAS_eff - e.CTAS_flat).abs().median()))])
gs, ls, gd, ld = [], [], [], []
for lab, t, f in TH:
    F = g.gene[g[f]]; fl = e.q_flat.reindex(F) < t; ef = e.q_eff.reindex(F) < t
    gs.append(dict(q_set=lab, final_genes=len(F), significant_in_flatN_sensitivity=int(fl.sum()), significant_in_effectiveN_sensitivity=int(ef.sum()),
                   final_genes_changing_status_flat_to_eff=int((fl != ef).sum()), changed_genes=';'.join(F[(fl != ef).values])))
    L = l.canonical_locus_id[l[f]]
    pf = {x: locus_present(e.q_flat, x, t) for x in L}; pe = {x: locus_present(e.q_eff, x, t) for x in L}
    ch = [x for x in L if pf[x] != pe[x]]
    ls.append(dict(q_set=lab, final_loci=len(L), present_in_flatN_sensitivity=sum(pf.values()), present_in_effectiveN_sensitivity=sum(pe.values()), final_loci_changing_status=len(ch), changed_loci=';'.join(ch)))
gs, ls = pd.DataFrame(gs), pd.DataFrame(ls)
F = g[g.in_q05]
gd = pd.DataFrame(dict(gene=F.gene, canonical_locus_id=F.canonical_locus_id_or_mapping_key, final_q=F.q_BH, final_q_tier=F.q_tier, CTAS_flatN_sens=e.CTAS_flat.reindex(F.gene).values,
    CTAS_effN_sens=e.CTAS_eff.reindex(F.gene).values, q_flatN_sens=e.q_flat.reindex(F.gene).values, q_effN_sens=e.q_eff.reindex(F.gene).values))
for lab, t, f in TH: gd[f'same_status_{lab}'] = (gd.q_flatN_sens < t) == (gd.q_effN_sens < t)
rep_rows = []
for lid in l.canonical_locus_id:
    m = mem[lid]; rf = max(m, key=lambda x: e.at[x, 'CTAS_flat']); re_ = max(m, key=lambda x: e.at[x, 'CTAS_eff'])
    rep_rows.append(dict(canonical_locus_id=lid, final_representative=reps[lid], top_member_flatN_sens=rf, top_member_effN_sens=re_, top_member_changes=rf != re_,
        **{f'present_flat_{lab}': locus_present(e.q_flat, lid, t) for lab, t, f in TH if l.set_index('canonical_locus_id').at[lid, f]},
        **{f'present_eff_{lab}': locus_present(e.q_eff, lid, t) for lab, t, f in TH if l.set_index('canonical_locus_id').at[lid, f]}))
ld = pd.DataFrame(rep_rows)
b = Book('ST6')
b.readme('Supplementary Table 6. Effective-sample-size sensitivity (MAGMA rerun with trait-specific N_eff), rebased to the final CTAS sets', headline_pairs() + [
    ('Method note', SENS_NOTE),
    ('Design', 'The 16 MAGMA analyses were rerun with trait-specific effective sample sizes N_eff = 4 / (1/N_cases + 1/N_controls) (same summary statistics, LD panel, annotation and SNPwise-mean model). Flat-N and effective-N gene matrices were scored by the same sensitivity engine; the comparison is therefore within that engine.'),
    ('Rebase', 'For each final CTAS gene and locus, status is compared between the flat-N and effective-N runs of the sensitivity engine at the matching threshold. A final locus is present in a run when at least one of its final member genes passes the threshold.'),
    ('Result', f"Genome-wide CTAS Spearman {conc.value[1]:.8f}. Final genes changing status between flat-N and effective-N: {', '.join(f'{r.q_set} {r.final_genes_changing_status_flat_to_eff}/{r.final_genes}' for r in gs.itertuples())}. Final loci changing status: {', '.join(f'{r.q_set} {r.final_loci_changing_status}/{r.final_loci}' for r in ls.itertuples())}."),
    ('Interpretation', 'Under the SNPwise-mean model the N argument does not enter the gene P value; residual movement is confined to hyper-significant cells and reflects rerun variability in the extreme tail.'),
    ('Sources', 'panel_v2_effective_n_sensitivity/v2_gene_comparison_all.csv (19,080 genes, both runs); per-trait N_eff from the signed effective-N workbook; final membership from the accepted primary result.')])
b.sheet('trait_effective_N', tr, 'Per-trait effective sample size and MAGMA gene-P movement between runs (signed values).')
b.sheet('genome_wide_concordance', conc, 'Within-engine concordance, flat-N vs effective-N (19,080 genes).')
b.sheet('final_gene_rebase', gs, 'Final CTAS genes: status in the flat-N and effective-N sensitivity runs.')
b.sheet('final_locus_rebase', ls, 'Final CTAS loci: presence in the flat-N and effective-N sensitivity runs.')
b.sheet('final_gene_detail', gd, '223 final CTAS genes with sensitivity-run CTAS and q.')
b.sheet('final_locus_detail', ld, '91 final CTAS loci: presence and top member in each sensitivity run.')
for r in gs.itertuples(): b.num('README / final_gene_rebase', f'final genes changing status {r.q_set}', f'{r.final_genes_changing_status_flat_to_eff}/{r.final_genes}', r.final_genes_changing_status_flat_to_eff, f'{r.final_genes} final genes', 'genes', EF)
for r in ls.itertuples(): b.num('README / final_locus_rebase', f'final loci changing status {r.q_set}', f'{r.final_loci_changing_status}/{r.final_loci}', r.final_loci_changing_status, f'{r.final_loci} final loci', 'loci', EF)
b.num('README', 'CTAS Spearman flat-N vs effective-N', f'{conc.value[1]:.8f}', conc.value[1], '19,080 genes', 'rho', EF)
b.sheet('validation', pd.DataFrame([dict(check='stored results cover all 19,080 genes', result=True), dict(check='final sets 223/99/38 and 91/46/21', result=list(gs.final_genes) == [223, 99, 38] and list(ls.final_loci) == [91, 46, 21])]), 'Checks.')
b.save('ST6_CTAS_effective_N_sensitivity.xlsx')

# ---------------- ST7 leave-one-trait-out
W = OUT + 'work/'
Q = pd.read_csv(W + 'loto_q_by_omitted_trait.tsv', sep='\t').set_index('gene'); CT_ = pd.read_csv(W + 'loto_ctas_by_omitted_trait.tsv', sep='\t').set_index('gene')
ref = pd.read_csv(W + 'loto_full16_reference.tsv', sep='\t').set_index('gene'); traits = list(Q.columns); assert len(Q) == 19080 and len(traits) == 16
gstab = pd.read_excel(PAPER + 'Supplementary Tables/ST7_V2_leave_one_trait_out.xlsx', sheet_name='genomewide_stability', header=0)
li = l.set_index('canonical_locus_id')
lr = []
for lid in l.canonical_locus_id:
    row = dict(canonical_locus_id=lid, chromosome=int(li.at[lid, 'chromosome']), representative_gene=reps[lid], locus_q_tier=li.at[lid, 'most_stringent_locus_tier'], n_final_members=len(mem[lid]),
               retained_full16_sensitivity_reference=locus_present(ref.q_full16_sens, lid, .05))
    for tname in traits: row[tname] = int(locus_present(Q[tname], lid, .05))
    row['n_runs_retained_q05'] = sum(row[t] for t in traits); row['retention_proportion'] = row['n_runs_retained_q05'] / 16
    row['rep_n_runs_q05'] = int((Q.loc[reps[lid]] < .05).sum())
    if li.at[lid, 'in_q01']: row['n_runs_retained_q01'] = sum(locus_present(Q[t], lid, .01) for t in traits)
    if li.at[lid, 'in_q001']: row['n_runs_retained_q001'] = sum(locus_present(Q[t], lid, .001) for t in traits)
    lr.append(row)
lr = pd.DataFrame(lr)
rr = pd.DataFrame([dict(canonical_locus_id=lid, representative_gene=reps[lid], final_CTAS=gi.at[reps[lid], 'CTAS'], final_q=gi.at[reps[lid], 'q_BH'], **{f'q_omit_{t}': Q.at[reps[lid], t] for t in traits}) for lid in l.canonical_locus_id])
F05 = set(g.gene[g.in_q05]); summ = []
for tname in traits:
    run05 = set(Q.index[Q[tname] < .05])
    summ.append(dict(omitted_trait=tname, loci_retained_q05_of_91=int(lr[tname].sum()), pct_q05=100 * lr[tname].mean(),
        loci_retained_q01_of_46=sum(locus_present(Q[tname], x, .01) for x in l.canonical_locus_id[l.in_q01]), loci_retained_q001_of_21=sum(locus_present(Q[tname], x, .001) for x in l.canonical_locus_id[l.in_q001]),
        representatives_retained_q05=int((Q.loc[list(reps), tname] < .05).sum()), final_genes_retained_q05_of_223=len(run05 & F05), gene_jaccard_vs_final_q05=len(run05 & F05) / len(run05 | F05),
        CTAS_spearman_vs_full16_sensitivity=spearmanr(CT_[tname], ref.CTAS_full16_sens.reindex(CT_.index))[0], n_run_genes_outside_final_q05=len(run05 - F05)))
summ = pd.DataFrame(summ).sort_values('loci_retained_q05_of_91')
stable = lr[lr.n_runs_retained_q05 == 16][['canonical_locus_id', 'representative_gene', 'locus_q_tier']]
newd = [dict(omitted_trait=t, gene=x, run_q=Q.at[x, t], locus_mapping_key=gi.at[x, 'canonical_locus_id_or_mapping_key'], final_q=gi.at[x, 'q_BH']) for t in traits for x in Q.index[Q[t] < .05] if x not in F05]
newd = pd.DataFrame(newd)
b = Book('ST7')
b.readme('Supplementary Table 7. Leave-one-trait-out stability, rebased to the final CTAS loci', headline_pairs() + [
    ('Method note', SENS_NOTE),
    ('Design', 'Each of the 16 runs is a fully recalibrated 15-trait analysis (one trait removed; working correlation matrix refitted; rewhitened; 10^6-draw null, seed 42; BH across 19,080 genes). The stored per-gene q values for all 19,080 genes and 16 runs are reused unchanged.'),
    ('Retention rule', 'A final CTAS locus is retained in an omitted-trait run at threshold t if at least one of its final member genes (q<0.05 in the final result) has run q < t. q<0.01 and q<0.001 retention are evaluated for loci in the final 46- and 21-locus sets.'),
    ('Result', f"Mean {summ.loci_retained_q05_of_91.mean():.1f}/91 final loci retained per omitted trait at q<0.05 (range {summ.loci_retained_q05_of_91.min()}-{summ.loci_retained_q05_of_91.max()}); {len(stable)} loci retained in all 16 runs. Lowest retention when omitting {summ.omitted_trait.iloc[0]}."),
    ('Diagnostic only', 'newly_sig_diagnostic lists genes passing q<0.05 in a run that are not final CTAS genes; they cannot alter the final locus set.'),
    ('Sources', 'panel_v2_leave_one_out/loo_results.rds (per-gene q and CTAS for 16 runs plus the 16-trait sensitivity-engine reference), exported by code/21_export_loto.R; run-level matrix diagnostics from the signed leave-one-trait-out workbook.')])
b.sheet('locus_retention_91', lr, '91 final loci x 16 omitted-trait runs (1 = retained at q<0.05), with q<0.01/q<0.001 retention counts for the nested sets.')
b.sheet('representative_q_by_run', rr, 'Representative-gene q in each omitted-trait run.')
b.sheet('omitted_trait_summary', summ, 'Per omitted trait: final loci and genes retained.')
b.sheet('stable_all_16_runs', stable, 'Final loci retained in all 16 runs.')
b.sheet('run_matrix_diagnostics', gstab, 'Run-level working-matrix diagnostics and genome-wide run counts (signed values; sensitivity engine).')
b.sheet('newly_sig_diagnostic', newd, 'Run-significant genes outside the final CTAS gene set (diagnostic only).')
b.num('README', 'mean final loci retained at q<0.05 per omitted trait', f'{summ.loci_retained_q05_of_91.mean():.1f}/91', summ.loci_retained_q05_of_91.mean(), '91 loci', 'loci', 'loo_results.rds')
b.num('README', 'loci retained in all 16 runs', len(stable), len(stable), '91 loci', 'loci', 'loo_results.rds')
b.num('README', 'retention range', f'{summ.loci_retained_q05_of_91.min()}-{summ.loci_retained_q05_of_91.max()}', '', '91 loci', 'loci', 'loo_results.rds')
b.sheet('validation', pd.DataFrame([dict(check='stored q cover 19,080 genes x 16 runs', result=Q.shape == (19080, 16)), dict(check='91 loci', result=len(lr) == 91 and lr.canonical_locus_id.is_unique),
    dict(check='final loci retained in the 16-trait sensitivity reference', result=f"{int(lr.retained_full16_sensitivity_reference.sum())}/91")]), 'Checks.')
b.save('ST7_CTAS_leave_one_trait_out.xlsx')

# ---------------- ST8 obstructive cluster
OS = SRC + 'obstructive_sensitivity/'
red = pd.read_csv(OS + 'reduced_gene_results.csv').set_index('gene'); newsig = pd.read_csv(OS + 'newly_significant.csv', dtype={'canonical_locus_id': str}).set_index('gene')
oldret = pd.read_csv(OS + 'locus_retention.csv').set_index('canonical_locus_id'); assert len(red) == 19080
def region(x):
    if x in newsig.index and isinstance(newsig.at[x, 'canonical_locus_id'], str) and newsig.at[x, 'canonical_locus_id']: return newsig.at[x, 'canonical_locus_id']
    k = gi.at[x, 'canonical_locus_id_or_mapping_key']; return k if str(k).startswith('V2L') else ('RED10:' + newsig.at[x, 'RED10_id'] if x in newsig.index and isinstance(newsig.at[x, 'RED10_id'], str) else k)
sig_red = red.index[red.q_BH < .05]; reg = {x: region(x) for x in sig_red}
lt = []
for lid in l.canonical_locus_id:
    row = dict(canonical_locus_id=lid, representative_gene=reps[lid], locus_q_tier=li.at[lid, 'most_stringent_locus_tier'], rep_reduced_q=red.at[reps[lid], 'q_BH'], rep_retained_q05=red.at[reps[lid], 'q_BH'] < .05)
    for lab, t, f in TH:
        if li.at[lid, f]:
            row[f'regional_retained_{lab}'] = any(reg[x] == lid and red.at[x, 'q_BH'] < t for x in sig_red)
            row[f'member_retained_{lab}'] = locus_present(red.q_BH, lid, t)
    row['signed_regional_retained_q05'] = bool(oldret.at[lid, 'regional_retained_q05'])
    lt.append(row)
lt = pd.DataFrame(lt)
agree = bool((lt['regional_retained_q<0.05'] == lt.signed_regional_retained_q05).all())
rs_ = []
for lab, t, f in TH:
    d = lt[lt[f'regional_retained_{lab}'].notna()]
    rs_.append(dict(q_set=lab, final_loci=len(d), regional_retained=int(d[f'regional_retained_{lab}'].astype(bool).sum()), member_retained=int(d[f'member_retained_{lab}'].astype(bool).sum()),
                    pct_regional=100 * d[f'regional_retained_{lab}'].astype(bool).mean(), reduced_genes_significant=int((red.q_BH < t).sum())))
rs_ = pd.DataFrame(rs_)
outside = [x for x in sig_red if x not in F05]
nd = pd.DataFrame([dict(gene=x, reduced_q=red.at[x, 'q_BH'], final_q=gi.at[x, 'q_BH'], region=reg[x], region_type=('final CTAS locus' if reg[x] in set(l.canonical_locus_id) else 'canonical locus outside the final set' if str(reg[x]).startswith('V2L') else 'reduced-panel-only region (diagnostic)')) for x in outside]).sort_values('reduced_q')
b = Book('ST8')
b.readme('Supplementary Table 8. Obstructive-cluster reduced-panel sensitivity, rebased to the final CTAS loci', headline_pairs() + [
    ('Method note', SENS_NOTE),
    ('Question', 'Does the CTAS architecture remain when the overlapping obstructive/atherosclerotic phenotype cluster is represented by CAD alone?'),
    ('Design', 'Fully recalibrated 10-trait run (removed: MI, angina, PAD, atherosclerosis, ischaemic stroke, TIA; retained: AF, aortic aneurysm, aortic valve, CAD, cardiomyopathy, DVT, HF, hypertension, mitral valve, VTE); fresh 10 x 10 working matrix, 10^6-draw null (seed 42), BH across 19,080 genes. Stored per-gene results reused unchanged.'),
    ('Retention', 'Regional: a final locus is retained at threshold t when any reduced-panel gene with q < t maps into its canonical neighbourhood (stored corrected mapping). Member: at least one final member gene has reduced-panel q < t.'),
    ('Result', '; '.join(f'{r.q_set}: regional {r.regional_retained}/{r.final_loci}, member {r.member_retained}/{r.final_loci}' for r in rs_.itertuples()) + '. No reduced-panel gene reaches q<0.001.'),
    ('Diagnostic only', 'Reduced-panel significant genes outside the final CTAS gene set, including reduced-panel-only regions, are listed for transparency; they cannot enter the final architecture.'),
    ('Sources', 'CTAS V2 paper/_source_assets/obstructive_sensitivity/ (reduced_gene_results.csv, newly_significant.csv, locus_retention.csv).')])
b.sheet('retention_summary', rs_, 'Final 91/46/21 loci retained in the reduced-panel run.')
b.sheet('locus_retention_91', lt, 'Per final locus: regional and member retention; the signed q<0.05 regional call is shown for reconciliation.')
b.sheet('reduced_gene_results', red.reset_index(), 'Reduced-panel CTAS, empirical P and BH q for 19,080 genes (signed values).')
b.sheet('outside_final_diagnostic', nd, 'Reduced-panel q<0.05 genes outside the final CTAS gene set (diagnostic only).')
for r in rs_.itertuples(): b.num('README / retention_summary', f'regional retention {r.q_set}', f'{r.regional_retained}/{r.final_loci}', r.regional_retained, f'{r.final_loci} final loci', 'loci', OS + 'reduced_gene_results.csv')
b.sheet('validation', pd.DataFrame([dict(check='q<0.05 regional calls reproduce the signed calls for all 91 loci', result=agree), dict(check='stored results cover 19,080 genes', result=len(red) == 19080)]), 'Checks.')
b.save('ST8_CTAS_obstructive_cluster_sensitivity.xlsx')
print(gs.drop(columns='changed_genes').to_string()); print(ls.to_string()); print(summ[['omitted_trait', 'loci_retained_q05_of_91', 'loci_retained_q01_of_46', 'loci_retained_q001_of_21']].to_string()); print('stable', len(stable)); print(rs_.to_string(), 'agree', agree)
