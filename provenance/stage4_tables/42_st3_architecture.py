#!/usr/bin/env python3
"""ST3 locus architecture: 91 final CTAS loci, 223 significant member genes; family composition, whitened
profiles and ordination from the final-whitener recomputation (code/20_architecture.R)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stlib import *
g = gene_table(); l = locus_table(); SD = OUT + 'source_data/'
A = pd.read_csv(SD + 'architecture_locus_ranked.tsv', sep='\t')
prof = pd.read_csv(SD + 'architecture_trait_profiles.tsv', sep='\t'); ordn = pd.read_csv(SD + 'architecture_ordination.tsv', sep='\t')
struct = pd.read_csv(SD + 'architecture_structure_diagnostics.tsv', sep='\t', keep_default_na=False); sens = pd.read_csv(SD + 'architecture_family_metric_sensitivity.tsv', sep='\t', keep_default_na=False)
cov = pd.read_csv(SD + 'architecture_raw_family_coverage.tsv', sep='\t'); covs = pd.read_csv(SD + 'architecture_raw_family_coverage_summary.tsv', sep='\t')
cmp = pd.read_csv(SD + 'architecture_vs_earlier_summaries.tsv', sep='\t'); mpr = pd.read_csv(SD + 'architecture_MPR_reproduction.tsv', sep='\t')
b = Book('ST3')
summ = A[['canonical_locus_id', 'ctas_rank', 'chromosome', 'locus_start', 'locus_end', 'representative_gene', 'representative_CTAS', 'magnitude_M', 'participation_ratio_PR',
          'representative_empirical_p', 'representative_q_BH', 'locus_q_tier', 'in_q05', 'in_q01', 'in_q001', 'n_q05_genes', 'all_q05_genes', 'Neff_family_breadth', 'largest_share_family', 'leading_family_display']].copy()
li = l.set_index('canonical_locus_id'); summ['n_q01_genes'] = summ.canonical_locus_id.map(li.n_q01_genes); summ['n_q001_genes'] = summ.canonical_locus_id.map(li.n_q001_genes)
fam = A[['canonical_locus_id', 'ctas_rank', 'representative_gene'] + [c for c in A.columns if c.startswith('A_') or c.startswith('S_')] + ['Neff_family_breadth', 'max_family_share', 'largest_share_family', 'leading_family_display']]
mem = g[g.in_q05][['canonical_locus_id_or_mapping_key', 'gene', 'CTAS', 'magnitude_M', 'participation_ratio_PR', 'empirical_p', 'q_BH', 'q_tier']].rename(columns={'canonical_locus_id_or_mapping_key': 'canonical_locus_id'})
mem['is_representative'] = mem.gene.isin(set(l.best_gene_by_final_C_CTAS)); mem = mem.sort_values(['canonical_locus_id', 'CTAS'], ascending=[True, False])
fsum = []
for k, f in [('q<0.05', 'in_q05'), ('q<0.01', 'in_q01'), ('q<0.001', 'in_q001')]:
    d = A[A[f]]; vc = d.largest_share_family.value_counts()
    fsum.append(dict(locus_set=k, n_loci=len(d), median_Neff=d.Neff_family_breadth.median(), min_Neff=d.Neff_family_breadth.min(), max_Neff=d.Neff_family_breadth.max(),
                     **{f'leading_{x}': int(vc.get(x, 0)) for x in ['coronary', 'pressure', 'arrhythmia', 'heart_failure', 'valve', 'aortic_vascular', 'cerebrovascular', 'thrombotic']}))
fsum = pd.DataFrame(fsum)
verdict = struct[struct.metric == 'verdict_hellinger_final'].key.iloc[0]
b.readme('Supplementary Table 3. CTAS locus architecture (91 loci, 223 significant member genes)', headline_pairs() + [
    ('Membership', 'Locus membership, representatives, CTAS, q and nested-set flags come from the final primary result. A gene counts toward a locus only if it is itself significant (q<0.05). Representative = highest-CTAS significant member; a statistical label, not a causal gene.'),
    ('Whitened profiles and families', 'Signed and positive whitened 16-trait profiles, family composition (size-adjusted mean positive evidence per family, shares S_f), Family N_eff = 1 / sum S_f^2 and the leading family were recomputed deterministically with the final ZCA-cor whitener from the unchanged 16-trait evidence. Recomputed M and PR reproduce the final values for all 223 significant genes (max abs difference < 1.1e-14).'),
    ('Ordination', f'Hellinger-transformed family shares (primary) and raw shares (sensitivity): PCA plus hierarchical, k-means and PAM clustering diagnostics (k = 2-10, 200 bootstrap resamples, seed 20260805). Verdict: {verdict.replace("_", " ").lower()} - no discrete cluster structure.'),
    ('Raw-trait coverage', 'Descriptive count of cardiovascular families with at least one raw gene-level P<0.01 for each representative (method-independent; identical to the earlier coverage for all 91 loci).'),
    ('Comparison with earlier summaries', f'Against the earlier summaries derived with the raw-Spearman whitener: representatives identical for 91/91 loci; leading family unchanged for {int((~cmp.leading_family_changed).sum())}/91; Family N_eff differs by at most {cmp.delta_Neff.abs().max():.3f} (median {cmp.delta_Neff.abs().median():.4f}).'),
    ('Units / precision', 'Full numeric precision in data sheets.'),
    ('Sources', 'Accepted primary-result bundle (membership, CTAS, q and the final ZCA-cor whitener); unchanged 16-trait matrix analytic_null_validation/inputs/UKBB_16trait_gene_matrix.tsv; recomputation code/20_architecture.R and source_data/architecture_*.tsv.')])
b.sheet('locus_summary', summ, '91 CTAS loci ranked by representative CTAS.')
b.sheet('gene_membership', mem, '223 significant member genes with locus assignment.')
b.sheet('family_composition', fam, 'Size-adjusted family evidence (A_f) and shares (S_f) from the final whitened profile of each representative.')
b.sheet('family_summary_by_qset', fsum, 'Family N_eff and leading-family counts for the nested 91/46/21 locus sets.')
b.sheet('trait_profiles', prof, 'Signed (w_) and positive-part (v_) whitened 16-trait profiles of each representative (final whitener).')
b.sheet('ordination', ordn, 'Hellinger (primary) and raw-share (sensitivity) principal-component coordinates.')
b.sheet('structure_diagnostics', struct, 'Variance explained, silhouette, ARI and bootstrap stability for k = 2-10.')
b.sheet('metric_sensitivity', sens, 'Mean-positive versus RMS-positive family metric.')
b.sheet('raw_family_coverage', cov, 'Raw P<0.01 family coverage per representative (descriptive).')
b.sheet('raw_coverage_summary', covs, 'Raw-coverage category counts for the nested sets.')
b.sheet('vs_earlier_summaries', cmp, 'Final-whitener family summaries versus the earlier raw-Spearman-whitener summaries (same representatives).')
b.sheet('MPR_reproduction', mpr, 'Recomputed M, PR and CTAS versus the final values (223 genes).')
for r in fsum.itertuples(): b.num('family_summary_by_qset', f'median Family N_eff {r.locus_set}', f'{r.median_Neff:.2f}', r.median_Neff, f'{r.n_loci} loci', 'N_eff', 'source_data/architecture_locus_ranked.tsv')
b.num('README', 'loci / significant member genes', '91 / 223', '91/223', '', 'loci / genes', 'final_C_locus_results_seed42.tsv')
b.num('README', 'leading family unchanged vs earlier', f'{int((~cmp.leading_family_changed).sum())}/91', int((~cmp.leading_family_changed).sum()), '91 loci', 'loci', 'source_data/architecture_vs_earlier_summaries.tsv')
b.num('README', 'max |delta Family N_eff| vs earlier', f'{cmp.delta_Neff.abs().max():.3f}', cmp.delta_Neff.abs().max(), '91 loci', '', 'source_data/architecture_vs_earlier_summaries.tsv')
v = pd.DataFrame([dict(check='91 unique loci', result=summ.canonical_locus_id.is_unique and len(summ) == 91), dict(check='223 unique significant member rows', result=len(mem) == 223 and mem.gene.is_unique),
    dict(check='nested 91/46/21', result=[int(summ.in_q05.sum()), int(summ.in_q01.sum()), int(summ.in_q001.sum())] == [91, 46, 21]),
    dict(check='members of each locus match the final locus table', result=bool((mem.groupby('canonical_locus_id').size() == li.n_q05_genes.reindex(mem.canonical_locus_id.unique())).all())),
    dict(check='recomputed M/PR/CTAS reproduce final values (max abs diff)', result=float(mpr.value.max())), dict(check='representatives identical to earlier summaries', result=bool(cmp.representative_matches_earlier.all()))])
b.sheet('validation', v, 'Checks.')
b.save('ST3_CTAS_locus_architecture.xlsx')
