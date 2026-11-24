#!/usr/bin/env python3
"""ST9 colocalisation (Task 3), ST10 prior-study convergence (signed four-study matrix subset to the final sets;
coordinate-Foley primary), ST11 FinnGen external support (Task 2)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stlib import *
g = gene_table(); l = locus_table(); li = l.set_index('canonical_locus_id')
rd = lambda p, **k: pd.read_csv(p, sep='\t', **k)

# ---------------- ST9
L3 = rd(T3 + 'final_locus_colocalisation.tsv'); cells = rd(T3 + 'final_cell_colocalisation.tsv'); sup = rd(T3 + 'supported_cell_detail.tsv')
cls = rd(T3 + 'nested_locus_class_summary.tsv'); cs = rd(T3 + 'nested_cell_summary.tsv'); tis = rd(T3 + 'tissue_summary_by_qset.tsv'); trt = rd(T3 + 'trait_summary_by_qset.tsv')
core = rd(T3 + 'strict_core_21_locus_colocalisation.tsv'); snp = rd(T3 + 'matched_snp_diagnostic.tsv'); val3 = rd(T3 + 'validation_checks.tsv')
c = cls.set_index('set')
b = Book('ST9')
b.readme('Supplementary Table 9. Colocalisation evidence for the final CTAS loci', headline_pairs() + [
    ('Design', '91 final CTAS loci x 1 statistical representative gene x 7 GTEx cis-eQTL datasets x 16 CVD traits = 10,192 cells; coloc.abf (p1 = p2 = 1e-4, p12 = 1e-5), one consistent EUR effect-size mode, hg19 to hg38 liftover. Values are the signed colocalisation evidence remapped to the final loci (representatives identical for 91/91 loci); no colocalisation was rerun.'),
    ('Classes', 'strong_PP4: max PP4 >= 0.8; moderate_PP4: 0.5 <= max PP4 < 0.8; PP3_only: no PP4 >= 0.5 and some PP3 >= 0.5 (evidence favouring distinct variants, not a failed colocalisation); neither; not_assessed (MHC exclusion or no cis-eQTL in any tissue).'),
    ('Headline', '; '.join(f"{k}: {c.at[k, 'PP4_supported_loci']}/{c.at[k, 'n_loci_total']} loci with PP4 >= 0.5 ({c.at[k, 'PP4_supported_loci']}/{c.at[k, 'n_loci_evaluable']} evaluable)" for k in ['q<0.05', 'q<0.01', 'q<0.001']) + '.'),
    ('Classes by set', '; '.join(f"{k} {c.at[k, 'strong_PP4']}/{c.at[k, 'moderate_PP4']}/{c.at[k, 'PP3_only']}/{c.at[k, 'neither']}/{c.at[k, 'not_assessed']}" for k in ['q<0.05', 'q<0.01', 'q<0.001']) + ' (strong/moderate/PP3-only/neither/not assessed).'),
    ('Non-evaluable handling', 'not_assessed loci are excluded only from evaluable denominators and are never counted or described as negative results.'),
    ('Interpretation', 'PP4 supports a shared-variant model for the tested representative-gene eQTL and CVD association; it does not prove the representative gene is causal.'),
    ('Precision', 'PP0-PP4 are carried at full precision from the signed long results.'),
    ('Source', 'outputs/final_ctas_colocalisation_20261021/ (accepted colocalisation remap and resummary).')])
b.sheet('locus_colocalisation', L3, '91 final loci: class, disposition, maxima, best tissue/trait, cell counts, q-set flags.')
b.sheet('class_summary', cls, 'Locus classes and both denominators by nested q-set.')
b.sheet('cell_summary', cs, 'Cell counts by nested q-set.')
b.sheet('supported_cells', sup, 'q<0.05 assessed cells with PP4 >= 0.5 (150), strength and q-set flags.')
b.sheet('strong_cells', sup[sup.PP4_strength == 'strong'], 'Subset with PP4 >= 0.8.')
b.sheet('by_tissue', tis, 'Supported/strong cells by tissue and q-set.')
b.sheet('by_trait', trt, 'Supported/strong cells by trait and q-set.')
b.sheet('strict_core_21', core, 'The 21-locus q<0.001 core.')
b.sheet('matched_snp_diagnostic', snp, 'Descriptive Spearman of log10(matched SNPs) vs PP4 among assessed cells.')
b.sheet('cell_results', cells, 'All 10,192 final cells (signed values, verbatim) with q-set flags.')
b.sheet('validation', val3, 'Accepted colocalisation validation gates.')
for k in ['q<0.05', 'q<0.01', 'q<0.001']:
    b.num('README / class_summary', f'PP4-supported loci {k}', f"{c.at[k, 'PP4_supported_loci']}/{c.at[k, 'n_loci_total']}", c.at[k, 'prop_supported_of_total'], f"{c.at[k, 'n_loci_total']} loci", 'loci', T3 + 'nested_locus_class_summary.tsv')
    b.num('README / class_summary', f'PP4-supported evaluable loci {k}', f"{c.at[k, 'PP4_supported_loci']}/{c.at[k, 'n_loci_evaluable']}", c.at[k, 'prop_supported_of_evaluable'], f"{c.at[k, 'n_loci_evaluable']} evaluable", 'loci', T3 + 'nested_locus_class_summary.tsv')
    b.num('README', f'classes {k}', f"{c.at[k, 'strong_PP4']}/{c.at[k, 'moderate_PP4']}/{c.at[k, 'PP3_only']}/{c.at[k, 'neither']}/{c.at[k, 'not_assessed']}", '', '', 'loci', T3 + 'nested_locus_class_summary.tsv')
st9 = dict(loci=L3.canonical_locus_id.nunique(), cells=len(cells), core=core.canonical_locus_id.nunique())
b.save('ST9_CTAS_colocalisation.xlsx')

# ---------------- ST10
FOUR = CT + 'outputs/four_comparator_harmonisation_20260919/'
m = rd(FOUR + 'ctas_104_locus_four_study_support_matrix.tsv', dtype=str, keep_default_na=False).set_index('canonical_locus_id').loc[l.canonical_locus_id].reset_index()
fm = g[g.in_q05].groupby('canonical_locus_id_or_mapping_key').gene.apply(lambda s: ';'.join(s))
m = m.drop(columns=['ctas_q_tier', 'n_q05_members', 'all_q05_genes', 'n_supporting_studies', 'support_combination']).rename(columns={'Foley_support': 'Foley_coordinate_support', 'Foley_gene_support': 'Foley_gene_symbol_legacy_provenance'})
m.insert(5, 'final_locus_q_tier', m.canonical_locus_id.map(li.most_stringent_locus_tier)); m.insert(6, 'final_q05_members', m.canonical_locus_id.map(fm))
m.insert(7, 'in_q01', m.canonical_locus_id.map(li.in_q01)); m.insert(8, 'in_q001', m.canonical_locus_id.map(li.in_q001))
m['statistical_representative_gene'] = m.canonical_locus_id.map(li.best_gene_by_final_C_CTAS)
Y = lambda col: m[col] == 'yes'
m['n_supporting_studies_of_4'] = Y('Qiao_support').astype(int) + Y('Liu_support').astype(int) + Y('Foley_coordinate_support').astype(int) + Y('PLEIO_support').astype(int)
summ = []
for k, mask in [('q<0.05', np.ones(len(m), bool)), ('q<0.01', m.in_q01.values), ('q<0.001', m.in_q001.values)]:
    d = m[mask]; y = lambda col: int((d[col] == 'yes').sum())
    summ.append(dict(q_set=k, n_loci=len(d), Qiao=y('Qiao_support'), Liu=y('Liu_support'), Foley_coordinate=y('Foley_coordinate_support'), PLEIO=y('PLEIO_support'),
                     four_study_union=int((d.n_supporting_studies_of_4 > 0).sum()), none_of_4=int((d.n_supporting_studies_of_4 == 0).sum()), Foley_gene_symbol_legacy_provenance=y('Foley_gene_symbol_legacy_provenance')))
summ = pd.DataFrame(summ)
exp10 = [[19, 31, 12, 48, 58], [15, 22, 11, 28, 33], [8, 13, 9, 17, 17]]
ok10 = summ[['Qiao', 'Liu', 'Foley_coordinate', 'PLEIO', 'four_study_union']].values.tolist() == exp10
qc = pd.read_csv(SRC + 'R5_review/R5_literature_coordinate_QC.csv', dtype=str); qc = qc[qc.canonical_locus_id.isin(set(l.canonical_locus_id))]
pl = pd.read_csv(SRC + 'R5_review/R5_pleio_support_by_locus.csv', dtype=str, keep_default_na=False); pl = pl[pl.canonical_locus_id.isin(set(l.canonical_locus_id))]
srcman = rd(FOUR + 'comparator_source_manifest.tsv', dtype=str, keep_default_na=False)
notes = pd.DataFrame([('Genome build', 'GRCh37/hg19 for all four studies and the canonical CTAS loci; no lift-over.'),
    ('Qiao (native unit: region)', 'Interval overlap of 38 reported MTAG cardiovascular loci with canonical locus boundaries.'),
    ('Liu (native unit: gene)', 'HGNC-normalised overlap of 180 PLACO pleiotropic genes with each locus\'s significant member genes; an independent 618-SNP coordinate cross-check is provided (Liu_snp_coord_support). All Liu-matched genes remain final member genes of their loci.'),
    ('Foley (native unit: candidate variant)', 'Primary: candidate-variant coordinate containment of 43 CHD colocalisation regions. The earlier gene-symbol call is retained only as labelled provenance; the two differ only at V2L068/APOA5 (rs964184), matched by coordinate containment and a final q<0.05 (not q<0.01) locus.'),
    ('PLEIO (native unit: index variant)', 'Index-variant coordinate containment of 625 reported hits; 13 novel variants checked separately (PLEIO_hits_by_locus).'),
    ('Alias handling', 'HGNC alias and previous symbols normalised to approved symbols on both sides.'),
    ('Not carried over', 'The earlier locus-set PLEIO permutation enrichment was computed on a superseded locus set and is not reproduced here.'),
    ('Framing', 'Convergence with prior cross-trait resources, not independent replication; it never defines CTAS membership. Statistical representatives label loci; matching never relied on representative identity.')], columns=['item', 'detail'])
b = Book('ST10')
b.readme('Supplementary Table 10. Convergence of the final CTAS loci with four prior cross-trait cardiovascular resources', headline_pairs() + [
    ('Scope', 'Signed per-locus support calls for Qiao (MTAG loci), Liu (PLACO genes), Foley (CHD colocalisation regions) and PLEIO (index variants), each matched on its native unit, subset to the final 91/46/21 CTAS loci. No matching or literature search was rerun.'),
    ('Headline', '; '.join(f'{r.q_set} ({r.n_loci}): Qiao {r.Qiao}, Liu {r.Liu}, Foley {r.Foley_coordinate}, PLEIO {r.PLEIO}, four-study union {r.four_study_union}' for r in summ.itertuples()) + '.'),
    ('Foley definition', 'Paper-facing Foley = candidate-variant coordinate containment (12 / 11 / 9). The legacy gene-symbol count (11 / 11 / 9) is shown only as labelled provenance.'),
    ('None of 4', '"None of 4" = comparator-unmatched; not evidence of novelty or absent biology.'),
    ('Source', 'outputs/four_comparator_harmonisation_20260919/ctas_104_locus_four_study_support_matrix.tsv (signed calls); coordinate QC R5_literature_coordinate_QC.csv; PLEIO per-locus R5_pleio_support_by_locus.csv.')])
b.sheet('four_study_matrix', m, 'Per final locus: support call, evidence unit and source identifier for each comparator.')
b.sheet('four_study_summary', summ, 'Supported-locus counts for the nested final sets.')
b.sheet('methods_notes', notes, 'Native-unit matching rules and framing.')
b.sheet('coordinate_QC', qc, 'Coordinate QC of the canonical locus intervals (hg19/hg38).')
b.sheet('PLEIO_hits_by_locus', pl, 'PLEIO index hits per final locus.')
b.sheet('comparator_sources', srcman, 'Comparator source manifest.')
b.sheet('validation', pd.DataFrame([dict(check='19/31/12/48/58, 15/22/11/28/33, 8/13/9/17/17 (Qiao/Liu/coordinate-Foley/PLEIO/union)', result=ok10), dict(check='91 unique loci', result=m.canonical_locus_id.is_unique and len(m) == 91)]), 'Checks.')
for r in summ.itertuples(): b.num('README / four_study_summary', f'Qiao/Liu/Foley(coord)/PLEIO/union {r.q_set}', f'{r.Qiao}/{r.Liu}/{r.Foley_coordinate}/{r.PLEIO}/{r.four_study_union}', '', f'{r.n_loci} loci', 'loci', FOUR + 'ctas_104_locus_four_study_support_matrix.tsv')
b.save('ST10_CTAS_literature_overlap.xlsx')

# ---------------- ST11
gs = rd(T2 + 'gene_nested_set_summary.tsv'); gd = rd(T2 + 'gene_support_distribution.tsv'); ge = rd(T2 + 'gene_enrichment_tests.tsv'); led = rd(T2 + 'ctas_gene_finngen_ledger.tsv', dtype={'endpoints_p01': str}, keep_default_na=False)
miss = rd(T2 + 'gene_missingness.tsv'); lr = rd(T2 + 'locus_representative_finngen.tsv', keep_default_na=False); ln = rd(T2 + 'locus_nested_set_summary.tsv', keep_default_na=False)
core = rd(T2 + 'strict_core_21_locus_finngen.tsv', keep_default_na=False); il = rd(T2 + 'interpretation_ledger.tsv'); v2 = rd(T2 + 'validation_checks.tsv')
s = gs.set_index('set')
head = [(k, int(s.at[k, 'supported_ge3']), int(s.at[k, 'n_evaluable']), s.at[k, 'pct_ge3_exact'], int(s.at[k, 'pct_ge3_rounded'])) for k in ['q<0.05', 'q<0.01', 'q<0.001', 'background_q>=0.05']]
headline = pd.DataFrame([dict(set=k, supported=a, evaluable=n, exact_pct=p, rounded_pct=r, statement=f'{a:,}/{n:,} = {p:.4f}% -> {r}%') for k, a, n, p, r in head])
ok11 = [(a, n) for _, a, n, _, _ in head] == [(91, 218), (50, 97), (19, 37), (1843, 18326)]
b = Book('ST11')
b.readme('Supplementary Table 11. FinnGen external support for the final CTAS genes and loci', headline_pairs() + [
    ('Definition', 'FinnGen external support: MAGMA gene-level P < 0.01 in at least 3 of 10 prespecified FinnGen R12 endpoints (AF, angina, CAD, cardiomyopathy, HF, hypertension, MI, stroke, TIA, VTE).'),
    ('Gene-level headline', '; '.join(f'{k}: {a:,}/{n:,} = {p:.4f}% -> {r}%' for k, a, n, p, r in head) + '.'),
    ('Trend', 'Support increases from q<0.05 to q<0.01 and then plateaus at q<0.001 (51.55% vs 51.35%); the three CTAS sets are nested.'),
    ('Missing genes', 'Genes absent from the FinnGen gene outputs are not evaluable and are never scored as P = 1 or zero endpoints (gene_missingness).'),
    ('Locus-representative analysis', 'A separate companion analysis (one highest-CTAS representative per locus): see the locus_* sheets; it does not replace the gene-level headline.'),
    ('Interpretation', 'Independent-cohort external support for UK Biobank-selected CTAS genes and loci. CTAS is not rerun in FinnGen, so this is not a replication of the CTAS statistic.'),
    ('Source', 'outputs/final_ctas_finngen_20261021/ (accepted FinnGen bundle).')])
b.sheet('gene_headline', headline, 'GENE-LEVEL HEADLINE: support in >=3 of 10 endpoints, nested CTAS gene sets vs background.')
b.sheet('gene_nested_set_summary', gs, 'Gene-level support at >=1, >=2 and >=3 endpoints with evaluable denominators.')
b.sheet('gene_enrichment_tests', ge, 'One-sided Fisher tests vs fixed background (BH over nine tests), disjoint-band checks.')
b.sheet('gene_support_distribution', gd, 'Distribution of endpoint counts 0-10 by set.')
b.sheet('gene_missingness', miss, 'Genes lacking FinnGen coverage (not evaluable).')
b.sheet('gene_ledger', led, 'All 19,080 genes: support fields and all ten endpoint P values (blank = not available).')
b.sheet('locus_rep_summary', ln, 'LOCUS-REPRESENTATIVE COMPANION ANALYSIS: nested 91/46/21 locus sets, evaluable denominators.')
b.sheet('locus_rep_detail', lr, 'Per final locus representative: availability, endpoint counts and P<0.01 endpoints.')
b.sheet('strict_core_21', core, 'Strict-core 21-locus FinnGen evidence.')
b.sheet('interpretation', il, 'Downstream-safe statements.')
b.sheet('validation', v2, 'Accepted FinnGen validation gates.')
for k, a, n, p, r in head: b.num('README / gene_headline', f'FinnGen support {k}', f'{a:,}/{n:,} ({r}%)', p, f'{n:,} evaluable genes', '% genes', T2 + 'gene_nested_set_summary.tsv')
b.save('ST11_CTAS_FinnGen_support.xlsx')
print('ST9', st9, '| ST10 ok', ok10, summ.values.tolist(), '| ST11 ok', ok11)
open(OUT + 'work/external_checks.tsv', 'w').write(f'st9_loci\t{st9["loci"]}\nst9_cells\t{st9["cells"]}\nst9_core\t{st9["core"]}\nst10_ok\t{ok10}\nst11_ok\t{ok11}\n')
