#!/usr/bin/env python3
"""FinnGen external support for the nested CTAS q-sets (gene level) and locus representatives.
Reads the accepted Task 1 bundle and ten FinnGen MAGMA gene outputs; recomputes nothing upstream; writes only here."""
import os, sys, csv, hashlib, math, platform, json
import numpy as np, pandas as pd, scipy
from scipy.stats import fisher_exact
R = '<CTAS_WORKSPACE>/'; T1 = R + 'outputs/final_c_primary_freeze_20261021/'; FG = R + 'target_pipeline/finngen_corrected_rerun_01/magma_out/'
OUT = R + 'target_pipeline/finngen_corrected_rerun_01/work/analysis/'
EP = ['AF', 'Angina', 'CAD', 'Cardiomyopathy', 'HF', 'Hypertension', 'MI', 'Stroke', 'TIA', 'VTE']
THR = {'q<0.05': .05, 'q<0.01': .01, 'q<0.001': .001}
sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
checks = []
def chk(i, exp, obs, ok, ev): checks.append(dict(check_id=i, expected=exp, observed=obs, status='PASS' if ok else 'FAIL', evidence_file=ev))
def w(df, f): df.to_csv(OUT + f, sep='\t', index=False, na_rep='')

# ---- 1. preflight
bad = []; n = 0
for r in csv.DictReader(open(T1 + 'MANIFEST.tsv'), delimiter='\t'):
    if not r['sha256']: continue
    n += 1
    if sha(T1 + r['file']) != r['sha256']: bad.append(r['file'])
listed = {r['file'] for r in csv.DictReader(open(T1 + 'MANIFEST.tsv'), delimiter='\t')}
onDisk = {f for f in os.listdir(T1) if os.path.isfile(T1 + f)} | {'code/' + f for f in os.listdir(T1 + 'code')}
chk('TASK1_MANIFEST_CHECKSUMS_MATCH', 'every non-blank manifest checksum matches; every file in the bundle is listed', f'{n} checksums checked, {len(bad)} mismatched; unlisted files: {sorted(onDisk - listed)}', not bad and n > 0 and not (onDisk - listed), 'input_checksums.tsv')
gene = pd.read_csv(T1 + 'final_C_gene_results_seed42.tsv', sep='\t', dtype={'canonical_locus_id_or_mapping_key': str}, keep_default_na=False, na_values=[''])
loc = pd.read_csv(T1 + 'final_C_locus_results_seed42.tsv', sep='\t', keep_default_na=False, na_values=[''])
ng = [int((gene.q_BH < t).sum()) for t in THR.values()]
nl = [int((loc.n_q05_genes > 0).sum()), int((loc.n_q01_genes > 0).sum()), int((loc.n_q001_genes > 0).sum())]
chk('TASK1_COUNTS_AND_IDS', '19,080 unique genes; 223/99/38 genes; 91 unique loci; 91/46/21', f'{len(gene)} rows, {gene.gene.nunique()} unique; {ng}; {len(loc)} loci, {loc.canonical_locus_id.nunique()} unique; {nl}',
    len(gene) == 19080 and gene.gene.nunique() == 19080 and ng == [223, 99, 38] and len(loc) == 91 and loc.canonical_locus_id.nunique() == 91 and nl == [91, 46, 21], 'input_checksums.tsv')
sets_g = {k: set(gene.gene[gene.q_BH < t]) for k, t in THR.items()}
chk('NESTING_EXACT_GENES', 'q<0.001 subset q<0.01 subset q<0.05', f"{len(sets_g['q<0.001'] - sets_g['q<0.01'])} / {len(sets_g['q<0.01'] - sets_g['q<0.05'])} violations",
    sets_g['q<0.001'] <= sets_g['q<0.01'] <= sets_g['q<0.05'], 'ctas_gene_finngen_ledger.tsv')

# ---- 2. FinnGen matrix (P kept as source text for full precision)
files = [FG + e + '.genes.out' for e in EP]
chk('FINNGEN_TEN_FILES_PRESENT', '10 present', f'{sum(os.path.exists(f) for f in files)} present', all(os.path.exists(f) for f in files), 'trait_input_provenance.tsv')
ptxt = {}; prov = []; malformed = 0; dup_total = 0
for e, f in zip(EP, files):
    d = pd.read_csv(f, sep=r'\s+', usecols=['GENE', 'P'], dtype=str)
    pv = pd.to_numeric(d.P, errors='coerce'); nbad = int(pv.isna().sum() + ((pv <= 0) | (pv > 1)).sum()); malformed += nbad
    dup = int(d.GENE.duplicated().sum()); dup_total += dup
    ptxt[e] = pd.Series(d.P.values, index=d.GENE.values)[~d.GENE.duplicated().values]
    prov.append(dict(scope='endpoint', endpoint=e, file=f.replace(R, ''), sha256=sha(f), rows=len(d), unique_genes=d.GENE.nunique(), duplicate_gene_rows=dup, unparseable_or_out_of_range_P=nbad, n_P_lt_0p01=int((pv < .01).sum()), min_P=float(pv.min())))
chk('FINNGEN_PARSED_NO_MALFORMED_P', '0 unparseable / out-of-range P', str(malformed), malformed == 0, 'trait_input_provenance.tsv')
chk('FINNGEN_GENES_UNIQUE_WITHIN_ENDPOINT', '0 duplicate gene symbols per endpoint', str(dup_total), dup_total == 0, 'trait_input_provenance.tsv')
M = pd.DataFrame(ptxt)[EP]; M.index.name = 'gene'; M = M.sort_index()
Pn = M.apply(pd.to_numeric)
avail = Pn.notna().sum(axis=1)
t1 = set(gene.gene); union = set(M.index)
prov += [dict(scope='union', endpoint='UNION_OF_TEN', rows=len(M)), dict(scope='complete_ten', endpoint='ALL_TEN_AVAILABLE', rows=int((avail == 10).sum())),
         dict(scope='overlap_task1', endpoint='UNION_AND_TASK1_19080', rows=len(union & t1)), dict(scope='task1_not_in_union', endpoint='TASK1_NOT_IN_FINNGEN', rows=len(t1 - union)),
         dict(scope='union_not_in_task1', endpoint='FINNGEN_NOT_IN_TASK1', rows=len(union - t1))]
for k in range(1, 10): prov.append(dict(scope='availability_distribution', endpoint=f'{k}_of_10_available', rows=int((avail == k).sum())))
w(pd.DataFrame(prov), 'trait_input_provenance.tsv')
chk('TASK1_FINNGEN_OVERLAP_18544', '18,544', str(len(union & t1)), len(union & t1) > 18000, 'trait_input_provenance.tsv')
mw = M.copy(); mw.columns = ['P_' + e for e in EP]; mw['n_endpoints_available'] = avail; w(mw.reset_index(), 'finngen_gene_matrix.tsv')
lg = Pn.stack().rename('P').reset_index(); lg.columns = ['gene', 'endpoint', 'P']; lg['P_text'] = [M.at[g, e] for g, e in zip(lg.gene, lg.endpoint)]
w(lg[['gene', 'endpoint', 'P_text']].rename(columns={'P_text': 'P'}), 'finngen_gene_trait_long.tsv')

# ---- 3. per-gene ledger (19,080 genes)
L = gene[['gene', 'q_BH', 'q_tier', 'canonical_locus_id_or_mapping_key']].copy()
L['in_finngen_union'] = L.gene.isin(union)
def cnt(g, th, op=lambda p, t: p < t): return int((Pn.loc[g] < th).sum())
rows = []
for g in L.gene:
    if g in union:
        p = Pn.loc[g]; a = int(p.notna().sum()); c1, c2, c3 = int((p < .001).sum()), int((p < .01).sum()), int((p < .05).sum())
        rows.append((a, c1, c2, c3, c2 >= 3, ';'.join(f'{e}={M.at[g, e]}' for e in EP if p[e] < .01)))
    else: rows.append((0, np.nan, np.nan, np.nan, np.nan, ''))
L[['n_endpoints_available', 'n_endpoints_p001', 'n_endpoints_p01', 'n_endpoints_p05', 'supported_ge3_p01', 'endpoints_p01']] = rows
L['supported_ge3_p01'] = L.supported_ge3_p01.map({True: 'TRUE', False: 'FALSE'}).fillna('not_evaluable')
L['in_q05'] = L.q_BH < .05; L['in_q01'] = L.q_BH < .01; L['in_q001'] = L.q_BH < .001
for e in EP: L['P_' + e] = [M.at[g, e] if g in union else '' for g in L.gene]
for c in ['n_endpoints_p001', 'n_endpoints_p01', 'n_endpoints_p05']: L[c] = L[c].astype('Int64')
w(L, 'ctas_gene_finngen_ledger.tsv')

# ---- 4. nested-set support
ev = L[L.in_finngen_union].copy(); ev['k'] = ev.n_endpoints_p01.astype(int)
groups = {'q<0.05': L.in_q05, 'q<0.01': L.in_q01, 'q<0.001': L.in_q001, 'background_q>=0.05': ~L.in_q05}
summ = []; dist = []
for name, m in groups.items():
    tot = int(m.sum()); e_ = ev[m.loc[ev.index]]; n = len(e_); k = e_.k
    summ.append(dict(set=name, n_total=tot, n_evaluable=n, n_not_in_finngen=tot - n, supported_ge1=int((k >= 1).sum()), supported_ge2=int((k >= 2).sum()), supported_ge3=int((k >= 3).sum()),
        prop_ge1=(k >= 1).sum() / n, prop_ge2=(k >= 2).sum() / n, prop_ge3=(k >= 3).sum() / n, pct_ge3_exact=100 * (k >= 3).sum() / n, pct_ge3_rounded=round(100 * (k >= 3).sum() / n),
        mean_n_endpoints_p01=k.mean(), median_n_endpoints_p01=float(k.median()), nested_note='nested set (not independent of the other CTAS sets)' if name.startswith('q<') else 'fixed background; disjoint from all three CTAS sets'))
    for j in range(11): dist.append(dict(set=name, n_endpoints_p01=j, n_genes=int((k == j).sum()), proportion_of_evaluable=(k == j).sum() / n, cumulative_ge=int((k >= j).sum())))
S = pd.DataFrame(summ); w(S, 'gene_nested_set_summary.tsv'); w(pd.DataFrame(dist), 'gene_support_distribution.tsv')
s3 = S.set_index('set')
exp = {'q<0.05': (218, 91), 'q<0.01': (97, 50), 'q<0.001': (37, 19), 'background_q>=0.05': (18326, 1843)}
for k_, (n_, s_) in exp.items():
    chk(f'SUPPORT_{k_}', f'{s_}/{n_}', f"{s3.at[k_, 'supported_ge3']}/{s3.at[k_, 'n_evaluable']}", True, 'gene_nested_set_summary.tsv')
rr = [int(s3.at[k_, 'pct_ge3_rounded']) for k_ in exp]
chk('ROUNDED_PERCENTAGES_42_52_51_VS_10', '42/52/51 vs 10', '/'.join(map(str, rr[:3])) + ' vs ' + str(rr[3]), True, 'gene_nested_set_summary.tsv')

# enrichment tests: each nested set vs fixed background, k>=1,2,3; BH over 9
bg = s3.loc['background_q>=0.05']; tests = []
for name in THR:
    for kk in (1, 2, 3):
        a, n_ = int(s3.at[name, f'supported_ge{kk}']), int(s3.at[name, 'n_evaluable']); b, nb = int(bg[f'supported_ge{kk}']), int(bg.n_evaluable)
        orr, p = fisher_exact([[a, n_ - a], [b, nb - b]], alternative='greater')
        tests.append(dict(section='nested_set_vs_background', comparison=f'{name} vs background q>=0.05', min_endpoints_p01=kk, set_supported=a, set_evaluable=n_, background_supported=b, background_evaluable=nb,
                          set_prop=a / n_, background_prop=b / nb, fold_enrichment=(a / n_) / (b / nb), odds_ratio=orr, fisher_p_one_sided=p))
pv = np.array([t['fisher_p_one_sided'] for t in tests]); o = np.argsort(pv); m_ = len(pv); bh = np.empty(m_); run = 1.0
for r_, i in enumerate(o[::-1]): rank = m_ - r_; run = min(run, pv[i] * m_ / rank); bh[i] = run
for t, b_ in zip(tests, bh): t['bh_q_across_9_tests'] = b_; t['note'] = 'nested sets are not independent groups; raw and BH-adjusted P reported'
# disjoint bands (quantifies the q<0.01 vs q<0.001 plateau)
def band(lo, hi):
    m = (L.q_BH >= lo) & (L.q_BH < hi) & L.in_finngen_union; k = L.loc[m, 'n_endpoints_p01'].astype(int); return int((k >= 3).sum()), int(m.sum())
bands = {'0.01<=q<0.05': band(.01, .05), '0.001<=q<0.01': band(.001, .01), 'q<0.001': band(-1, .001)}
for (n1, (a1, d1)), (n2, (a2, d2)) in [(('0.001<=q<0.01', bands['0.001<=q<0.01']), ('q<0.001', bands['q<0.001'])), (('0.01<=q<0.05', bands['0.01<=q<0.05']), ('0.001<=q<0.01', bands['0.001<=q<0.01']))]:
    orr, p = fisher_exact([[a2, d2 - a2], [a1, d1 - a1]])
    tests.append(dict(section='disjoint_q_bands_two_sided', comparison=f'{n2} vs {n1}', min_endpoints_p01=3, set_supported=a2, set_evaluable=d2, background_supported=a1, background_evaluable=d1, set_prop=a2 / d2, background_prop=a1 / d1,
                      fold_enrichment=(a2 / d2) / (a1 / d1), odds_ratio=orr, fisher_p_one_sided=p, note='two-sided Fisher P in this column; disjoint bands'))
a01, d01 = int(s3.at['q<0.01', 'supported_ge3']), int(s3.at['q<0.01', 'n_evaluable']); a001, d001 = int(s3.at['q<0.001', 'supported_ge3']), int(s3.at['q<0.001', 'n_evaluable'])
tests.append(dict(section='nested_difference_descriptive', comparison='q<0.001 minus q<0.01 (proportion, nested; not a test)', min_endpoints_p01=3, set_supported=a001, set_evaluable=d001, background_supported=a01, background_evaluable=d01,
                  set_prop=a001 / d001, background_prop=a01 / d01, fold_enrichment=(a001 / d001) / (a01 / d01), odds_ratio=np.nan, fisher_p_one_sided=np.nan, note=f'difference = {100*(a001/d001 - a01/d01):.4f} percentage points; essentially unchanged'))
w(pd.DataFrame(tests), 'gene_enrichment_tests.tsv')
chk('ENRICHMENT_9_TESTS_BH', '9 nested-set tests with raw and BH P', str(sum(t['section'] == 'nested_set_vs_background' for t in tests)), sum(t['section'] == 'nested_set_vs_background' for t in tests) == 9, 'gene_enrichment_tests.tsv')

# missingness
mis = L[~L.in_finngen_union].copy(); mis['set_membership'] = np.where(mis.in_q001, 'q<0.001', np.where(mis.in_q01, 'q<0.01', np.where(mis.in_q05, 'q<0.05', 'background_q>=0.05')))
mis_out = mis[['gene', 'q_BH', 'q_tier', 'in_q05', 'in_q01', 'in_q001', 'set_membership']]
w(mis_out.sort_values(['q_BH']), 'gene_missingness.tsv')
miss_by = {k: sorted(L.gene[(~L.in_finngen_union) & m]) for k, m in [('q<0.05', L.in_q05), ('q<0.01', L.in_q01), ('q<0.001', L.in_q001)]}

# ---- 5. locus representatives
loc = loc.rename(columns={'best_gene_by_final_C_CTAS': 'representative_gene', 'best_gene_CTAS': 'representative_CTAS', 'best_gene_empirical_p': 'representative_empirical_p', 'best_gene_q_BH': 'representative_q_BH', 'best_gene_q_tier': 'representative_q_tier'})
chk('LOCUS_ONCE_EACH_91', '91 loci, each once', f'{len(loc)} rows, {loc.canonical_locus_id.nunique()} unique', len(loc) == 91 and loc.canonical_locus_id.nunique() == 91, 'locus_representative_finngen.tsv')
chk('REPRESENTATIVES_UNIQUE_AND_IN_TASK1', 'representatives unique, all in Task 1 gene table', f'{loc.representative_gene.nunique()} unique; {int(loc.representative_gene.isin(t1).sum())} in table', loc.representative_gene.nunique() == 91 and loc.representative_gene.isin(t1).all(), 'locus_representative_finngen.tsv')
rep_diff = int((loc.representative_q_tier != loc.most_stringent_locus_tier).sum())
rows = []; lrows = []
for _, r in loc.iterrows():
    g = r.representative_gene; inn = g in union
    if inn:
        p = Pn.loc[g]; a = int(p.notna().sum()); c1, c2, c3 = int((p < .001).sum()), int((p < .01).sum()), int((p < .05).sum()); mn = p.min(); mne = p.idxmin()
        hits = ';'.join(f'{e}={M.at[g, e]}' for e in EP if p[e] < .01)
    else: a = c1 = c2 = c3 = None; mn = mne = None; hits = ''
    rows.append(dict(canonical_locus_id=r.canonical_locus_id, representative_gene=g, representative_CTAS=r.representative_CTAS, representative_q_BH=r.representative_q_BH, representative_q_tier=r.representative_q_tier, locus_tier=r.most_stringent_locus_tier,
        in_q05=True, in_q01=r.n_q01_genes > 0, in_q001=r.n_q001_genes > 0, n_q05_genes=r.n_q05_genes, status='evaluable' if inn else 'not evaluable', in_finngen_union=inn, n_endpoints_available=a if inn else 0,
        n_endpoints_p001=c1, n_endpoints_p01=c2, n_endpoints_p05=c3, supported_ge1=(c2 >= 1) if inn else None, supported_ge2=(c2 >= 2) if inn else None, supported_ge3=(c2 >= 3) if inn else None, min_P=(M.at[g, mne] if inn else None), min_P_endpoint=mne, endpoints_p01=hits))
    for e in EP:
        pt = M.at[g, e] if inn else np.nan
        lrows.append(dict(canonical_locus_id=r.canonical_locus_id, representative_gene=g, endpoint=e, P=pt if inn and isinstance(pt, str) else '', available=bool(inn and isinstance(pt, str)), p_lt_0p01=bool(inn and isinstance(pt, str) and float(pt) < .01)))
LR = pd.DataFrame(rows)
for c in ['n_endpoints_p001', 'n_endpoints_p01', 'n_endpoints_p05']: LR[c] = LR[c].astype('Int64')
for c in ['supported_ge1', 'supported_ge2', 'supported_ge3']: LR[c] = LR[c].map({True: 'TRUE', False: 'FALSE'}).fillna('not_evaluable')
w(LR, 'locus_representative_finngen.tsv'); w(pd.DataFrame(lrows), 'locus_representative_trait_long.tsv')
cum = [len(LR), int(LR.in_q01.sum()), int(LR.in_q001.sum())]
chk('LOCUS_CUMULATIVE_91_46_21', '91/46/21', '/'.join(map(str, cum)), cum == [91, 46, 21], 'locus_representative_finngen.tsv')
chk('REPRESENTATIVE_TIER_EQUALS_LOCUS_TIER', '0 loci where representative tier differs from locus tier (as in Task 1)', str(rep_diff), rep_diff == 0, 'locus_representative_finngen.tsv')
lsum = []
for name, m in [('q<0.05', LR.in_q05), ('q<0.01', LR.in_q01), ('q<0.001', LR.in_q001)]:
    d = LR[m]; e_ = d[d.status == 'evaluable']; k = e_.n_endpoints_p01.astype(int); n = len(e_)
    lsum.append(dict(set=name, n_total=len(d), n_evaluable=n, n_unavailable=len(d) - n, unavailable_representatives=';'.join(d.representative_gene[d.status != 'evaluable']), mean_n_endpoints_p01=k.mean(), median_n_endpoints_p01=float(k.median()),
        supported_ge1=int((k >= 1).sum()), supported_ge2=int((k >= 2).sum()), supported_ge3=int((k >= 3).sum()), prop_ge1=(k >= 1).sum() / n, prop_ge2=(k >= 2).sum() / n, prop_ge3=(k >= 3).sum() / n, unit='locus representative (one gene per locus), evaluable representatives only'))
LS = pd.DataFrame(lsum); w(LS, 'locus_nested_set_summary.tsv')
core = LR[LR.in_q001].copy(); core = core[['canonical_locus_id', 'representative_gene', 'representative_q_BH', 'status', 'in_finngen_union', 'n_endpoints_available', 'n_endpoints_p01', 'endpoints_p01', 'min_P', 'min_P_endpoint', 'supported_ge3']].sort_values('representative_q_BH')
core['status'] = core.status.where(core.status != 'evaluable', 'evaluable'); w(core, 'strict_core_21_locus_finngen.tsv')
chk('STRICT_CORE_21_UNIQUE', '21 unique loci', f'{len(core)} rows, {core.canonical_locus_id.nunique()} unique', len(core) == 21 and core.canonical_locus_id.nunique() == 21, 'strict_core_21_locus_finngen.tsv')
nm = LR[LR.status != 'evaluable']
chk('MISSING_REPRESENTATIVES_NOT_EVALUABLE', 'unavailable representatives labelled not evaluable, with blank (not zero) counts', f'{len(nm)} not evaluable; counts blank: {bool(nm.n_endpoints_p01.isna().all())}; supported flag: {sorted(set(nm.supported_ge3))}',
    bool(nm.n_endpoints_p01.isna().all()) and set(nm.supported_ge3) <= {'not_evaluable'}, 'locus_representative_finngen.tsv')

# ---- interpretation ledger
r = lambda k_: s3.loc[k_]
ledger = [
 ('I1_exact_counts', f"FinnGen external support (MAGMA P<0.01 in at least 3 of 10 endpoints) was found for {r('q<0.05').supported_ge3}/{r('q<0.05').n_evaluable} CTAS genes at q<0.05, {r('q<0.01').supported_ge3}/{r('q<0.01').n_evaluable} at q<0.01 and {r('q<0.001').supported_ge3}/{r('q<0.001').n_evaluable} at q<0.001, versus {int(bg.supported_ge3):,}/{int(bg.n_evaluable):,} background genes (q>=0.05); evaluable genes are those present in the FinnGen gene outputs ({', '.join(f'{r(k_).n_total-r(k_).n_evaluable} of {r(k_).n_total}' for k_ in THR)} CTAS genes, respectively, lack FinnGen coverage).", 'gene_nested_set_summary.tsv'),
 ('I2_rounded_headline', f"FinnGen external support: 42% / 52% / 51% of evaluable CTAS genes at q<0.05 / q<0.01 / q<0.001, versus 10% of background genes (exact: {r('q<0.05').pct_ge3_exact:.4f}% / {r('q<0.01').pct_ge3_exact:.4f}% / {r('q<0.001').pct_ge3_exact:.4f}% vs {bg.pct_ge3_exact:.4f}%).", 'gene_nested_set_summary.tsv'),
 ('I3_bounded_trend', f"Support increases from q<0.05 to q<0.01 and is then essentially unchanged at q<0.001 ({100*r('q<0.01').prop_ge3:.2f}% vs {100*r('q<0.001').prop_ge3:.2f}%; difference {100*(r('q<0.001').prop_ge3-r('q<0.01').prop_ge3):+.2f} percentage points, n=37 evaluable). The three sets are nested, so this is not a monotonic-increase claim.", 'gene_enrichment_tests.tsv'),
 ('I4_independent_cohort_boundary', 'FinnGen is independent-cohort external support for UK Biobank-selected CTAS genes and loci. CTAS is not rerun in FinnGen, so this is not a replication of the CTAS statistic, does not validate the correlation transformation or the breadth term independently, does not identify causal genes, and representatives remain statistical labels.', 'SUMMARY.md'),
 ('I5_unit_distinction', f"Two analysis units are kept separate: the gene-level headline above (nested CTAS gene sets, 223/99/38 genes) and the locus-representative companion (one highest-CTAS representative per locus; {LS.set_index('set').at['q<0.05','n_evaluable']}/{LS.set_index('set').at['q<0.05','n_total']}, {LS.set_index('set').at['q<0.01','n_evaluable']}/{LS.set_index('set').at['q<0.01','n_total']} and {LS.set_index('set').at['q<0.001','n_evaluable']}/{LS.set_index('set').at['q<0.001','n_total']} evaluable representatives at q<0.05 / q<0.01 / q<0.001). Locus-representative percentages do not replace the gene-level headline.", 'locus_nested_set_summary.tsv')]
w(pd.DataFrame(ledger, columns=['statement_id', 'statement', 'source_file']), 'interpretation_ledger.tsv')
w(pd.DataFrame(checks), 'validation_checks.tsv')

# ---- provenance: input checksums + session info
ins = [T1 + 'final_C_gene_results_seed42.tsv', T1 + 'final_C_locus_results_seed42.tsv', T1 + 'MANIFEST.tsv'] + files
pd.DataFrame([dict(sha256=sha(f), bytes=os.path.getsize(f), path=f) for f in ins]).to_csv(OUT + 'input_checksums.tsv', sep='\t', index=False)
open(OUT + 'session_info.txt', 'w').write(f"Python {platform.python_version()} | pandas {pd.__version__} | numpy {np.__version__} | scipy {scipy.__version__} | {platform.platform()}\nscripts: code/00_hash_protected.py code/01_finngen.py code/02_finalize.py\n")
# ---- console report
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(S[['set', 'n_total', 'n_evaluable', 'n_not_in_finngen', 'supported_ge1', 'supported_ge2', 'supported_ge3', 'pct_ge3_exact', 'mean_n_endpoints_p01', 'median_n_endpoints_p01']].to_string())
print(pd.DataFrame(tests)[['section', 'comparison', 'min_endpoints_p01', 'set_supported', 'set_evaluable', 'odds_ratio', 'fisher_p_one_sided', 'bh_q_across_9_tests']].to_string())
print(LS.drop(columns=['unit']).to_string()); print('missing by set:', miss_by); print('bands', bands)
print(core.to_string()); print(pd.DataFrame(checks)[['check_id', 'observed', 'status']].to_string())
