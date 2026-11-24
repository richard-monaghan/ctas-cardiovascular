#!/usr/bin/env python3
"""Deterministic remap and resummary of the signed colocalisation evidence onto the final nested CTAS locus sets.
No coloc is rerun; PP values are carried as their signed source text. Writes only to this folder."""
import os, csv, hashlib, platform
import numpy as np, pandas as pd, scipy
from scipy.stats import spearmanr
R = '<CTAS_WORKSPACE>/'; T1 = R + 'outputs/final_c_primary_freeze_20261021/'; SRC = '<CTAS_V2_PAPER_WORKSPACE>/_source_assets/coloc_review/'
OUT = R + 'outputs/final_ctas_colocalisation_20261021/'
LOCF = T1 + 'final_C_locus_results_seed42.tsv'
SIGNED = ['V2_COLOC_LOCUS_ABSTRACTION.tsv', 'V2_COLOC_LONG_RESULTS.tsv', 'V2_COLOC_METHODS.md', 'V2_COLOC_SIGNOFF.md', 'V2_COLOC_VALIDATION.md', 'V2_COLOC_AUTHORITY_CHECKSUMS.txt']
SETS = {'q<0.05': 'in_q05', 'q<0.01': 'in_q01', 'q<0.001': 'in_q001'}
CLASSES = ['strong_PP4', 'moderate_PP4', 'PP3_only', 'neither', 'not_assessed']
sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
checks = []; auth = []
def chk(i, exp, obs, ok, ev): checks.append(dict(check_id=i, expected=exp, observed=str(obs), status='PASS' if ok else 'FAIL', evidence_file=ev))
def w(df, f): df.to_csv(OUT + f, sep='\t', index=False, na_rep='')

# ---- provenance
for f in SIGNED:
    if not os.path.exists(SRC + f): raise SystemExit(f'STOP: signed input missing: {f}')
man = {r['file']: r['sha256'] for r in csv.DictReader(open(T1 + 'MANIFEST.tsv'), delimiter='\t')}
h_loc = sha(LOCF)
chk('FINAL_LOCUS_CHECKSUM_MATCHES_TASK1_MANIFEST', man.get('final_C_locus_results_seed42.tsv'), h_loc, h_loc == man.get('final_C_locus_results_seed42.tsv'), 'input_checksums.tsv')
ins = [LOCF, T1 + 'MANIFEST.tsv'] + [SRC + f for f in SIGNED]
w(pd.DataFrame([dict(sha256=sha(f), bytes=os.path.getsize(f), path=f) for f in ins]), 'input_checksums.tsv')

# ---- final loci
loc = pd.read_csv(LOCF, sep='\t', keep_default_na=False, na_values=[''])
loc['in_q05'] = loc.n_q05_genes > 0; loc['in_q01'] = loc.n_q01_genes > 0; loc['in_q001'] = loc.n_q001_genes > 0
sz = [int(loc[c].sum()) for c in SETS.values()]
nest = bool((~loc.in_q001 | loc.in_q01).all() and (~loc.in_q01 | loc.in_q05).all())
chk('FINAL_QSETS_91_46_21_NESTED', '91/46/21, nested, 91 unique ids', f"{'/'.join(map(str, sz))}; nested={nest}; unique={loc.canonical_locus_id.nunique()}", sz == [91, 46, 21] and nest and loc.canonical_locus_id.nunique() == 91 == len(loc), 'final_locus_colocalisation.tsv')

# ---- signed sources (strings preserved)
A = pd.read_csv(SRC + 'V2_COLOC_LOCUS_ABSTRACTION.tsv', sep='\t', dtype=str, keep_default_na=False)
raw_lines = open(SRC + 'V2_COLOC_LONG_RESULTS.tsv').read().split('\n')
hdr = raw_lines[0].split('\t'); body = [l for l in raw_lines[1:] if l != '']
Lg = pd.DataFrame([l.split('\t') for l in body], columns=hdr); Lg['_line'] = body
num = lambda s: pd.to_numeric(s, errors='coerce')
for c in ['n_snps', 'PP0', 'PP1', 'PP2', 'PP3', 'PP4']: Lg[c + '_num'] = num(Lg[c])
key_dup = int(Lg.duplicated(['canonical_locus_id', 'tissue', 'trait']).sum())
per_locus = Lg.groupby('canonical_locus_id').size()
auth += [dict(check='abstraction_rows_unique_loci', expected='104 / 104', observed=f'{len(A)} / {A.canonical_locus_id.nunique()}', ok=len(A) == 104 == A.canonical_locus_id.nunique()),
         dict(check='long_rows', expected='11648', observed=len(Lg), ok=len(Lg) == 11648),
         dict(check='long_tissues_traits', expected='7 tissues, 16 traits', observed=f'{Lg.tissue.nunique()} tissues, {Lg.trait.nunique()} traits', ok=Lg.tissue.nunique() == 7 and Lg.trait.nunique() == 16),
         dict(check='long_unique_locus_tissue_trait', expected='0 duplicate keys', observed=key_dup, ok=key_dup == 0),
         dict(check='long_112_cells_per_locus', expected='104 loci x 112', observed=f'{per_locus.size} loci; min {per_locus.min()} max {per_locus.max()}', ok=per_locus.size == 104 and per_locus.min() == per_locus.max() == 112),
         dict(check='long_gene_matches_abstraction_gene', expected='0 mismatches', observed=int((Lg.gene != Lg.canonical_locus_id.map(A.set_index('canonical_locus_id').gene)).sum()), ok=int((Lg.gene != Lg.canonical_locus_id.map(A.set_index('canonical_locus_id').gene)).sum()) == 0),
         dict(check='assessed_cells_have_numeric_PP', expected='all assessed rows have PP0-PP4 and n_snps', observed=int(Lg[Lg.disposition == 'assessed'][['PP0_num', 'PP1_num', 'PP2_num', 'PP3_num', 'PP4_num', 'n_snps_num']].isna().any(axis=1).sum()), ok=int(Lg[Lg.disposition == 'assessed'][['PP0_num', 'PP1_num', 'PP2_num', 'PP3_num', 'PP4_num', 'n_snps_num']].isna().any(axis=1).sum()) == 0)]
chk('SOURCE_ABSTRACTION_104_UNIQUE', '104 unique loci', f'{len(A)} rows, {A.canonical_locus_id.nunique()} unique', len(A) == 104 == A.canonical_locus_id.nunique(), 'source_authority_checks.tsv')
chk('SOURCE_LONG_11648_ROWS', '11,648', len(Lg), len(Lg) == 11648, 'source_authority_checks.tsv')

# ---- locus summary from cells (used for both signed-rule check on 104 and final 91)
def summarise(cells):
    rows = []
    for k, d in cells.groupby('canonical_locus_id', sort=True):
        a = d[d.disposition == 'assessed']; n_ass = len(a)
        if n_ass:
            pp4, pp3 = a.PP4_num, a.PP3_num; b = a.loc[pp4.idxmax()]
            mx4, mx3 = pp4.max(), pp3.max()
            cls = 'strong_PP4' if mx4 >= .8 else 'moderate_PP4' if mx4 >= .5 else 'PP3_only' if mx3 >= .5 else 'neither'
            disp = 'assessed'; best = (b.tissue, b.trait, b.PP4, a.loc[pp3.idxmax()].PP3)
        else:
            cls = 'not_assessed'; disp = 'not_assessed_MHC' if (d.disposition == 'not_assessed_MHC').any() else 'not_assessed_no_eQTL'; best = ('', '', '', '')
        rows.append(dict(canonical_locus_id=k, signed_gene=d.gene.iloc[0], coloc_disposition=disp, coloc_class=cls, n_cells=len(d), n_assessed_cells=n_ass,
            n_not_evaluable_cells=len(d) - n_ass, n_PP4_supported_cells=int((a.PP4_num >= .5).sum()) if n_ass else 0, n_PP4_strong_cells=int((a.PP4_num >= .8).sum()) if n_ass else 0,
            n_PP3_supported_cells=int((a.PP3_num >= .5).sum()) if n_ass else 0, max_PP4=best[2], max_PP3=best[3], best_tissue=best[0], best_trait=best[1]))
    return pd.DataFrame(rows)
S104 = summarise(Lg)
cmp = S104.merge(A[['canonical_locus_id', 'coloc_class', 'n_assessed', 'n_PP4_sup', 'n_PP4_strong', 'max_PP4', 'best_tissue', 'best_trait']], on='canonical_locus_id', suffixes=('', '_signed'))
mism = cmp[(cmp.coloc_class != cmp.coloc_class_signed) | (cmp.n_assessed_cells.astype(str) != cmp.n_assessed) | (cmp.n_PP4_supported_cells.astype(str) != cmp.n_PP4_sup) | (cmp.n_PP4_strong_cells.astype(str) != cmp.n_PP4_strong)]
mx_diff = (num(cmp.max_PP4) - num(cmp.max_PP4_signed)).abs().max()
bt_m = int(((cmp.best_tissue != cmp.best_tissue_signed) | (cmp.best_trait != cmp.best_trait_signed))[cmp.coloc_class != 'not_assessed'].sum())
auth += [dict(check='classification_rule_reproduces_signed_class_104', expected='0 loci differ (class, n_assessed, n_PP4_sup, n_PP4_strong)', observed=f'{len(mism)} differ: {list(mism.canonical_locus_id)}', ok=len(mism) == 0),
         dict(check='signed_headline_104_classes', expected='19/14/28/38/5 (strong/moderate/PP3_only/neither/not_assessed)', observed='/'.join(str(int((S104.coloc_class == c).sum())) for c in CLASSES), ok=[int((S104.coloc_class == c).sum()) for c in CLASSES] == [19, 14, 28, 38, 5]),
         dict(check='signed_max_PP4_vs_recomputed_from_cells', expected='max abs diff <= 5e-5 (signed abstraction prints 4 dp)', observed=f'{mx_diff:.2e}', ok=mx_diff <= 5e-5),
         dict(check='signed_best_tissue_trait_vs_recomputed', expected='0 assessed loci differ', observed=bt_m, ok=bt_m == 0)]
w(pd.DataFrame(auth).rename(columns={'ok': 'pass'}), 'source_authority_checks.tsv')
chk('SOURCE_AUTHORITY_ALL', 'all source authority checks pass', f"{sum(a['ok'] for a in auth)}/{len(auth)}", all(a['ok'] for a in auth), 'source_authority_checks.tsv')

# ---- representative alignment
al = loc[['canonical_locus_id', 'best_gene_by_final_C_CTAS']].rename(columns={'best_gene_by_final_C_CTAS': 'final_representative_gene'})
cnt = A.canonical_locus_id.value_counts(); al['n_signed_rows'] = al.canonical_locus_id.map(cnt).fillna(0).astype(int)
al['signed_coloc_gene'] = al.canonical_locus_id.map(A.set_index('canonical_locus_id').gene); al['signed_ensg'] = al.canonical_locus_id.map(A.set_index('canonical_locus_id').ensg)
al['representative_match'] = al.final_representative_gene == al.signed_coloc_gene
w(al, 'representative_alignment.tsv')
chk('ALL_91_MAP_ONCE_TO_SIGNED', '91/91 found exactly once', f'{int((al.n_signed_rows == 1).sum())}/91', bool((al.n_signed_rows == 1).all()), 'representative_alignment.tsv')
chk('REPRESENTATIVE_ALIGNMENT_ZERO_CHANGES', '91/91 match, 0 changes', f"{int(al.representative_match.sum())}/91 match; {int((~al.representative_match).sum())} changes", bool(al.representative_match.all()), 'representative_alignment.tsv')
if not (al.n_signed_rows == 1).all() or not al.representative_match.all():
    w(pd.DataFrame(checks), 'validation_checks.tsv'); raise SystemExit('STOP: representative alignment failed')

# ---- final cell table: signed source rows verbatim + q-set flags
flags = loc.set_index('canonical_locus_id')[['in_q05', 'in_q01', 'in_q001']]
C = Lg[Lg.canonical_locus_id.isin(loc.canonical_locus_id)].copy()
with open(OUT + 'final_cell_colocalisation.tsv', 'w') as o:
    o.write('\t'.join(hdr + ['in_q05', 'in_q01', 'in_q001']) + '\n')
    for _, r in C.iterrows():
        f = flags.loc[r.canonical_locus_id]; o.write(r._line + '\t' + '\t'.join(str(bool(x)).upper() for x in f) + '\n')
src_set = {l for l in body}
back = [l.rsplit('\t', 3)[0] for l in open(OUT + 'final_cell_colocalisation.tsv').read().split('\n')[1:] if l]
chk('FINAL_CELLS_10192', '91 x 7 x 16 = 10,192', len(back), len(back) == 10192 and len(set(back)) == 10192, 'final_cell_colocalisation.tsv')
chk('FINAL_CELLS_BYTE_IDENTICAL_TO_SIGNED', 'every final cell row (minus flags) is a verbatim signed source line', f'{sum(l in src_set for l in back)}/{len(back)} verbatim', all(l in src_set for l in back), 'final_cell_colocalisation.tsv')
C = C.join(flags, on='canonical_locus_id')

# ---- final locus table
S = summarise(C).merge(loc[['canonical_locus_id', 'chromosome', 'locus_start', 'locus_end', 'best_gene_by_final_C_CTAS', 'best_gene_q_BH', 'best_gene_q_tier', 'most_stringent_locus_tier', 'in_q05', 'in_q01', 'in_q001']], on='canonical_locus_id')
S = S.rename(columns={'best_gene_by_final_C_CTAS': 'representative_gene', 'best_gene_q_BH': 'representative_q_BH', 'best_gene_q_tier': 'representative_q_tier', 'most_stringent_locus_tier': 'locus_q_tier'})
S['PP4_supported_locus'] = S.coloc_class.isin(['strong_PP4', 'moderate_PP4'])
S['evaluable'] = S.coloc_class != 'not_assessed'
S['PP4_supported_locus'] = np.where(S.evaluable, S.PP4_supported_locus.map({True: 'TRUE', False: 'FALSE'}), 'not_evaluable')
cols = ['canonical_locus_id', 'chromosome', 'locus_start', 'locus_end', 'representative_gene', 'representative_q_BH', 'representative_q_tier', 'locus_q_tier', 'in_q05', 'in_q01', 'in_q001', 'coloc_disposition', 'coloc_class', 'evaluable', 'PP4_supported_locus',
        'n_cells', 'n_assessed_cells', 'n_not_evaluable_cells', 'n_PP4_supported_cells', 'n_PP4_strong_cells', 'n_PP3_supported_cells', 'max_PP4', 'max_PP3', 'best_tissue', 'best_trait']
S = S[cols].sort_values('canonical_locus_id'); w(S, 'final_locus_colocalisation.tsv')
chk('FINAL_LOCUS_TABLE_91_UNIQUE_WITH_FLAGS', '91 unique loci, three q-set flags', f'{len(S)} rows, {S.canonical_locus_id.nunique()} unique', len(S) == 91 == S.canonical_locus_id.nunique(), 'final_locus_colocalisation.tsv')

# ---- nested summaries (regenerated from the final tables)
ls_rows = []; cs_rows = []; t_rows = []; tr_rows = []; snp_rows = []
for name, fl in SETS.items():
    d = S[S[fl]]; ev = d[d.evaluable]; sup = int(d.coloc_class.isin(['strong_PP4', 'moderate_PP4']).sum())
    r = dict(set=name, n_loci_total=len(d), n_loci_evaluable=len(ev), n_not_assessed=int((~d.evaluable).sum()), not_assessed_loci=';'.join(d.canonical_locus_id[~d.evaluable] + '/' + d.representative_gene[~d.evaluable]))
    for c in CLASSES: r[c] = int((d.coloc_class == c).sum())
    r.update(PP4_supported_loci=sup, prop_supported_of_total=sup / len(d), prop_supported_of_evaluable=sup / len(ev), PP3_only_prop_of_evaluable=r['PP3_only'] / len(ev)); ls_rows.append(r)
    c = C[C[fl]]; a = c[c.disposition == 'assessed']
    cs_rows.append(dict(set=name, n_cells=len(c), n_assessed=len(a), n_not_evaluable=len(c) - len(a), n_no_eqtl=int(c.disposition.str.startswith('no_eqtl').sum()), n_not_assessed_MHC=int((c.disposition == 'not_assessed_MHC').sum()),
                        n_PP4_ge_0p5=int((a.PP4_num >= .5).sum()), n_PP4_ge_0p8=int((a.PP4_num >= .8).sum()), n_PP4_0p5_to_0p8=int(((a.PP4_num >= .5) & (a.PP4_num < .8)).sum()), n_PP3_ge_0p5=int((a.PP3_num >= .5).sum())))
    for grp, store in (('tissue', t_rows), ('trait', tr_rows)):
        for g, cc in c.groupby(grp):
            aa = cc[cc.disposition == 'assessed']
            store.append({'set': name, grp: g, 'n_cells': len(cc), 'n_assessed': len(aa), 'n_PP4_supported_cells': int((aa.PP4_num >= .5).sum()), 'n_PP4_strong_cells': int((aa.PP4_num >= .8).sum()),
                          'n_loci_with_supported_cell': aa.loc[aa.PP4_num >= .5, 'canonical_locus_id'].nunique(), 'n_loci_with_strong_cell': aa.loc[aa.PP4_num >= .8, 'canonical_locus_id'].nunique()})
    rho = spearmanr(np.log10(a.n_snps_num), a.PP4_num).statistic
    snp_rows.append(dict(set=name, n_assessed_cells=len(a), spearman_rho_log10_nsnps_vs_PP4=rho, min_matched_snps=int(a.n_snps_num.min()), max_matched_snps=int(a.n_snps_num.max()), median_matched_snps=float(a.n_snps_num.median()),
                         note='descriptive only; not a claim that PP4 is independent of SNP count'))
a104 = Lg[Lg.disposition == 'assessed']
snp_rows.append(dict(set='reference: all 104 signed loci', n_assessed_cells=len(a104), spearman_rho_log10_nsnps_vs_PP4=spearmanr(np.log10(a104.n_snps_num), a104.PP4_num).statistic, min_matched_snps=int(a104.n_snps_num.min()), max_matched_snps=int(a104.n_snps_num.max()),
                     median_matched_snps=float(a104.n_snps_num.median()), note='signed-source reference matching the current supplementary-figure population'))
LS = pd.DataFrame(ls_rows); CS = pd.DataFrame(cs_rows)
w(LS, 'nested_locus_class_summary.tsv'); w(CS, 'nested_cell_summary.tsv')
w(pd.DataFrame(t_rows).sort_values(['set', 'n_PP4_supported_cells'], ascending=[True, False]), 'tissue_summary_by_qset.tsv')
w(pd.DataFrame(tr_rows).sort_values(['set', 'n_PP4_supported_cells'], ascending=[True, False]), 'trait_summary_by_qset.tsv')
w(pd.DataFrame(snp_rows), 'matched_snp_diagnostic.tsv')
ls = LS.set_index('set'); cs = CS.set_index('set')
expc = {'q<0.05': [17, 14, 25, 31, 4], 'q<0.01': [12, 10, 9, 13, 2], 'q<0.001': [5, 7, 4, 4, 1]}
for k, e in expc.items():
    o = [int(ls.at[k, c]) for c in CLASSES]; chk(f'CLASSES_{k}', '/'.join(map(str, e)), '/'.join(map(str, o)), o == e, 'nested_locus_class_summary.tsv')
o = [f"{ls.at[k, 'PP4_supported_loci']}/{ls.at[k, 'n_loci_total']}" for k in SETS]; chk('PP4_SUPPORTED_LOCI', '31/91, 22/46, 12/21', ', '.join(o), o == ['31/91', '22/46', '12/21'], 'nested_locus_class_summary.tsv')
o = [int(ls.at[k, 'n_loci_evaluable']) for k in SETS]; chk('EVALUABLE_DENOMINATORS', '87, 44, 20', o, o == [87, 44, 20], 'nested_locus_class_summary.tsv')
o = [f"{cs.at[k, 'n_PP4_ge_0p5']}/{cs.at[k, 'n_PP4_ge_0p8']}" for k in SETS]; chk('SUPPORTED_STRONG_CELLS', '150/83, 111/59, 47/25', ', '.join(o), o == ['150/83', '111/59', '47/25'], 'nested_cell_summary.tsv')
o = [f"{cs.at[k, 'n_cells']}/{cs.at[k, 'n_assessed']}" for k in SETS]; chk('ALL_AND_ASSESSED_CELLS', '10192/9056, 5152/4528, 2352/2000', ', '.join(o), o == ['10192/9056', '5152/4528', '2352/2000'], 'nested_cell_summary.tsv')

# supported-cell detail (q<0.05 assessed, PP4 >= 0.5)
sd = C[(C.disposition == 'assessed') & (C.PP4_num >= .5)].copy()
sd['PP4_strength'] = np.where(sd.PP4_num >= .8, 'strong', 'moderate')
sd = sd[['canonical_locus_id', 'gene', 'ensg', 'tissue', 'trait', 'n_snps', 'PP0', 'PP1', 'PP2', 'PP3', 'PP4', 'PP4_strength', 'in_q05', 'in_q01', 'in_q001']].assign(_s=sd.PP4_num).sort_values(['canonical_locus_id', '_s'], ascending=[True, False]).drop(columns='_s')
w(sd, 'supported_cell_detail.tsv')
chk('SUPPORTED_CELL_DETAIL_150', '150 q<0.05 assessed cells with PP4 >= 0.5', len(sd), len(sd) == 150, 'supported_cell_detail.tsv')

# strict core
core = S[S.in_q001][['canonical_locus_id', 'representative_gene', 'representative_q_BH', 'locus_q_tier', 'coloc_disposition', 'coloc_class', 'evaluable', 'PP4_supported_locus', 'n_assessed_cells', 'n_PP4_supported_cells', 'n_PP4_strong_cells', 'max_PP4', 'max_PP3', 'best_tissue', 'best_trait']].sort_values('representative_q_BH')
w(core, 'strict_core_21_locus_colocalisation.tsv')
na = core[~core.evaluable]
chk('STRICT_CORE_21_UNIQUE', '21 unique; sole not-assessed = V2L036/ATP6V1G2 (MHC)', f"{len(core)} rows, {core.canonical_locus_id.nunique()} unique; not assessed: {list(na.canonical_locus_id + '/' + na.representative_gene + '/' + na.coloc_disposition)}",
    len(core) == 21 == core.canonical_locus_id.nunique() and list(na.canonical_locus_id) == ['V2L036'] and list(na.coloc_disposition) == ['not_assessed_MHC'], 'strict_core_21_locus_colocalisation.tsv')
nas = S[~S.evaluable]
chk('NOT_ASSESSED_EXCLUDED_ONLY_FROM_EVALUABLE', 'not_assessed loci: in totals, out of evaluable denominators, labelled not_evaluable (never scored 0/negative)',
    f"{len(nas)} not assessed; flags {sorted(set(nas.PP4_supported_locus))}; total-evaluable diffs {[int(ls.at[k,'n_loci_total']-ls.at[k,'n_loci_evaluable']) for k in SETS]}",
    set(nas.PP4_supported_locus) == {'not_evaluable'} and all(int(ls.at[k, 'n_loci_total'] - ls.at[k, 'n_loci_evaluable']) == int(ls.at[k, 'not_assessed']) for k in SETS), 'final_locus_colocalisation.tsv')

# interpretation ledger
def line(k): return f"{ls.at[k,'PP4_supported_loci']}/{ls.at[k,'n_loci_total']} CTAS loci at {k} ({100*ls.at[k,'prop_supported_of_total']:.1f}%), or {ls.at[k,'PP4_supported_loci']}/{ls.at[k,'n_loci_evaluable']} evaluable ({100*ls.at[k,'prop_supported_of_evaluable']:.1f}%)"
led = [('C1_q05', f"Colocalisation evidence (max PP4 >= 0.5 for the representative-gene eQTL and a CVD trait) at {line('q<0.05')}: {ls.at['q<0.05','strong_PP4']} strong, {ls.at['q<0.05','moderate_PP4']} moderate, {ls.at['q<0.05','PP3_only']} PP3-only, {ls.at['q<0.05','neither']} neither, {ls.at['q<0.05','not_assessed']} not assessed.", 'nested_locus_class_summary.tsv'),
       ('C2_q01', f"At q<0.01: {line('q<0.01')}; {ls.at['q<0.01','strong_PP4']} strong, {ls.at['q<0.01','moderate_PP4']} moderate, {ls.at['q<0.01','PP3_only']} PP3-only, {ls.at['q<0.01','neither']} neither, {ls.at['q<0.01','not_assessed']} not assessed.", 'nested_locus_class_summary.tsv'),
       ('C3_q001', f"At q<0.001 (strict core): {line('q<0.001')}; {ls.at['q<0.001','strong_PP4']} strong, {ls.at['q<0.001','moderate_PP4']} moderate, {ls.at['q<0.001','PP3_only']} PP3-only, {ls.at['q<0.001','neither']} neither, 1 not assessed (V2L036/ATP6V1G2, MHC).", 'strict_core_21_locus_colocalisation.tsv'),
       ('C4_cells', f"Supported / strong cells: {cs.at['q<0.05','n_PP4_ge_0p5']}/{cs.at['q<0.05','n_PP4_ge_0p8']} of {cs.at['q<0.05','n_assessed']:,} assessed cells (q<0.05), {cs.at['q<0.01','n_PP4_ge_0p5']}/{cs.at['q<0.01','n_PP4_ge_0p8']} of {cs.at['q<0.01','n_assessed']:,} (q<0.01), {cs.at['q<0.001','n_PP4_ge_0p5']}/{cs.at['q<0.001','n_PP4_ge_0p8']} of {cs.at['q<0.001','n_assessed']:,} (q<0.001). The three sets are nested.", 'nested_cell_summary.tsv'),
       ('C5_denominators', 'Report both denominators: all loci in the set and evaluable loci. Loci classed not_assessed (MHC exclusion, or no cis-eQTL for the representative gene in any of the 7 tissues) are excluded only from the evaluable denominator and are never counted or described as negative results.', 'nested_locus_class_summary.tsv'),
       ('C6_PP3', 'PP3-only means the evidence favours distinct causal variants for the eQTL and the CVD association; it is not a failed colocalisation.', 'final_locus_colocalisation.tsv'),
       ('C7_causal_boundary', 'PP4 supports a shared-variant model for the tested representative-gene eQTL and CVD association; it does not prove the representative gene is causal. Representatives are statistical locus labels, not causal-gene assignments.', 'final_locus_colocalisation.tsv'),
       ('C8_snp_diagnostic', 'The matched-SNP Spearman correlation is descriptive only and is not evidence that PP4 is independent of SNP count.', 'matched_snp_diagnostic.tsv')]
w(pd.DataFrame(led, columns=['statement_id', 'statement', 'source_file']), 'interpretation_ledger.tsv')
w(pd.DataFrame(checks), 'validation_checks.tsv')
open(OUT + 'session_info.txt', 'w').write(f"Python {platform.python_version()} | pandas {pd.__version__} | numpy {np.__version__} | scipy {scipy.__version__} | {platform.platform()}\nscripts: code/00_hash_protected.py code/01_coloc_resummary.py code/02_finalize.py\nNo coloc.abf rerun; no data fetched; PP values carried verbatim from the signed long results.\n")
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
print(LS.to_string()); print(CS.to_string()); print(pd.DataFrame(snp_rows).drop(columns='note').to_string()); print(core.drop(columns=['locus_q_tier']).to_string())
print(pd.DataFrame(auth).to_string()); print(pd.DataFrame(checks)[['check_id', 'observed', 'status']].to_string())
for name in SETS: print(name, pd.DataFrame(t_rows).query('set==@name').nlargest(3, 'n_PP4_supported_cells')[['tissue', 'n_PP4_supported_cells', 'n_PP4_strong_cells']].values.tolist(), pd.DataFrame(tr_rows).query('set==@name').nlargest(4, 'n_PP4_supported_cells')[['trait', 'n_PP4_supported_cells', 'n_PP4_strong_cells']].values.tolist())
