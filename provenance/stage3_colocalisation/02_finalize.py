#!/usr/bin/env python3
"""Compare protected hashes, scan new prose/filenames for banned terminology, append gates, write MANIFEST.tsv."""
import os, re, csv, hashlib, pandas as pd
O = '<CTAS_WORKSPACE>/outputs/final_ctas_colocalisation_20261021/'
b = pd.read_csv(O + 'protected_hashes_before.tsv', sep='\t'); a = pd.read_csv(O + 'protected_hashes_after.tsv', sep='\t')
b = b.drop(columns='bytes'); a = a.drop(columns='bytes'); m = b.merge(a, on='path', how='outer', suffixes=('_before', '_after'))
m['status'] = ['not readable' if r.sha256_before == 'not readable' and r.sha256_after == 'not readable' else 'UNCHANGED' if (r.sha256_before == r.sha256_after and isinstance(r.sha256_before, str)) else 'CHANGED' for r in m.itertuples()]
m.to_csv(O + 'protected_hash_comparison.tsv', sep='\t', index=False)
ok_h = (m.status != 'CHANGED').all() and len(b) == len(a) == 5
v = pd.read_csv(O + 'validation_checks.tsv', sep='\t', keep_default_na=False)
v = v[~v.check_id.isin(['PROTECTED_HASHES_UNCHANGED', 'NO_BANNED_TERMINOLOGY_IN_NEW_PROSE_OR_FILENAMES', 'NO_RISING_WITH_STRINGENCY_PHRASE'])]
bad = re.compile(r'(?i:\bV[12]\b|estimator|\bmethod [AC]\b|\b(?:CTAS|score) [AC]\b|\b[AC] score\b|Score_C|whitening)|final_C_')  # internal label matched case-sensitively so final_cell_* is not a hit
hits = []
for f in sorted(os.listdir(O)):
    p = O + f
    if os.path.isdir(p): continue
    if bad.search(f): hits.append(f'filename:{f}')
    if f in ('validation_checks.tsv', 'input_checksums.tsv', 'protected_hashes_before.tsv', 'protected_hashes_after.tsv', 'protected_hash_comparison.tsv', 'MANIFEST.tsv', 'final_cell_colocalisation.tsv'): continue
    for i, line in enumerate(open(p, encoding='utf-8', errors='replace'), 1):
        if bad.search(line): hits.append(f'{f}:{i}')
for f in os.listdir(O + 'code'):
    if bad.search(f): hits.append('filename:code/' + f)
phr = [f for f in os.listdir(O) if os.path.isfile(O + f) and f not in ('validation_checks.tsv',) and re.search('rising with stringency', open(O + f, errors='replace').read(), re.I)]
extra = [dict(check_id='PROTECTED_HASHES_UNCHANGED', expected='all protected files identical before/after', observed=f"{(m.status == 'UNCHANGED').sum()} of {len(m)} unchanged", status='PASS' if ok_h else 'FAIL', evidence_file='protected_hash_comparison.tsv'),
         dict(check_id='NO_BANNED_TERMINOLOGY_IN_NEW_PROSE_OR_FILENAMES', expected='no letter-coded method/score labels, version labels or internal audit label in new prose/filenames', observed=f'{len(hits)} hits: {hits[:5]}', status='PASS' if not hits else 'FAIL', evidence_file='SUMMARY.md; interpretation_ledger.tsv; summary tables'),
         dict(check_id='NO_RISING_WITH_STRINGENCY_PHRASE', expected='phrase absent from deliverables', observed=f'{len(phr)} files: {phr}', status='PASS' if not phr else 'FAIL', evidence_file='SUMMARY.md; interpretation_ledger.tsv')]
v = pd.concat([v, pd.DataFrame(extra)]); v.to_csv(O + 'validation_checks.tsv', sep='\t', index=False)
print(v[['check_id', 'status']].to_string()); print('hits', hits, 'phrase', phr)
purpose = {'SUMMARY.md': 'Summary and output report', 'MANIFEST.tsv': 'This manifest (no self checksum)', 'input_checksums.tsv': 'SHA-256 of the final locus table, Task 1 manifest and six signed colocalisation inputs',
 'protected_hashes_before.tsv': 'Five named downstream files, before', 'protected_hashes_after.tsv': 'Five named downstream files, after', 'protected_hash_comparison.tsv': 'Before/after comparison',
 'source_authority_checks.tsv': 'Signed-source structure and classification-rule reproduction', 'representative_alignment.tsv': 'Final representative vs signed colocalisation gene (91 loci)',
 'final_locus_colocalisation.tsv': '91 final loci: class, disposition, maxima, best tissue/trait, cell counts, q-set flags', 'final_cell_colocalisation.tsv': '10,192 signed cell rows verbatim plus q-set flags',
 'supported_cell_detail.tsv': 'q<0.05 assessed cells with PP4 >= 0.5', 'nested_locus_class_summary.tsv': 'Locus classes and both denominators by nested q-set', 'nested_cell_summary.tsv': 'Cell counts by nested q-set',
 'tissue_summary_by_qset.tsv': 'Supported/strong cells by tissue', 'trait_summary_by_qset.tsv': 'Supported/strong cells by trait', 'strict_core_21_locus_colocalisation.tsv': 'Strict-core 21-locus colocalisation table',
 'matched_snp_diagnostic.tsv': 'Descriptive Spearman of log10 matched SNPs vs PP4', 'interpretation_ledger.tsv': 'Downstream-safe statements', 'validation_checks.tsv': 'Gate results', 'session_info.txt': 'Python/package/session information'}
rows = []
for f in sorted(x for x in os.listdir(O) if os.path.isfile(O + x)) + ['code/' + x for x in sorted(os.listdir(O + 'code'))]:
    h = '' if f == 'MANIFEST.tsv' else hashlib.sha256(open(O + f, 'rb').read()).hexdigest(); dim = ''
    if f.endswith('.tsv') and f != 'MANIFEST.tsv':
        d = pd.read_csv(O + f, sep='\t', dtype=str, keep_default_na=False); dim = f'{len(d)} rows x {d.shape[1]} cols'
    rows.append([f, purpose.get(f, 'script' if f.startswith('code/') else ''), dim, h])
w = csv.writer(open(O + 'MANIFEST.tsv', 'w'), delimiter='\t', lineterminator='\n'); w.writerow(['file', 'purpose', 'dimensions', 'sha256']); w.writerows(rows); print(len(rows), 'manifest rows')
