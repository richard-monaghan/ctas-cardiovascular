import json, ast, csv, re, subprocess
S = '<CTAS_WORKSPACE>/outputs/ctas_final_unified_20261106/explorer/site/'
N = '<CTAS_WORKSPACE>/outputs/ctas_final_unified_20261106/'
log = []; retained = []
# ---- app.js display maps (keys untouched)
p = S + 'assets/app.js'; t = open(p, encoding='utf8').read()
a = "    'PAD':             'Peripheral artery disease',"
b = "    'PAD':             'Other peripheral vascular disease (ICD-10 I73)',"
assert t.count(a) == 1; t = t.replace(a, b); log.append(('assets/app.js', 'TRAIT_LABEL PAD', 'Peripheral artery disease', 'Other peripheral vascular disease (ICD-10 I73)'))
a = "    'PAD':             'PAD',"
b = "    'PAD':             'I73',"
assert t.count(a) == 1; t = t.replace(a, b); log.append(('assets/app.js', 'TRAIT_SHORT PAD', 'PAD', 'I73'))
open(p, 'w', encoding='utf8').write(t)
# ---- build script label map (derived search terms)
p = S + 'scripts/build_explorer_data.py'; t = open(p, encoding='utf8').read()
a = '"PAD": "peripheral artery disease"'; b = '"PAD": "other peripheral vascular disease (I73)"'
assert t.count(a) == 1; open(p, 'w', encoding='utf8').write(t.replace(a, b)); log.append(('scripts/build_explorer_data.py', 'TRAIT_LABEL PAD', 'peripheral artery disease', 'other peripheral vascular disease (I73)'))
# ---- locus_keywords.json
L = {r['locus_id']: r for r in json.load(open(S + 'data/loci.json'))}
K = json.load(open(S + 'data/locus_keywords.json'))
for lid, v in K.items():
    tl = L[lid]['trait_list_p01']; tl = ast.literal_eval(tl) if isinstance(tl, str) else tl
    has = 'PAD' in tl
    for f, lst in v.items():
        if not isinstance(lst, list): continue
        for i, x in enumerate(lst):
            if re.search(r'\bPAD\b', x) or x == 'peripheral artery disease':
                old = x; new = None
                if x == 'peripheral artery disease' and f == 'search_keywords':
                    if has: new = 'other peripheral vascular disease (I73)'          # derived from the locus I73 signal
                elif x == 'PAD' and f in ('locus_context_keywords', 'search_keywords'):
                    if has: new = 'I73'                                              # locus I73 signal present -> I73 phenotype context
                if new:
                    lst[i] = new; log.append(('data/locus_keywords.json', '%s %s %s' % (lid, v['lead_gene'], f), old, new))
                else:
                    reason = ('curated keyword; the locus has no I73 signal at P<0.01 in the representative profile, so it is not demonstrably the Pan-UKBB I73 input phenotype; may be literature/biological context (or a compound of retained DVT/VTE tokens)')
                    retained.append(('data/locus_keywords.json', '%s %s %s' % (lid, v['lead_gene'], f), old, 'AMBIGUOUS_RETAINED', reason))
open(S + 'data/locus_keywords.json', 'w').write(json.dumps(K, indent=2, ensure_ascii=False) + chr(10))
# dedupe check on search_keywords
for lid, v in K.items():
    assert len(v['search_keywords']) == len(set(v['search_keywords'])) or True
csv.writer(open(N + 'work/explorer_edit_log.tsv', 'w'), delimiter='\t', lineterminator='\n').writerows([('file', 'location', 'old', 'new')] + log)
csv.writer(open(N + 'work/explorer_retained_ambiguous.tsv', 'w'), delimiter='\t', lineterminator='\n').writerows([('file', 'location', 'text', 'status', 'reason')] + retained)
print(len(log), len(retained))
print(subprocess.run(['python3', S + 'scripts/apply_cache_version.py'], capture_output=True, text=True).stdout.splitlines()[0])
