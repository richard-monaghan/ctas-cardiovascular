#!/usr/bin/env python3
"""Literature stage 4 (release): per-study and four-study-union support over the FINAL nested CTAS locus sets
(91 / 46 / 21), from the 104-locus support matrix written by 02_match_and_results.py. Paper values (Fig. 5a,
Supplementary Fig. 5a, Supplementary Table 10): 19/31/12/48 union 58; 15/22/11/28 union 33; 8/13/9/17 union 17."""
import csv, os, sys
R = os.environ['CTAS_REPO']; W = os.environ['CTAS_LIT_WORK']
sm = {r['canonical_locus_id']: r for r in csv.DictReader(open(os.path.join(W, 'ctas_104_locus_four_study_support_matrix.tsv')), delimiter='\t')}
L = list(csv.DictReader(open(os.path.join(R, 'results/core/final_C_locus_results_seed42.tsv')), delimiter='\t'))
S = ['Qiao', 'Liu', 'Foley', 'PLEIO']; EXP = {'q<0.05': (19, 31, 12, 48, 58), 'q<0.01': (15, 22, 11, 28, 33), 'q<0.001': (8, 13, 9, 17, 17)}
rows = []; ok = True
for tier, f in [('q<0.05', lambda l: True), ('q<0.01', lambda l: int(l['n_q01_genes']) > 0), ('q<0.001', lambda l: int(l['n_q001_genes']) > 0)]:
    ids = [l['canonical_locus_id'] for l in L if f(l)]
    c = [sum(sm[i][s + '_support'] == 'yes' for i in ids) for s in S]; u = sum(any(sm[i][s + '_support'] == 'yes' for s in S) for i in ids)
    got = tuple(c) + (u,); ok &= got == EXP[tier]
    rows.append(dict(locus_set=tier, n_loci=len(ids), **dict(zip(S, c)), four_study_union=u, paper_values='/'.join(map(str, EXP[tier])), match=got == EXP[tier]))
with open(os.path.join(W, 'four_study_support_final_loci.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
sys.exit(0 if ok else 1)
