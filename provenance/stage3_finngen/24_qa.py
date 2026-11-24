import os,hashlib,subprocess,pandas as pd,csv
os.chdir('<CTAS_WORKSPACE>/target_pipeline/finngen_corrected_rerun_01'); R='<CTAS_WORKSPACE>/'
EP=['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']
res=[]; chk=lambda n,ok,d='': res.append((n,'PASS' if ok else 'FAIL',d))
chk('Q1 prespecification hash unchanged',hashlib.sha256(open('RERUN_PRESPECIFICATION.md','rb').read()).hexdigest() in open('FROZEN_PRESPEC_HASH.txt').read())
l=pd.read_csv('LIFTOVER_QC.tsv',sep='\t'); chk('Q2 G1 all ten >=95%% (min %.3f)'%l.pct_agree_GRCh37.min(),(l.gate_G1_ge95pct=='PASS').all() and len(l)==10)
s=pd.read_csv('SENTINEL_CHECKS.tsv',sep='\t'); chk('Q3 sentinels 120/120 PASS',(s.status=='PASS').all() and len(s)==120)
a=pd.read_csv('ANNOTATION_QC.tsv',sep='\t'); chk('Q4 G3 and G4 pass all ten',(a.G3_gate=='PASS').all() and (a.G4_gate=='PASS').all() and (a.N_values==500348).all())
# independent recount from raw genes.out
G={e:pd.read_csv('magma_out/%s.genes.out'%e,sep=r'\s+').set_index('GENE').P for e in EP}
M=pd.DataFrame(G); k=(M<0.01).sum(axis=1)
gt=pd.read_csv(R+'outputs/final_c_primary_freeze_20261021/final_C_gene_results_seed42.tsv',sep='\t')
ok=True; out=[]
for t,th in [('q05',.05),('q01',.01),('q001',.001)]:
    gs=set(gt.gene[gt.q_BH<th]); e=[g for g in gs if g in k.index]; out.append((len(gs),len(e),int((k[e]>=3).sum())))
bg=set(gt.gene[gt.q_BH>=.05]); e=[g for g in bg if g in k.index]; out.append((len(bg),len(e),int((k[e]>=3).sum())))
t=pd.read_csv('CORRECTED_TIER_SUMMARY.tsv',sep='\t',nrows=4)
chk('Q5 independent recount of tier support equals table',[tuple(x) for x in out]==[(int(r.n_total),int(r.n_evaluable),int(r.supported_ge3)) for r in t.itertuples()],str(out))
chk('Q6 gene sets 223/99/38',[o[0] for o in out[:3]]==[223,99,38])
ap=pd.read_csv('work/anchors_long.tsv',sep='\t'); ap=ap[ap.gene=='APOE'].set_index('endpoint')
chk('Q7 APOE table equals genes.out and count 7',all(abs(ap.at[e,'corrected_P']-M.at['APOE',e])<1e-12 for e in EP) and int((M.loc['APOE']<0.01).sum())==7)
chk('Q8 no old gene result used as input',all('magma_out' in open('code/20_analysis.py').read() for _ in [0]) and 'CTAS_FinnGen/results/magma_results' not in open('code/20_analysis.py').read())
chk('Q9 required outputs exist',all(os.path.exists(x) for x in 'START_HERE.md RERUN_PRESPECIFICATION.md INPUT_AND_BUILD_AUDIT.tsv LIFTOVER_QC.tsv SENTINEL_CHECKS.tsv ANNOTATION_QC.tsv CORRECTED_FINNGEN_GENE_TRAIT_LEDGER.tsv OLD_VS_CORRECTED_SUPPORT.tsv CORRECTED_TIER_SUMMARY.tsv CORRECTED_LOCUS_SUPPORT.tsv CORRECTED_TABLE3_FINNGEN.tsv APOE_CORRECTED_RESULT.md MANUSCRIPT_CHANGE_MAP.md EXPLORER_CHANGE_MAP.md BERNARD_RESPONSE_FINAL_DRAFT.md DEVIATIONS.tsv SOURCE_LEDGER.tsv code logs'.split()))
v=pd.read_csv('work/analysis/validation_checks.tsv',sep='\t'); chk('Q10 analysis internal checks %d/%d'%((v.status=='PASS').sum(),len(v)),(v.status=='PASS').all())
def rd(p): return dict(pd.read_csv(p,sep='\t')[['path','sha256']].values)
pre,post=rd('work/PROTECTED_HASHES_pre.tsv'),rd('work/PROTECTED_HASHES_post.tsv')
diff=[p for p in pre if 'apoe_finngen_forensic_01' not in p and pre[p]!=post.get(p)]
chk('Q11 protected inputs/outputs unchanged (%d files)'%len(pre),not diff,';'.join(diff[:3]))
o=subprocess.run(['find',R,'-newer','RERUN_PRESPECIFICATION.md','-type','f','-not','-path','*/finngen_corrected_rerun_01/*','-not','-name','.DS_Store'],capture_output=True,text=True).stdout.split('\n'); o=[x for x in o if x]
chk('Q12 nothing outside folder modified',not o,';'.join(o[:4]))
w=open('QA_REPORT.md','w'); w.write('# QA_REPORT\n\n| check | result | detail |\n|---|---|---|\n')
for r in res: w.write('| %s | %s | %s |\n'%r)
w.write('\n%d/%d PASS\n'%(sum(r[1]=='PASS' for r in res),len(res))); w.close(); print(open('QA_REPORT.md').read())
