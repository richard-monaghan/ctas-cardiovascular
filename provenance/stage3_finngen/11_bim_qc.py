import pandas as pd,numpy as np,json,csv,sys
OUT='<CTAS_WORKSPACE>/target_pipeline/finngen_corrected_rerun_01/'
EPS=['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']
bim=pd.read_csv('<REFERENCE_DIR>/reference/g1000_eur.bim',sep=r'\s+',header=None,usecols=[0,1,3,4,5],names=['bchr','SNP','bpos','a1','a2'],dtype={'bchr':str,'SNP':str,'a1':str,'a2':str})
nb=len(bim); dup=bim.SNP.duplicated(keep=False); nbdup=int(dup.sum())
bu=bim[~dup].set_index('SNP')
print('bim rows',nb,'dup-id rows',nbdup,flush=True)
comp={'A':'T','T':'A','C':'G','G':'C'}
rows=[]
for e in EPS:
    st=json.load(open(OUT+'work/lift_stats_%s.json'%e))
    d=pd.read_csv(OUT+'input_lifted/%s.lifted.tsv.gz'%e,sep='\t',dtype={'SNP':str,'CHR':str,'REF':str,'ALT':str},usecols=['SNP','CHR','BP37','BP38','REF','ALT'])
    fd=d.SNP.duplicated(keep=False); nfd=int(fd.sum())
    u=d[~fd]; j=u.join(bu,on='SNP',how='inner')
    n=len(j); same_chr=(j.CHR==j.bchr)
    m37=same_chr&(j.BP37==j.bpos); m38=same_chr&(j.BP38==j.bpos)
    pm=j[m37]
    ok=((pm.REF==pm.a1)&(pm.ALT==pm.a2))|((pm.REF==pm.a2)&(pm.ALT==pm.a1))
    cm=((pm.REF.map(comp)==pm.a1)&(pm.ALT.map(comp)==pm.a2))|((pm.REF.map(comp)==pm.a2)&(pm.ALT.map(comp)==pm.a1))
    rows.append(dict(endpoint=e,raw_rows=st['rows'],empty_rsid_excluded=st['empty_rsid'],retained=st['retained'],lift_ok=st['lift_ok'],lift_failed=st['lift_failed'],lift_ok_pct=round(100*st['lift_ok']/st['retained'],3),multi_rsid_strings=st['multi_rsid'],
      rows_in_lifted_layer=len(d),rsid_duplicated_in_layer_rows=nfd,bim_rows_with_duplicated_id=nbdup,n_unique_comparable_rsids=n,
      n_chr_differs=int((~same_chr).sum()),n_agree_GRCh37=int(m37.sum()),pct_agree_GRCh37=round(100*m37.sum()/n,3),n_agree_if_NOT_lifted=int(m38.sum()),pct_agree_if_not_lifted=round(100*m38.sum()/n,3),
      n_pos_disagree_same_chr=int((same_chr&~m37).sum()),allele_set_match=int(ok.sum()),allele_match_via_complement_only=int((~ok&cm).sum()),allele_mismatch=int((~ok&~cm).sum()),
      gate_G1_ge95pct='PASS' if 100*m37.sum()/n>=95 else 'FAIL'))
    print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(OUT+'LIFTOVER_QC.tsv',sep='\t',index=False)
