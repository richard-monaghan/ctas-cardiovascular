import pandas as pd,re,json,csv,numpy as np
R='<CTAS_WORKSPACE>/'; N=R+'target_pipeline/finngen_corrected_rerun_01/'; O=R+'CTAS_FinnGen/results/magma_results/'
EP=['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']
code=dict(AF='I9_AF',Angina='I9_ANGINA',CAD='I9_CORATHER',Cardiomyopathy='I9_CARDMYO',HF='I9_HEARTFAIL',Hypertension='I9_HYPTENS',MI='I9_MI_STRICT',Stroke='I9_STR',TIA='I9_TIA',VTE='I9_VTE')
def g(pat,t,cast=int):
    m=re.search(pat,t); return cast(m.group(1).replace(',','')) if m else None
g3=pd.read_csv(N+'work/annot_G3.tsv',sep='\t').set_index('endpoint')
rows=[]
for e in EP:
    a=open(N+'annot/%s_annot.log'%e).read(); ao=open(O+'%s_annot.log'%e).read()
    lg=open(N+'magma_out/%s.log'%e).read(); lo=open(O+'%s.log'%e).read()
    ge=pd.read_csv(N+'magma_out/%s.genes.out'%e,sep=r'\s+'); go=pd.read_csv(O+'%s.genes.out'%e,sep=r'\s+')
    ns=ge.NSNPS; nso=go.NSNPS
    ap=ge[ge.GENE=='APOE']
    rows.append(dict(endpoint=e,finngen_file=code[e],
      annot_snps_read=g(r'(\d+) SNP locations read',a),annot_pct_snps_mapped=g(r'\((\d+\.\d+)%\) mapped',a,float),annot_genes_with_snp=g(r'total of (\d+) genes',a),
      old_annot_pct_snps_mapped_INVALID=g(r'\((\d+\.\d+)%\) mapped',ao,float),
      G3_sample_pct_inside_GRCh37_gene=g3.at[e,'pct_inside'],G3_gate=g3.at[e,'gate_G3_ge95pct'],
      magma_snps_with_valid_P=g(r'valid SNP p-values for (\d+) SNPs',lg),magma_dup_snps_dropped=g(r'file contained (\d+) SNPs',lg),magma_ended='End time' in lg,
      genes_analysed=len(ge),old_genes_analysed_INVALID=len(go),genes_in_old_not_new=len(set(go.GENE)-set(ge.GENE)),genes_in_new_not_old=len(set(ge.GENE)-set(go.GENE)),
      NSNPS_median=float(ns.median()),NSNPS_q1=float(ns.quantile(.25)),NSNPS_q3=float(ns.quantile(.75)),old_NSNPS_median_INVALID=float(nso.median()),
      N_values=','.join(map(str,sorted(ge.N.unique()))),APOE_NSNPS=int(ap.NSNPS.iloc[0]) if len(ap) else None,APOE_P=float(ap.P.iloc[0]) if len(ap) else None,
      G4_gate='PASS' if ('End time' in lg and len(ge)>=18000 and len(ap) and int(ap.NSNPS.iloc[0])>=10) else 'FAIL'))
pd.DataFrame(rows).to_csv(N+'ANNOTATION_QC.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).drop(columns=['finngen_file']).T.to_string())
