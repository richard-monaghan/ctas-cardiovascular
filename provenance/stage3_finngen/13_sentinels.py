import random,csv,re,subprocess,collections
import pandas as pd
OUT='<CTAS_WORKSPACE>/target_pipeline/finngen_corrected_rerun_01/'
REF='<REFERENCE_DIR>/reference/'
EPS=['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']
loc={}
for l in open(REF+'NCBI37.3.gene.loc'):
    f=l.split()
    if len(f)>=4: loc[f[0]]=(f[1],int(f[2]),int(f[3]))
sent=['rs429358','rs7412','rs7903146','rs10455872','rs12740374','rs11591147','rs1333049','rs3184504','rs6025','rs1801133','rs9939609','rs699']
random.seed(20261104)
samp={}; assign={}; nann={}
for e in EPS:
    a=collections.defaultdict(list); s=[]; ng=0; ns=0
    for l in open(OUT+'annot/%s_annot.genes.annot'%e):
        if l.startswith('#'): continue
        f=l.rstrip('\n').split('\t')
        if len(f)<3: continue
        ng+=1
        for t in f[2:]:
            for x in t.split(','):
                ns+=1
                if x in sent: a[x].append(f[0])
                if random.random()<0.0125: s.append((f[0],x))
    samp[e]=s; assign[e]=dict(a); nann[e]=(ng,ns)
    print(e,ng,ns,len(s),flush=True)
bim=pd.read_csv(REF+'g1000_eur.bim',sep=r'\s+',header=None,usecols=[0,1,3],names=['bchr','SNP','bpos'],dtype={'bchr':str,'SNP':str})
bi=bim.set_index('SNP')
rows=[]; qc=[]
def geneat(c,p): return [g for g,(gc,s,e_) in loc.items() if gc==c and s<=p<=e_]
exp={}
for rs in sent:
    if rs in bi.index: exp[rs]=(bi.at[rs,'bchr'],int(bi.at[rs,'bpos']),geneat(bi.at[rs,'bchr'],int(bi.at[rs,'bpos'])))
for e in EPS:
    lay=subprocess.run(['awk','-F\t','BEGIN{s="%s"; n=split(s,a," "); for(i=1;i<=n;i++)h[a[i]]=1} ($1 in h){print}'%' '.join(sent),OUT+'input_lifted/%s.magma37.txt'%e],capture_output=True,text=True).stdout.strip().split('\n')
    lay={l.split('\t')[0]:l.split('\t') for l in lay if l}
    for rs in sent:
        want=({'rs429358':('19',45411941,['APOE']),'rs7412':('19',45412079,['APOE'])}.get(rs) or exp.get(rs))
        present=rs in lay
        got=(lay[rs][1],int(lay[rs][2])) if present else None
        ag=assign[e].get(rs,[])
        if want is None: status='not_testable_no_BIM_entry'
        elif not present: status='not_testable_absent_from_endpoint' if rs not in('rs429358','rs7412') else 'FAIL_mandatory_variant_absent'
        else:
            posok=(got[0]==want[0] and got[1]==want[1]); geneok=any(g in ag for g in want[2]) if want[2] else None
            status='PASS' if posok and (geneok or geneok is None) else 'FAIL'
        rows.append(dict(endpoint=e,variant=rs,expected_chr=want[0] if want else '',expected_GRCh37_pos=want[1] if want else '',expected_gene='/'.join(want[2]) if want else '',in_MAGMA_input=present,input_chr=got[0] if got else '',input_pos=got[1] if got else '',assigned_gene_in_annot='/'.join(ag) or 'NONE',status=status))
    s=pd.DataFrame(samp[e],columns=['gene','SNP']).join(bi,on='SNP',how='inner')
    s['inside']=[ (g in loc and loc[g][0]==c and loc[g][1]<=p<=loc[g][2]) for g,c,p in zip(s.gene,s.bchr,s.bpos)]
    qc.append(dict(endpoint=e,annot_gene_lines=nann[e][0],annot_snp_tokens=nann[e][1],sample_assignments_with_bim_pos=len(s),sample_inside_GRCh37_gene_interval=int(s.inside.sum()),pct_inside=round(100*s.inside.mean(),3),gate_G3_ge95pct='PASS' if 100*s.inside.mean()>=95 else 'FAIL'))
    print(qc[-1],flush=True)
pd.DataFrame(rows).to_csv(OUT+'SENTINEL_CHECKS.tsv',sep='\t',index=False)
pd.DataFrame(qc).to_csv(OUT+'work/annot_G3.tsv',sep='\t',index=False)
d=pd.DataFrame(rows); print(d.status.value_counts()); print(d[d.status.str.startswith('FAIL')])
