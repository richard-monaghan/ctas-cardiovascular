import pandas as pd,numpy as np,math,json
from scipy.stats import fisher_exact
R='<CTAS_WORKSPACE>/'; OLD=R+'outputs/final_ctas_finngen_20261021/'; NEW=R+'target_pipeline/finngen_corrected_rerun_01/'; A=NEW+'work/analysis/'
EP=['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']
AB={'AF':'AF','Angina':'Ang','CAD':'CAD','Cardiomyopathy':'CM','HF':'HF','Hypertension':'HTN','MI':'MI','Stroke':'Stroke','TIA':'TIA','VTE':'VTE'}
rd=lambda p: pd.read_csv(p,sep='\t',keep_default_na=False,na_values=[''],dtype={'canonical_locus_id_or_mapping_key':str})
INV='INVALID_UNCORRECTED_BUILD'
def wilson(k,n,z=1.959964):
    if n==0: return (np.nan,np.nan)
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return (c-h,c+h)
def woolf(a,n1,b,n2,z=1.959964):
    c,d=n1-a,n2-b
    if min(a,b,c,d)==0: a,b,c,d=a+.5,b+.5,c+.5,d+.5
    o=(a*d)/(b*c); se=math.sqrt(1/a+1/b+1/c+1/d); return o,math.exp(math.log(o)-z*se),math.exp(math.log(o)+z*se)
# ---- sets
so,sn=rd(OLD+'gene_nested_set_summary.tsv').set_index('set'),rd(A+'gene_nested_set_summary.tsv').set_index('set')
to,tn=rd(OLD+'gene_enrichment_tests.tsv'),rd(A+'gene_enrichment_tests.tsv')
rows=[]
for s in ['q<0.05','q<0.01','q<0.001','background_q>=0.05']:
    for k in (1,2,3):
        on,od=int(so.at[s,f'supported_ge{k}']),int(so.at[s,'n_evaluable']); cn,cd=int(sn.at[s,f'supported_ge{k}']),int(sn.at[s,'n_evaluable'])
        lo,hi=wilson(cn,cd); olo,ohi=wilson(on,od)
        rows.append(dict(unit='gene (nested set)',set=s,min_endpoints_p01=k,old_supported=on,old_evaluable=od,old_pct=100*on/od,old_status=INV,corrected_supported=cn,corrected_evaluable=cd,corrected_pct=100*cn/cd,corrected_pct_CI95_low=100*lo,corrected_pct_CI95_high=100*hi,diff_pct_points=100*(cn/cd-on/od),old_n_not_in_finngen=int(so.at[s,'n_not_in_finngen']),corrected_n_not_in_finngen=int(sn.at[s,'n_not_in_finngen'])))
OS=pd.DataFrame(rows)
# tests with CIs
tt=[]
for _,r in tn[tn.section=='nested_set_vs_background'].iterrows():
    o=to[(to.comparison==r.comparison)&(to.min_endpoints_p01==r.min_endpoints_p01)].iloc[0]
    orr,lo,hi=woolf(int(r.set_supported),int(r.set_evaluable),int(r.background_supported),int(r.background_evaluable))
    tt.append(dict(comparison=r.comparison,min_endpoints_p01=int(r.min_endpoints_p01),corrected_set_prop=r.set_prop,corrected_background_prop=r.background_prop,corrected_OR=r.odds_ratio,corrected_OR_CI95_low=lo,corrected_OR_CI95_high=hi,corrected_fisher_p_one_sided=r.fisher_p_one_sided,corrected_BH_q=r.bh_q_across_9_tests,old_OR=o.odds_ratio,old_fisher_p=o.fisher_p_one_sided,old_BH_q=o.bh_q_across_9_tests,old_status=INV,fold_enrichment_corrected=r.fold_enrichment,fold_enrichment_old=o.fold_enrichment))
TT=pd.DataFrame(tt)
for _,r in tn[tn.section!='nested_set_vs_background'].iterrows():
    tt.append(dict(comparison=r.section+': '+r.comparison,min_endpoints_p01=int(r.min_endpoints_p01),corrected_set_prop=r.set_prop,corrected_background_prop=r.background_prop,corrected_OR=r.odds_ratio,corrected_fisher_p_one_sided=r.fisher_p_one_sided,old_status=INV,fold_enrichment_corrected=r.fold_enrichment))
TT=pd.DataFrame(tt)
# dist
dn=rd(A+'gene_support_distribution.tsv'); dold=rd(OLD+'gene_support_distribution.tsv')
# ---- loci
lo_,ln=rd(OLD+'locus_representative_finngen.tsv'),rd(A+'locus_representative_finngen.tsv')
m=ln.merge(lo_[['canonical_locus_id','n_endpoints_p01','min_P','min_P_endpoint','endpoints_p01','supported_ge3','status']],on='canonical_locus_id',suffixes=('_corrected','_old'))
m=m.rename(columns={'n_endpoints_p01_old':'old_n_endpoints_p01_INVALID','min_P_old':'old_min_P_INVALID','min_P_endpoint_old':'old_min_P_endpoint_INVALID','endpoints_p01_old':'old_endpoints_INVALID','supported_ge3_old':'old_supported_ge3_INVALID','status_old':'old_status_INVALID','n_endpoints_p01_corrected':'corrected_n_endpoints_p01','min_P_corrected':'corrected_min_P','min_P_endpoint_corrected':'corrected_min_P_endpoint','endpoints_p01_corrected':'corrected_endpoints','supported_ge3_corrected':'corrected_supported_ge3','status_corrected':'corrected_status'})
keep=['canonical_locus_id','representative_gene','representative_q_BH','locus_tier','corrected_status','corrected_n_endpoints_p01','corrected_min_P','corrected_min_P_endpoint','corrected_endpoints','corrected_supported_ge3','old_status_INVALID','old_n_endpoints_p01_INVALID','old_min_P_INVALID','old_min_P_endpoint_INVALID','old_endpoints_INVALID','old_supported_ge3_INVALID']
m[keep].to_csv(NEW+'CORRECTED_LOCUS_SUPPORT.tsv',sep='\t',index=False)
lso,lsn=rd(OLD+'locus_nested_set_summary.tsv'),rd(A+'locus_nested_set_summary.tsv')
for _,r in lsn.iterrows():
    o=lso[lso.set==r.set].iloc[0]
    for k in (1,2,3):
        cn,cd=int(r[f'supported_ge{k}']),int(r.n_evaluable); lo,hi=wilson(cn,cd); on,od=int(o[f'supported_ge{k}']),int(o.n_evaluable)
        OS.loc[len(OS)]=dict(unit='locus representative',set=r.set,min_endpoints_p01=k,old_supported=on,old_evaluable=od,old_pct=100*on/od,old_status=INV,corrected_supported=cn,corrected_evaluable=cd,corrected_pct=100*cn/cd,corrected_pct_CI95_low=100*lo,corrected_pct_CI95_high=100*hi,diff_pct_points=100*(cn/cd-on/od))
OS.to_csv(NEW+'OLD_VS_CORRECTED_SUPPORT.tsv',sep='\t',index=False)
# ---- tier summary
ts=[]
for s in ['q<0.05','q<0.01','q<0.001','background_q>=0.05']:
    r=sn.loc[s]; lo,hi=wilson(int(r.supported_ge3),int(r.n_evaluable))
    ts.append(dict(set=s,n_total=int(r.n_total),n_evaluable=int(r.n_evaluable),n_not_in_finngen=int(r.n_not_in_finngen),supported_ge1=int(r.supported_ge1),supported_ge2=int(r.supported_ge2),supported_ge3=int(r.supported_ge3),pct_ge3=100*r.prop_ge3,pct_ge3_CI95_low=100*lo,pct_ge3_CI95_high=100*hi,mean_endpoints=r.mean_n_endpoints_p01,median_endpoints=r.median_n_endpoints_p01,old_pct_ge3_INVALID=so.at[s,'pct_ge3_exact']))
TS=pd.DataFrame(ts); TS.to_csv(NEW+'CORRECTED_TIER_SUMMARY.tsv',sep='\t',index=False)
TT.to_csv(NEW+'work/enrichment_with_CI.tsv',sep='\t',index=False)
# append enrichment to tier summary file as second section
with open(NEW+'CORRECTED_TIER_SUMMARY.tsv','a') as f:
    f.write('\n# section: enrichment vs background (one-sided Fisher; BH across 9; Woolf 95% CI)\n'); TT.to_csv(f,sep='\t',index=False)
# ---- Table 3
core=rd(A+'strict_core_21_locus_finngen.tsv'); oc=rd(OLD+'strict_core_21_locus_finngen.tsv').set_index('canonical_locus_id')
def fmt(r):
    if r.status!='evaluable': return 'Not evaluable'
    k=int(r.n_endpoints_p01)
    return f'{k}/10' if k==0 else f'{k}/10: '+', '.join(AB[x.split('=')[0]] for x in r.endpoints_p01.split(';'))
rows=[]
for _,r in core.iterrows():
    o=oc.loc[r.canonical_locus_id]; rows.append(dict(canonical_locus_id=r.canonical_locus_id,representative=r.representative_gene,corrected_table3_cell=fmt(r),corrected_n=r.n_endpoints_p01,corrected_min_P=r.min_P,corrected_min_P_endpoint=r.min_P_endpoint,corrected_ge3=r.supported_ge3,old_table3_cell_INVALID=fmt(o.rename(None)) if True else '',old_n_INVALID=o.n_endpoints_p01,old_min_P_INVALID=o.min_P,old_min_P_endpoint_INVALID=o.min_P_endpoint,old_ge3_INVALID=o.supported_ge3))
T3=pd.DataFrame(rows); T3['cell_changed']=T3.corrected_table3_cell!=T3.old_table3_cell_INVALID
T3.to_csv(NEW+'CORRECTED_TABLE3_FINNGEN.tsv',sep='\t',index=False)
# ---- ledger deliverable
L=rd(A+'ctas_gene_finngen_ledger.tsv'); Lo=rd(OLD+'ctas_gene_finngen_ledger.tsv').set_index('gene')
L['old_n_endpoints_p01_INVALID']=L.gene.map(Lo.n_endpoints_p01); L['old_supported_ge3_INVALID']=L.gene.map(Lo.supported_ge3_p01)
L.to_csv(NEW+'CORRECTED_FINNGEN_GENE_TRAIT_LEDGER.tsv',sep='\t',index=False)
# ---- anchors & QC
G={e:pd.read_csv(NEW+'magma_out/%s.genes.out'%e,sep=r'\s+').set_index('GENE') for e in EP}
an=[]
for g in ['APOE','APOC1','TOMM40','PVRL2','LPA','CELSR2','CDKN2B','LDLR']:
    for e in EP:
        r=G[e].loc[g] if g in G[e].index else None
        an.append(dict(gene=g,endpoint=e,corrected_P=r.P if r is not None else np.nan,corrected_NSNPS=int(r.NSNPS) if r is not None else np.nan,old_P_INVALID=float(Lo.at[g,'P_'+e]) if g in Lo.index and Lo.at[g,'P_'+e]!='' else np.nan))
AN=pd.DataFrame(an); AN.to_csv(NEW+'work/anchors_long.tsv',sep='\t',index=False)
pc=sum(1 for g in ['APOE','LPA','CDKN2B','CELSR2','LDLR'] if min(AN[(AN.gene==g)&AN.endpoint.isin(['CAD','MI'])].corrected_P)<1e-5)
print('positive controls with P<1e-5 in CAD/MI:',pc)
print(TS.to_string()); print(TT[['comparison','min_endpoints_p01','corrected_OR','corrected_OR_CI95_low','corrected_OR_CI95_high','corrected_BH_q','old_OR']].head(9).to_string())
print(OS[OS.min_endpoints_p01==3].to_string())
print(T3[['representative','corrected_table3_cell','old_table3_cell_INVALID','cell_changed']].to_string())
print(AN[AN.gene=='APOE'].to_string())
json.dump(dict(positive_controls_p1e5=pc),open(NEW+'work/poscontrol.json','w'))
