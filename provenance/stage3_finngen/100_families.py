import pandas as pd,numpy as np,math
from scipy.stats import fisher_exact
R='<CTAS_WORKSPACE>/'; N=R+'target_pipeline/finngen_corrected_rerun_01/'; O=R+'outputs/finngen_corrected_assets_20261104/families/'
EP=['AF','Angina','CAD','Cardiomyopathy','HF','Hypertension','MI','Stroke','TIA','VTE']
FAM={'Coronary':['Angina','CAD','MI'],'Arrhythmia':['AF'],'Myocardial':['Cardiomyopathy','HF'],'Pressure':['Hypertension'],'Cerebrovascular':['Stroke','TIA'],'Venous':['VTE']}
FN=list(FAM)
L=pd.read_csv(N+'CORRECTED_FINNGEN_GENE_TRAIT_LEDGER.tsv',sep='\t',dtype=str,keep_default_na=False)
ev=L.in_finngen_union=='True'
P=L[['P_'+e for e in EP]].replace('',np.nan).astype(float); P.columns=EP
F=pd.DataFrame({f:(P[m]<0.01).any(axis=1) for f,m in FAM.items()})
G=pd.DataFrame({'gene':L.gene,'q_BH':L.q_BH.astype(float),'in_finngen_union':ev})
G['n_endpoints_p01']=(P<0.01).sum(axis=1)
for f in FN: G['fam_'+f]=F[f]
G['n_families_p01']=F.sum(axis=1); G['n_noncoronary_families']=F.drop(columns='Coronary').sum(axis=1)
G['families']=[';'.join(f for f in FN if F.at[i,f]) for i in F.index]
G.loc[~ev,['n_endpoints_p01','n_families_p01','n_noncoronary_families']]=np.nan
G['supported_ge3_families']=np.where(~ev,'not_evaluable',np.where(G.n_families_p01>=3,'TRUE','FALSE'))
G['in_q05']=G.q_BH<.05; G['in_q01']=G.q_BH<.01; G['in_q001']=G.q_BH<.001
G['canonical_locus_id_or_mapping_key']=L.canonical_locus_id_or_mapping_key
G.to_csv(O+'FAMILY_GENE_LEDGER.tsv',sep='\t',index=False)
def wilson(k,n,z=1.959964):
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
def woolf(a,n1,b,n2,z=1.959964):
    c,d=n1-a,n2-b
    if min(a,b,c,d)==0: a,b,c,d=a+.5,b+.5,c+.5,d+.5
    o=a*d/(b*c); se=math.sqrt(1/a+1/b+1/c+1/d); return o,math.exp(math.log(o)-z*se),math.exp(math.log(o)+z*se)
E=G[G.in_finngen_union].copy(); E['k']=E.n_families_p01.astype(int); E['nc']=E.n_noncoronary_families.astype(int)
sets={'q<0.05':E.in_q05,'q<0.01':E.in_q01,'q<0.001':E.in_q001,'background_q>=0.05':~E.in_q05}
tot={'q<0.05':G.in_q05,'q<0.01':G.in_q01,'q<0.001':G.in_q001,'background_q>=0.05':~G.in_q05}
rows=[]
for s,m in sets.items():
    k=E.k[m]; n=len(k); nc=E.nc[m]; lo,hi=wilson((k>=3).sum(),n)
    rows.append(dict(set=s,n_total=int(tot[s].sum()),n_evaluable=n,n_not_in_finngen=int(tot[s].sum())-n,supported_ge1=int((k>=1).sum()),supported_ge2=int((k>=2).sum()),supported_ge3=int((k>=3).sum()),
      pct_ge1=100*(k>=1).mean(),pct_ge2=100*(k>=2).mean(),pct_ge3=100*(k>=3).mean(),pct_ge3_CI95_low=100*lo,pct_ge3_CI95_high=100*hi,mean_families=k.mean(),median_families=float(k.median()),
      noncoronary_ge2=int((nc>=2).sum()),pct_noncoronary_ge2=100*(nc>=2).mean()))
TS=pd.DataFrame(rows)
for kk in (1,2):
    TS[f'pct_ge{kk}_CI95_low']=[100*wilson(int(r[f'supported_ge{kk}']),int(r.n_evaluable))[0] for _,r in TS.iterrows()]; TS[f'pct_ge{kk}_CI95_high']=[100*wilson(int(r[f'supported_ge{kk}']),int(r.n_evaluable))[1] for _,r in TS.iterrows()]
TS.to_csv(O+'FAMILY_TIER_SUMMARY.tsv',sep='\t',index=False)
bg=TS.set_index('set').loc['background_q>=0.05']; tests=[]
for s in ['q<0.05','q<0.01','q<0.001']:
    r=TS.set_index('set').loc[s]
    for kk in (1,2,3):
        a,n=int(r[f'supported_ge{kk}']),int(r.n_evaluable); b,nb=int(bg[f'supported_ge{kk}']),int(bg.n_evaluable)
        o,p=fisher_exact([[a,n-a],[b,nb-b]],alternative='greater'); oo,lo,hi=woolf(a,n,b,nb)
        tests.append(dict(comparison=f'{s} vs background q>=0.05',min_families=kk,set_supported=a,set_evaluable=n,background_supported=b,background_evaluable=nb,set_prop=a/n,background_prop=b/nb,fold_enrichment=(a/n)/(b/nb),odds_ratio=oo,OR_CI95_low=lo,OR_CI95_high=hi,fisher_p_one_sided=p))
pv=np.array([t['fisher_p_one_sided'] for t in tests]); order=np.argsort(pv); m_=len(pv); bh=np.empty(m_); run=1.0
for r_,i in enumerate(order[::-1]): rank=m_-r_; run=min(run,pv[i]*m_/rank); bh[i]=run
for t,b_ in zip(tests,bh): t['bh_q_across_9_tests']=b_
TT=pd.DataFrame(tests); TT.to_csv(O+'FAMILY_ENRICHMENT_TESTS.tsv',sep='\t',index=False)
def band(lo,hi):
    m=(E.q_BH>=lo)&(E.q_BH<hi); k=E.k[m]; return int((k>=3).sum()),int(m.sum())
bands={'0.01<=q<0.05':band(.01,.05),'0.001<=q<0.01':band(.001,.01),'q<0.001':band(-1,.001)}
br=[dict(band=k,supported_ge3=a,n=n,pct=100*a/n,CI95_low=100*wilson(a,n)[0],CI95_high=100*wilson(a,n)[1]) for k,(a,n) in bands.items()]
for (n1,n2) in [('q<0.001','0.001<=q<0.01'),('0.001<=q<0.01','0.01<=q<0.05')]:
    (a1,d1),(a2,d2)=bands[n1],bands[n2]; o,p=fisher_exact([[a1,d1-a1],[a2,d2-a2]]); br.append(dict(band=f'{n1} vs {n2} (two-sided Fisher)',supported_ge3=None,n=None,pct=None,OR=o,fisher_p_two_sided=p))
pd.DataFrame(br).to_csv(O+'FAMILY_BANDS.tsv',sep='\t',index=False)
D=[]
for s,m in sets.items():
    k=E.k[m]; n=len(k)
    for j in range(7): D.append(dict(set=s,n_families=j,n_genes=int((k==j).sum()),proportion=(k==j).mean(),cumulative_ge=int((k>=j).sum())))
pd.DataFrame(D).to_csv(O+'FAMILY_DISTRIBUTION.tsv',sep='\t',index=False)
# fam-level rates
FR=[]
for s,m in sets.items():
    for f in FN: FR.append(dict(set=s,family=f,supported=int(E.loc[m,'fam_'+f].sum()),n=int(m.sum()),pct=100*E.loc[m,'fam_'+f].mean()))
pd.DataFrame(FR).to_csv(O+'FAMILY_PER_FAMILY_RATES.tsv',sep='\t',index=False)
# loci
loc=pd.read_csv(R+'outputs/final_c_primary_freeze_20261021/final_C_locus_results_seed42.tsv',sep='\t')
Gi=G.set_index('gene'); rr=[]
for _,r in loc.iterrows():
    g=r.best_gene_by_final_C_CTAS; x=Gi.loc[g]
    rr.append(dict(canonical_locus_id=r.canonical_locus_id,representative_gene=g,locus_tier=r.most_stringent_locus_tier,in_q01=r.n_q01_genes>0,in_q001=r.n_q001_genes>0,evaluable=bool(x.in_finngen_union),
       n_endpoints_p01=x.n_endpoints_p01,n_families_p01=x.n_families_p01,families=x.families,noncoronary_families=x.n_noncoronary_families,supported_ge3_families=x.supported_ge3_families,representative_q_BH=r.best_gene_q_BH))
LR=pd.DataFrame(rr); LR.to_csv(O+'FAMILY_LOCUS_SUPPORT.tsv',sep='\t',index=False)
ls=[]
for s,m in [('q<0.05',LR.locus_tier.notna()),('q<0.01',LR.in_q01),('q<0.001',LR.in_q001)]:
    d=LR[m&LR.evaluable]; k=d.n_families_p01.astype(int); n=len(d)
    ls.append(dict(set=s,n_total=int(m.sum()),n_evaluable=n,supported_ge1=int((k>=1).sum()),supported_ge2=int((k>=2).sum()),supported_ge3=int((k>=3).sum()),pct_ge3=100*(k>=3).mean(),pct_ge3_CI95_low=100*wilson((k>=3).sum(),n)[0],pct_ge3_CI95_high=100*wilson((k>=3).sum(),n)[1],mean_families=k.mean()))
pd.DataFrame(ls).to_csv(O+'FAMILY_LOCUS_SUMMARY.tsv',sep='\t',index=False)
print(TS[['set','n_evaluable','supported_ge1','supported_ge2','supported_ge3','pct_ge3','pct_ge3_CI95_low','pct_ge3_CI95_high','mean_families','noncoronary_ge2','pct_noncoronary_ge2']].round(2).to_string())
print(TT[['comparison','min_families','odds_ratio','OR_CI95_low','OR_CI95_high','fold_enrichment','bh_q_across_9_tests']].round(3).to_string())
print(pd.DataFrame(br).round(3).to_string()); print(pd.DataFrame(ls).round(2).to_string())
print(LR[LR.in_q001][['representative_gene','n_endpoints_p01','n_families_p01','families']].to_string())
