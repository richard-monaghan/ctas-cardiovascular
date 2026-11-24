import sys,os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from figlib import *
FD='<CTAS_WORKSPACE>/outputs/finngen_corrected_assets_20261104/families/'
ts=pd.read_csv(FD+'FAMILY_TIER_SUMMARY.tsv',sep='\t').set_index('set'); dd=pd.read_csv(FD+'FAMILY_DISTRIBUTION.tsv',sep='\t')
items=[('q<0.05','q<0.05 CTAS genes',SKY),('q<0.01','q<0.01 CTAS genes',BLUE),('q<0.001','q<0.001 CTAS genes',DARK),('background_q>=0.05','Background genes (q ≥ 0.05)',LGREY)]
EXPF={'q<0.05':(109,223,49),'q<0.01':(56,99,57),'q<0.001':(23,38,61),'background_q>=0.05':(1154,18701,6)}
# F5b
fig,ax=plt.subplots(figsize=(7.087,2.9)); fig.subplots_adjust(left=0.10,right=0.99,bottom=0.16,top=0.93)
for i,(k,lab,c_) in enumerate(items):
    a_,n_=int(ts.at[k,'supported_ge3']),int(ts.at[k,'n_evaluable']); v=pct(a_,n_); lo,hi=ts.at[k,'pct_ge3_CI95_low'],ts.at[k,'pct_ge3_CI95_high']
    ax.bar(i,v,0.6,color=c_); ax.errorbar(i,v,yerr=[[v-lo],[hi-v]],color=INK,capsize=3,lw=0.8)
    ax.text(i+0.34,v,f'{v:.0f}%\n{a_:,}/{n_:,}',ha='left',va='center',fontsize=7.5)
    num('F5','b',f'FinnGen family support {k}',f'{v:.0f}% ({a_:,}/{n_:,})','FAMILY_TIER_SUMMARY.tsv','supported_ge3/n_evaluable',f'{a_}/{n_}={round(v)}',f'{EXPF[k][0]}/{EXPF[k][1]}={EXPF[k][2]}')
ax.set_xticks(range(4)); ax.set_xticklabels([l for _,l,_ in items],fontsize=7.5); ax.set_ylabel('Genes with FinnGen support (%)'); ax.set_ylim(0,100); ax.set_xlim(-0.5,3.85)
ax.text(0.01,0.97,'Support: gene-level P<0.01 in ≥3 of 6 clinical families of FinnGen endpoints\n(error bars: 95% Wilson CI; increases across the nested CTAS sets)',transform=ax.transAxes,ha='left',va='top',fontsize=7)
src('F5','b','gene-level FinnGen family support (>=3 of 6 families) for nested CTAS gene sets and background','FAMILY_TIER_SUMMARY.tsv','supported_ge3 / n_evaluable, Wilson CI','rebuilt (family measure)')
save(fig,'F5b',ts.reset_index())
# SF5b
fig,ax=plt.subplots(figsize=(3.45,3.0)); fig.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.83)
for j,(k,lab,c) in enumerate([(a,b,c) for a,b,c in items]):
    for i,thr in enumerate([1,2,3]):
        v=ts.at[k,f'pct_ge{thr}']; ax.bar(i+(j-1.5)*0.2,v,0.19,color=c,label=(k.replace('background_q>=0.05','background (q ≥ 0.05)')+f' (n = {int(ts.at[k,"n_evaluable"]):,})') if i==0 else None)
        if thr==3: ax.text(i+(j-1.5)*0.2,v+1,f'{v:.0f}',ha='center',fontsize=6.3)
        num('SF5','b',f'{k} support >= {thr} families',f'{v:.0f}%' if thr==3 else '(bar)','FAMILY_TIER_SUMMARY.tsv',f'supported_ge{thr}/n_evaluable',f'{int(ts.at[k,f"supported_ge{thr}"])}/{int(ts.at[k,"n_evaluable"])}',{3:{'q<0.05':'109/223','q<0.01':'56/99','q<0.001':'23/38','background_q>=0.05':'1154/18701'}[k]}.get(thr))
ax.set_xticks(range(3)); ax.set_xticklabels(['≥1','≥2','≥3']); ax.set_xlabel('FinnGen clinical families (of 6) with gene-level P<0.01'); ax.set_ylabel('CTAS genes supported (%)'); ax.set_ylim(0,108); ax.set_yticks(range(0,101,20))
ax.legend(frameon=False,fontsize=6.3,loc='lower left',bbox_to_anchor=(-0.02,1.0),ncol=2,handlelength=1,columnspacing=0.8)
src('SF5','b','gene-level FinnGen family support at >=1/2/3 families','FAMILY_TIER_SUMMARY.tsv','pct_ge1/2/3','rebuilt (family measure)')
save(fig,'SF5b',ts.reset_index())
# SF5c
d5=dd[dd.set=='q<0.05']
fig,ax=plt.subplots(figsize=(3.45,3.0)); fig.subplots_adjust(left=0.15,bottom=0.15,right=0.98,top=0.86)
ax.bar(d5.n_families,d5.n_genes,color=[DARK if k>=3 else SKY for k in d5.n_families],width=0.75)
for r in d5.itertuples():
    if r.n_genes: ax.text(r.n_families,r.n_genes+1,str(r.n_genes),ha='center',fontsize=6.5)
ax.set_xticks(range(7)); ax.set_xlabel('Supported FinnGen clinical families (of 6)'); ax.set_ylabel('q<0.05 CTAS genes')
n_ev=int(d5.n_genes.sum()); ge3=int(d5[d5.n_families>=3].n_genes.sum())
ax.set_title(f'{n_ev} evaluable q<0.05 genes (all present in FinnGen)\n{ge3} with ≥3 families (dark)',fontsize=7,pad=3)
num('SF5','c','evaluable q<0.05 genes',n_ev,'FAMILY_DISTRIBUTION.tsv','sum n_genes (q<0.05)',n_ev,223); num('SF5','c','genes with >=3 families',ge3,'FAMILY_DISTRIBUTION.tsv','sum n_genes k>=3',ge3,109)
for r in d5.itertuples(): num('SF5','c',f'genes with {r.n_families} families',r.n_genes,'FAMILY_DISTRIBUTION.tsv','n_genes',int(r.n_genes))
src('SF5','c','distribution of supported-family counts among evaluable q<0.05 genes','FAMILY_DISTRIBUTION.tsv','set == q<0.05','rebuilt (family measure)')
save(fig,'SF5c',d5)
dump('f5sf5_families'); print('done',len(LEDGER),sum(r['status']=='FAIL' for r in LEDGER))
