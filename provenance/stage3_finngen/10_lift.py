# Strict, tab-explicit GRCh38->GRCh37 conversion of one FinnGen R12 endpoint. Usage: 10_lift.py <EP> <FILECODE>
import sys,subprocess,bisect,collections,gzip,json,time
EP,CODE=sys.argv[1],sys.argv[2]
R='<CTAS_WORKSPACE>/'; OUT=R+'target_pipeline/finngen_corrected_rerun_01/'
RAW=R+'CTAS_FinnGen/data/finngen_gwas/finngen_R12_%s.gz'%CODE
HDR=['#chrom','pos','ref','alt','rsids','nearest_genes','pval','mlogp','beta','sebeta','af_alt','af_alt_cases','af_alt_controls']
blocks=collections.defaultdict(list); cur=None
for l in open('<REFERENCE_DIR>/reference/hg38ToHg19.over.chain'):
    f=l.split()
    if not f: continue
    if f[0]=='chain':
        score=int(f[1]); tn,qn=f[2],f[7]
        cur=None if tn!=qn or f[9]!='+' else [tn,int(f[5]),int(f[10]),score]; continue
    if cur is None: continue
    sz=int(f[0]); tn,t,q,sc=cur; blocks[tn].append((t,t+sz,q,sc))
    if len(f)==3: cur[1]=t+sz+int(f[1]); cur[2]=q+sz+int(f[2])
starts={}
for k in blocks: blocks[k].sort(); starts[k]=[b[0] for b in blocks[k]]
cn={'23':'chrX','24':'chrY','25':'chrM','MT':'chrM'}
def lift(c,p):
    k=cn.get(c,'chr'+c)
    if k not in starts: return None
    i=bisect.bisect_right(starts[k],p-1)-1; best=None
    for j in range(i,max(i-8,-1),-1):
        b=blocks[k][j]
        if b[0]<=p-1<b[1] and (best is None or b[3]>best[3]): best=b
    return None if best is None else best[2]+(p-1-best[0])+1
st=collections.Counter(); t0=time.time()
p=subprocess.Popen(['gzip','-dc',RAW],stdout=subprocess.PIPE,bufsize=1<<20)
fo=gzip.open(OUT+'input_lifted/%s.lifted.tsv.gz'%EP,'wt',compresslevel=1)
fm=open(OUT+'input_lifted/%s.magma37.txt'%EP,'w')
first=True; fo.write('SNP\tCHR\tBP37\tBP38\tREF\tALT\tP\n')
for raw in p.stdout:
    l=raw.decode().rstrip('\n'); f=l.split('\t')
    if first:
        first=False
        if f!=HDR: print('STOP: header mismatch',f); sys.exit(2)
        continue
    st['rows']+=1
    if len(f)!=13: print('STOP: row with %d fields at row %d'%(len(f),st['rows'])); sys.exit(3)
    rs=f[4]
    if rs in('','.'): st['empty_rsid']+=1; continue
    try: pv=float(f[6])
    except: st['bad_P']+=1; continue
    if not(pv>0 and pv<=1): st['bad_P']+=1; continue
    if pv<1e-300: pv=1e-300
    st['retained']+=1
    if ',' in rs: st['multi_rsid']+=1
    try: pos=int(f[1])
    except: st['bad_pos']+=1; continue
    q=lift(f[0],pos)
    if q is None: st['lift_failed']+=1; continue
    st['lift_ok']+=1
    fo.write('%s\t%s\t%d\t%d\t%s\t%s\t%.10e\n'%(rs,f[0],q,pos,f[2],f[3],pv))
    fm.write('%s\t%s\t%d\t%.10e\n'%(rs,f[0],q,pv))
fo.close(); fm.close()
st['seconds']=int(time.time()-t0); st['endpoint']=EP; st['file']=CODE
json.dump(dict(st),open(OUT+'work/lift_stats_%s.json'%EP,'w'))
print(EP,dict(st))
