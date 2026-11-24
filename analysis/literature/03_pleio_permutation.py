# CTAS V2 public release copy of four_comparator_harmonisation_20260919/scripts/03_pleio_permutation.py. Changes: paths only.
#!/usr/bin/env python3
"""Phase 5: PLEIO enrichment null at BOTH hit level (625 index variants) and the fairer
region level (1Mb-merged PLEIO regions). Chromosome- and size-preserving. Fixed seed."""
import csv, os, random
ROOT=os.environ["CTAS_REPO"]; HERE=os.environ["CTAS_LIT_WORK"]   # release: paths from run_literature.sh
LOCKM=os.path.join(ROOT,"inputs/historical_prefinal/prefinal_locus_master.csv"); LOCKG=os.path.join(ROOT,"inputs/historical_prefinal/prefinal_gene_to_locus.csv")
SEED=20260919; NPERM=50000
HG19={1:249250621,2:243199373,3:198022430,4:191154276,5:180915260,6:171115067,7:159138663,
8:146364022,9:141213431,10:135534747,11:135006516,12:133851895,13:115169878,14:107349540,
15:102531392,16:90354753,17:81195210,18:78077248,19:59128983,20:63025520,21:48129895,22:51304566}
loci=[]
with open(LOCKM) as f:
    for r in csv.DictReader(f):
        loci.append((int(r["chromosome"]),int(r["locus_start"]),int(r["locus_end"])))
by_chr={}
for c,a,b in loci: by_chr.setdefault(c,[]).append((a,b))
def hit_in(c,p):
    for a,b in by_chr.get(c,[]):
        if a<=p<=b: return True
    return False
def region_hits(c,a,b):
    for x,y in by_chr.get(c,[]):
        if a<=y and x<=b: return True
    return False

pleio=list(csv.DictReader(open(os.path.join(HERE,"pleio_clean_hg19.tsv")),delimiter="\t"))
hits=[(int(p["chr"]),int(p["pos_hg19"])) for p in pleio]
# regions: 1Mb merge
MERGE=1_000_000; regions=[]
for c in sorted({h[0] for h in hits}):
    ps=sorted(p for cc,p in hits if cc==c); s=e=None
    for p in ps:
        if s is None: s=e=p
        elif p-e<=MERGE: e=p
        else: regions.append((c,s,e)); s=e=p
    if s is not None: regions.append((c,s,e))

obs_hit=sum(1 for c,p in hits if hit_in(c,p))
obs_reg=sum(1 for c,a,b in regions if region_hits(c,a,b))
# per-chr locus coverage for analytic expectation
cover={c:sum(b-a for a,b in by_chr.get(c,[])) for c in HG19}
exp_hit=sum(cover[c]/HG19[c] for c,p in hits if c in HG19)

random.seed(SEED)
ge_hit=0; hmean=0
for _ in range(NPERM):
    cnt=0
    for c,p in hits:
        if c in HG19 and hit_in(c,random.randint(1,HG19[c])): cnt+=1
    hmean+=cnt
    if cnt>=obs_hit: ge_hit+=1
hmean/=NPERM; p_hit=(ge_hit+1)/(NPERM+1)

random.seed(SEED)
ge_reg=0; rmean=0
for _ in range(NPERM):
    cnt=0
    for c,a,b in regions:
        if c not in HG19: continue
        L=b-a; start=random.randint(1,max(1,HG19[c]-L));
        if region_hits(c,start,start+L): cnt+=1
    rmean+=cnt
    if cnt>=obs_reg: ge_reg+=1
rmean/=NPERM; p_reg=(ge_reg+1)/(NPERM+1)

rows=[["level","unit","n_units","observed_overlaps","expected_mean_null","empirical_p_ge_obs","enrichment_obs_over_exp","seed","n_perm","note"],
["hit","PLEIO index variant",len(hits),obs_hit,round(hmean,2),p_hit,round(obs_hit/hmean,2),SEED,NPERM,
 "matches prior 2026-09-19 audit; hits cluster at genes so uniform null is coarse/conservative"],
["region","1Mb-merged PLEIO region",len(regions),obs_reg,round(rmean,2),p_reg,round(obs_reg/rmean,2),SEED,NPERM,
 "fairer locus/region-unit null; preserves region size and chromosome"]]
with open(os.path.join(HERE,"pleio_permutation_results.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(rows)
for r in rows: print(r)
