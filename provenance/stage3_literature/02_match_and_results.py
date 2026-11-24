#!/usr/bin/env python3
"""Phase 3-4: native-unit matching, four-study support matrix, tier summaries,
3-vs-4-study union, and reconciliation vs the live ST10 workbook. hg19 throughout."""
import csv, os, hashlib
ROOT="<CTAS_WORKSPACE>"; HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT=os.path.join(ROOT,"manuscript/literature_overlap")
ST10="<CTAS_V2_PAPER_WORKSPACE>/Supplementary Tables/ST10_V2_literature_overlap.xlsx"
csv.field_size_limit(10**7)
def rd(p): return list(csv.DictReader(open(p),delimiter="\t"))

# ---- HGNC norm (rebuild map) ----
approved=set(); alias2app={}
with open(os.path.join(HERE,"data_raw/hgnc_complete_set.txt"),encoding="utf-8") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        sym=r["symbol"].strip()
        if not sym: continue
        approved.add(sym.upper())
        for col in ("alias_symbol","prev_symbol"):
            for a in (r.get(col) or "").replace('"','').split("|"):
                a=a.strip().upper()
                if a and a not in alias2app: alias2app[a]=sym
def norm(g):
    u=(g or "").strip().upper()
    if not u: return ""
    if u in approved: return u
    return alias2app.get(u,u).upper() if u in alias2app else u

# ---- CTAS loci ----
loci=[]
with open(os.path.join(ROOT,"panel_v2_locus_lock/v2_locus_master.csv")) as f:
    for row in csv.DictReader(f):
        members=[m for m in row["all_q05_genes"].split(";") if m]
        loci.append(dict(id=row["canonical_locus_id"],chr=int(row["chromosome"]),
            start=int(row["locus_start"]),end=int(row["locus_end"]),rep=row["statistical_representative_gene"],
            q05=True,q01=row["survives_q01"]=="TRUE",q001=row["survives_q001"]=="TRUE",
            members=members,members_norm=set(norm(m) for m in members)))
by_chr={}
for L in loci: by_chr.setdefault(L["chr"],[]).append(L)
def contain(chrom,pos):
    for L in by_chr.get(chrom,[]):
        if L["start"]<=pos<=L["end"]: return L
    return None

# ---- comparator inputs ----
qiao=rd(os.path.join(HERE,"qiao_clean_hg19.tsv"))
liu=rd(os.path.join(HERE,"liu_clean_hg19.tsv"))
foley=rd(os.path.join(HERE,"foley_clean_hg19.tsv"))
pleio=rd(os.path.join(HERE,"pleio_clean_hg19.tsv"))
liu_snps=rd(os.path.join(HERE,"data_raw/liu_snps_clean_hg19.tsv"))
liu_genes_norm={r["hgnc_normalized_symbol"] for r in liu}
liu_genes_raw={r["original_symbol"].upper() for r in liu}

# ---- per-locus support ----
for L in loci:
    # Qiao: region interval overlap
    qhit=[q for q in qiao if int(q["chr"])==L["chr"] and int(q["start_hg19"])<=L["end"] and L["start"]<=int(q["stop_hg19"])]
    L["qiao"]=qhit
    # Liu gene (raw ST10-style: raw member vs raw liu symbol)
    L["liu_gene_raw"]=sorted({m for m in L["members"] if m.upper() in liu_genes_raw})
    # Liu gene (HGNC-normalised both sides)
    L["liu_gene_norm"]=sorted({m for m in L["members"] if norm(m) in liu_genes_norm})
    # Liu SNP coordinate containment
    L["liu_snp"]=[s for s in liu_snps if int(s["chr_hg19"])==L["chr"] and L["start"]<=int(s["pos_hg19"])<=L["end"]]
    # Foley gene (ST10-style: normalised member vs normalised foley region genes)
    L["foley_gene"]=[]
    for fo in foley:
        fg={x for x in fo["all_region_genes_hgnc_norm"].split(";") if x}
        if L["members_norm"] & fg: L["foley_gene"].append(fo)
    # Foley coordinate (candidate SNP containment)
    L["foley_coord"]=[fo for fo in foley if fo["snp_chr_hg19"] and int(fo["snp_chr_hg19"])==L["chr"] and L["start"]<=int(fo["snp_pos_hg19"])<=L["end"]]
    # PLEIO coordinate
    L["pleio"]=[p for p in pleio if int(p["chr"])==L["chr"] and L["start"]<=int(p["pos_hg19"])<=L["end"]]

def cnt(pred, tier=None):
    return sum(1 for L in loci if pred(L) and (tier is None or L[tier]))

# reproduce ST10-arm calls: qiao(coord), liu(raw gene), foley(gene)
def st10_qiao(L): return len(L["qiao"])>0
def st10_liu(L):  return len(L["liu_gene_raw"])>0
def st10_foley(L):return len(L["foley_gene"])>0
# refined-arm: qiao(coord), liu(norm gene), foley(coord), pleio(coord)
def r_qiao(L): return len(L["qiao"])>0
def r_liu(L):  return len(L["liu_gene_norm"])>0
def r_foley(L):return len(L["foley_coord"])>0
def r_pleio(L):return len(L["pleio"])>0

# ---- load live ST10 per-locus calls for reconciliation ----
import openpyxl
wb=openpyxl.load_workbook(ST10,data_only=True); ws=wb["locus_overlap_104"]
st10={}
hdr=[c.value for c in ws[2]]
idx={h:i for i,h in enumerate(hdr)}
for row in ws.iter_rows(min_row=3,values_only=True):
    if not row[0]: continue
    st10[row[0]]=dict(in_qiao=row[idx["in_qiao"]],in_liu=row[idx["in_liu"]],in_foley=row[idx["in_foley"]],
        support_combination=row[idx["support_combination"]],n=row[idx["n_supporting_studies"]])

# ============ OUTPUT 1: four-study support matrix ============
mat_rows=[]
for L in loci:
    tier="q<0.001" if L["q001"] else ("q<0.01" if L["q01"] else "q<0.05")
    calls=dict(Qiao=r_qiao(L),Liu=r_liu(L),Foley=r_foley(L),PLEIO=r_pleio(L))
    supp=[k for k in ["Qiao","Liu","Foley","PLEIO"] if calls[k]]
    mat_rows.append(dict(
        canonical_locus_id=L["id"],chromosome=L["chr"],locus_start_hg19=L["start"],locus_end_hg19=L["end"],
        statistical_representative_gene=L["rep"],ctas_q_tier=tier,n_q05_members=len(L["members"]),
        all_q05_genes=";".join(L["members"]),
        Qiao_support=("yes" if calls["Qiao"] else "no"),
        Qiao_evidence=("region_interval_overlap" if calls["Qiao"] else ""),
        Qiao_source=";".join(q["qiao_locus_id"] for q in L["qiao"]),
        Liu_support=("yes" if calls["Liu"] else "no"),
        Liu_evidence=("gene_symbol_HGNC_normalised" if calls["Liu"] else ""),
        Liu_source=";".join(L["liu_gene_norm"]),
        Liu_snp_coord_support=("yes" if L["liu_snp"] else "no"),
        Liu_snp_source=";".join(s["rsID"] for s in L["liu_snp"]),
        Foley_support=("yes" if calls["Foley"] else "no"),
        Foley_evidence=("candidate_SNP_coordinate_containment" if calls["Foley"] else ""),
        Foley_source=";".join(f["candidate_snp"] for f in L["foley_coord"]),
        Foley_gene_support=("yes" if st10_foley(L) else "no"),
        PLEIO_support=("yes" if calls["PLEIO"] else "no"),
        PLEIO_evidence=("index_variant_coordinate_containment" if calls["PLEIO"] else ""),
        PLEIO_n_hits=len(L["pleio"]),
        PLEIO_source=";".join(p["rsID"] for p in L["pleio"]),
        PLEIO_novel_here=";".join(p["rsID"] for p in L["pleio"] if p["is_novel_S4"]=="yes"),
        n_supporting_studies=len(supp),support_combination=("+".join(supp) if supp else "none")))
with open(os.path.join(HERE,"ctas_104_locus_four_study_support_matrix.tsv"),"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(mat_rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(mat_rows)
assert len(mat_rows)==104

# ============ OUTPUT 2: native-unit results ============
def recip(predname):
    return (cnt(predname),cnt(predname,"q01"),cnt(predname,"q001"))
def pct(a,b): return "%.1f%%"%(100*a/b) if b else "0%"
native=[["comparator","native_unit","source_units","harmonised_units","overlap_metric",
         "loci_supported_q05_104","loci_supported_q01_47","loci_supported_q001_16",
         "pct_q05","pct_q01","pct_q001","supported_locus_ids"]]
def ids(pred): return ";".join(L["id"] for L in loci if pred(L))
# Qiao
q05,q01,q001=recip(r_qiao)
native.append(["Qiao","region (38)",38,38,"region∩locus interval overlap (hg19)",q05,q01,q001,
    pct(q05,104),pct(q01,47),pct(q001,16),ids(r_qiao)])
# Liu (norm gene primary)
q05,q01,q001=recip(r_liu)
native.append(["Liu","gene (180)",180,180,"HGNC-normalised gene∈locus membership",q05,q01,q001,
    pct(q05,104),pct(q01,47),pct(q001,16),ids(r_liu)])
# Liu SNP coord (secondary)
q05s,q01s,q001s=recip(lambda L: len(L["liu_snp"])>0)
native.append(["Liu (SNP coord x-check)","index SNP (618)",618,618,"SNP coordinate containment (hg19)",q05s,q01s,q001s,
    pct(q05s,104),pct(q01s,47),pct(q001s,16),ids(lambda L: len(L["liu_snp"])>0)])
# Foley coord primary
q05,q01,q001=recip(r_foley)
native.append(["Foley","region (43)",43,43,"candidate-SNP coordinate containment (hg19)",q05,q01,q001,
    pct(q05,104),pct(q01,47),pct(q001,16),ids(r_foley)])
# Foley gene secondary
q05g,q01g,q001g=recip(st10_foley)
native.append(["Foley (gene x-check)","region genes (43)",43,43,"HGNC-normalised gene overlap",q05g,q01g,q001g,
    pct(q05g,104),pct(q01g,47),pct(q001g,16),ids(st10_foley)])
# PLEIO
q05,q01,q001=recip(r_pleio)
native.append(["PLEIO","index variant (625)",625,625,"index-variant coordinate containment (hg19)",q05,q01,q001,
    pct(q05,104),pct(q01,47),pct(q001,16),ids(r_pleio)])
with open(os.path.join(HERE,"comparator_native_unit_results.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(native)

# ============ OUTPUT 3: tier summary by comparator ============
tiers=[("q05_104",104,lambda L:True),("q01_47",47,lambda L:L["q01"]),("q001_16",16,lambda L:L["q001"])]
tsum=[["qset","n_loci","Qiao","Liu","Foley","PLEIO","any_of_4","none_of_4",
       "any_of_3_QLF","Liu_raw_ST10","Foley_gene_ST10"]]
for name,n,tf in tiers:
    a4=cnt(lambda L: (r_qiao(L) or r_liu(L) or r_foley(L) or r_pleio(L)),None) if False else sum(1 for L in loci if tf(L) and (r_qiao(L) or r_liu(L) or r_foley(L) or r_pleio(L)))
    a3=sum(1 for L in loci if tf(L) and (r_qiao(L) or r_liu(L) or r_foley(L)))
    row=[name,n,
        sum(1 for L in loci if tf(L) and r_qiao(L)),
        sum(1 for L in loci if tf(L) and r_liu(L)),
        sum(1 for L in loci if tf(L) and r_foley(L)),
        sum(1 for L in loci if tf(L) and r_pleio(L)),
        a4, n-a4, a3,
        sum(1 for L in loci if tf(L) and st10_liu(L)),
        sum(1 for L in loci if tf(L) and st10_foley(L))]
    tsum.append(row)
with open(os.path.join(HERE,"tier_summary_by_comparator.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(tsum)

# ============ OUTPUT 4: 3-study vs 4-study union ============
u=[["metric","q05_104","q01_47","q001_16"]]
def union3(L): return r_qiao(L) or r_liu(L) or r_foley(L)
def union4(L): return union3(L) or r_pleio(L)
for label,pred in [("supported_by_3study_union",union3),("supported_by_4study_union",union4),
                   ("gained_by_adding_PLEIO",lambda L: union4(L) and not union3(L))]:
    u.append([label]+[sum(1 for L in loci if tf(L) and pred(L)) for _,_,tf in tiers])
# by number of supporting studies (of 4)
for k in range(5):
    u.append([f"supported_by_exactly_{k}_of_4"]+[sum(1 for L in loci if tf(L) and sum([r_qiao(L),r_liu(L),r_foley(L),r_pleio(L)])==k) for _,_,tf in tiers])
with open(os.path.join(HERE,"three_vs_four_study_union.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(u)
# also list exact loci by n-of-4 at q05
with open(os.path.join(HERE,"data_raw","union_membership_q05.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["n_of_4","locus_id","rep","combination"])
    for L in loci:
        k=sum([r_qiao(L),r_liu(L),r_foley(L),r_pleio(L)])
        combo="+".join([s for s,c in [("Qiao",r_qiao(L)),("Liu",r_liu(L)),("Foley",r_foley(L)),("PLEIO",r_pleio(L))] if c]) or "none"
        w.writerow([k,L["id"],L["rep"],combo])

# ============ OUTPUT 5: old(ST10) vs rebuilt reconciliation + changed calls ============
st10_counts=dict(Qiao=sum(1 for v in st10.values() if v["in_qiao"]=="yes"),
                 Liu=sum(1 for v in st10.values() if v["in_liu"]=="yes"),
                 Foley=sum(1 for v in st10.values() if v["in_foley"]=="yes"))
rebuilt_reproduction=dict(Qiao=cnt(st10_qiao),Liu=cnt(st10_liu),Foley=cnt(st10_foley))  # SAME methods as ST10
rebuilt_refined=dict(Qiao=cnt(r_qiao),Liu=cnt(r_liu),Foley=cnt(r_foley),PLEIO=cnt(r_pleio))
recon=[["study","ST10_reported","rebuilt_same_method","rebuilt_refined_method","reason_for_any_change"]]
recon.append(["Qiao",st10_counts["Qiao"],rebuilt_reproduction["Qiao"],rebuilt_refined["Qiao"],
    "coordinate method unchanged; canonical 104/hg19 already used by ST10"])
recon.append(["Liu",st10_counts["Liu"],rebuilt_reproduction["Liu"],rebuilt_refined["Liu"],
    "refined = HGNC alias normalisation of both sides (may recover alias-mismatched genes)"])
recon.append(["Foley",st10_counts["Foley"],rebuilt_reproduction["Foley"],rebuilt_refined["Foley"],
    "refined = candidate-SNP COORDINATE containment instead of gene-symbol overlap"])
recon.append(["PLEIO","(absent)","(n/a)",rebuilt_refined["PLEIO"],"new fourth comparator (coordinate)"])
with open(os.path.join(HERE,"old_vs_rebuilt_reconciliation.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(recon)

# changed support calls (per locus, ST10 vs refined) for Qiao/Liu/Foley
changed=[["canonical_locus_id","representative","study","ST10_call","rebuilt_same_method","refined_call","reason"]]
for L in loci:
    sid=L["id"]; s=st10.get(sid,{})
    for study,st10call,samem,refm,reason in [
        ("Qiao",s.get("in_qiao"),("yes" if st10_qiao(L) else "no"),("yes" if r_qiao(L) else "no"),"coordinate (unchanged)"),
        ("Liu",s.get("in_liu"),("yes" if st10_liu(L) else "no"),("yes" if r_liu(L) else "no"),"HGNC alias normalisation"),
        ("Foley",s.get("in_foley"),("yes" if st10_foley(L) else "no"),("yes" if r_foley(L) else "no"),"coordinate vs gene-symbol")]:
        if not (st10call==samem==refm):
            changed.append([sid,L["rep"],study,st10call,samem,refm,reason])
with open(os.path.join(HERE,"changed_support_calls.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(changed)

# ---- print key numbers ----
print("ST10 reported (live):",st10_counts)
print("rebuilt SAME method  :",rebuilt_reproduction, "  <-- should reproduce ST10")
print("rebuilt REFINED      :",rebuilt_refined)
print("tier summary:")
for r in tsum: print("  ",r)
print("changed support-call rows (excl header):",len(changed)-1)
print("PLEIO reciprocal q05/q01/q001:",recip(r_pleio))
print("PLEIO loci NOT recovered in q001:",[L["id"]+"/"+L["rep"] for L in loci if L["q001"] and not r_pleio(L)])
