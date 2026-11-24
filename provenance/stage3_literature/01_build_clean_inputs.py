#!/usr/bin/env python3
"""Phase 1-2: build canonical CTAS universes, HGNC alias map, and clean hg19 comparator inputs
for Qiao / Liu / Foley / PLEIO. All four are GRCh37/hg19 (verified). AUDIT+ST update task.
Nothing canonical is modified; only files under this output folder are written."""
import csv, os, hashlib

ROOT="<CTAS_WORKSPACE>"
HERE=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT=os.path.join(ROOT,"manuscript/literature_overlap")
AUDIT=os.path.join(ROOT,"outputs/pleio_ctas_overlap_audit_20260919")
csv.field_size_limit(10**7)

def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(65536),b''): h.update(b)
    return h.hexdigest()

# ---------------- HGNC alias/prev -> approved map ----------------
hgnc=os.path.join(HERE,"data_raw/hgnc_complete_set.txt")
approved=set(); alias2app={}
with open(hgnc,encoding="utf-8") as f:
    rd=csv.DictReader(f,delimiter="\t")
    for r in rd:
        sym=r["symbol"].strip()
        if not sym: continue
        approved.add(sym.upper())
        for col in ("alias_symbol","prev_symbol"):
            raw=r.get(col) or ""
            for a in raw.replace('"','').split("|"):
                a=a.strip().upper()
                if a and a not in alias2app:
                    alias2app[a]=sym  # first-wins; approved symbols take precedence below
def norm(g):
    if g is None: return ""
    u=g.strip().upper()
    if not u: return ""
    if u in approved: return u                 # already an approved symbol
    if u in alias2app: return alias2app[u].upper()
    return u                                    # unknown -> keep as-is (uppercased)
print("HGNC: approved symbols=%d ; alias/prev mappings=%d"%(len(approved),len(alias2app)))

# ---------------- CTAS canonical universes ----------------
loci=[]
with open(os.path.join(ROOT,"panel_v2_locus_lock/v2_locus_master.csv")) as f:
    for row in csv.DictReader(f):
        members=[m for m in row["all_q05_genes"].split(";") if m]
        loci.append(dict(id=row["canonical_locus_id"],chr=int(row["chromosome"]),
            start=int(row["locus_start"]),end=int(row["locus_end"]),
            rep=row["statistical_representative_gene"],
            q05=row["survives_q05"]=="TRUE",q01=row["survives_q01"]=="TRUE",q001=row["survives_q001"]=="TRUE",
            members=members, members_norm=set(norm(m) for m in members)))
assert len(loci)==104
genes={}
with open(os.path.join(ROOT,"panel_v2_locus_lock/v2_gene_to_locus.csv")) as f:
    for row in csv.DictReader(f):
        genes[row["gene"]]=row["q_tier"]
g05={norm(g) for g in genes}
g01={norm(g) for g,t in genes.items() if t in ("q<0.01","q<0.001")}
g001={norm(g) for g,t in genes.items() if t=="q<0.001"}
assert (len(genes),len(g01),len(g001))==(250,102,32), (len(genes),len(g01),len(g001))
print("CTAS genes 250/102/32 OK ; loci 104/47/16 OK")

# save normalized CTAS gene universe for provenance
with open(os.path.join(HERE,"data_raw/ctas_gene_universe_normalized.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["original_symbol","hgnc_normalized","q_tier"])
    for g,t in genes.items(): w.writerow([g,norm(g),t])

# ---------------- Qiao clean (38 regions, hg19) ----------------
qiao=[]
with open(os.path.join(LIT,"processed/qiao_2025_sd17a_loci.tsv")) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        qiao.append(dict(locus_id=r["locus_id"],chr=int(r["chr"]),start=int(r["start"]),stop=int(r["stop"]),
            candidate_snp=r.get("candidate_snp",""),cytoband=r.get("region_annotation",""),
            gene_annot=r.get("variant_annotation",""),traits=r.get("hyprcoloc_traits",""),
            intersections=r.get("intersections","")))
assert len(qiao)==38
with open(os.path.join(HERE,"qiao_clean_hg19.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["qiao_locus_id","chr","start_hg19","stop_hg19","candidate_snp","cytoband",
                "annotated_gene","mtag_traits","native_unit","source_build","liftover_status"])
    for q in qiao:
        w.writerow([q["locus_id"],q["chr"],q["start"],q["stop"],q["candidate_snp"],q["cytoband"],
                    q["gene_annot"],q["traits"] or q["intersections"],"region","hg19","native_hg19_no_liftover"])

# ---------------- Liu clean (180 genes native; 618 SNPs hg19 for coord check) ----------------
liu_genes=[]
with open(os.path.join(LIT,"processed/liu_2023_pleiotropic_genes.tsv")) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        liu_genes.append(dict(gene=r["gene"],n_traits=r["n_traits"],traits=r["traits"]))
assert len(liu_genes)==180
with open(os.path.join(HERE,"liu_clean_hg19.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["original_symbol","hgnc_normalized_symbol","symbol_changed","n_traits","traits","native_unit","source_build"])
    for g in liu_genes:
        n=norm(g["gene"]); w.writerow([g["gene"],n,("yes" if n!=g["gene"].upper() else "no"),g["n_traits"],g["traits"],"gene","hg19"])
# Liu SNPs (unique)
liu_snps={}
with open(os.path.join(LIT,"processed/liu_2023_pleiotropic_snps.tsv")) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        liu_snps[r["snp"]]=(int(r["chr"]),int(r["pos"]))
with open(os.path.join(HERE,"data_raw/liu_snps_clean_hg19.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["rsID","chr_hg19","pos_hg19"])
    for rs,(c,p) in sorted(liu_snps.items()): w.writerow([rs,c,p])
print("Liu unique SNPs (hg19):",len(liu_snps))

# ---------------- Foley clean (43 regions, coords via candidate SNP rsID -> hg19) ----------------
foley_pos={}
with open(os.path.join(HERE,"data_raw/foley_snp_hg19.tsv")) as f:
    for r in csv.DictReader(f,delimiter="\t"): foley_pos[r["rsID"]]=(int(r["chr_hg19"]),int(r["pos_hg19"]))
foley=[]
with open(os.path.join(LIT,"processed/foley_2021_chd43_regions.tsv")) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        rs=r["candidate_snp"]; c,p=foley_pos.get(rs,("",""))
        genes_here=[x.strip() for x in (r["locus"].replace(";",",").split(",")) if x.strip()]
        genes_here+=[r["nearest_gene"].strip()] if r.get("nearest_gene") else []
        foley.append(dict(chr=r["chr"],locus=r["locus"],traits=r["traits"],candidate_snp=rs,
            snp_chr=c,snp_pos=p,nearest_gene=r.get("nearest_gene",""),
            genes=";".join(sorted(set(genes_here))),genes_norm=";".join(sorted(set(norm(g) for g in genes_here))),
            known=r.get("known_chd_locus",""),ppfc=r.get("ppfc",""),ppe=r.get("ppe","")))
assert len(foley)==43
with open(os.path.join(HERE,"foley_clean_hg19.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["foley_region_genes","traits","candidate_snp","snp_chr_hg19","snp_pos_hg19",
                "nearest_gene","all_region_genes","all_region_genes_hgnc_norm","known_chd_locus","ppfc","ppe",
                "native_unit","coord_source"])
    for x in foley:
        w.writerow([x["locus"],x["traits"],x["candidate_snp"],x["snp_chr"],x["snp_pos"],x["nearest_gene"],
                    x["genes"],x["genes_norm"],x["known"],x["ppfc"],x["ppe"],"region","candidate_snp_rsID_to_hg19_Ensembl_GRCh37"])

# ---------------- PLEIO clean (625 hits + novel flag) ----------------
s3=list(csv.DictReader(open(os.path.join(AUDIT,"pleio_s3_625_hits_clean.tsv")),delimiter="\t"))
s4=list(csv.DictReader(open(os.path.join(AUDIT,"pleio_s4_13_novel_clean.tsv")),delimiter="\t"))
assert len(s3)==625 and len(s4)==13
novel_rs={r["rsID"] for r in s4}
with open(os.path.join(HERE,"pleio_clean_hg19.tsv"),"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["rsID","chr","pos_hg19","A1","A2","pleio_LRT","pleio_P","is_novel_S4","native_unit","source_build"])
    for r in s3:
        w.writerow([r["rsID"],r["chr"],r["pos_hg19"],r["A1"],r["A2"],r["pleio_LRT"],r["pleio_P"],
                    ("yes" if r["rsID"] in novel_rs else "no"),"index_variant","hg19"])

# ---------------- source manifest ----------------
man=[["study","publication","doi_pmid","source_supplement","native_unit","source_build","n_source_units",
      "prior_filtering","prior_ctas_match_rule","this_build_match_rule","source_md5"]]
man.append(["Qiao 2025","Nat Commun 2025 (MTAG+HyPrColoc CVD)","10.1038/s41467-025-62419-0",
  "MOESM19 SD17 section a (38 CVD MTAG loci)","region (chr,start,stop)","GRCh37/hg19 (verified rs880315/rs12129500)","38",
  "SD17a primary CVD-only section (b/c/d cross-disease excluded)","coordinate interval intersection (hg19)",
  "coordinate interval intersection (hg19) [unchanged]", md5(os.path.join(LIT,"qiao_2025/41467_2025_62419_MOESM19_ESM.xlsx"))])
man.append(["Liu 2023","iScience 2023 (PLACO ten CVD)","10.1016/j.isci.2023.106804",
  "mmc2 Table S2 (genes); mmc3 Table S3 (SNPs)","gene (primary) + index SNP","GRCh37/hg19 (verified rs10077410/rs10212121)","180 genes / 618 SNPs",
  "264 pairwise gene records -> 180 unique; SNP 690 records -> 618 unique","gene-symbol overlap vs all_q05_genes (ST10)",
  "HGNC-normalised gene overlap (native) + independent SNP coordinate containment cross-check",
  md5(os.path.join(LIT,"liu_2023_ten_cvd/mmc2.xlsx"))])
man.append(["Foley 2021","Nat Commun 2021 (HyPrColoc CHD)","10.1038/s41467-020-20885-8",
  "MOESM tables -> 43 CHD colocalisation regions (candidate SNP + genes)","region (candidate SNP + gene)","GRCh37/hg19 (rsID-resolved)","43",
  "43 CHD regions (Table1/2); no positions in local extract","gene-symbol overlap vs all_q05_genes (ST10)",
  "candidate-SNP rsID -> hg19 coordinate containment (primary) + HGNC-normalised gene overlap (secondary)",
  md5(os.path.join(LIT,"processed/foley_2021_chd43_regions.tsv"))])
man.append(["Lee PLEIO 2021","Am J Hum Genet 2021","10.1016/j.ajhg.2020.11.017 / PMID 33352115",
  "mmc2.xlsx Table S3 (625 hits), Table S4 (13 novel)","index variant","GRCh37/hg19 (VEP GRCh37; rs662799 verified)","625 hits / 13 novel",
  "S3 approx-independent LD-clumped hits; S4 novel","(new comparator; prior PLEIO audit 2026-09-19)",
  "coordinate containment (hg19) + 1Mb-merged regions + novel-13 coordinate check",
  md5(os.path.join(AUDIT,"data_raw/mmc2.xlsx"))])
with open(os.path.join(HERE,"comparator_source_manifest.tsv"),"w",newline="") as f:
    csv.writer(f,delimiter="\t").writerows(man)

print("clean inputs + manifest written.")
