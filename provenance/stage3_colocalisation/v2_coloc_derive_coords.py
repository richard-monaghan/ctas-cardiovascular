#!/usr/bin/env python3
"""
Derive hg38 gene-body coords + ENSG for all 104 V2 targets.
 - 48 cached (unchanged) targets: reuse V1 all_101_gene_coords.tsv (same ENSG that
   produced the cached slices).
 - 56 uncached targets: fresh Ensembl GRCh38 symbol lookup.
Writes v2_target_coords.tsv and logs any lookup failures. Resumable: existing rows
with a valid ENSG are kept; only missing genes are (re)looked up.
"""
import csv, json, os, sys, time, urllib.request, urllib.error

BASE   = os.path.expanduser("<CTAS_WORKSPACE>")
PAPER  = "<CTAS_V2_PAPER_WORKSPACE>"
COLOC  = os.path.join(PAPER, "_source_assets", "coloc_review")
MANI   = os.path.join(COLOC, "V2_104_COLOC_MANIFEST.csv")
V1CO   = os.path.join(BASE, "eqtl_vignette", "coords", "all_101_gene_coords.tsv")
OUT    = os.path.join(COLOC, "v2_target_coords.tsv")
LOG    = os.path.join(BASE, "eqtl_vignette", "logs_v2", "coord_lookup_failures.log")
PAD    = 1_000_000
MHC_LO, MHC_HI = 28_500_000, 33_500_000

def ensembl_lookup(gene, retries=5):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?expand=0"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for a in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.load(r)
            return d.get("id",""), str(d.get("seq_region_name","")), d.get("start"), d.get("end")
        except Exception as e:
            if a == retries-1:
                return None
            time.sleep(2*(a+1))
    return None

def main():
    rows = list(csv.DictReader(open(MANI)))
    v1 = {}
    # V1 coords columns: gene(0) ensg(1) chr(2) start(3) end(4) wstart(5) wend(6)
    for r in csv.reader(open(V1CO), delimiter="\t"):
        if r and r[0] != "gene":
            v1[r[0]] = dict(ensg=r[1], chr=r[2], gstart=int(r[3]), gend=int(r[4]))
    # resume: keep already-derived rows
    have = {}
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT), delimiter="\t"):
            if r.get("ensg","").startswith("ENSG"):
                have[r["canonical_locus_id"]] = r
    out_rows, fails = [], []
    for r in rows:
        lid  = r["canonical_locus_id"]
        gene = r["statistical_representative_gene"]
        cached = r["n_cached_valid_slices"] == "7"
        if lid in have:
            out_rows.append(have[lid]); continue
        if cached and gene in v1:
            c = v1[gene]; ensg, chrom, gs, ge, src = c["ensg"], c["chr"], c["gstart"], c["gend"], "v1_coords"
        else:
            res = ensembl_lookup(gene)
            if res is None or not res[0] or not res[2] or not res[3] or not res[1].isdigit():
                fails.append(f"{lid}\t{gene}\tENSEMBL_LOOKUP_FAILED")
                print(f"  FAIL {lid} {gene}", file=sys.stderr); continue
            ensg, chrom, gs, ge, src = res[0], res[1], int(res[2]), int(res[3]), "ensembl_hg38"
        ws = max(1, gs - PAD); we = ge + PAD
        is_mhc = (chrom == "6" and gs <= MHC_HI and ge >= MHC_LO)
        out_rows.append(dict(canonical_locus_id=lid, gene=gene, ensg=ensg, chr=chrom,
                             gstart=gs, gend=ge, wstart=ws, wend=we,
                             cached=str(cached), coord_source=src, is_mhc=str(is_mhc)))
        print(f"  {lid} {gene:12s} chr{chrom}:{ws}-{we} {ensg} [{src}]", flush=True)
    cols = ["canonical_locus_id","gene","ensg","chr","gstart","gend","wstart","wend","cached","coord_source","is_mhc"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader()
        for r in out_rows: w.writerow({k: r[k] for k in cols})
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    open(LOG, "w").write("\n".join(fails) + ("\n" if fails else ""))
    print(f"\nDerived {len(out_rows)}/104 targets; {len(fails)} failures -> {OUT}")
    if fails:
        print("FAILURES:\n" + "\n".join(fails)); sys.exit(2)

if __name__ == "__main__":
    main()
