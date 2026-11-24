#!/usr/bin/env python3
"""
Fetch the 392 missing eQTL slices for the 56 uncached V2 targets (7 datasets each)
into eqtl_vignette/gtex_slices_v2/. Cached V1 slices are never touched.

Resumable: a slice with a non-empty file is skipped. Each slice is validated by
checking whether the target ENSG appears in the fetched cis-window; an empty ENSG
result is recorded as 'no_ensg_in_window' (NOT a coloc negative), while a tabix
transport failure is recorded as 'fetch_failed', its partial file removed, so a
rerun retries only genuine failures. Per-slice QC is appended to
logs_v2/v2_slice_fetch_qc.tsv so progress is monitorable live.
"""
import csv, os, subprocess, sys, time

BASE   = os.path.expanduser("<CTAS_WORKSPACE>")
VIG    = os.path.join(BASE, "eqtl_vignette")
OUTDIR = os.path.join(VIG, "gtex_slices_v2")
LOGDIR = os.path.join(VIG, "logs_v2")
COORDS = "<CTAS_V2_PAPER_WORKSPACE>/_source_assets/coloc_review/v2_target_coords.tsv"
QC     = os.path.join(LOGDIR, "v2_slice_fetch_qc.tsv")

TISSUES = [("heart_LV","QTD000256"),("heart_AA","QTD000251"),("aorta","QTD000131"),
           ("coronary","QTD000136"),("blood","QTD000356"),("liver","QTD000266"),
           ("kidney_cortex","QTD000261")]
URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{q}/{q}.all.tsv.gz"
RETRIES = 5

def already_logged():
    done = set()
    if os.path.exists(QC):
        for r in csv.DictReader(open(QC), delimiter="\t"):
            if r["status"] in ("ok","no_ensg_in_window"):
                done.add((r["gene"], r["tissue"]))
    return done

def fetch_one(url, region, out):
    for a in range(RETRIES):
        with open(out, "w") as fh:
            rc = subprocess.run(["tabix", url, region], stdout=fh,
                                stderr=subprocess.DEVNULL).returncode
        if rc == 0 and os.path.getsize(out) >= 0:
            return True
        time.sleep(3*(a+1))
    return False

def main():
    os.makedirs(OUTDIR, exist_ok=True); os.makedirs(LOGDIR, exist_ok=True)
    os.chdir(VIG)  # use local .tbi if present
    rows = [r for r in csv.DictReader(open(COORDS), delimiter="\t") if r["cached"] == "False"]
    done = already_logged()
    new_qc = not os.path.exists(QC)
    qcf = open(QC, "a", newline=""); qw = csv.writer(qcf, delimiter="\t")
    if new_qc:
        qw.writerow(["gene","canonical_locus_id","tissue","qtd","ensg","region",
                     "status","n_rows","n_ensg_rows"]); qcf.flush()
    n_ok=n_noensg=n_fail=n_skip=0
    total = len(rows)*len(TISSUES)
    i=0
    for r in rows:
        gene, lid, ensg = r["gene"], r["canonical_locus_id"], r["ensg"]
        region = f'{r["chr"]}:{r["wstart"]}-{r["wend"]}'
        for label, qtd in TISSUES:
            i+=1
            out = os.path.join(OUTDIR, f"{gene}__{label}.tsv")
            if (gene,label) in done or (os.path.exists(out) and os.path.getsize(out) > 0
                    and (gene,label) in done):
                n_skip+=1; continue
            if os.path.exists(out) and os.path.getsize(out) > 0 and (gene,label) not in done:
                pass  # file exists but unvalidated -> re-validate below without refetch
                ok = True
            else:
                ok = fetch_one(URL.format(q=qtd), region, out)
            if not ok:
                if os.path.exists(out): os.remove(out)
                qw.writerow([gene,lid,label,qtd,ensg,region,"fetch_failed",0,0]); qcf.flush()
                n_fail+=1
                print(f"[{i}/{total}] FAIL {gene}__{label}", flush=True); continue
            n_rows=0; n_ensg=0
            with open(out) as fh:
                for line in fh:
                    n_rows+=1
                    f=line.split("\t")
                    if len(f)>=17 and f[16]==ensg: n_ensg+=1
            status = "ok" if n_ensg>0 else "no_ensg_in_window"
            qw.writerow([gene,lid,label,qtd,ensg,region,status,n_rows,n_ensg]); qcf.flush()
            if status=="ok": n_ok+=1
            else: n_noensg+=1
            time.sleep(1.5)
        print(f"[{i}/{total}] {gene} done (ok={n_ok} noensg={n_noensg} fail={n_fail})", flush=True)
    qcf.close()
    print(f"\nFETCH SUMMARY: ok={n_ok} no_ensg_in_window={n_noensg} fetch_failed={n_fail} skipped={n_skip} / {total}")
    if n_fail>0: sys.exit(3)

if __name__ == "__main__":
    main()
