# CTAS V2 release copy of the 16 Nov 2026 refetch script (paths configurable). To fetch ALL 728 slices, run it once and then with
# EXTRA_PAIRS for the remaining targets, or adapt the target filter; expected checksums are in inputs/colocalisation/EQTL_SLICE_MANIFEST.tsv.
#!/usr/bin/env python3
"""P3(b) refetch (16 Nov 2026): fetch the public eQTL Catalogue (GTEx v8, study QTS000015) cis-window slices needed to
re-execute the CTAS V2 coloc run. Same URLs, regions and tabix call as v2_coloc_fetch.py / fetch_data_all_101.sh.
 - REQUIRED: the 49 'cached' targets whose original slices (eqtl_vignette/gtex_slices/) no longer exist (7 datasets each)
 - CHECK:    the 7 slices of one target whose 7 Aug 2026 slices still exist, to detect changes in the external data
Writes only into this folder; records URL, dataset id, region, UTC access time, tabix return code, rows and SHA-256."""
import csv, os, subprocess, hashlib, datetime, sys
HERE = os.environ.get("CTAS_EQTL_FETCH_DIR", "work/eqtl_fetch"); OUT = os.path.join(HERE, "slices"); LOG = os.path.join(HERE, "logs", "fetch_manifest.tsv")
os.makedirs(OUT, exist_ok=True); os.makedirs(os.path.dirname(LOG), exist_ok=True)
COORDS = os.path.join(os.environ.get("CTAS_REPO", "."), "inputs/colocalisation/v2_target_coords.tsv")
TISSUES = [("heart_LV","QTD000256"),("heart_AA","QTD000251"),("aorta","QTD000131"),("coronary","QTD000136"),("blood","QTD000356"),("liver","QTD000266"),("kidney_cortex","QTD000261")]
URL = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{q}/{q}.all.tsv.gz"
CHECK_TARGET = sys.argv[1] if len(sys.argv) > 1 else "ACAD10"
rows = list(csv.DictReader(open(COORDS), delimiter="\t"))
todo = [(r, "required_missing_cached_slice") for r in rows if r["cached"] == "True"] + [(r, "external_change_check") for r in rows if r["gene"] == CHECK_TARGET]
# extra target x tissue pairs missing from the surviving Aug 2026 folder (env EXTRA_PAIRS="GENE:tissue,..."; recorded as purpose required_missing_aug_slice)
EXTRA = [x.split(':') for x in os.environ.get("EXTRA_PAIRS", "").split(',') if x]
if EXTRA: todo = [(r, "required_missing_aug_slice:" + t) for g, t in EXTRA for r in rows if r["gene"] == g]
done = set()
if os.path.exists(LOG):
    for r in csv.DictReader(open(LOG), delimiter="\t"):
        if r["status"] == "ok": done.add(r["file"])
new = not os.path.exists(LOG); lf = open(LOG, "a", newline=""); w = csv.writer(lf, delimiter="\t")
if new: w.writerow(["file","purpose","gene","canonical_locus_id","tissue","dataset_id","url","region","accessed_utc","tabix_rc","status","n_rows","n_target_ensg_rows","sha256","bytes"])
import threading, concurrent.futures as cf
lock = threading.Lock(); jobs = []
for r, purpose in todo:
    region = f'{r["chr"]}:{r["wstart"]}-{r["wend"]}'
    for label, q in TISSUES:
        if purpose.startswith("required_missing_aug_slice:") and label != purpose.split(":")[1]: continue
        fn = f'{r["gene"]}__{label}.tsv' if purpose.startswith("required") else f'CHECK_{r["gene"]}__{label}.tsv'
        if fn not in done: jobs.append((r, purpose.split(":")[0], region, label, q, fn))
def work(j):
    r, purpose, region, label, q, fn = j; out = os.path.join(OUT, fn); url = URL.format(q=q)
    for attempt in range(3):
        t = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        with open(out, "wb") as fh: rc = subprocess.run(["tabix", url, region], stdout=fh, stderr=subprocess.DEVNULL).returncode
        if rc == 0: break
    b = open(out, "rb").read(); n = b.count(b"\n"); ne = sum(1 for l in b.split(b"\n") if l and l.split(b"\t")[16:17] == [r["ensg"].encode()])
    st = "ok" if rc == 0 else "fetch_failed"
    if rc != 0: os.remove(out)
    with lock:
        w.writerow([fn, purpose, r["gene"], r["canonical_locus_id"], label, q, url, region, t, rc, st, n, ne, hashlib.sha256(b).hexdigest() if rc == 0 else "", len(b) if rc == 0 else 0]); lf.flush()
        print(st, fn, n, ne, flush=True)
with cf.ThreadPoolExecutor(max_workers=int(os.environ.get("FETCH_WORKERS", "6"))) as ex: list(ex.map(work, jobs))
print("FETCH COMPLETE")
