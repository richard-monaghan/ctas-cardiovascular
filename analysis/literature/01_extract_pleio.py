# CTAS V2 public release copy of pleio_ctas_overlap_audit_20260919/scripts/01_extract_pleio.py. Changes: paths only.
#!/usr/bin/env python3
"""Extract PLEIO Supplementary Table S3 (625 top hits) and S4 (13 novel variants)
from the authoritative Elsevier supplement mmc2.xlsx (AJHG 2021, PMID 33352115).
Coordinates are GRCh37/hg19 (stated in the paper; landmark-verified). AUDIT ONLY."""
import openpyxl, csv, os, hashlib

HERE = os.environ["CTAS_LIT_PLEIO"]   # release: output dir from run_literature.sh
RAW  = os.path.join(HERE, "data_raw", "mmc2.xlsx")

# task-provided annotation list (to be checked against the ACTUAL supplement)
TASK_S4_GENES = ["ACVR2A","ZNF827","LRBA","TNKS1BP1","ERGIC2","HPN","MTTP","RASA3",
                 "FGFBP1","RP11-175E9.1","RPL35AP19","RP11-466I1.1","AC106729.1"]

def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(65536),b''): h.update(b)
    return h.hexdigest()

wb = openpyxl.load_workbook(RAW, read_only=True, data_only=True)

# ---- Table S3: 625 top hits ----
ws = wb["Table S3"]
rows = [list(r) for r in ws.iter_rows(values_only=True)]
# row0 title, row1 header, rest data
data = [r for r in rows[2:] if r[0] not in (None, "")]
s3_path = os.path.join(HERE, "pleio_s3_625_hits_clean.tsv")
with open(s3_path,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["rsID","chr","pos_hg19","A1","A2","pleio_LRT","pleio_P"])
    for r in data:
        w.writerow([r[0], int(r[1]), int(r[2]), r[3], r[4], r[5], r[6]])
n3 = len(data)
uniq_rs = len({r[0] for r in data})
uniq_pos = len({(int(r[1]),int(r[2])) for r in data})

# ---- Table S4: 13 novel ----
ws = wb["Table S4"]
rows = [list(r) for r in ws.iter_rows(values_only=True)]
data4 = [r for r in rows[2:] if r[0] not in (None,"")]
s4_path = os.path.join(HERE, "pleio_s4_13_novel_clean.tsv")
with open(s4_path,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    # S4_HGNC_symbol is the AUTHORITATIVE annotation from the supplement itself.
    # The task-provided gene list is reconciled as an unordered SET (not row-aligned) in step 03.
    w.writerow(["rsID","chr","pos_hg19","A1","A2","pleio_LRT","pleio_P","S4_HGNC_symbol","novelty_status"])
    for r in data4:
        sym = r[7] if len(r)>7 and r[7] not in (None,"") else ""
        w.writerow([r[0], int(r[1]), int(r[2]), r[3], r[4], r[5], r[6], sym, "novel_PLEIO_S4"])
n4 = len(data4)

# reconcile task-provided gene set vs actual S4 symbols (set membership, NOT row order)
s4_syms = {(r[7] if len(r)>7 and r[7] not in (None,"") else None) for r in data4} - {None}
task_set = set(TASK_S4_GENES)
print("S4 non-null symbols (%d):"%len(s4_syms), sorted(s4_syms))
print("task genes in S4 symbols (%d):"%len(task_set & s4_syms), sorted(task_set & s4_syms))
print("task genes NOT in S4 symbol col:", sorted(task_set - s4_syms))
print("S4 symbols NOT in task list:", sorted(s4_syms - task_set))

print("SOURCE mmc2.xlsx md5:", md5(RAW))
print("Table S3: %d data rows | unique rsID=%d | unique (chr,pos)=%d" % (n3, uniq_rs, uniq_pos))
print("Table S4: %d data rows" % n4)
print("wrote:", os.path.basename(s3_path), os.path.basename(s4_path))
# landmark build check
for r in data:
    if r[0]=="rs662799":
        print("landmark rs662799 -> chr%s:%s (expect hg19 chr11:116663707)"%(r[1],r[2]))
