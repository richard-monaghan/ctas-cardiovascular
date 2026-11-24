#!/usr/bin/env python3
"""
Build effective_n_table.csv for the 16 CVD-CTAS traits.

Design ruling #3: N_eff = 4 / (1/N_cases + 1/N_controls), using the case/control
counts for the SAME ancestry/definition the GWAS used.

Established empirically in Step 0: the p-values fed to MAGMA are Pan-UKBB
`neglog10_pval_meta` (the all-ancestry meta-analysis), NOT a single-ancestry
column. Therefore the matching sample size is the META sample size = the sum of
per-ancestry case/control counts over the ancestries included in the meta
(`pops` column of phenotype_manifest.tsv).

We also emit EUR-only N_eff as a transparency/sensitivity column, since the LD
reference panel is g1000_eur.
"""
import csv, sys

MANIFEST = "<CTAS_WORKSPACE>/phenotype_manifest.tsv"
OUT      = "<CTAS_WORKSPACE>/magma_effective_n/effective_n_table.csv"

# trait_name -> (trait_type, phenocode)  [pheno_sex = both_sexes for all]
# Mapping taken verbatim from the original driver scripts:
#   scripts/run_magma_cvd.sh (12 icd10 traits)
#   scripts/rerun_4_existing_cvd_traits.sh (4 "existing" traits)
TRAITS = [
    ("CAD_I25",          "icd10",  "I25"),
    ("AF_I48",           "icd10",  "I48"),
    ("HF_I50",           "icd10",  "I50"),
    ("Hypertension_401", "phecode","401"),
    ("MI",               "icd10",  "I21"),
    ("Angina",           "icd10",  "I20"),
    ("Ischemic_Stroke",  "icd10",  "I63"),
    ("TIA",              "icd10",  "G45"),
    ("VTE",              "icd10",  "I82"),
    ("Cardiomyopathy",   "icd10",  "I42"),
    ("Aortic_Valve",     "icd10",  "I35"),
    ("Mitral_Valve",     "icd10",  "I34"),
    ("PAD",              "icd10",  "I73"),
    ("Aortic_Aneurysm",  "icd10",  "I71"),
    ("Atherosclerosis",  "icd10",  "I70"),
    ("DVT",              "icd10",  "I80"),
]

POPS_ALL = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]

def fnum(x):
    if x is None or x == "" or x == "NA":
        return None
    try:
        return float(x)
    except ValueError:
        return None

# Index manifest rows by (trait_type, phenocode, pheno_sex)
rows = {}
with open(MANIFEST, newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        key = (row["trait_type"], row["phenocode"], row["pheno_sex"])
        rows[key] = row

def neff(nc, nk):
    if not nc or not nk:
        return None
    return 4.0 / (1.0/nc + 1.0/nk)

out_rows = []
flags = []
for name, ttype, pheno in TRAITS:
    key = (ttype, pheno, "both_sexes")
    m = rows.get(key)
    if m is None:
        flags.append(f"[FLAG] {name}: no manifest row for {key}")
        out_rows.append([name, ttype, pheno, "", "", "", "", 500000, "", "", ""])
        continue

    pops = (m.get("pops") or "").split(",")
    pops = [p for p in pops if p]

    # META counts = sum over the pops actually included in the meta-analysis
    nc_meta = 0.0; nk_meta = 0.0; missing = []
    for p in pops:
        c = fnum(m.get(f"n_cases_{p}"))
        k = fnum(m.get(f"n_controls_{p}"))
        if c is None or k is None:
            missing.append(p)
            continue
        nc_meta += c; nk_meta += k
    if missing:
        flags.append(f"[FLAG] {name}: missing per-pop counts for pops {missing} (pops={pops})")

    # EUR-only (transparency / LD-panel-matched alternative)
    nc_eur = fnum(m.get("n_cases_EUR"))
    nk_eur = fnum(m.get("n_controls_EUR"))

    neff_meta = neff(nc_meta, nk_meta)
    neff_eur  = neff(nc_eur, nk_eur)
    ratio = (neff_meta/500000.0) if neff_meta else ""

    out_rows.append([
        name, ttype, pheno, ",".join(pops),
        int(nc_meta) if nc_meta else "", int(nk_meta) if nk_meta else "",
        round(neff_meta,1) if neff_meta else "",
        500000,
        round(ratio,4) if ratio != "" else "",
        int(nc_eur) if nc_eur else "", int(nk_eur) if nk_eur else "",
        round(neff_eur,1) if neff_eur else "",
    ])

header = ["trait","trait_type","phenocode","pops_meta",
          "N_cases_meta","N_controls_meta","N_eff_meta",
          "flat_N","ratio_meta_over_flat",
          "N_cases_EUR","N_controls_EUR","N_eff_EUR"]

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(header)
    w.writerows(out_rows)

# Console summary
print(f"{'trait':18s} {'pops':16s} {'Ncase':>8s} {'Nctrl':>9s} {'N_eff':>10s} {'ratio':>7s} {'N_eff_EUR':>10s}")
for row in out_rows:
    print(f"{row[0]:18s} {row[3]:16s} {str(row[4]):>8s} {str(row[5]):>9s} {str(row[6]):>10s} {str(row[8]):>7s} {str(row[11]):>10s}")
print()
if flags:
    print("FLAGS:")
    for f in flags: print(" ", f)
else:
    print("No flags: all 16 traits resolved unambiguously against the manifest.")
print(f"\nWritten: {OUT}")
