#!/bin/bash
###############################################################################
# CVD-CTAS effective-N sensitivity pipeline  (bash 3.2 compatible)
#
# Changes ONLY the MAGMA sample-size argument N. Same sumstats, same reference
# panel, same gene annotation, same MAGMA model (snpwise-mean default).
# Everything is written under magma_effective_n/. The circulated outputs in
# magma_cvd_results/ are read-only and never modified.
#
# Original per-trait invocation (from magma_cvd_results/<trait>.log), verbatim:
#   MAGMA_Mac/magma --bfile <g1000_eur> \
#     --pval <trait>_rsids.magma use=1,4 N=500000 \
#     --gene-annot <trait>_annot.genes.annot --out <trait>
###############################################################################
set -u

BASE=${CTAS_WORKSPACE:?set CTAS_WORKSPACE}
WORK=$BASE/magma_effective_n
MAGMA=$BASE/MAGMA_Mac/magma
REF=${REFERENCE_DIR:?set REFERENCE_DIR}/reference/g1000_eur
BIM=$REF.bim
ORIG=$BASE/magma_cvd_results          # immutable, read-only
FMT=$BASE/gwas_magma_formatted        # icd10-*.magma.gz (SNP CHR POS PVAL)

REGEN=$WORK/rsids_regenerated
N500=$WORK/reruns_n500k
NEFF=$WORK/reruns_neff
LOG=$WORK/logs
mkdir -p "$REGEN" "$N500" "$NEFF" "$LOG"

# trait  N_eff(int)  source     icd10file(or - for existing)
# source: EXIST = use ORIG/<trait>_rsids.magma ; REGEN = regenerate from icd10 file
RECORDS='
CAD_I25            121878  EXIST  -
AF_I48             74709   EXIST  -
HF_I50             30642   EXIST  -
Hypertension_401   309156  EXIST  -
MI                 41350   REGEN  icd10-I21-both_sexes
Angina             85369   REGEN  icd10-I20-both_sexes
Ischemic_Stroke    17567   REGEN  icd10-I63-both_sexes
TIA                12736   REGEN  icd10-G45-both_sexes
VTE                2054    REGEN  icd10-I82-both_sexes
Cardiomyopathy     6622    REGEN  icd10-I42-both_sexes
Aortic_Valve       14091   REGEN  icd10-I35-both_sexes
Mitral_Valve       12768   REGEN  icd10-I34-both_sexes
PAD                19497   REGEN  icd10-I73-both_sexes
Aortic_Aneurysm    8085    REGEN  icd10-I71-both_sexes
Atherosclerosis    3724    REGEN  icd10-I70-both_sexes
DVT                19325   REGEN  icd10-I80-both_sexes
'

ts () { date '+%Y-%m-%d %H:%M:%S'; }

rsids_path () {  # trait source
  if [ "$2" = "EXIST" ]; then echo "$ORIG/${1}_rsids.magma"
  else echo "$REGEN/${1}_rsids.magma"; fi
}

run_magma () {   # trait  N  outdir  source
  local t=$1 N=$2 od=$3 src=$4
  local pv annot
  pv=$(rsids_path "$t" "$src")
  annot=$ORIG/${t}_annot.genes.annot
  "$MAGMA" --bfile "$REF" \
    --pval "$pv" use=1,4 N="$N" \
    --gene-annot "$annot" \
    --out "$od/${t}" > "$LOG/${t}_N${N}.stdout" 2>&1
}

echo "=== [$(ts)] PHASE 1: build chr:pos -> rsID lookup from reference bim (once) ==="
LOOKUP=$WORK/tmp_pos_to_rsid.txt
if [ ! -s "$LOOKUP" ]; then
  awk '{print $1":"$4, $2}' "$BIM" > "$LOOKUP"
fi
echo "    lookup lines: $(wc -l < "$LOOKUP")"

echo "=== [$(ts)] PHASE 2: regenerate rsids for the 12 icd10 traits ==="
echo "$RECORDS" | while read -r t N src icd; do
  [ -z "$t" ] && continue
  [ "$src" = "EXIST" ] && continue
  out=$REGEN/${t}_rsids.magma
  if [ -s "$out" ]; then echo "    $t: exists, skip"; continue; fi
  echo "    [$(ts)] $t <- $icd"
  gunzip -c "$FMT/${icd}.magma.gz" | awk 'BEGIN{OFS="\t"}
    NR==FNR { lookup[$1]=$2; next }
    { pk=$2":"$3; if (pk in lookup) print lookup[pk], $2, $3, $4 }' "$LOOKUP" - > "$out"
  echo "        -> $(wc -l < "$out") SNPs with rsIDs"
done

echo "=== [$(ts)] PHASE 3: SANITY GATE — reproduce CAD & TIA at N=500000 ==="
run_magma CAD_I25 500000 "$N500" EXIST
run_magma TIA     500000 "$N500" REGEN
echo "    [$(ts)] sanity reruns written to $N500 (compared by separate script)."

echo "=== [$(ts)] PHASE 4: run all 16 traits at effective N (batches of 6) ==="
n=0
echo "$RECORDS" | { while read -r t N src icd; do
  [ -z "$t" ] && continue
  if [ -s "$NEFF/${t}.genes.out" ]; then echo "    $t: exists, skip"; continue; fi
  echo "    [$(ts)] launch $t (N=$N, $src)"
  run_magma "$t" "$N" "$NEFF" "$src" &
  n=$((n+1))
  if [ $((n % 6)) -eq 0 ]; then wait; fi
done; wait; }

echo "=== [$(ts)] PIPELINE COMPLETE ==="
echo "neff genes.out files: $(ls -1 "$NEFF"/*.genes.out 2>/dev/null | wc -l)"
