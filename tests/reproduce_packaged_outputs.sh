#!/usr/bin/env bash
# Reruns every analysis that is reproducible from the packaged inputs alone (other than the 10^9 null,
# see run_core_pipeline.sh) and compares the results with the released outputs. ~2 minutes.
# Optional: CTAS_NCBI37_GENE_LOC=<NCBI37.3.gene.loc> also runs the obstructive-panel sensitivity analysis.
# The file must match the SHA-256 in inputs/EXTERNAL_INPUT_CHECKSUMS.tsv; any other copy is a FAIL, never a silent run.
#   - behaviour simulations          -> Supplementary Table 4 (Fig. 2a-b, Suppl. Fig. 2b-c)
#   - family architecture            -> source data of Suppl. Table 3, Fig. 3, Suppl. Fig. 3
#   - primary supplementary tables   -> Supplementary Tables 1, 2, 5 (cell-by-cell)
#   - trait-matrix missingness       -> Supplementary Table 12 per-trait counts and headline
# Needs R (digest, cluster) and Python 3 (pandas, numpy, openpyxl).
set -uo pipefail
cd "$(dirname "$0")/.."; export CTAS_REPO="$PWD"; W=work; fail=0
ok() { if [ "$1" = 0 ]; then echo "PASS  $2"; else echo "FAIL  $2"; fail=1; fi; }
Rscript analysis/simulations/run_behaviour_simulations.R > /dev/null 2>&1; ok $? "behaviour simulations ran"
python3 - <<'PY'; ok $? "ST4 sheets equal regenerated simulations (|diff| < 1e-12)"
import pandas as pd, numpy as np, sys
x='results/supplementary_tables/CTAS_Supplementary_Table_4_behaviour_simulations_101126.xlsx'; bad=0
for sh,c in [('calibration_typeI','sim1_calibration_typeI.csv'),('fixed_total_signal','sim2_fixed_total_signal.csv'),('breadth_by_strength','sim3_breadth_by_strength.csv')]:
    a=pd.read_excel(x,sheet_name=sh,header=1); b=pd.read_csv('work/simulations/'+c)
    d=max(np.nanmax(np.abs(a[k].astype(float).values-b[k].astype(float).values)) for k in b.columns)
    print(f'      {sh}: max |diff| {d:.2e}'); bad|= d>1e-12
sys.exit(1 if bad else 0)
PY
Rscript analysis/architecture/20_architecture.R > /dev/null 2>&1; ok $? "family architecture ran"
python3 - <<'PY'; ok $? "architecture locus table equals released figure source data (Fig. 3a)"
import pandas as pd, sys
a=pd.read_csv('work/architecture/architecture_locus_ranked.tsv',sep='\t'); f=pd.read_csv('results/figure_source_data/panel_F3a.tsv',sep='\t')
cols=[c for c in f.columns if c in a.columns]; m=a.set_index('canonical_locus_id').loc[f['canonical_locus_id'] if 'canonical_locus_id' in f else a.canonical_locus_id]
d=max((m[c].values.astype(float)-f[c].values.astype(float)).__abs__().max() for c in cols if c!='canonical_locus_id' and pd.api.types.is_numeric_dtype(f[c]))
print(f'      {len(cols)} shared columns, max |diff| {d:.2e}'); sys.exit(0 if d<1e-12 else 1)
PY
python3 - <<'PY'; ok $? "ST12 per-trait missingness and headline recomputed from the packaged gene matrix"
import pandas as pd, numpy as np, sys
m=pd.read_csv('inputs/UKBB_16trait_gene_matrix.tsv',sep='\t',index_col=0)
x='results/supplementary_tables/CTAS_Supplementary_Table_12_trait_missingness_101126.xlsx'
t=pd.read_excel(x,sheet_name='per_trait',header=1).set_index('trait')
v=m.notna().sum(); ok1=all(int(t.loc[k,'valid_genes'])==int(v[k]) and int(t.loc[k,'missing_count'])==int(len(m)-v[k]) for k in m.columns)
comp=int(m.notna().all(axis=1).sum()); miss=int(m.isna().sum().sum())
print(f'      16 traits per-trait counts match: {ok1}; complete genes {comp} (paper 18,303); missing cells {miss} (paper 3,144)')
sys.exit(0 if ok1 and comp==18303 and miss==3144 else 1)
PY
mkdir -p "$W/tables/Supplementary Tables" "$W/tables/work" "$W/tables/source_data"
python3 analysis/tables/41_st_primary.py > /dev/null 2>&1; ok $? "Supplementary Tables 1, 2, 5 built"
for n in 1_gene_architecture 2_null_calibration 5_correlation_sensitivity; do k=${n%%_*}
  python3 tests/compare_xlsx.py "$(ls "$W/tables/Supplementary Tables/ST${k}_"*.xlsx)" "results/supplementary_tables/CTAS_Supplementary_Table_${n}_101126.xlsx"; ok $? "ST${k} cell-identical to released"; done

# ---- completion stages (16 Nov 2026): locus map, sensitivity, FinnGen, colocalisation re-summary
Rscript analysis/locus_map/build_locus_map.R > /dev/null 2>&1; ok $? "locus map rebuilt from discovery-stage inputs"
CTAS_LOCUSMAP_OVERRIDE="$W/locus_map" Rscript -e 'source("config/ctas_paths.R"); LOCUSMAP <- Sys.getenv("CTAS_LOCUSMAP_OVERRIDE"); source("R/common.R"); LOCUSMAP <- Sys.getenv("CTAS_LOCUSMAP_OVERRIDE")
  G <- read.delim("results/core/final_C_gene_results_seed42.tsv", stringsAsFactors = FALSE); LM <- load_locus_map(); lk <- locus_key(G$gene, LM)
  n <- sapply(c(.05,.01,.001), function(a) length(unique(lk[G$q_BH < a]))); cat(sprintf("      rebuilt map -> keys identical: %s; loci %s\n", identical(lk, G$canonical_locus_id_or_mapping_key), paste(n, collapse="/")))
  quit(status = if (identical(lk, G$canonical_locus_id_or_mapping_key) && identical(as.integer(n), c(91L,46L,21L))) 0 else 1)' ; ok $? "91/46/21 loci reproduced from the rebuilt locus map"
Rscript analysis/sensitivity/00_reference_null_1e6.R > /dev/null 2>&1; ok $? "historical 10^6 reference null regenerated"
Rscript analysis/sensitivity/run_all_loo.R > /dev/null 2>&1; ok $? "leave-one-trait-out (16 reruns)"
Rscript analysis/sensitivity/run_v2_effective_n.R > /dev/null 2>&1; ok $? "effective-N sensitivity"
NCBI37_SHA=$(awk -F'\t' '$1=="NCBI37.3.gene.loc"{print $2}' inputs/EXTERNAL_INPUT_CHECKSUMS.tsv)
sha256_of() { if command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1; else sha256sum "$1" | cut -d' ' -f1; fi; }
if [ -n "${CTAS_NCBI37_GENE_LOC:-}" ] && { [ ! -f "$CTAS_NCBI37_GENE_LOC" ] || [ "$(sha256_of "$CTAS_NCBI37_GENE_LOC")" != "$NCBI37_SHA" ]; }; then
  echo "FAIL  obstructive panel not run: CTAS_NCBI37_GENE_LOC=$CTAS_NCBI37_GENE_LOC is missing or does not match SHA-256 $NCBI37_SHA (inputs/EXTERNAL_INPUT_CHECKSUMS.tsv)"; fail=1
elif [ -n "${CTAS_NCBI37_GENE_LOC:-}" ]; then
  echo "      NCBI37.3.gene.loc verified (SHA-256 $NCBI37_SHA)"
  Rscript analysis/sensitivity/run_obstructive_10trait.R > /dev/null 2>&1; ok $? "obstructive panel (first run)"
  mkdir -p "$W/sensitivity/obstructive/superseded_pre_mapping_correction_20260806"; cp "$W"/sensitivity/obstructive/{newly_significant,locus_retention,summary}.csv "$W/sensitivity/obstructive/superseded_pre_mapping_correction_20260806/"
  Rscript analysis/sensitivity/correct_obstructive_regional_mapping.R > /dev/null 2>&1; ok $? "obstructive panel (regional-mapping correction)"
else echo "SKIP  obstructive panel (set CTAS_NCBI37_GENE_LOC to the public MAGMA NCBI37.3.gene.loc file)"; fi
python3 analysis/finngen/20_analysis.py > /dev/null 2>&1; ok $? "FinnGen gene-level analysis from packaged MAGMA outputs"
python3 analysis/finngen/21_ledger.py > /dev/null 2>&1 && python3 analysis/finngen/100_families.py > /dev/null 2>&1; ok $? "FinnGen six-family measure"
python3 analysis/colocalisation/01_coloc_resummary.py > /dev/null 2>&1; ok $? "colocalisation re-summary from packaged cell-level results"
python3 - <<'PY'; ok $? "regenerated intermediates byte-identical to the original outputs (tests/expected_intermediate_sha256.tsv)"
import csv, hashlib, os, sys
rows = list(csv.DictReader(open('tests/expected_intermediate_sha256.tsv'), delimiter='\t')); bad = []; n = 0; skip = 0
for r in rows:
    p = r['work_path']
    if not os.path.exists(p):
        if '/obstructive/' in p: skip += 1; continue
        bad.append(p + ' (missing)'); continue
    n += 1
    if hashlib.sha256(open(p, 'rb').read()).hexdigest() != r['sha256']: bad.append(p)
print(f'      {n} files checked, {len(bad)} differ, {skip} skipped (obstructive needs NCBI37.3)'); [print('      DIFF', b) for b in bad]
sys.exit(1 if bad else 0)
PY
echo; [ $fail = 0 ] && echo "PACKAGED-INPUT REPRODUCTION: PASS" || echo "PACKAGED-INPUT REPRODUCTION: FAIL"; exit $fail
