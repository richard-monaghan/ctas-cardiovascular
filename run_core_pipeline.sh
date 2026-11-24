#!/usr/bin/env bash
# Runs the CTAS V2 core analysis from the packaged derived inputs.
#   step 01  Sigma_A check and Sigma_C                         (seconds)
#   step 02  10^9-draw nulls, estimators A and C, seeds 42/1/2  (~10 min on 10 cores; ~500 MB of .rds in work/null_runs)
#   step 03  primary freeze: 223/99/38 genes, 91/46/21 loci     (~1 min)
# Then compare work/freeze/ with results/core/ (tests/compare_freeze.R).
set -euo pipefail
cd "$(dirname "$0")"; export CTAS_REPO="$PWD"
Rscript analysis/core/01_sigma_matrices.R
Rscript analysis/core/02_run_large_null.R 1e9
Rscript analysis/core/03_freeze_primary.R
Rscript tests/compare_freeze.R
