#!/usr/bin/env Rscript
# After a full rerun (run_core_pipeline.sh): compare regenerated work/freeze/ tables with the
# published results/core/ tables (numeric tolerance 1e-12 relative; strings exact).
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R"))
fs <- c("final_C_gene_results_seed42.tsv", "final_C_locus_results_seed42.tsv", "seed_stability_C.tsv", "eigenvalues_C.tsv")
ok <- TRUE
for (f in fs) { a <- read.delim(file.path(ROOT, "results/core", f), stringsAsFactors = FALSE); b <- read.delim(file.path(OUT_DIR, f), stringsAsFactors = FALSE)
  same <- identical(dim(a), dim(b)) && all(mapply(function(x, y) if (is.numeric(x)) isTRUE(all.equal(x, y, tolerance = 1e-12)) else identical(as.character(x), as.character(y)), a, b))
  cat(sprintf("%-38s %s\n", f, if (same) "MATCH" else "DIFFERS")); ok <- ok && same }
quit(status = if (ok) 0 else 1)
