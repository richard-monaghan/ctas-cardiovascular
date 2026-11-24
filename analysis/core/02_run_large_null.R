# CTAS V2 public release copy of outputs/large_null_rerun_20261005/code/02_run_large_null.R (SHA-256 see SCRIPT_PROVENANCE.tsv)
# Only path handling was changed (see SCRIPT_PROVENANCE.tsv); all statistical code is verbatim.
# Usage: Rscript analysis/core/02_run_large_null.R [M]   (default 1e9; writes work/null_runs/run_<A|C>_seed<42|1|2>_M1e+09.rds)
#!/usr/bin/env Rscript
# STEP 3 -- M=1e8, estimators A and C, seeds 42/1/2. Saves one rds per run.
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R")); source(file.path(ROOT, "R", "engine.R"))
args <- commandArgs(TRUE); M <- if (length(args)) as.numeric(args[1]) else 1e9   # release default: the published M = 10^9
I <- load_inputs(); log <- list()
for (est in c("A", "C")) {
  set.seed(7); pil <- ctas_only(matrix(rnorm(2e6 * 16), ncol = 16) %*% (chol(I$Sig[[est]]) %*% I$Sinv[[est]])); t0 <- unname(quantile(pil, 0.985))
  for (seed in c(42, 1, 2)) {
    tt <- system.time(run <- run_null(I$Sig[[est]], I$Sinv[[est]], I$real[[est]], M, seed, chunk = 2e6, t0 = t0))
    saveRDS(run, file.path(W2, sprintf("run_%s_seed%d_M%.0e.rds", est, seed, M)))
    cq <- counts_to_q(I$real[[est]], run)
    cat(sprintf("%s seed %d M=%.0e: %.1f min | genes q<.05/.01/.001 = %d/%d/%d | count0=%d extrap=%d | tail kept %d\n", est, seed, M, tt[["elapsed"]] / 60,
      sum(cq$q < .05), sum(cq$q < .01), sum(cq$q < .001), cq$n_count0, cq$n_extrap, length(run$top)))
    log[[length(log) + 1]] <- data.frame(est = est, seed = seed, M = M, minutes = tt[["elapsed"]] / 60)
  }
}
write.table(do.call(rbind, log), file.path(W2, sprintf("runtime_M%.0e.tsv", M)), sep = "\t", quote = FALSE, row.names = FALSE)
