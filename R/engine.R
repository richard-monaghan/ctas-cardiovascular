# CTAS V2 public release copy of outputs/large_null_rerun_20261005/code/engine.R (SHA-256 80c2e8da...)
# Only path handling was changed (see SCRIPT_PROVENANCE.tsv); all statistical code is verbatim.
################################################################################
# engine.R -- chunked large-null CTAS engine (audit only). Reuses the canonical
# functions in whitening_scale_alignment_audit_20261002/code/common.R (read-only).
################################################################################
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R")); source(file.path(ROOT, "R", "common.R"))
suppressPackageStartupMessages(library(parallel))
OUT2 <- WORK_DIR; W2 <- file.path(WORK_DIR, "null_runs"); dir.create(W2, showWarnings = FALSE, recursive = TRUE)

ctas_only <- function(ZW) { ZW[ZW < 0] <- 0; W <- ZW^2; s <- rowSums(W); ss <- rowSums(W^2)
  ifelse(s > 0, sqrt(s) * log2(1 + ifelse(s > 0, s^2 / ss, 0)), 0) }

# Chunked null. Each chunk gets its own L'Ecuyer stream derived from `seed` (deterministic, independent of worker count).
# Returns per-gene count of null scores > real (canonical convention), the retained upper tail (> t0), and chunk bookkeeping.
run_null <- function(Sigma, Sinv, real, M, seed, chunk = 2e6, t0, cores = CORES) {
  nch <- ceiling(M / chunk); sizes <- rep(chunk, nch); sizes[nch] <- M - chunk * (nch - 1)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed); s <- .Random.seed; seeds <- vector("list", nch)
  for (i in seq_len(nch)) { seeds[[i]] <- s; s <- parallel::nextRNGStream(s) }
  Q <- chol(Sigma) %*% Sinv      # (z chol(Sigma)) Sinv == z Q ; Q orthogonal
  res <- mclapply(seq_len(nch), function(i) {
    assign(".Random.seed", seeds[[i]], envir = .GlobalEnv)
    sc <- ctas_only(matrix(rnorm(sizes[i] * 16), ncol = 16) %*% Q); ss <- sort(sc)
    list(cnt = sizes[i] - findInterval(real, ss), top = sc[sc > t0], n = sizes[i], mx = max(sc))
  }, mc.cores = cores, mc.preschedule = FALSE)
  stopifnot(!any(sapply(res, inherits, "try-error")))
  list(count = Reduce(`+`, lapply(res, `[[`, "cnt")), top = sort(unlist(lapply(res, `[[`, "top")), decreasing = TRUE),
       M = sum(sapply(res, `[[`, "n")), max = max(sapply(res, `[[`, "mx")), nchunks = nch)
}
# type-7 quantile of the full null from its descending upper tail T (length K, all values > t0)
q7_from_top <- function(T, M, p) { h <- (M - 1) * p + 1; lo <- floor(h); hi <- ceiling(h)
  at <- function(i) T[M - i + 1]                     # ascending index i -> descending position M-i+1
  at(lo) + (h - lo) * (at(hi) - at(lo)) }
# canonical p / tail / BH from counts (tail model of stringency_sweep.R, evaluated on the retained upper tail)
counts_to_q <- function(real, run) {
  M <- run$M; T <- run$top; floor_p <- 1 / (M + 1); p <- (1 + run$count) / (1 + M); at_floor <- run$count == 0
  lo <- q7_from_top(T, M, 0.99); hi <- q7_from_top(T, M, 0.99999)
  stopifnot(length(T) > 0.0101 * M)   # tail set must reach below the 0.99 quantile
  xg <- seq(lo, hi, length.out = 30); Pg <- sapply(xg, function(x) sum(T >= x) / M); k <- Pg > 0
  fit <- lm(log(Pg) ~ xg, data = data.frame(xg = xg[k], Pg = Pg[k]))
  extrap <- rep(FALSE, length(p))
  if (any(at_floor)) { p[at_floor] <- pmin(floor_p, exp(predict(fit, newdata = data.frame(xg = real[at_floor])))); extrap <- at_floor & p < floor_p - 1e-15 }
  list(p = p, q = p.adjust(p, "BH"), n_count0 = sum(at_floor), n_extrap = sum(extrap), slope = unname(coef(fit)[2]), band = c(lo, hi), floor_p = floor_p)
}
load_inputs <- function() {
  P <- read_matrix(); SA <- as.matrix(read.csv(file.path(INPUTS, "sigma_A_raw_spearman.csv"), row.names = 1))
  SC <- sigma_C(SA); dimnames(SC) <- dimnames(SA); Z <- na0(p_to_z_na(P))
  list(P = P, Sig = list(A = SA, C = SC), Sinv = list(A = inv_sqrt(SA), C = inv_sqrt(SC)), Z = Z, genes = rownames(P),
       real = list(A = ctas_parts(Z, inv_sqrt(SA))$ctas, C = ctas_parts(Z, inv_sqrt(SC))$ctas))
}
