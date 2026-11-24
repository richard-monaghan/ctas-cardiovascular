# CTAS V2 public release copy of CTAS V2 paper/_source_assets/R2_review/run_R2_simulations.R (SHA-256 see RECOVERED_SCRIPT_MANIFEST.tsv)
# Generates Supplementary Table 4 (behaviour simulations) and the data for Fig. 2a-b and Supplementary Fig. 2b-c.
# Changes: paths only; the sim0 null-reproduction self-check runs only if the superseded cached null is supplied.
# Usage: Rscript analysis/simulations/run_behaviour_simulations.R   (writes work/simulations/; needs R package digest)
#!/usr/bin/env Rscript
################################################################################
# Results Section 2 — CTAS validation simulations, run on the IMPLEMENTED CTAS
# statistic (Sigma whitening -> positive-part squared weights -> magnitude x
# participation-ratio breadth). NOT a V1 recurrence-count simulation.
#
# Key identity: under the Gaussian-copula null z ~ N(0, Sigma), the whitened
# vector z_white = Sigma^(-1/2) z ~ N(0, I_16). CTAS depends on z only through
# z_white, so controlled decorrelated evidence is specified directly in whitened
# space (identity covariance). Sigma enters only when mapping real genes.
#
# Records: seeds, replicate counts, parameter grid, Sigma checksum, score
# implementation, Monte-Carlo uncertainty.
################################################################################
suppressPackageStartupMessages({ library(digest) })
source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R"))
OUT <- file.path(WORK_DIR, "simulations"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REF_NULL_CACHE <- Sys.getenv("CTAS_REFERENCE_NULL_CACHE", "")   # optional superseded 1e6 cached null (not redistributed); only for the sim0 check

Sigma <- as.matrix(read.csv(file.path(INPUTS, "sigma_A_raw_spearman.csv"), row.names = 1))
SIGMA_MD5 <- digest(Sigma)
NT <- nrow(Sigma)
cat(sprintf("Sigma checksum (md5 digest): %s ; n_traits=%d\n", SIGMA_MD5, NT))

# vectorised CTAS on whitened evidence (rows = genes/reps, cols = 16 axes)
ctas_from_white <- function(ZW) {
  ZW[ZW < 0] <- 0            # positive part (preserves matrix dims, unlike pmax)
  W <- ZW^2
  s <- rowSums(W)
  ss <- rowSums(W^2)
  pr <- ifelse(s > 0, s^2 / ss, 0)
  ifelse(s > 0, sqrt(s) * log2(1 + pr), 0)
}
SCORE_MD5 <- digest(deparse(ctas_from_white))

################################################################################
# 0. VALIDATE the simulator reproduces the cached 1e6 null (seed 42, N(0,I))
################################################################################
M_NULL <- 1e6; SEED_NULL <- 42
set.seed(SEED_NULL)
null_white <- matrix(rnorm(M_NULL * NT), ncol = NT)
ctas_null <- ctas_from_white(null_white)
if (nzchar(REF_NULL_CACHE) && file.exists(REF_NULL_CACHE)) {
cache <- readRDS(REF_NULL_CACHE)
qs <- c(0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999)
val <- data.frame(quantile = qs,
                  sim_Nident = quantile(ctas_null, qs),
                  cached_NSigma = quantile(cache$Score_C_null, qs))
val$abs_diff <- abs(val$sim_Nident - val$cached_NSigma)
cat("\n[VALIDATION] whitened-space null vs cached N(0,Sigma) null quantiles:\n"); print(val, row.names = FALSE)
# assert on well-resolved quantiles (<=0.9999); the 0.99999 quantile rests on
# ~10 of 1e6 draws and its cross-sample MC scatter is expected to be larger.
stopifnot(max(val$abs_diff[val$quantile <= 0.9999]) < 0.05)
write.csv(val, file.path(OUT, "sim0_null_reproduction.csv"), row.names = FALSE)
} else cat("[VALIDATION] sim0 skipped: optional reference null cache not supplied (CTAS_REFERENCE_NULL_CACHE); released values are in results/supplementary_tables ST4 sheet validation_null_reproduction
")

# empirical p helper against the reference null
ref_sorted <- sort(ctas_null)
emp_p <- function(x) 1 - (findInterval(x, ref_sorted) - 0) / M_NULL
# CTAS cutoffs at the three empirical-p tiers (for selection-probability criteria)
c_alpha <- function(a) as.numeric(quantile(ctas_null, 1 - a))
CUT <- c(a05 = c_alpha(0.05), a01 = c_alpha(0.01), a001 = c_alpha(0.001))
cat(sprintf("\nCTAS cutoffs (null survival): a=0.05 -> %.3f ; a=0.01 -> %.3f ; a=0.001 -> %.3f\n",
            CUT["a05"], CUT["a01"], CUT["a001"]))

################################################################################
# 1. EMPIRICAL-NULL CALIBRATION (F2a): fresh independent null genes -> empirical
#    p should be ~Uniform; nominal vs observed type-I across alpha.
################################################################################
SEED_CAL <- 20260805; N_CAL <- 200000
set.seed(SEED_CAL)
cal_white <- matrix(rnorm(N_CAL * NT), ncol = NT)
cal_ctas  <- ctas_from_white(cal_white)
cal_p     <- emp_p(cal_ctas)
alphas <- c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001)
calib <- data.frame(alpha = alphas,
  observed_typeI = sapply(alphas, function(a) mean(cal_p <= a)),
  n_reject = sapply(alphas, function(a) sum(cal_p <= a)),
  n_test = N_CAL)
# Wilson-ish binomial 95% CI
calib$se <- sqrt(calib$observed_typeI * (1 - calib$observed_typeI) / N_CAL)
calib$ci_lo <- pmax(0, calib$observed_typeI - 1.96 * calib$se)
calib$ci_hi <- calib$observed_typeI + 1.96 * calib$se
calib$within_MC <- (calib$alpha >= calib$ci_lo & calib$alpha <= calib$ci_hi)
cat("\n[CALIBRATION] nominal vs observed type-I (independent null genes):\n"); print(calib, row.names = FALSE)
write.csv(calib, file.path(OUT, "sim1_calibration_typeI.csv"), row.names = FALSE)
# QQ data (observed empirical-p quantiles vs uniform), thinned
probs <- sort(cal_p); idx <- unique(round(seq(1, N_CAL, length.out = 2000)))
qq <- data.frame(expected_uniform = idx / N_CAL, observed_p = probs[idx])
write.csv(qq, file.path(OUT, "sim1_calibration_qq.csv"), row.names = FALSE)

################################################################################
# 2. FIXED-TOTAL-SIGNAL breadth selectivity (ST4): hold total non-centrality
#    lambda = sum(mu^2) constant; spread over k active axes (mu = sqrt(lambda/k)).
################################################################################
SEED_FTS <- 77; R_FTS <- 30000
lambdas <- c(20, 40, 80)
fts <- data.frame()
set.seed(SEED_FTS)
for (lam in lambdas) for (k in 1:NT) {
  mu <- sqrt(lam / k)
  M <- matrix(0, R_FTS, NT); M[, 1:k] <- mu
  ZW <- matrix(rnorm(R_FTS * NT), ncol = NT) + M
  sc <- ctas_from_white(ZW)
  fts <- rbind(fts, data.frame(lambda = lam, k_axes = k, per_axis_mu = mu,
     mean_CTAS = mean(sc), sd_CTAS = sd(sc),
     se_mean = sd(sc)/sqrt(R_FTS),
     sel_prob_a01 = mean(sc >= CUT["a01"]),
     sel_prob_a001 = mean(sc >= CUT["a001"]),
     n_rep = R_FTS))
}
cat("\n[FIXED-TOTAL-SIGNAL] mean CTAS by breadth k (lambda=40):\n")
print(fts[fts$lambda==40, c("k_axes","per_axis_mu","mean_CTAS","sel_prob_a01")], row.names = FALSE)
write.csv(fts, file.path(OUT, "sim2_fixed_total_signal.csv"), row.names = FALSE)

################################################################################
# 3. BREADTH x STRENGTH (ST5): grid of k active axes x per-axis strength mu.
################################################################################
SEED_BS <- 131; R_BS <- 20000
mus <- c(0.5, 1.0, 1.5, 2.0, 3.0)
ks  <- c(1, 2, 4, 8, 12, 16)
bs <- data.frame()
set.seed(SEED_BS)
for (mu in mus) for (k in ks) {
  M <- matrix(0, R_BS, NT); M[, 1:k] <- mu
  ZW <- matrix(rnorm(R_BS * NT), ncol = NT) + M
  sc <- ctas_from_white(ZW)
  bs <- rbind(bs, data.frame(per_axis_mu = mu, k_axes = k, total_lambda = k*mu^2,
     mean_CTAS = mean(sc), se_mean = sd(sc)/sqrt(R_BS),
     sel_prob_a05 = mean(sc >= CUT["a05"]),
     sel_prob_a01 = mean(sc >= CUT["a01"]),
     sel_prob_a001 = mean(sc >= CUT["a001"]),
     n_rep = R_BS))
}
cat("\n[BREADTH x STRENGTH] mean CTAS grid (rows mu, cols k):\n")
print(round(tapply(bs$mean_CTAS, list(bs$per_axis_mu, bs$k_axes), identity), 2))
write.csv(bs, file.path(OUT, "sim3_breadth_by_strength.csv"), row.names = FALSE)

# provenance
prov <- data.frame(
  item = c("sigma_md5","score_impl_md5","M_null","seed_null","seed_calibration","N_calibration",
           "seed_fixed_total_signal","R_fixed_total_signal","lambdas",
           "seed_breadth_strength","R_breadth_strength","mu_grid","k_grid",
           "cut_a05","cut_a01","cut_a001","null_model","signal_model"),
  value = c(SIGMA_MD5, SCORE_MD5, M_NULL, SEED_NULL, SEED_CAL, N_CAL,
            SEED_FTS, R_FTS, paste(lambdas,collapse=";"),
            SEED_BS, R_BS, paste(mus,collapse=";"), paste(ks,collapse=";"),
            round(CUT["a05"],4), round(CUT["a01"],4), round(CUT["a001"],4),
            "z_white ~ N(0, I_16) (== N(0,Sigma) after whitening)",
            "active axes shifted by mu in whitened space; noise N(0,1)"))
write.csv(prov, file.path(OUT, "sim_provenance.csv"), row.names = FALSE)
cat("\nAll R2 simulation review CSVs written to _source_assets/R2_review/.\n")
