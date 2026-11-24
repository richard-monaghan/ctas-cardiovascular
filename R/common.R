# CTAS V2 public release copy of outputs/whitening_scale_alignment_audit_20261002/code/common.R (SHA-256 ddcf4283...)
# Only path handling was changed (see SCRIPT_PROVENANCE.tsv); all statistical code is verbatim.
################################################################################
# common.R -- shared functions for the whitening scale-alignment audit (2026-10-02)
# AUDIT ONLY. Reads canonical inputs read-only; writes only inside this output dir.
# CTAS, clipping, NA handling, tail model and BH are reproduced verbatim from the
# locked scripts (compute_magnitude_score.R / stringency_sweep.R / loo_engine.R).
################################################################################
suppressPackageStartupMessages({ library(Matrix) })
if (!exists("ROOT")) source(file.path(Sys.getenv("CTAS_REPO", "."), "config", "ctas_paths.R"))
M_NULL <- 1e6; SEED_NULL <- 42; SEED_CAL <- 20260805; N_CAL <- 200000

read_matrix <- function() {
  pm <- read.delim(file.path(INPUTS, "UKBB_16trait_gene_matrix.tsv"), check.names = FALSE)
  P <- as.matrix(pm[, -1]); rownames(P) <- pm$GENE
  P
}
clip_P <- function(P) pmin(pmax(P, 1e-300), 1 - 1e-10)           # canonical clipping
p_to_z_na <- function(P) qnorm(clip_P(P), lower.tail = FALSE)      # NA kept
na0 <- function(Z) { Z[is.na(Z)] <- 0; Z }                          # canonical missing rule: z=0

# canonical pseudo-observations (same denominator convention as fit_gaussian_copula)
pseudo_obs <- function(P) apply(clip_P(P), 2, function(x)
  rank(x, ties.method = "average", na.last = "keep") / (sum(!is.na(x)) + 1))
# Gaussian-rank normal scores, ORIENTED like z = Phi^-1(1-p): small P -> large positive score
rank_scores_na <- function(P) qnorm(pseudo_obs(P), lower.tail = FALSE)

# canonical PD safeguard (threshold 0.01, nearPD corr=TRUE); returns list(S, invoked, eig_before, eig_after)
pd_safeguard <- function(S) {
  eb <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  inv <- min(eb) < 0.01
  if (inv) S <- as.matrix(Matrix::nearPD(S, corr = TRUE)$mat)
  ea <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  list(S = S, invoked = inv, eig_before = eb, eig_after = ea)
}
inv_sqrt <- function(S) { eg <- eigen(S, symmetric = TRUE)
  eg$vectors %*% diag(1 / sqrt(pmax(eg$values, 1e-10))) %*% t(eg$vectors) }   # canonical formula

# --- estimators (A-D). P = raw gene x trait P matrix (NA allowed) -------------
sigma_A <- function(P) cor(qnorm(pseudo_obs(P)), use = "pairwise.complete.obs", method = "spearman")  # VERBATIM canonical
sigma_B <- function(P) cor(rank_scores_na(P),   use = "pairwise.complete.obs", method = "pearson")
sigma_C <- function(SA) { S <- 2 * sin(pi * SA / 6); diag(S) <- 1; S }
sigma_D <- function(P) cor(p_to_z_na(P),        use = "pairwise.complete.obs", method = "pearson")

# --- CTAS (verbatim from loo_engine.R ctas_rows, == score_C) -------------------
ctas_parts <- function(Z, Sinv) {
  ZW <- Z %*% Sinv; ZW[ZW < 0] <- 0; W <- ZW^2
  s <- rowSums(W); ss <- rowSums(W^2)
  pr <- ifelse(s > 0, s^2 / ss, 0)
  mag <- sqrt(s)
  list(ctas = ifelse(s > 0, mag * log2(1 + pr), 0), mag = mag, pr = pr)
}
null_ctas <- function(Sigma, Sinv, m = M_NULL, seed = SEED_NULL) {   # N(0,Sigma), chol, seed 42 (canonical)
  set.seed(seed)
  draws <- matrix(rnorm(m * ncol(Sigma)), ncol = ncol(Sigma)) %*% chol(Sigma)
  ctas_parts(draws, Sinv)$ctas
}

# --- empirical p, exponential tail extrapolation, BH: VERBATIM stringency_sweep.R
score_to_q <- function(real, null, extrapolate = TRUE) {
  M <- length(null); floor_p <- 1 / (M + 1)
  ns <- sort(null)
  rank_ge <- length(ns) - findInterval(real, ns)
  p <- (1 + rank_ge) / (1 + M)
  at_floor <- p <= floor_p
  decades <- rep(0, length(p))
  if (extrapolate && any(at_floor)) {
    band <- quantile(null, c(0.99, 0.99999))
    xg <- seq(band[1], band[2], length.out = 30)
    Pg <- sapply(xg, function(x) mean(null >= x))
    keep <- Pg > 0
    fit <- lm(log(Pg) ~ xg, data = data.frame(xg = xg[keep], Pg = Pg[keep]))
    p[at_floor] <- pmin(floor_p, exp(predict(fit, newdata = data.frame(xg = real[at_floor]))))
    decades[at_floor] <- pmax(0, log10(floor_p) - log10(p[at_floor]))
  }
  list(p = p, q = p.adjust(p, method = "BH"), at_floor = at_floor, decades = decades,
       n_above_nominal_1e6 = sum(real > quantile(null, 1 - 1e-6)))
}

# --- locus mapping: locked machinery (as in run_v2_effective_n.R) ------------
load_locus_map <- function() {
  gla <- read.csv(file.path(LOCUSMAP, "gene_neighbourhood_map.csv"), stringsAsFactors = FALSE)
  g2l <- read.csv(file.path(LOCUSMAP, "stable_to_canonical_locus.csv"), stringsAsFactors = FALSE)
  master <- read.csv(file.path(LOCUSMAP, "canonical_locus_master.csv"), stringsAsFactors = FALSE)
  s2c <- unique(g2l[, c("stable_map_source_id", "canonical_locus_id")])
  list(gla = gla, s2c = s2c, master = master, g2l = g2l)
}
# key = canonical V2L id; else stable-map id (not in locked 104); else UNMAPPED:gene
locus_key <- function(genes, LM) {
  sid <- LM$gla$locus_id[match(genes, LM$gla$gene)]
  cid <- LM$s2c$canonical_locus_id[match(sid, LM$s2c$stable_map_source_id)]
  ifelse(!is.na(cid), cid, ifelse(!is.na(sid), paste0("STABLE:", sid), paste0("UNMAPPED:", genes)))
}
sessioninfo_line <- function() sprintf("R %s | Matrix %s | %s", R.version.string, packageVersion("Matrix"), Sys.info()[["sysname"]])
