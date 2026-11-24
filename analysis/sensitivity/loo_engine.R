# CTAS V2 public release copy of panel_v2_leave_one_out/loo_engine.R. Changes: paths only.
#!/usr/bin/env Rscript
################################################################################
# CTAS V2 leave-one-trait-out engine. Each run is a FULLY RECALIBRATED analysis
# on a subset of traits: refit Sigma via the SAME fit_gaussian_copula(), rewhiten
# with Sigma^(-1/2), recompute the SAME locked CTAS (magnitude x breadth), draw a
# correctly-dimensioned N(0,Sigma) null (M draws, seed 42) and compute genome-wide
# empirical p and BH q. Not an approximation by dropping a component.
################################################################################
suppressPackageStartupMessages({ library(Matrix) })
REPO <- normalizePath(Sys.getenv("CTAS_REPO", "."))
source(file.path(REPO,"analysis/sensitivity/lib/00_utility_functions_REFERENCE_COPY.R"))  # fit_gaussian_copula
M <- as.numeric(Sys.getenv("LOO_M","1e6")); SEED <- 42

pm <- read.delim(file.path(REPO,"inputs/UKBB_16trait_gene_matrix.tsv"), check.names=FALSE)
genes <- pm$GENE; Pfull <- as.matrix(pm[,-1]); trait_names <- colnames(Pfull)
stopifnot(ncol(Pfull)==16, nrow(Pfull)==19080)

# CTAS on a matrix of z (genes x k), given symmetric Sigma^(-1/2)
ctas_rows <- function(Z, Sinv){
  ZW <- Z %*% Sinv                 # whiten (Sinv symmetric): row_i -> Sinv %*% z_i
  ZW[ZW < 0] <- 0; W <- ZW^2
  s <- rowSums(W); ss <- rowSums(W^2)
  pr <- ifelse(s>0, s^2/ss, 0)
  ifelse(s>0, sqrt(s)*log2(1+pr), 0)
}
p_to_z <- function(P){ Pc <- pmin(pmax(P,1e-300),1-1e-10); z <- qnorm(Pc, lower.tail=FALSE); z[is.na(z)] <- 0; z }

run_ctas <- function(Psub, seed=SEED, m=M){
  fit <- fit_gaussian_copula(Psub, verbose=FALSE); Sigma <- fit$cor_matrix
  rownames(Sigma) <- colnames(Sigma) <- colnames(Psub)
  stopifnot(isSymmetric(unname(Sigma)), min(eigen(Sigma,symmetric=TRUE)$values) > 0)
  eg <- eigen(Sigma, symmetric=TRUE); Sinv <- eg$vectors %*% diag(1/sqrt(pmax(eg$values,1e-10))) %*% t(eg$vectors)
  ctas_real <- ctas_rows(p_to_z(Psub), Sinv)
  set.seed(seed)
  draws <- matrix(rnorm(m*ncol(Psub)), ncol=ncol(Psub)) %*% chol(Sigma)
  ctas_null <- ctas_rows(draws, Sinv)
  ref <- sort(ctas_null)
  # empirical p = P(null >= real); +1/(M+1) convention avoids zero (retention uses q<0.05/0.01, well above floor)
  emp_p <- (m - findInterval(ctas_real - .Machine$double.eps, ref) + 1) / (m + 1)
  q <- p.adjust(emp_p, method="BH")
  list(sigma=Sigma, ctas=ctas_real, p=emp_p, q=q,
       eig_min=min(eg$values), cor_rng=range(Sigma[lower.tri(Sigma)]))
}
