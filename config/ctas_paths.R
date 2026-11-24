################################################################################
# ctas_paths.R -- single place where CTAS V2 paths are defined.
# Every path is relative to the repository root unless overridden by environment
# variables, so nothing depends on the original author's machine.
#   CTAS_REPO      repository root (default: auto-detected from this file / getwd())
#   CTAS_WORK      scratch directory for regenerated null runs (default: <repo>/work)
#   CTAS_OUT       directory for regenerated frozen results (default: <repo>/work/freeze)
#   CTAS_CORES     cores for the chunked null (default: parallel::detectCores() - 1, min 1)
################################################################################
ctas_repo_root <- function() {
  env <- Sys.getenv("CTAS_REPO", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  d <- normalizePath(getwd())
  repeat { if (file.exists(file.path(d, "config", "ctas_paths.R"))) return(d)
    p <- dirname(d); if (p == d) stop("Cannot find repository root; set CTAS_REPO or run from inside the repository."); d <- p }
}
ROOT    <- ctas_repo_root()
INPUTS  <- file.path(ROOT, "inputs")
LOCUSMAP <- file.path(INPUTS, "locus_map")
WORK_DIR <- Sys.getenv("CTAS_WORK", file.path(ROOT, "work"))
OUT_DIR  <- Sys.getenv("CTAS_OUT",  file.path(WORK_DIR, "freeze"))
CORES    <- as.integer(Sys.getenv("CTAS_CORES", max(1L, parallel::detectCores() - 1L)))
dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)
