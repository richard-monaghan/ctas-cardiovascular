# Resolve the public MAGMA NCBI37.3.gene.loc file for the obstructive-panel sensitivity analysis.
# Added 17 Nov 2026 (release hygiene). The file must be given explicitly via CTAS_NCBI37_GENE_LOC and must
# match the SHA-256 recorded in inputs/EXTERNAL_INPUT_CHECKSUMS.tsv. There is no fallback path, so a missing,
# truncated or different copy stops the analysis instead of being selected silently.
ctas_ncbi37_gene_loc <- function(repo) {
  ref <- read.delim(file.path(repo, "inputs", "EXTERNAL_INPUT_CHECKSUMS.tsv"), stringsAsFactors = FALSE)
  ref <- ref[ref$file == "NCBI37.3.gene.loc", ]
  stopifnot(nrow(ref) == 1)
  p <- Sys.getenv("CTAS_NCBI37_GENE_LOC", "")
  if (!nzchar(p)) stop("CTAS_NCBI37_GENE_LOC is not set: point it at the public MAGMA NCBI37.3.gene.loc file (SHA-256 ", ref$sha256, ").")
  if (!file.exists(p)) stop("CTAS_NCBI37_GENE_LOC does not exist: ", p)
  h <- digest::digest(file = p, algo = "sha256")
  if (!identical(h, ref$sha256))
    stop(sprintf("CTAS_NCBI37_GENE_LOC checksum %s does not match the expected %s (%d lines); truncated or different copies are rejected.",
                 h, ref$sha256, ref$lines))
  normalizePath(p)
}
