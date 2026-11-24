# CTAS V2 public release copy of <CTAS_V2_PAPER_WORKSPACE>/_source_assets/scripts/v2_coloc_run.R. Changes: config path only.
#!/usr/bin/env Rscript
################################################################################
# CTAS V2 colocalisation — main resumable run engine.
#
# Design: 104 loci x 1 representative gene x 7 datasets x 16 traits, beta_EUR mode.
# Per-locus checkpoints (one TSV per locus in coloc_review/checkpoints/) => fully
# resumable; a locus already checkpointed is skipped. MHC loci (hg38 gene body in
# chr6:28.5-33.5 Mb) are recorded as not_assessed_MHC and skipped from coloc.
#
# Usage:
#   Rscript v2_coloc_run.R              # full 104-locus run (resumable)
#   Rscript v2_coloc_run.R V2L004,...   # run only listed loci (validation subset)
################################################################################
suppressPackageStartupMessages({ library(data.table); library(coloc) })
args <- commandArgs(trailingOnly = TRUE)
source(file.path(normalizePath(Sys.getenv("CTAS_REPO", ".")), "analysis/colocalisation/coloc_engine/v2_coloc_config.R"))
V2$load_trait_N()
dir.create(V2$CKPT_DIR, showWarnings = FALSE, recursive = TRUE)

coords <- fread(V2$COORDS)
setkey(coords, canonical_locus_id)
mani   <- fread(V2$MANIFEST)
if (length(args) >= 1 && nzchar(args[1])) {
  want <- strsplit(args[1], ",")[[1]]
  coords <- coords[canonical_locus_id %in% want]
  cat(sprintf("Subset run: %d loci (%s)\n", nrow(coords), paste(want, collapse=",")))
}

run_locus <- function(cr) {
  lid <- cr$canonical_locus_id; gene <- cr$gene; ensg <- cr$ensg
  ckpt <- file.path(V2$CKPT_DIR, sprintf("%s.tsv", lid))
  if (file.exists(ckpt) && file.size(ckpt) > 0) { cat(sprintf("  [skip] %s cached\n", lid)); return(invisible()) }
  rows <- list()
  add <- function(...) rows[[length(rows)+1]] <<- data.table(...)
  is_mhc <- as.logical(cr$is_mhc) || V2$is_mhc(cr$chr, cr$gstart, cr$gend)

  if (is_mhc) {
    for (ti in seq_len(nrow(V2$tissues))) for (tri in seq_len(nrow(V2$traits)))
      add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
          trait=V2$traits$trait[tri], disposition="not_assessed_MHC", n_snps=NA_integer_,
          PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,
          n_gwas=NA_integer_,n_eqtl=NA_integer_,n_lift_in=NA_integer_,n_lift_out=NA_integer_,
          n_match=NA_integer_,n_flip=NA_integer_,n_comp=NA_integer_)
    fwrite(rbindlist(rows), ckpt, sep="\t"); cat(sprintf("  [MHC] %s %s not_assessed\n", lid, gene)); return(invisible())
  }

  # Load 7 eQTL slices (cached V1 dir or new V2 dir), filter to target ENSG.
  eqtl <- vector("list", nrow(V2$tissues)); estat <- character(nrow(V2$tissues))
  for (ti in seq_len(nrow(V2$tissues))) {
    s <- V2$load_slice(gene, V2$tissues$label[ti], ensg)
    # NB: assigning NULL via [[<- would DELETE the element and shrink the list,
    # misaligning tissue indices. Leave the pre-allocated NULL in place instead.
    if (!is.null(s$dt)) eqtl[[ti]] <- s$dt
    estat[ti] <- s$status
  }

  # Per trait: extract GWAS window, liftover hg19->hg38, harmonise, coloc vs each tissue.
  for (tri in seq_len(nrow(V2$traits))) {
    tr <- V2$traits[tri]
    ext <- V2$extract_gwas(tr, as.integer(cr$chr), as.integer(cr$wstart), as.integer(cr$wend))
    if (is.null(ext$dt)) {
      for (ti in seq_len(nrow(V2$tissues)))
        add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
            trait=tr$trait, disposition=paste0("no_gwas:",ext$reason), n_snps=NA_integer_,
            PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,
            n_gwas=0L,n_eqtl=NA_integer_,n_lift_in=0L,n_lift_out=0L,
            n_match=NA_integer_,n_flip=NA_integer_,n_comp=NA_integer_)
      next
    }
    lift <- V2$liftover_dt(ext$dt)
    g38 <- lift$dt
    if (!is.null(g38)) g38 <- g38[chr==as.integer(cr$chr) & pos>=cr$wstart & pos<=cr$wend]
    for (ti in seq_len(nrow(V2$tissues))) {
      edt <- eqtl[[ti]]
      if (is.null(edt)) {
        add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
            trait=tr$trait, disposition=paste0("no_eqtl:",estat[ti]), n_snps=0L,
            PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,
            n_gwas=nrow(ext$dt),n_eqtl=0L,n_lift_in=lift$n_in,n_lift_out=lift$n_out,
            n_match=NA_integer_,n_flip=NA_integer_,n_comp=NA_integer_); next
      }
      if (is.null(g38) || nrow(g38)==0) {
        add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
            trait=tr$trait, disposition="no_gwas_in_window_post_liftover", n_snps=0L,
            PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,
            n_gwas=nrow(ext$dt),n_eqtl=nrow(edt),n_lift_in=lift$n_in,n_lift_out=lift$n_out,
            n_match=0L,n_flip=0L,n_comp=0L); next
      }
      hm <- V2$harmonise_beta(g38, edt); m <- hm$merged
      if (is.null(m) || nrow(m) < V2$MIN_SNPS) {
        add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
            trait=tr$trait, disposition="too_few_matched_snps", n_snps=if(is.null(m))0L else nrow(m),
            PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,
            n_gwas=nrow(g38),n_eqtl=nrow(edt),n_lift_in=lift$n_in,n_lift_out=lift$n_out,
            n_match=hm$n_match,n_flip=hm$n_flip,n_comp=hm$n_comp); next
      }
      res <- V2$run_coloc(m, V2$tissues$n_eqtl[ti], tr$n_cases, tr$n_controls)
      if (is.null(res)) {
        add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
            trait=tr$trait, disposition="coloc_error", n_snps=nrow(m),
            PP0=NA_real_,PP1=NA_real_,PP2=NA_real_,PP3=NA_real_,PP4=NA_real_,
            n_gwas=nrow(g38),n_eqtl=nrow(edt),n_lift_in=lift$n_in,n_lift_out=lift$n_out,
            n_match=hm$n_match,n_flip=hm$n_flip,n_comp=hm$n_comp); next
      }
      add(canonical_locus_id=lid, gene=gene, ensg=ensg, tissue=V2$tissues$label[ti],
          trait=tr$trait, disposition="assessed", n_snps=nrow(m),
          PP0=unname(res["PP.H0.abf"]),PP1=unname(res["PP.H1.abf"]),PP2=unname(res["PP.H2.abf"]),
          PP3=unname(res["PP.H3.abf"]),PP4=unname(res["PP.H4.abf"]),
          n_gwas=nrow(g38),n_eqtl=nrow(edt),n_lift_in=lift$n_in,n_lift_out=lift$n_out,
          n_match=hm$n_match,n_flip=hm$n_flip,n_comp=hm$n_comp)
    }
  }
  out <- rbindlist(rows, fill=TRUE)
  fwrite(out, ckpt, sep="\t")
  nass <- sum(out$disposition=="assessed"); mx <- suppressWarnings(max(out$PP4, na.rm=TRUE))
  cat(sprintf("  [done] %s %-10s assessed=%d/%d maxPP4=%s\n", lid, gene, nass, nrow(out),
              ifelse(is.finite(mx), sprintf("%.3f",mx), "NA")))
}

for (i in seq_len(nrow(coords))) {
  cr <- coords[i]
  cat(sprintf("[%d/%d] %s %s chr%s:%s-%s\n", i, nrow(coords),
              cr$canonical_locus_id, cr$gene, cr$chr, cr$wstart, cr$wend))
  tryCatch(run_locus(cr), error=function(e) cat(sprintf("  [ERROR] %s: %s\n", cr$canonical_locus_id, conditionMessage(e))))
}
cat("RUN COMPLETE\n")
