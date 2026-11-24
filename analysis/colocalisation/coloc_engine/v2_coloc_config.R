# CTAS V2 public release copy of <CTAS_V2_PAPER_WORKSPACE>/_source_assets/scripts/v2_coloc_config.R. Changes: path block only (environment variables).
#!/usr/bin/env Rscript
################################################################################
# CTAS V2 colocalisation — shared config & method library (SINGLE SOURCE)
#
# Encodes the corrected, consistent V2 design:
#   * 104 frozen loci V2L001-V2L104; target = statistical_representative_gene.
#   * 7 valid GTEx cis-eQTL datasets ONLY. QTD000439/Schmiedel Tfh_memory EXCLUDED.
#   * ALL 16 Pan-UKBB traits via ONE consistent EUR beta/SE mode (no p-value path,
#     no beta_meta, no af_cases_meta shortcut). EUR case/control N from the
#     effective-N authority; EUR *weighted* AF from af_cases_EUR + af_controls_EUR.
#   * Column indices derived PER FILE by header name (layouts differ: e.g. MI=28
#     cols, CAD=52 cols) — never hard-coded positions.
#   * cis-window = gene body +/- 1 Mb (hg38); coloc.abf priors p1=p2=1e-4 p12=1e-5.
#   * GWAS hg19 -> hg38 liftOver; harmonise on chr:pos:ref:alt (match/flip/complement).
#   * MHC (chr6:28.5-33.5 Mb hg38) genes flagged; assessed explicitly, never called
#     "no colocalisation" if not evaluable.
#
# This file defines objects/functions only; it runs nothing. Sourced by the
# fetch, run and validation scripts.
################################################################################

suppressPackageStartupMessages({ library(data.table) })

## ---- Paths -----------------------------------------------------------------
V2 <- new.env()
## ---- release: every location comes from the environment (see analysis/colocalisation/README_COLOC.md) ----
V2$REPO        <- normalizePath(Sys.getenv("CTAS_REPO", "."))
V2$SLICES_V1   <- Sys.getenv("CTAS_EQTL_SLICES_V1")          # slices of the 49 originally cached targets
V2$SLICES_V2   <- Sys.getenv("CTAS_EQTL_SLICES_V2")          # slices of the 55 targets fetched in Aug 2026
V2$GWAS_BGZ    <- Sys.getenv("CTAS_PANUKBB_BGZ_DIR")         # 12 local Pan-UKBB .bgz (+ .tbi); 4 traits are read remotely
V2$CHAIN_FILE  <- Sys.getenv("CTAS_CHAIN_HG19_TO_HG38")      # UCSC hg19ToHg38.over.chain
V2$LIFTOVER    <- Sys.getenv("CTAS_LIFTOVER_BIN")            # UCSC liftOver binary
V2$COLOC_DIR   <- Sys.getenv("CTAS_COLOC_OUT", file.path(V2$REPO, "work", "coloc_engine"))
V2$CKPT_DIR    <- file.path(V2$COLOC_DIR, "checkpoints")
V2$LOG_DIR     <- file.path(V2$COLOC_DIR, "logs")
V2$MANIFEST    <- file.path(V2$REPO, "inputs/colocalisation/V2_104_COLOC_MANIFEST.csv")
V2$COORDS      <- file.path(V2$REPO, "inputs/colocalisation/v2_target_coords.tsv")
V2$EFFN        <- file.path(V2$REPO, "inputs/effective_n/effective_n_table.csv")
for (d in c(V2$CKPT_DIR, V2$LOG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
V2$REMOTE_BASE <- "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"

V2$BUFFER  <- 1000000L
V2$MHC_CHR <- 6L; V2$MHC_LO <- 28500000L; V2$MHC_HI <- 33500000L
V2$MIN_SNPS <- 50L   # min matched SNPs per coloc test (V1-consistent)

## ---- 7 valid eQTL datasets (Tfh EXCLUDED) ----------------------------------
V2$tissues <- data.table(
  label  = c("heart_LV","heart_AA","aorta","coronary","blood","liver","kidney_cortex"),
  qtd    = c("QTD000256","QTD000251","QTD000131","QTD000136","QTD000356","QTD000266","QTD000261"),
  n_eqtl = c(382L,       372L,       387L,       213L,       670L,       208L,       73L)
)
V2$EQT_COLS <- c("molecular_trait_id","chr","pos","ref","alt","variant_id",
                 "ma_samples","maf","pvalue","beta","se","type","ac","an",
                 "r2","molecular_trait_object_id","gene_id","median_tpm","rsid")

## ---- 16 traits: ONE consistent beta_EUR mode -------------------------------
# prefix -> {prefix}-both_sexes.tsv.bgz (local for 12, remote S3 for 4).
# EUR case/control N sourced from effective_n_table.csv at load time.
V2$traits <- data.table(
  trait  = c("CAD_I25","AF_I48","HF_I50","Hypertension_401","MI","Angina",
             "Ischemic_Stroke","TIA","VTE","Cardiomyopathy","Aortic_Valve",
             "Mitral_Valve","PAD","Aortic_Aneurysm","Atherosclerosis","DVT"),
  prefix = c("icd10-I25","icd10-I48","icd10-I50","phecode-401","icd10-I21","icd10-I20",
             "icd10-I63","icd10-G45","icd10-I82","icd10-I42","icd10-I35",
             "icd10-I34","icd10-I73","icd10-I71","icd10-I70","icd10-I80"),
  remote = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
             FALSE, FALSE, FALSE, FALSE, FALSE,
             FALSE, FALSE, FALSE, FALSE, FALSE)
)

V2$load_trait_N <- function() {
  effn <- fread(V2$EFFN)
  m <- match(V2$traits$trait, effn$trait)
  stopifnot(!anyNA(m))
  V2$traits[, n_cases    := effn$N_cases_EUR[m]]
  V2$traits[, n_controls := effn$N_controls_EUR[m]]
  invisible(V2$traits)
}

## ---- GWAS source resolver --------------------------------------------------
V2$gwas_src <- function(prefix, remote) {
  if (remote) sprintf("%s/%s-both_sexes.tsv.bgz", V2$REMOTE_BASE, prefix)
  else        file.path(V2$GWAS_BGZ, sprintf("%s-both_sexes.tsv.bgz", prefix))
}

# Header of a Pan-UKBB flat file (map columns by NAME, never by position).
# Local: bgzip -dc | head -1.  Remote: byte-range of first block -> bgzip -dc.
V2$gwas_header <- function(src, remote) {
  if (remote) {
    cmd <- sprintf("curl -sL --max-time 60 --range 0-131071 '%s' | bgzip -dc 2>/dev/null | head -1", src)
  } else {
    cmd <- sprintf("bgzip -dc '%s' 2>/dev/null | head -1", src)
  }
  h <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  if (length(h) == 0) return(NULL)
  strsplit(h[1], "\t", fixed = TRUE)[[1]]
}

# Extract a hg19 window from a Pan-UKBB file, returning a data.table with
# canonical columns: chr,pos,ref,alt,beta_gwas,se_gwas,maf_gwas,nlp,lowconf.
# Column positions are resolved from THIS file's header by name.
V2$extract_gwas <- function(trait_row, chr, wstart, wend, retries = 4L) {
  src <- V2$gwas_src(trait_row$prefix, trait_row$remote)
  hdr <- NULL
  for (a in seq_len(retries)) { hdr <- V2$gwas_header(src, trait_row$remote); if (!is.null(hdr)) break; Sys.sleep(2*a) }
  if (is.null(hdr)) return(list(dt = NULL, reason = "header_fetch_failed"))
  idx <- function(nm) { i <- match(nm, hdr); if (is.na(i)) NA_integer_ else i }
  need <- c(chr="chr", pos="pos", ref="ref", alt="alt",
            beta="beta_EUR", se="se_EUR", nlp="neglog10_pval_EUR",
            afca="af_cases_EUR", afco="af_controls_EUR", low="low_confidence_EUR")
  ix <- vapply(need, idx, integer(1))
  if (anyNA(ix[c("chr","pos","ref","alt","beta","se")]))
    return(list(dt = NULL, reason = "required_columns_missing"))

  s19 <- max(1L, wstart - V2$BUFFER); e19 <- wend + V2$BUFFER
  region <- sprintf("%d:%d-%d", chr, s19, e19)
  raw <- NULL
  for (a in seq_len(retries)) {
    raw <- tryCatch(
      fread(cmd = sprintf("tabix '%s' '%s'", src, region),
            header = FALSE, showProgress = FALSE, sep = "\t"),
      error = function(e) NULL)
    if (!is.null(raw) && nrow(raw) > 0) break
    Sys.sleep(2*a)
  }
  if (is.null(raw) || nrow(raw) == 0) return(list(dt = NULL, reason = "no_gwas_rows"))

  getc <- function(k) if (is.na(ix[[k]]) || ix[[k]] > ncol(raw)) NA else raw[[ix[[k]]]]
  afca <- suppressWarnings(as.numeric(getc("afca")))
  afco <- suppressWarnings(as.numeric(getc("afco")))
  Nca <- trait_row$n_cases; Nco <- trait_row$n_controls
  # EUR weighted allele frequency (documented; not af_cases_meta)
  afw <- (Nca*afca + Nco*afco) / (Nca + Nco)
  afw[is.na(afw) & !is.na(afco)] <- afco[is.na(afw) & !is.na(afco)]
  afw[is.na(afw) & !is.na(afca)] <- afca[is.na(afw) & !is.na(afca)]
  dt <- data.table(
    chr = suppressWarnings(as.integer(getc("chr"))),
    pos = suppressWarnings(as.integer(getc("pos"))),
    ref = as.character(getc("ref")),
    alt = as.character(getc("alt")),
    beta_gwas = suppressWarnings(as.numeric(getc("beta"))),
    se_gwas   = suppressWarnings(as.numeric(getc("se"))),
    maf_gwas  = pmax(0.001, pmin(0.499, afw)),
    nlp       = suppressWarnings(as.numeric(getc("nlp"))),
    lowconf   = if (is.na(ix[["low"]])) NA else as.character(getc("low"))
  )
  if (!all(is.na(dt$lowconf)))
    dt <- dt[!(lowconf %in% c("true","TRUE","True"))]
  dt <- dt[!is.na(chr) & !is.na(pos) & !is.na(beta_gwas) & !is.na(se_gwas) & se_gwas > 0]
  if (nrow(dt) == 0) return(list(dt = NULL, reason = "no_valid_gwas_rows"))
  list(dt = dt, reason = "ok")
}

## ---- liftover hg19 -> hg38 --------------------------------------------------
V2$liftover_dt <- function(dt, chr_col = "chr", pos_col = "pos") {
  bed_in <- tempfile(fileext=".bed"); bed_out <- tempfile(fileext=".bed"); bed_un <- tempfile(fileext=".bed")
  on.exit(unlink(c(bed_in, bed_out, bed_un)), add = TRUE)
  bed <- data.table(chrom = paste0("chr", dt[[chr_col]]),
                    start = dt[[pos_col]] - 1L, end = dt[[pos_col]], name = seq_len(nrow(dt)))
  fwrite(bed, bed_in, sep = "\t", col.names = FALSE)
  rc <- system2(V2$LIFTOVER, c(bed_in, V2$CHAIN_FILE, bed_out, bed_un), stdout = FALSE, stderr = FALSE)
  if (rc != 0 || !file.exists(bed_out) || file.size(bed_out) == 0) return(list(dt = NULL, n_in = nrow(dt), n_out = 0L, n_unmapped = nrow(dt)))
  lifted <- fread(bed_out, col.names = c("chrom","start38","end38","ridx"), showProgress = FALSE)
  out <- dt[lifted$ridx]
  out[[chr_col]] <- as.integer(sub("chr","",lifted$chrom))
  out[[pos_col]] <- as.integer(lifted$end38)
  list(dt = out, n_in = nrow(dt), n_out = nrow(out), n_unmapped = nrow(dt) - nrow(out))
}

## ---- allele harmonisation (match / flip / complement) ----------------------
V2$COMP <- c(A="T",T="A",C="G",G="C")
V2$flip_allele <- function(a) vapply(strsplit(toupper(a),""),
  function(x){ y<-V2$COMP[x]; if(anyNA(y)) NA_character_ else paste(y,collapse="") }, character(1))

# Returns list(merged, n_gwas, n_eqtl, n_match, n_flip, n_comp).
V2$harmonise_beta <- function(gwas_dt, eqtl_dt) {
  e <- copy(eqtl_dt)[, .(chr, pos, ref, alt, beta_eqtl=beta, se_eqtl=se, maf_eqtl=maf, pvalue=pvalue)]
  e[, key38 := paste(chr,pos,ref,alt,sep=":")]
  g <- copy(gwas_dt)[, .(chr,pos,ref_g=ref,alt_g=alt,beta_gwas,se_gwas,maf_gwas)]
  g[, key38      := paste(chr,pos,ref_g,alt_g,sep=":")]
  g[, key38_flip := paste(chr,pos,alt_g,ref_g,sep=":")]
  m1 <- merge(e, g[, .(key38,beta_gwas,se_gwas,maf_gwas)], by="key38")
  m2 <- merge(e[!key38 %in% m1$key38],
              g[, .(key38=key38_flip, beta_gwas=-beta_gwas, se_gwas, maf_gwas=1-maf_gwas)], by="key38")
  e_comp <- copy(e); e_comp[, key38_comp := paste(chr,pos,V2$flip_allele(ref),V2$flip_allele(alt),sep=":")]
  already <- c(m1$key38, m2$key38)
  m3 <- merge(e_comp[!key38 %in% already],
              g[, .(key38_comp=key38, beta_gwas, se_gwas, maf_gwas)], by="key38_comp")
  if (nrow(m3)) m3[, key38 := key38_comp]
  out <- rbindlist(list(m1,m2,m3), fill = TRUE, use.names = TRUE)
  n_match <- nrow(m1); n_flip <- nrow(m2); n_comp <- nrow(m3)
  if (nrow(out) == 0) return(list(merged=NULL, n_gwas=nrow(g), n_eqtl=nrow(e), n_match=0L, n_flip=0L, n_comp=0L))
  out[, snp_id := key38]
  out <- out[!is.na(beta_gwas) & !is.na(se_gwas) & !is.na(beta_eqtl) & !is.na(se_eqtl) & se_gwas>0 & se_eqtl>0]
  setorder(out, snp_id, pvalue); out <- out[!duplicated(snp_id)]
  list(merged=out, n_gwas=nrow(g), n_eqtl=nrow(e), n_match=n_match, n_flip=n_flip, n_comp=n_comp)
}

## ---- coloc.abf (GWAS cc x eQTL quant) --------------------------------------
V2$run_coloc <- function(merged, n_eqtl, n_cases, n_controls) {
  if (is.null(merged) || nrow(merged) < V2$MIN_SNPS) return(NULL)
  n_tot <- n_cases + n_controls; s <- n_cases / n_tot
  d_eqtl <- list(beta=merged$beta_eqtl, varbeta=merged$se_eqtl^2, N=n_eqtl,
                 MAF=merged$maf_eqtl, type="quant", snp=merged$snp_id)
  d_gwas <- list(beta=merged$beta_gwas, varbeta=merged$se_gwas^2, N=n_tot,
                 MAF=pmax(0.001, pmin(0.499, merged$maf_gwas)), s=s, type="cc", snp=merged$snp_id)
  tryCatch(suppressMessages(suppressWarnings(
    coloc::coloc.abf(dataset1=d_gwas, dataset2=d_eqtl, p1=1e-4, p2=1e-4, p12=1e-5)$summary)),
    error = function(e) NULL)
}

## ---- eQTL slice loader (cached V1 dir first, then V2 dir) -------------------
V2$load_slice <- function(gene, tissue_label, ensg) {
  f1 <- file.path(V2$SLICES_V1, sprintf("%s__%s.tsv", gene, tissue_label))
  f2 <- file.path(V2$SLICES_V2, sprintf("%s__%s.tsv", gene, tissue_label))
  f  <- if (file.exists(f1) && file.size(f1) > 0) f1 else if (file.exists(f2) && file.size(f2) > 0) f2 else NA
  if (is.na(f)) return(list(dt=NULL, status="missing", src=NA))
  dt <- tryCatch(fread(f, header=FALSE, col.names=V2$EQT_COLS, showProgress=FALSE), error=function(e) NULL)
  if (is.null(dt) || nrow(dt)==0) return(list(dt=NULL, status="empty_file", src=f))
  dt <- dt[gene_id == ensg & type == "SNP" & maf > 0.01]
  if (nrow(dt)==0) return(list(dt=NULL, status="zero_snps_after_filter", src=f))
  list(dt=dt, status="ok", src=f)
}

V2$is_mhc <- function(chr, gstart, gend)
  (as.integer(chr)==V2$MHC_CHR & gstart <= V2$MHC_HI & gend >= V2$MHC_LO)

invisible(NULL)
