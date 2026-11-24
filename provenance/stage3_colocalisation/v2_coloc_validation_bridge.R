#!/usr/bin/env Rscript
################################################################################
# CTAS V2 colocalisation — Phase 3 deterministic validation bridge.
#
# Step A (legacy reproduction): re-implement the V1 method faithfully (12 traits
#   beta_meta/se_meta/af_cases_meta; 4 traits CAD/AF/HF/HTN via MAGMA p-values) on
#   a fixed 8-locus subset of the 48 unchanged targets, and compare cell-by-cell
#   to the STORED V1 results (results_all_101/coloc_pp4_all_101.tsv). This proves
#   the rebuilt liftover/harmonise/coloc plumbing reproduces V1 before we trust it.
#
# Step B (corrected mode): read the corrected beta_EUR engine output (v2_coloc_run.R
#   checkpoints for the same 8 loci) and quantify the per-cell change legacy->corrected,
#   attributing it to the 3 corrections (Tfh removal / beta_meta->beta_EUR /
#   4 p-value->beta,se / EUR weighted MAF).
#
# Writes: coloc_review/validation/bridge_legacy_vs_v1.tsv,
#         coloc_review/validation/bridge_legacy_vs_corrected.tsv,
#         coloc_review/validation/BRIDGE_SUMMARY.txt
################################################################################
suppressPackageStartupMessages({ library(data.table); library(coloc) })
source("<CTAS_V2_PAPER_WORKSPACE>/_source_assets/scripts/v2_coloc_config.R")
V2$load_trait_N()

SUBSET <- c("V2L034","V2L017","V2L004","V2L033","V2L069","V2L040","V2L095","V2L055")
VDIR   <- file.path(V2$COLOC_DIR, "validation"); dir.create(VDIR, showWarnings=FALSE, recursive=TRUE)
coords <- fread(V2$COORDS)[canonical_locus_id %in% SUBSET]
v1res  <- fread(file.path(V2$VIGNETTE, "results_all_101", "coloc_pp4_all_101.tsv"))

## ---- Legacy V1 GWAS extraction (faithful to run_coloc_all_101.R) -----------
GWAS_MAG <- file.path(V2$BASE, "magma_cvd_results")
legacy_bgz <- data.table(
  trait  = c("Angina","Aortic_Aneurysm","Aortic_Valve","Atherosclerosis","Cardiomyopathy",
             "DVT","Ischemic_Stroke","MI","Mitral_Valve","PAD","TIA","VTE"),
  prefix = c("icd10-I20","icd10-I71","icd10-I35","icd10-I70","icd10-I42","icd10-I80",
             "icd10-I63","icd10-I21","icd10-I34","icd10-I73","icd10-G45","icd10-I82"))
legacy_mag <- data.table(
  trait = c("CAD_I25","AF_I48","HF_I50","Hypertension_401"),
  stem  = c("CAD_I25_rsids","AF_I48_rsids","HF_I50_rsids","Hypertension_401_rsids"))
Nmap <- V2$traits[, .(trait, n_cases, n_controls)]

extract_bgz_legacy <- function(prefix, chr, wstart, wend) {
  f <- file.path(V2$GWAS_BGZ, paste0(prefix, "-both_sexes.tsv.bgz")); if (!file.exists(f)) return(NULL)
  s19 <- max(1L, wstart - V2$BUFFER); e19 <- wend + V2$BUFFER
  cmd <- sprintf("(bgzip -d -c '%s' | head -1; tabix '%s' '%d:%d-%d') | cat", f, f, chr, s19, e19)
  raw <- tryCatch(fread(cmd = cmd, showProgress = FALSE), error = function(e) NULL)
  if (is.null(raw) || nrow(raw)==0) return(NULL)
  raw
}
extract_mag_legacy <- function(stem, chr, wstart, wend) {
  f <- file.path(GWAS_MAG, paste0(stem, ".magma")); if (!file.exists(f)) return(NULL)
  s19 <- max(1L, wstart - V2$BUFFER); e19 <- wend + V2$BUFFER
  cmd <- sprintf("awk -F'\\t' '$2==%d && $3>=%d && $3<=%d' '%s'", chr, s19, e19, f)
  tryCatch(fread(cmd=cmd, col.names=c("rsid","chr","pos","pvalue"), showProgress=FALSE,
                 colClasses=list(character="rsid")), error=function(e) NULL)
}
harmonise_pval_legacy <- function(gwas_dt, eqtl_dt) {
  e <- copy(eqtl_dt)[, .(rsid, chr, pos, ref, alt, beta_eqtl=beta, se_eqtl=se, maf_eqtl=maf)]
  g <- copy(gwas_dt)[, .(rsid, pvalue_gwas=pvalue)]
  out <- merge(e, g, by="rsid"); if (nrow(out)==0) return(NULL)
  out[, snp_id := paste(chr,pos,ref,alt,sep=":")]
  out <- out[!is.na(pvalue_gwas) & pvalue_gwas>0 & !is.na(beta_eqtl) & !is.na(se_eqtl) & se_eqtl>0]
  setorder(out, snp_id, pvalue_gwas); out[!duplicated(snp_id)]
}
run_coloc_legacy <- function(merged, n_eqtl, n_cases, n_controls, mode) {
  if (is.null(merged) || nrow(merged) < V2$MIN_SNPS) return(NULL)
  n_tot <- n_cases+n_controls; s <- n_cases/n_tot
  d_eqtl <- list(beta=merged$beta_eqtl, varbeta=merged$se_eqtl^2, N=n_eqtl, MAF=merged$maf_eqtl, type="quant", snp=merged$snp_id)
  if (mode=="beta") {
    d_gwas <- list(beta=merged$beta_gwas, varbeta=merged$se_gwas^2, N=n_tot,
                   MAF=pmax(0.001,pmin(0.499,merged$maf_gwas)), s=s, type="cc", snp=merged$snp_id)
  } else {
    d_gwas <- list(pvalues=merged$pvalue_gwas, N=n_tot, MAF=merged$maf_eqtl, s=s, type="cc", snp=merged$snp_id)
  }
  tryCatch(suppressMessages(suppressWarnings(
    coloc.abf(dataset1=d_gwas, dataset2=d_eqtl, p1=1e-4, p2=1e-4, p12=1e-5)$summary)), error=function(e) NULL)
}

## ---- Step A: legacy reproduction -------------------------------------------
leg <- list()
for (i in seq_len(nrow(coords))) {
  cr <- coords[i]; gene <- cr$gene; ensg <- cr$ensg
  cat(sprintf("[legacy %d/%d] %s %s\n", i, nrow(coords), cr$canonical_locus_id, gene))
  eqtl <- lapply(seq_len(nrow(V2$tissues)), function(ti) V2$load_slice(gene, V2$tissues$label[ti], ensg)$dt)
  for (tri in seq_len(nrow(legacy_bgz))) {
    tr <- legacy_bgz[tri]; nn <- Nmap[trait==tr$trait]
    raw <- extract_bgz_legacy(tr$prefix, as.integer(cr$chr), as.integer(cr$wstart), as.integer(cr$wend)); if (is.null(raw)) next
    l <- V2$liftover_dt(raw); g38 <- l$dt; if (is.null(g38)) next
    g38 <- g38[chr==as.integer(cr$chr) & pos>=cr$wstart & pos<=cr$wend]
    if ("low_confidence_EUR" %in% names(g38)) g38 <- g38[!(low_confidence_EUR %in% c("true",TRUE))]
    if (nrow(g38)==0) next
    bv <- if ("beta_meta" %in% names(g38)) g38$beta_meta else g38$beta_EUR
    sv <- if ("se_meta" %in% names(g38)) g38$se_meta else g38$se_EUR
    av <- if ("af_cases_meta" %in% names(g38)) g38$af_cases_meta else g38$af_cases_EUR
    g38[, `:=`(beta_gwas=bv, se_gwas=sv, maf_gwas=pmax(0.001,pmin(0.499,as.numeric(av))))]
    for (ti in seq_len(nrow(V2$tissues))) {
      if (is.null(eqtl[[ti]])) next
      hm <- V2$harmonise_beta(g38, eqtl[[ti]]); m <- hm$merged
      res <- run_coloc_legacy(m, V2$tissues$n_eqtl[ti], nn$n_cases, nn$n_controls, "beta"); if (is.null(res)) next
      leg[[length(leg)+1]] <- data.table(canonical_locus_id=cr$canonical_locus_id, gene=gene,
        tissue=V2$tissues$label[ti], trait=tr$trait, mode="beta_meta",
        n_snps=nrow(m), PP3_leg=unname(res["PP.H3.abf"]), PP4_leg=unname(res["PP.H4.abf"]))
    }
  }
  for (tri in seq_len(nrow(legacy_mag))) {
    tr <- legacy_mag[tri]; nn <- Nmap[trait==tr$trait]
    raw <- extract_mag_legacy(tr$stem, as.integer(cr$chr), as.integer(cr$wstart), as.integer(cr$wend)); if (is.null(raw)) next
    l <- V2$liftover_dt(raw); g38 <- l$dt; if (is.null(g38)) next
    g38 <- g38[chr==as.integer(cr$chr) & pos>=cr$wstart & pos<=cr$wend]; if (nrow(g38)==0) next
    for (ti in seq_len(nrow(V2$tissues))) {
      if (is.null(eqtl[[ti]])) next
      m <- harmonise_pval_legacy(g38, eqtl[[ti]])
      res <- run_coloc_legacy(m, V2$tissues$n_eqtl[ti], nn$n_cases, nn$n_controls, "pval"); if (is.null(res)) next
      leg[[length(leg)+1]] <- data.table(canonical_locus_id=cr$canonical_locus_id, gene=gene,
        tissue=V2$tissues$label[ti], trait=tr$trait, mode="magma_pval",
        n_snps=nrow(m), PP3_leg=unname(res["PP.H3.abf"]), PP4_leg=unname(res["PP.H4.abf"]))
    }
  }
}
legdt <- rbindlist(leg, fill=TRUE)
cmp <- merge(legdt, v1res[, .(gene, tissue, trait, PP4_v1=PP4, PP3_v1=PP3, gwas_source)],
             by=c("gene","tissue","trait"), all.x=TRUE)
cmp[, dPP4 := PP4_leg - PP4_v1]
fwrite(cmp, file.path(VDIR, "bridge_legacy_vs_v1.tsv"), sep="\t")

## ---- Step B: corrected beta_EUR (read checkpoints) -------------------------
corr <- rbindlist(lapply(SUBSET, function(l) {
  f <- file.path(V2$CKPT_DIR, paste0(l, ".tsv")); if (file.exists(f)) fread(f) else NULL }), fill=TRUE)
corr_ass <- corr[disposition=="assessed", .(gene, tissue, trait, PP4_cor=PP4, PP3_cor=PP3, n_cor=n_snps)]
bothcmp <- merge(legdt[, .(gene,tissue,trait,PP4_leg,PP3_leg,n_snps)], corr_ass, by=c("gene","tissue","trait"), all=TRUE)
bothcmp[, dPP4_corr := PP4_cor - PP4_leg]
fwrite(bothcmp, file.path(VDIR, "bridge_legacy_vs_corrected.tsv"), sep="\t")

## ---- Summary ---------------------------------------------------------------
sink(file.path(VDIR, "BRIDGE_SUMMARY.txt"))
cat("CTAS V2 coloc — Phase 3 validation bridge\n")
cat("Subset (8 unchanged targets):", paste(SUBSET, collapse=", "), "\n\n")
cat("=== Step A: legacy reproduction vs stored V1 ===\n")
ok <- cmp[!is.na(PP4_v1)]
cat(sprintf("cells compared: %d\n", nrow(ok)))
cat(sprintf("PP4 correlation (legacy vs V1): %.5f\n", cor(ok$PP4_leg, ok$PP4_v1)))
cat(sprintf("max |dPP4|: %.4f ; mean |dPP4|: %.4f ; median |dPP4|: %.4f\n",
            max(abs(ok$dPP4)), mean(abs(ok$dPP4)), median(abs(ok$dPP4))))
cat(sprintf("cells |dPP4|<0.01: %d (%.1f%%) ; <0.05: %d (%.1f%%)\n",
            sum(abs(ok$dPP4)<0.01), 100*mean(abs(ok$dPP4)<0.01),
            sum(abs(ok$dPP4)<0.05), 100*mean(abs(ok$dPP4)<0.05)))
cat("\nlargest legacy-vs-V1 discrepancies:\n")
print(ok[order(-abs(dPP4))][1:min(10,nrow(ok)), .(gene,tissue,trait,gwas_source,PP4_leg=round(PP4_leg,3),PP4_v1=round(PP4_v1,3),dPP4=round(dPP4,3))])
cat("\nby gwas_source:\n")
print(ok[, .(n=.N, max_absdPP4=round(max(abs(dPP4)),4), mean_absdPP4=round(mean(abs(dPP4)),4)), by=gwas_source])
cat("\n=== Step B: corrected beta_EUR vs legacy ===\n")
bb <- bothcmp[!is.na(PP4_cor) & !is.na(PP4_leg)]
if (nrow(bb)>0) {
  cat(sprintf("cells in both: %d ; PP4 correlation: %.4f\n", nrow(bb), cor(bb$PP4_cor, bb$PP4_leg)))
  cat(sprintf("max |dPP4 corr|: %.4f ; mean: %.4f\n", max(abs(bb$dPP4_corr)), mean(abs(bb$dPP4_corr))))
  clsf <- function(p) fifelse(p>=0.8,"strong",fifelse(p>=0.5,"moderate","below"))
  bb[, `:=`(cl_leg=clsf(PP4_leg), cl_cor=clsf(PP4_cor))]
  cat("\nclass-change table (legacy -> corrected):\n"); print(bb[, .N, by=.(cl_leg,cl_cor)][order(-N)])
  cat("\nlargest corrected-vs-legacy shifts:\n")
  print(bb[order(-abs(dPP4_corr))][1:min(12,nrow(bb)), .(gene,tissue,trait,PP4_leg=round(PP4_leg,3),PP4_cor=round(PP4_cor,3),dPP4=round(dPP4_corr,3))])
} else cat("corrected checkpoints not present yet.\n")
sink()
cat("BRIDGE DONE. See", file.path(VDIR,"BRIDGE_SUMMARY.txt"), "\n")
