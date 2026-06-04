#' Allele Information for WGS Data (No GC Correction)
#'
#' Detects per-SNP duplication (and optionally deletion) signals in WGS data
#' using a joint genotype + sequencing-depth likelihood-ratio framework.
#' Depth is modelled per sample with a negative-binomial distribution; no
#' per-bin GC correction is applied. Use this variant when GC bias is
#' negligible (e.g. PCR-free libraries on small target regions) or when a
#' reference FASTA is not available.
#'
#' @section Pipeline overview:
#' 1. Optional region filter (\code{chroms} or \code{bed}).
#' 2. Global inbreeding coefficient (\eqn{F_{is}}), per-sample
#'    negative-binomial depth parameters, and beta-binomial
#'    allelic-ratio parameters are fit on the full filtered set.
#' 3. Per-chromosome loop: depth-based likelihoods are computed with
#'    \code{stats::dnbinom()} using the global per-sample NB parameters;
#'    genotype likelihoods are computed under beta-binomial (or binomial
#'    when \eqn{\rho < 10^{-4}}).
#' 4. Likelihood ratios for duplication (and deletion) are computed
#'    per SNP; if \code{perm_test = TRUE}, significance thresholds are
#'    estimated from a permutation pool drawn across all chromosomes.
#' 5. All per-chromosome results are merged into one \code{data.frame}.
#'
#' @param vcf VCF input. Either (a) a file path (character) that will be read
#'   via \code{readVCF()}, (b) a list returned by \code{readVCF()} with a
#'   \code{$vcf} element, or (c) a \code{data.frame} whose first column holds
#'   chromosome names.
#' @param fis Numeric scalar, global inbreeding coefficient. If \code{NULL}
#'   (default), it is estimated internally with \code{h.zygosity()}.
#' @param parallel Logical. If \code{TRUE}, parallelise the per-sample and
#'   permutation loops with a \code{parallel} cluster.
#' @param numCores Number of cores to use when \code{parallel = TRUE}.
#'   Defaults to \code{parallel::detectCores() - 1}.
#' @param nb_para Optional pre-computed negative-binomial parameter table
#'   (one row per sample, with columns \code{mean}, \code{size}, \code{mean2},
#'   \code{size2}). If \code{NULL}, the table is fit internally.
#' @param sig Numeric in \code{(0, 1)}, quantile threshold for the
#'   permutation-based significance cut-off (default \code{0.95}).
#' @param deletion Logical. If \code{TRUE}, also test for deletions
#'   (N = 1 vs. N = 2) in addition to duplications (N = 4 vs. N = 2).
#' @param perm_test Logical. If \code{TRUE} (default), estimate per-Nhet
#'   significance thresholds by permutation; if \code{FALSE}, the function
#'   returns raw likelihood ratios only (no \code{*.sig} or \code{*.stat}
#'   columns).
#' @param chroms Optional character vector of chromosome names to process
#'   (e.g. \code{c("chr1", "chr2")}). Mutually exclusive with \code{bed}.
#' @param bed Optional path to a BED file. Only SNPs overlapping at least one
#'   BED interval are processed. Mutually exclusive with \code{chroms}.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @importFrom methods is
#' @importFrom utils head
#'
#' @return A single `data.frame` containing one row per processed SNP
#'   (merged across all chromosomes / regions). Columns include:
#'
#'   - **First 4 columns**: Inherited from the input VCF (typically
#'     `CHROM`, `POS`, `REF`, `ALT` or similar metadata).
#'   - **`mdep, Nhet, Nmiss, p_alt, delta`**: Mean depth across samples,
#'     heterozygous-genotype count, missing-genotype count, alt-allele
#'     frequency, and the deviation from expected Nhet given \eqn{F_{is}}.
#'   - **`lhr.dup, dup.sig, dup.stat`**: Duplication likelihood ratio,
#'     permutation threshold for the SNP's Nhet, and classification
#'     (`"duplicated"` / `"non-duplicated"`).
#'   - **`lhr.del, del.sig, del.stat`**: (only when `deletion = TRUE`)
#'     deletion analogues of the above.
#'
#'   When `perm_test = FALSE`, the `*.sig` and `*.stat` columns are omitted
#'   and only the raw `lhr.*` columns are returned.
#'
#' @section Region-filter behaviour:
#'   \code{chroms} and \code{bed} cannot be supplied together. With no
#'   filter, every SNP in the VCF is processed. With \code{chroms}, SNPs whose
#'   \code{CHROM} is in the vector are kept (missing chromosomes generate a
#'   warning). With \code{bed}, SNPs overlapping any interval are kept.
#'   In all cases, the per-chromosome loop iterates over the chromosomes that
#'   actually appear in the retained SNPs, and outputs are merged before
#'   return.
#'
#' @section When to use GC correction instead:
#'   The negative-binomial depth model used here treats all SNPs from one
#'   sample as draws from a single distribution, so any systematic GC-driven
#'   over- or under-coverage will inflate the false-positive rate at the
#'   extremes of the GC spectrum. If your sequencing protocol shows visible
#'   GC bias (e.g. PCR-amplified libraries, hybrid-capture targets), use
#'   \code{\link{allele.info.WGS.GCcor}} instead.
#'
#' @section Memory:
#'   The dominant memory term is the per-chromosome likelihood matrix
#'   (\code{nSNP_chr} rows by \code{4 * nsample} or \code{6 * nsample}
#'   columns). It is allocated and freed inside each iteration of the
#'   chromosome loop, so peak memory scales with the largest chromosome, not
#'   the full genome.
#'
#' @seealso \code{\link{allele.info.WGS.GCcor}} (GC-bias-corrected variant),
#'   \code{parse_region_filter} (internal region filter helper).
#'
#' @examples
#' \dontrun{
#' # Whole genome
#' res <- allele.info.WGS3(vcf = "sample.vcf.gz")
#'
#' # Subset of chromosomes, with deletion testing, in parallel
#' res <- allele.info.WGS3(
#'   vcf      = "sample.vcf.gz",
#'   chroms   = c("chr1", "chr2"),
#'   deletion = TRUE,
#'   parallel = TRUE,
#'   numCores = 8
#' )
#'
#' # Restrict to a BED file, raw LHR only (no permutation thresholds)
#' res <- allele.info.WGS3(
#'   vcf       = "sample.vcf.gz",
#'   bed       = "candidate_regions.bed",
#'   perm_test = FALSE
#' )
#' }
#'
#' @export
allele.info.WGS3 <- function(vcf       = NULL,
                             fis       = NULL,
                             parallel  = FALSE,
                             numCores  = NULL,
                             nb_para   = NULL,
                             sig       = 0.95,
                             deletion  = FALSE,
                             perm_test = TRUE,
                             chroms    = NULL,   # NEW: character vector of chr names
                             bed       = NULL,   # NEW: path to BED file
                             ...) {

  ## -- 0. Parallel setup ------------------------------------------------------
  cl <- NULL
  if (parallel) {
    if (is.null(numCores)) numCores <- parallel::detectCores() - 1L
    cl <- parallel::makeCluster(numCores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, {
      requireNamespace("stringi",    quietly = TRUE)
      requireNamespace("extraDistr", quietly = TRUE)
    })
  }

  pb_apply <- if (parallel) {
    function(X, MARGIN, FUN, ...) pbapply::pbapply(X, MARGIN, FUN, ..., cl = cl)
  } else {
    function(X, MARGIN, FUN, ...) pbapply::pbapply(X, MARGIN, FUN, ...)
  }
  pb_lapply <- if (parallel) {
    function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl)
  } else {
    function(X, FUN, ...) pbapply::pblapply(X, FUN, ...)
  }

  eps      <- 1e-300
  safe_log <- function(z) log(pmax(z, eps))

  ## -- 1. Parse VCF ------------------------------------------------------------
  if (is.character(vcf)) {
    vcf <- readVCF(vcf)
    vcf <- vcf$vcf
  } else if (!is.data.frame(vcf)) {
    vcf <- vcf$vcf
  }
  colnames(vcf)[1] <- "CHROM"

  gt <- hetTgen(vcf, "GT")
  ad <- hetTgen(vcf, "AD")
  ad[ad == ""] <- "0,0"

  ## -- 2. Region filter --------------------------------------------------------
  row_mask <- parse_region_filter(ad[, 1:2], chroms = chroms, bed = bed)

  vcf <- vcf[row_mask, , drop = FALSE]
  gt  <- gt[row_mask,  , drop = FALSE]
  ad  <- ad[row_mask,  , drop = FALSE]
  rm(row_mask); gc()

  # Per-chromosome row chunks (indices into the now-contiguous filtered rows)
  chr_col    <- as.character(ad[, 1])
  chr_chunks <- split(seq_len(nrow(ad)), chr_col)
  chr_names  <- names(chr_chunks)
  n_chr      <- length(chr_names)
  message(sprintf("Will process %d chromosome(s): %s",
                  n_chr, paste(chr_names, collapse = ", ")))

  ## -- 3. Global Fis ----------------------------------------------------------
  if (is.null(fis)) {
    fis_est <- h.zygosity(vcf, cl = if (parallel) cl else NULL)
    fis     <- round(mean(fis_est$Fis, na.rm = TRUE), 2)
  }
  message(sprintf("Estimated global inbreeding coefficient: Fis = %.2f", fis))
  rm(vcf); gc()

  # Keep SNP metadata for per-chr slicing into ad.info
  meta    <- ad[, 1:4, drop = FALSE]
  AD      <- ad[, -(1:4), drop = FALSE]
  GT_raw  <- gt[, -(1:4), drop = FALSE]
  rm(ad, gt); gc()

  nsample    <- ncol(AD)
  nrow_filt  <- nrow(AD)
  idx_sample <- seq_len(nsample)

  GT_full <- matrix(-9L, nrow_filt, nsample)
  GT_full[GT_raw == "0/0" | GT_raw == "0|0"]                       <- 0L
  GT_full[GT_raw == "0/1" | GT_raw == "0|1" | GT_raw == "1|0"]    <- 1L
  GT_full[GT_raw == "1/1" | GT_raw == "1|1"]                       <- 2L
  rm(GT_raw); gc()

  ## -- 4. Global NB + beta-binomial parameter estimation ---------------------
  message("step 1/4: global parameter estimation")
  message("  Extracting depth of coverage (all filtered SNPs)...")

  tmp_g      <- pb_apply(AD, 2, ad_parse, mode = "ref_alt")
  ref_g      <- matrix(0L, nrow_filt, nsample)
  alt_g      <- matrix(0L, nrow_filt, nsample)
  for (j in seq_along(tmp_g)) {
    ref_g[, j] <- tmp_g[[j]]$ref
    alt_g[, j] <- tmp_g[[j]]$alt
  }
  rm(tmp_g); gc()

  tot_g <- ref_g + alt_g
  rm(alt_g); gc()

  # NB depth parameters (per sample, global)
  if (is.null(nb_para)) {
    message("  Fitting NB depth model...")
    colstat       <- pb_apply(tot_g, 2, db_mod)
    colstat       <- do.call(rbind, colstat)
    colstat$mean2 <- 2 * colstat$mean
    colstat$size2 <- colstat$size
    nb_stat       <- colstat
  } else {
    nb_stat <- nb_para
  }

  # Beta-binomial parameters from heterozygous inlier sites
  message("  Fitting beta-binomial model...")
  H_full      <- (GT_full == 1L)
  tot.rmOutl  <- pb_apply(AD, 2, ad_parse, mode = "tot_inlier", outlier_sd = 1.5)
  het_inl_idx <- which(H_full & tot.rmOutl, arr.ind = TRUE)
  rm(tot.rmOutl); gc()

  ad_het_vec <- AD[het_inl_idx]
  if (length(ad_het_vec) > 1e5)
    ad_het_vec <- ad_het_vec[sample(seq_along(ad_het_vec), 1e5)]

  alle_ad_het <- do.call(rbind, pb_lapply(ad_het_vec, ad_parse, mode = "ref_tot"))
  rm(ad_het_vec); gc()

  fit_bb         <- betabinom_mle(alle_ad_het$ref, alle_ad_het$tot)
  alpha          <- round(fit_bb$estimates[1], 2)
  bb_beta        <- round(fit_bb$estimates[2], 2)
  rho            <- 1 / (alpha + bb_beta + 1)
  use_beta_binom <- (rho > 1e-4)

  message(sprintf(
    "  Beta-binomial: alpha = %.2f, beta = %.2f%s",
    alpha, bb_beta,
    if (!use_beta_binom) "  [rho < 1e-4 -> switching to binomial]" else ""
  ))
  p.cor <- if (!use_beta_binom)
    mean(alle_ad_het$ref / alle_ad_het$tot, na.rm = TRUE) else NULL
  rm(alle_ad_het); gc()

  # Free global depth matrices - per-chromosome data will be re-derived in the loop
  rm(ref_g, tot_g); gc()

  ## -- 5. Helper: build (lhr, gp, dp) from a single chromosome's lld block ---
  build_lhr <- function(lld_chr, g_off, d_off, base_off) {
    gp <- safe_log(lld_chr[, idx_sample + g_off    * nsample]) -
      safe_log(lld_chr[, idx_sample            ])
    dp <- safe_log(lld_chr[, idx_sample + d_off    * nsample]) -
      safe_log(lld_chr[, idx_sample + base_off * nsample])
    list(lhr = 2 * rowSums(pmax(gp + dp, 0)), gp = gp, dp = dp)
  }

  ## -- 6. Per-chromosome loop -------------------------------------------------
  message(sprintf("step 2/4: per-chromosome likelihood computation (%d chromosome(s))",
                  n_chr))

  results_list      <- vector("list", n_chr)
  names(results_list) <- chr_names

  # Permutation sample accumulators (5k per chr; pooled then resampled at the end)
  perm_acc <- list(geno_dup = numeric(0), dep_dup = numeric(0))
  if (deletion) {
    perm_acc$geno_del <- numeric(0)
    perm_acc$dep_del  <- numeric(0)
  }

  for (ci in seq_along(chr_names)) {
    chr       <- chr_names[ci]
    idx       <- chr_chunks[[chr]]
    n_snp_chr <- length(idx)

    message(sprintf("  [%d/%d] %s  (%d SNPs)", ci, n_chr, chr, n_snp_chr))

    AD_chr   <- AD[idx,      , drop = FALSE]
    GT_chr   <- GT_full[idx, , drop = FALSE]
    meta_chr <- meta[idx,    , drop = FALSE]

    H_chr     <- (GT_chr == 1L)
    Nhet_chr  <- rowSums(H_chr)
    Nmiss_chr <- rowSums(GT_chr == -9L)

    # ref / alt / tot for this chromosome
    tmp_chr <- pb_apply(AD_chr, 2, ad_parse, mode = "ref_alt")
    ref_chr <- matrix(0L, n_snp_chr, nsample)
    alt_chr <- matrix(0L, n_snp_chr, nsample)
    for (j in seq_along(tmp_chr)) {
      ref_chr[, j] <- tmp_chr[[j]]$ref
      alt_chr[, j] <- tmp_chr[[j]]$alt
    }
    rm(tmp_chr); gc()

    tot_chr <- ref_chr + alt_chr
    rm(alt_chr); gc()

    # info table for this chromosome
    p_alt_chr  <- round(
      rowMeans(replace(GT_chr, GT_chr == -9L, NA_integer_), na.rm = TRUE) / 2, 2
    )
    E.Nhet_chr <- 2 * p_alt_chr * (1 - p_alt_chr) * (1 - fis) * (nsample - Nmiss_chr)

    ad_info_chr <- data.frame(
      meta_chr,
      mdep  = rowMeans(tot_chr),
      Nhet  = Nhet_chr,
      Nmiss = Nmiss_chr,
      p_alt = p_alt_chr,
      delta = round(E.Nhet_chr - Nhet_chr, 2),
      stringsAsFactors = FALSE
    )

    # Depth-based likelihoods (NB, per-sample global params)
    dep_lld2 <- stats::dnbinom(t(tot_chr), size = nb_stat$size,  mu = nb_stat$mean)
    dep_lld4 <- stats::dnbinom(t(tot_chr), size = nb_stat$size2, mu = nb_stat$mean2)
    if (deletion)
      dep_lld1 <- stats::dnbinom(t(tot_chr), size = nb_stat$size, mu = nb_stat$mean / 2)

    # Capture zero-mask BEFORE modifying dep_lld2 (silent-no-op fix carried over)
    zero_mask <- (dep_lld2 == 0)
    dep_lld4[zero_mask] <- 1e-50
    dep_lld2[zero_mask] <- 1e-100
    if (deletion) dep_lld1[zero_mask] <- 1e-150
    rm(zero_mask)

    # Genotype + depth likelihoods
    lld_chr <- if (use_beta_binom) {
      pb_apply(cbind(ref_chr, tot_chr, GT_chr), 1, cal_lld_v2.1,
               inb = fis, nb_stat = nb_stat,
               alpha = alpha, beta = bb_beta, deletion = deletion)
    } else {
      pb_apply(cbind(ref_chr, tot_chr, GT_chr), 1, cal_lld_v2.1,
               inb = fis, nb_stat = nb_stat,
               p.cor = p.cor, deletion = deletion)
    }
    rm(ref_chr, tot_chr, GT_chr); gc()

    lld_chr <- if (deletion) {
      t(rbind(lld_chr, dep_lld2, dep_lld1, dep_lld4))
    } else {
      t(rbind(lld_chr, dep_lld2, dep_lld4))
    }
    rm(dep_lld2, dep_lld4)
    if (deletion) rm(dep_lld1)
    gc()

    # LHR + permutation sample accumulation
    if (deletion) {
      # layout: geno = (N2, N1, N4); depth = (N2, N1, N4)
      dup_chr <- build_lhr(lld_chr, g_off = 2L, d_off = 5L, base_off = 3L)
      del_chr <- build_lhr(lld_chr, g_off = 1L, d_off = 4L, base_off = 3L)

      if (perm_test) {
        n_acc <- min(5000L, sum(!is.na(dup_chr$gp[H_chr])))
        perm_acc$geno_dup <- c(perm_acc$geno_dup,
                               sample(na.omit(dup_chr$gp[H_chr]), n_acc))
        perm_acc$dep_dup  <- c(perm_acc$dep_dup,
                               sample_from_mat(dup_chr$dp,
                                               min(5000L, length(dup_chr$dp))))
        n_acc <- min(5000L, sum(!is.na(del_chr$gp[H_chr])))
        perm_acc$geno_del <- c(perm_acc$geno_del,
                               sample(na.omit(del_chr$gp[H_chr]), n_acc))
        perm_acc$dep_del  <- c(perm_acc$dep_del,
                               sample_from_mat(del_chr$dp,
                                               min(5000L, length(del_chr$dp))))
      }

      results_list[[chr]] <- list(
        ad_info = ad_info_chr,
        lhr_dup = dup_chr$lhr,
        lhr_del = del_chr$lhr,
        Nhet    = Nhet_chr
      )

    } else {
      # layout: geno = (N2, N4); depth = (N2, N4)
      dup_chr <- build_lhr(lld_chr, g_off = 1L, d_off = 3L, base_off = 2L)

      if (perm_test) {
        n_acc <- min(5000L, sum(!is.na(dup_chr$gp[H_chr])))
        perm_acc$geno_dup <- c(perm_acc$geno_dup,
                               sample(na.omit(dup_chr$gp[H_chr]), n_acc))
        perm_acc$dep_dup  <- c(perm_acc$dep_dup,
                               sample_from_mat(dup_chr$dp,
                                               min(5000L, length(dup_chr$dp))))
      }

      results_list[[chr]] <- list(
        ad_info = ad_info_chr,
        lhr_dup = dup_chr$lhr,
        Nhet    = Nhet_chr
      )
    }

    rm(lld_chr, H_chr, ad_info_chr); gc()
    message(sprintf("    %s complete.", chr))
  }

  rm(AD, GT_full, H_full, meta); gc()

  ## -- 7. Merge per-chromosome results ----------------------------------------
  message("step 4/4: merging results")

  info_all    <- do.call(rbind, lapply(results_list, `[[`, "ad_info"))
  lhr_dup_all <- do.call(c,     lapply(results_list, `[[`, "lhr_dup"))
  Nhet_all    <- do.call(c,     lapply(results_list, `[[`, "Nhet"))
  if (deletion)
    lhr_del_all <- do.call(c,   lapply(results_list, `[[`, "lhr_del"))
  rm(results_list); gc()

  ## -- 8. Early return (no permutation test) -----------------------------------
  if (!perm_test) {
    if (deletion) {
      return(data.frame(info_all, lhr.dup = lhr_dup_all, lhr.del = lhr_del_all))
    } else {
      return(data.frame(info_all, lhr.dup = lhr_dup_all))
    }
  }

  ## -- 9. Global permutation thresholds (pooled across chromosomes) ------------
  message("step 3/4: estimating permutation thresholds (global pool)")

  run_perm_thresh <- function(geno_acc, dep_acc, label) {
    message(sprintf("  %s...", label))
    geno_s <- sample(geno_acc, min(1e5L, length(geno_acc)))
    dep_s  <- sample(dep_acc,  min(1e5L, length(dep_acc)))
    thre <- unlist(pb_lapply(
      unique(Nhet_all), sample_and_quantile2,
      nsample = nsample, geno = geno_s, dep = dep_s,
      nrep = 10000L, quan = sig
    ))
    names(thre) <- unique(Nhet_all)
    thre
  }

  thre_dup <- run_perm_thresh(perm_acc$geno_dup, perm_acc$dep_dup,
                              "duplication thresholds")
  if (deletion)
    thre_del <- run_perm_thresh(perm_acc$geno_del, perm_acc$dep_del,
                                "deletion thresholds")
  rm(perm_acc); gc()

  ## -- 10. Assemble final output -----------------------------------------------
  if (deletion) {
    info <- data.frame(
      info_all,
      lhr.dup = lhr_dup_all, dup.sig = thre_dup[as.character(Nhet_all)],
      lhr.del = lhr_del_all, del.sig = thre_del[as.character(Nhet_all)]
    )
    info$dup.stat <- ifelse(info$lhr.dup > info$dup.sig, "duplicated", "non-duplicated")
    info$del.stat <- ifelse(info$lhr.del > info$del.sig, "deleted",    "non-deleted")
  } else {
    info <- data.frame(
      info_all,
      lhr.dup = lhr_dup_all,
      dup.sig = thre_dup[as.character(Nhet_all)]
    )
    info$dup.stat <- ifelse(info$lhr.dup > info$dup.sig, "duplicated", "non-duplicated")
  }

  info
}

#' Allele Information with GC-Bias Correction for WGS Data
#'
#' Detects per-SNP duplication (and optionally deletion) signals in
#' whole-genome sequencing (WGS) data using a joint genotype + sequencing-depth
#' likelihood-ratio framework, with GC-bias correction applied at the
#' single-SNP level. Sequencing depth is modelled separately within each
#' 1-% GC bin using SNP-centred windows, so the GC and depth signals share the
#' same spatial footprint (avoiding the geometry mismatch that arises when GC
#' is estimated from sliding windows but depth from SNP positions).
#'
#' @section Pipeline overview:
#' 1. Optional region filter (\code{chroms} or \code{bed}).
#' 2. GC content estimated for a \code{gc_win_size}-bp window centred on
#'    each retained SNP, then binned at 1-\% resolution.
#' 3. Global inbreeding coefficient (\eqn{F_{is}}), negative-binomial
#'    depth parameters, and beta-binomial allelic-ratio parameters are
#'    fit on the full filtered set so per-chromosome likelihoods remain
#'    directly comparable.
#' 4. Per-chromosome loop: depth mean/SD is re-fit within each 1-\% GC
#'    bin, depth-based and allelic-ratio likelihoods are computed for
#'    copy numbers N = 2, N = 4 (and optionally N = 1).
#' 5. Likelihood ratios for duplication (and deletion) are computed
#'    per SNP; if \code{perm_test = TRUE}, significance thresholds are
#'    estimated from a permutation pool drawn across all chromosomes.
#' 6. All per-chromosome results are merged into one \code{data.frame}.
#'
#' @param vcf VCF input. Either (a) a file path (character) that will be read
#'   via \code{readVCF()}, (b) a list returned by \code{readVCF()} with a
#'   \code{$vcf} element, or (c) a \code{data.frame} whose first column holds
#'   chromosome names.
#' @param fis Numeric scalar, global inbreeding coefficient. If \code{NULL}
#'   (default), it is estimated internally with \code{h.zygosity()} and
#'   reported via message.
#' @param parallel Logical. If \code{TRUE}, parallelise the per-sample and
#'   permutation loops with a \code{parallel} cluster.
#' @param ref Character path to the reference FASTA (required). Must be
#'   indexable by \code{Rsamtools::indexFa()}; the index will be created on
#'   the fly if missing.
#' @param gc_win_size Integer window size in bp used to compute GC content
#'   around each SNP (default \code{1000}). The window is centred on
#'   \code{POS}: \code{[POS - floor((win-1)/2), POS + floor(win/2)]}.
#' @param numCores Number of cores to use when \code{parallel = TRUE}.
#'   Defaults to \code{parallel::detectCores() - 1}.
#' @param nb_para Optional pre-computed negative-binomial parameter table
#'   (one row per sample, with columns \code{mean}, \code{size}, \code{mean2},
#'   \code{size2}). If \code{NULL}, the table is fit internally.
#' @param sig Numeric in \code{(0, 1)}, quantile threshold for the
#'   permutation-based significance cut-off (default \code{0.95}).
#' @param deletion Logical. If \code{TRUE}, also test for deletions
#'   (N = 1 vs. N = 2) in addition to duplications (N = 4 vs. N = 2).
#' @param perm_test Logical. If \code{TRUE} (default), estimate per-Nhet
#'   significance thresholds by permutation; if \code{FALSE}, the function
#'   returns raw likelihood ratios only (no \code{*.sig} or \code{*.stat}
#'   columns).
#' @param chroms Optional character vector of chromosome names to process
#'   (e.g. \code{c("chr1", "chr2")}). Mutually exclusive with \code{bed}.
#' @param bed Optional path to a BED file. Only SNPs overlapping at least one
#'   BED interval are processed. Mutually exclusive with \code{chroms}.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return A single `data.frame` containing one row per processed SNP
#'   (merged across all chromosomes / regions). Columns include:
#'
#'   - **`CHROM, POS, win_start, win_end, win_width, win_clipped`**: SNP
#'     coordinates and the GC window actually used (clipped at chromosome
#'     boundaries).
#'   - **`gc_percent, nonN_fraction, gc_idx`**: GC content of the window,
#'     fraction of non-`N` bases, and the 1-% GC bin.
#'   - **`mdep, Nhet, Nmiss, p_alt, delta`**: Mean depth across samples,
#'     heterozygous-genotype count, missing-genotype count, alt-allele
#'     frequency, and the deviation from expected Nhet given \eqn{F_{is}}.
#'   - **`lhr.dup, dup.sig, dup.stat`**: Duplication likelihood ratio,
#'     permutation threshold for the SNP's Nhet, and classification
#'     (`"duplicated"` / `"non-duplicated"`).
#'   - **`lhr.del, del.sig, del.stat`**: (only when `deletion = TRUE`)
#'     deletion analogues of the above.
#'
#'   When `perm_test = FALSE`, the `*.sig` and `*.stat` columns are omitted
#'   and only the raw `lhr.*` columns are returned.
#'
#' @section Region-filter behaviour:
#'   \code{chroms} and \code{bed} cannot be supplied together. With no
#'   filter, every SNP in the VCF is processed. With \code{chroms}, SNPs whose
#'   \code{CHROM} is in the vector are kept (missing chromosomes generate a
#'   warning). With \code{bed}, SNPs overlapping any interval are kept.
#'   In all cases, the per-chromosome loop iterates over the chromosomes that
#'   actually appear in the retained SNPs, and outputs are merged before
#'   return.
#'
#' @section Statistical caveat:
#'   The per-chromosome GC-depth model uses only SNPs on that chromosome. If a
#'   chromosome carries fewer than ~500 SNPs, the 1-\% GC bins become sparse
#'   and the fitted depth mean/SD may be unstable; the function emits a
#'   warning in that case but still returns the result.
#'
#' @section Memory:
#'   The dominant memory term is the per-chromosome likelihood matrix
#'   (\code{nSNP_chr} rows by \code{4 * nsample} or \code{6 * nsample}
#'   columns). It is allocated and freed inside each iteration of the
#'   chromosome loop, so peak memory scales with the largest chromosome, not
#'   the full genome.
#'
#' @seealso \code{\link{allele.info.WGS3}} (non-GC-corrected variant using
#'   the negative-binomial depth model directly),
#'   \code{run_gc_estimation}, \code{parse_region_filter} (internal helpers).
#'
#' @examples
#' \dontrun{
#' # Whole genome with default 1 kb GC window
#' res <- allele.info.WGS.GCcor(
#'   vcf = "sample.vcf.gz",
#'   ref = "genome.fa"
#' )
#'
#' # Restrict to a few chromosomes, also test deletions
#' res <- allele.info.WGS.GCcor(
#'   vcf      = "sample.vcf.gz",
#'   ref      = "genome.fa",
#'   chroms   = c("chr1", "chr2", "chr3"),
#'   deletion = TRUE
#' )
#'
#' # Restrict to a BED file of candidate regions, in parallel
#' res <- allele.info.WGS.GCcor(
#'   vcf       = "sample.vcf.gz",
#'   ref       = "genome.fa",
#'   bed       = "candidate_regions.bed",
#'   parallel  = TRUE,
#'   numCores  = 8,
#'   perm_test = TRUE
#' )
#' }
#'
#' @export
allele.info.WGS.GCcor <- function(vcf        = NULL,
                                  fis        = NULL,
                                  parallel   = FALSE,
                                  ref        = NULL,
                                  gc_win_size = 1000,
                                  numCores   = NULL,
                                  nb_para    = NULL,
                                  sig        = 0.95,
                                  deletion   = FALSE,
                                  perm_test  = TRUE,
                                  chroms     = NULL,  # character vector of chr names
                                  bed        = NULL,  # path to BED file
                                  ...) {

  # -- 0. Parallel setup -------------------------------------------------------
  cl <- NULL
  if (parallel) {
    if (is.null(numCores)) numCores <- parallel::detectCores() - 1L
    cl <- parallel::makeCluster(numCores)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterEvalQ(cl, {
      requireNamespace("stringi",    quietly = TRUE)
      requireNamespace("extraDistr", quietly = TRUE)
    })
  }

  pb_apply <- if (parallel) {
    function(X, MARGIN, FUN, ...) pbapply::pbapply(X, MARGIN, FUN, ..., cl = cl)
  } else {
    function(X, MARGIN, FUN, ...) pbapply::pbapply(X, MARGIN, FUN, ...)
  }
  pb_lapply <- if (parallel) {
    function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl)
  } else {
    function(X, FUN, ...) pbapply::pblapply(X, FUN, ...)
  }

  # -- 1. Parse VCF ------------------------------------------------------------
  if (is.null(ref) || !file.exists(ref)) {
    stop("Please provide a valid path to the reference FASTA.")
  }

  if (is.character(vcf)) {
    vcf <- readVCF(vcf)
    vcf <- vcf$vcf
  }else {
    if (!is.data.frame(vcf)) {
      vcf <- vcf$vcf
    }
  }

  gt <- hetTgen(vcf, "GT")
  ad <- hetTgen(vcf, "AD")
  ad[ad == ""] <- "0,0"

  # -- 2. Region filter --------------------------------------------------------
  # Returns a logical mask over the original VCF rows
  row_mask <- parse_region_filter(ad[, 1:2], chroms = chroms, bed = bed)

  # Subset everything to the filtered rows
  vcf <- vcf[row_mask,  , drop = FALSE]
  gt  <- gt[row_mask,   , drop = FALSE]
  ad  <- ad[row_mask,   , drop = FALSE]
  rm(row_mask); gc()

  # Build per-chromosome chunks as indices into the now-contiguous filtered rows
  chr_col    <- as.character(ad[, 1])
  chr_chunks <- split(seq_len(nrow(ad)), chr_col)
  chr_names  <- names(chr_chunks)
  n_chr      <- length(chr_names)
  message(sprintf("Will process %d chromosome(s): %s",
                  n_chr, paste(chr_names, collapse = ", ")))

  # -- 3. GC annotation (all filtered SNPs at once) ----------------------------
  message("Annotating SNPs with GC content (SNP-centred windows)...")
  snp_annotated <- run_gc_estimation(
    snps_df  = ad[, 1:2],
    ref      = ref,
    win_size = gc_win_size,
    ignore_N = TRUE
  )

  # Bin GC into 1-% factors for the per-chromosome GC-depth model.
  # Done here (rather than inside run_gc_estimation) so that GC estimation
  # stays a pure measurement step and binning lives next to its consumer.
  if (all(is.na(snp_annotated$gc_percent))) {
    warning("All GC values are NA; gc_idx will be entirely NA.")
    snp_annotated$gc_idx <- factor(rep(NA_character_, nrow(snp_annotated)))
  } else {
    breaks <- seq(
      floor(min(snp_annotated$gc_percent,   na.rm = TRUE) - 0.5),
      ceiling(max(snp_annotated$gc_percent, na.rm = TRUE) + 0.5),
      by = 1
    )
    snp_annotated$gc_idx <- cut(snp_annotated$gc_percent,
                                breaks = breaks, include.lowest = TRUE)
  }

  # -- 4. Global Fis -----------------------------------------------------------
  if (is.null(fis)) {
    fis_est <- h.zygosity(vcf, cl = if (parallel) cl else NULL)
    fis     <- round(mean(fis_est$Fis, na.rm = TRUE), 2)
  }
  message(sprintf("Estimated global inbreeding coefficient: Fis = %.2f", fis))
  rm(vcf); gc()

  AD      <- ad[, -(1:4), drop = FALSE]
  GT_raw  <- gt[, -(1:4), drop = FALSE]
  rm(ad, gt); gc()

  nsample    <- ncol(AD)
  nrow_filt  <- nrow(AD)
  idx_sample <- seq_len(nsample)

  # Recode GT once globally
  GT_full <- matrix(-9L, nrow_filt, nsample)
  GT_full[GT_raw == "0/0" | GT_raw == "0|0"]                       <- 0L
  GT_full[GT_raw == "0/1" | GT_raw == "0|1" | GT_raw == "1|0"]    <- 1L
  GT_full[GT_raw == "1/1" | GT_raw == "1|1"]                       <- 2L
  rm(GT_raw); gc()

  # -- 5. Global NB + beta-binomial parameter estimation -----------------------
  #   These must be estimated from all filtered data to stay consistent across
  #   chromosomes.  The full tot matrix is freed immediately after.
  message("step 1/4: global parameter estimation")
  message("  Extracting depth of coverage (all filtered SNPs)...")

  tmp_g      <- pb_apply(AD, 2, ad_parse, mode = "ref_alt")
  ref_g      <- matrix(0L, nrow_filt, nsample)
  alt_g      <- matrix(0L, nrow_filt, nsample)
  for (j in seq_along(tmp_g)) {
    ref_g[, j] <- tmp_g[[j]]$ref
    alt_g[, j] <- tmp_g[[j]]$alt
  }
  rm(tmp_g); gc()

  tot_g <- ref_g + alt_g
  colnames(tot_g) <- colnames(AD)
  rm(alt_g); gc()

  # NB parameters
  if (is.null(nb_para)) {
    message("  Fitting NB depth model...")
    colstat       <- pb_apply(tot_g, 2, db_mod)
    colstat       <- do.call(rbind, colstat)
    colstat$mean2 <- 2 * colstat$mean
    colstat$size2 <- colstat$size
    nb_stat       <- colstat
  } else {
    nb_stat <- nb_para
  }

  # Beta-binomial parameters
  message("  Fitting beta-binomial model...")
  H_full      <- (GT_full == 1L)
  tot.rmOutl  <- pb_apply(AD, 2, ad_parse, mode = "tot_inlier", outlier_sd = 1.5)
  het_inl_idx <- which(H_full & tot.rmOutl, arr.ind = TRUE)
  rm(tot.rmOutl); gc()

  ad_het_vec <- AD[het_inl_idx]
  if (length(ad_het_vec) > 1e5)
    ad_het_vec <- ad_het_vec[sample(seq_along(ad_het_vec), 1e5)]

  alle_ad_het <- do.call(rbind, pb_lapply(ad_het_vec, ad_parse, mode = "ref_tot"))
  rm(ad_het_vec); gc()

  fit_bb  <- betabinom_mle(alle_ad_het$ref, alle_ad_het$tot)
  bb_alpha <- round(fit_bb$estimates[1], 2)
  bb_beta  <- round(fit_bb$estimates[2], 2)
  rho      <- 1 / (bb_alpha + bb_beta + 1)

  use_bb <- (rho > 1e-4)
  message(sprintf(
    "  Beta-binomial: alpha = %.2f, beta = %.2f%s",
    bb_alpha, bb_beta,
    if (!use_bb) "  [rho < 1e-4 -> switching to binomial]" else ""
  ))
  p.cor <- if (!use_bb) mean(alle_ad_het$ref / alle_ad_het$tot, na.rm = TRUE) else NULL
  rm(alle_ad_het); gc()

  # Free global depth matrices - no longer needed
  rm(ref_g, tot_g); gc()

  # -- 6. Per-chromosome loop ---------------------------------------------------
  message(sprintf("step 2/4: per-chromosome likelihood computation (%d chromosome(s))",
                  n_chr))

  results_list      <- vector("list", n_chr)
  names(results_list) <- chr_names

  # Accumulators for permutation sampling (capped at 5k per chromosome)
  perm_acc <- list(
    geno_dup = numeric(0),
    dep_dup  = numeric(0)
  )
  if (deletion) {
    perm_acc$geno_del <- numeric(0)
    perm_acc$dep_del  <- numeric(0)
  }

  for (ci in seq_along(chr_names)) {
    chr       <- chr_names[ci]
    idx       <- chr_chunks[[chr]]
    n_snp_chr <- length(idx)

    message(sprintf("  [%d/%d] %s  (%d SNPs)", ci, n_chr, chr, n_snp_chr))

    # Warn if this chromosome is too sparse for reliable per-bin GC modelling
    if (n_snp_chr < 500L)
      warning(sprintf(
        "Chromosome %s has only %d SNPs; GC-bin model may be unreliable.",
        chr, n_snp_chr
      ))

    # Slice chromosome-specific matrices
    AD_chr      <- AD[idx,      , drop = FALSE]
    GT_chr      <- GT_full[idx, , drop = FALSE]
    snp_ann_chr <- snp_annotated[idx, , drop = FALSE]

    H_chr     <- (GT_chr == 1L)
    Nhet_chr  <- rowSums(H_chr)
    Nmiss_chr <- rowSums(GT_chr == -9L)

    # Extract ref / alt counts
    tmp_chr <- pb_apply(AD_chr, 2, ad_parse, mode = "ref_alt")
    ref_chr <- matrix(0L, n_snp_chr, nsample)
    alt_chr <- matrix(0L, n_snp_chr, nsample)
    for (j in seq_along(tmp_chr)) {
      ref_chr[, j] <- tmp_chr[[j]]$ref
      alt_chr[, j] <- tmp_chr[[j]]$alt
    }
    rm(tmp_chr); gc()

    tot_chr <- ref_chr + alt_chr
    colnames(tot_chr) <- colnames(AD_chr)
    rm(alt_chr); gc()

    # Per-chromosome GC-depth model
    # (uses SNP-centred GC bins for this chromosome only)
    gc_fit_chr <- lapply(seq_len(nsample), function(j) {
      dep_vec <- tot_chr[, j]
      m <- do.call(rbind, tapply(dep_vec, snp_ann_chr$gc_idx, function(x) {
        q05 <- quantile(x, 0.05, na.rm = TRUE)
        q95 <- quantile(x, 0.95, na.rm = TRUE)
        inl <- x[x > q05 & x < q95]
        data.frame(mean = mean(inl, na.rm = TRUE),
                   sd   = sd(inl,   na.rm = TRUE))
      }))
      m$gc_idx <- rownames(m)
      m
    })

    # ad.info for this chromosome
    p_alt_chr  <- round(
      rowMeans(replace(GT_chr, GT_chr == -9L, NA_integer_), na.rm = TRUE) / 2, 2
    )
    E.Nhet_chr <- 2 * p_alt_chr * (1 - p_alt_chr) * (1 - fis) * (nsample - Nmiss_chr)

    ad_info_chr <- data.frame(
      snp_ann_chr,
      mdep  = rowMeans(tot_chr),
      Nhet  = Nhet_chr,
      Nmiss = Nmiss_chr,
      p_alt = p_alt_chr,
      delta = round(E.Nhet_chr - Nhet_chr, 2),
      stringsAsFactors = FALSE
    )

    # Depth-based likelihoods (pre-allocated, single loop over samples)
    dep_lld2 <- matrix(NA_real_, nsample, n_snp_chr)
    dep_lld4 <- matrix(NA_real_, nsample, n_snp_chr)
    if (deletion) dep_lld1 <- matrix(NA_real_, nsample, n_snp_chr)

    for (j in seq_len(nsample)) {
      gf      <- gc_fit_chr[[j]]
      mi      <- match(snp_ann_chr$gc_idx, gf$gc_idx)
      gc_mean <- gf$mean[mi]
      gc_sd   <- gf$sd[mi]
      tot_j   <- tot_chr[, j]

      dep_lld2[j, ] <- pnorm(tot_j + 0.5, gc_mean,     gc_sd) -
        pnorm(tot_j - 0.5, gc_mean,     gc_sd)
      dep_lld4[j, ] <- pnorm(tot_j + 0.5, gc_mean * 2, gc_sd) -
        pnorm(tot_j - 0.5, gc_mean * 2, gc_sd)
      if (deletion)
        dep_lld1[j, ] <- pnorm(tot_j + 0.5, gc_mean / 2, gc_sd) -
        pnorm(tot_j - 0.5, gc_mean / 2, gc_sd)
    }
    rm(gc_fit_chr); gc()

    dep_lld4[dep_lld2 == 0] <- 1e-50
    dep_lld2[dep_lld2 == 0] <- 1e-100
    if (deletion) dep_lld1[dep_lld2 == 0] <- 1e-150

    # Genotype + depth likelihoods
    lld_chr <- if (use_bb) {
      pb_apply(cbind(ref_chr, tot_chr, GT_chr), 1, cal_lld_v2.1,
               inb = fis, nb_stat = nb_stat,
               alpha = bb_alpha, beta = bb_beta, deletion = deletion)
    } else {
      pb_apply(cbind(ref_chr, tot_chr, GT_chr), 1, cal_lld_v2.1,
               inb = fis, nb_stat = nb_stat,
               p.cor = p.cor, deletion = deletion)
    }
    rm(ref_chr, tot_chr, GT_chr); gc()

    lld_chr <- if (deletion) {
      t(rbind(lld_chr, dep_lld2, dep_lld1, dep_lld4))
    } else {
      t(rbind(lld_chr, dep_lld2, dep_lld4))
    }
    rm(dep_lld2, dep_lld4)
    if (deletion) rm(dep_lld1)
    gc()

    # -- Compute LHR + accumulate permutation samples -------------------------
    eps <- 1e-300

    if (deletion) {
      gp_dup <- log(pmax(lld_chr[, idx_sample + 2*nsample], eps)) -
        log(pmax(lld_chr[, idx_sample            ], eps))
      dp_dup <- log(pmax(lld_chr[, idx_sample + 5*nsample], eps)) -
        log(pmax(lld_chr[, idx_sample + 3*nsample], eps))
      lhr_dup <- 2 * rowSums(pmax(gp_dup + dp_dup, 0))

      gp_del <- log(pmax(lld_chr[, idx_sample + 1*nsample], eps)) -
        log(pmax(lld_chr[, idx_sample            ], eps))
      dp_del <- log(pmax(lld_chr[, idx_sample + 4*nsample], eps)) -
        log(pmax(lld_chr[, idx_sample + 3*nsample], eps))
      lhr_del <- 2 * rowSums(pmax(gp_del + dp_del, 0))

      if (perm_test) {
        n_acc <- min(5000L, sum(!is.na(gp_dup[H_chr])))
        perm_acc$geno_dup <- c(perm_acc$geno_dup,
                               sample(na.omit(gp_dup[H_chr]), n_acc))
        perm_acc$dep_dup  <- c(perm_acc$dep_dup,
                               sample_from_mat(dp_dup, min(5000L, length(dp_dup))))
        n_acc <- min(5000L, sum(!is.na(gp_del[H_chr])))
        perm_acc$geno_del <- c(perm_acc$geno_del,
                               sample(na.omit(gp_del[H_chr]), n_acc))
        perm_acc$dep_del  <- c(perm_acc$dep_del,
                               sample_from_mat(dp_del, min(5000L, length(dp_del))))
      }

      results_list[[chr]] <- list(
        ad_info = ad_info_chr,
        lhr_dup = lhr_dup,
        lhr_del = lhr_del,
        Nhet    = Nhet_chr
      )

    } else {
      gp_dup <- log(pmax(lld_chr[, idx_sample + 1*nsample], eps)) -
        log(pmax(lld_chr[, idx_sample            ], eps))
      dp_dup <- log(pmax(lld_chr[, idx_sample + 3*nsample], eps)) -
        log(pmax(lld_chr[, idx_sample + 2*nsample], eps))
      lhr_dup <- 2 * rowSums(pmax(gp_dup + dp_dup, 0))

      if (perm_test) {
        n_acc <- min(5000L, sum(!is.na(gp_dup[H_chr])))
        perm_acc$geno_dup <- c(perm_acc$geno_dup,
                               sample(na.omit(gp_dup[H_chr]), n_acc))
        perm_acc$dep_dup  <- c(perm_acc$dep_dup,
                               sample_from_mat(dp_dup, min(5000L, length(dp_dup))))
      }

      results_list[[chr]] <- list(
        ad_info = ad_info_chr,
        lhr_dup = lhr_dup,
        Nhet    = Nhet_chr
      )
    }

    rm(lld_chr, H_chr, ad_info_chr); gc()
    message(sprintf("    %s complete.", chr))
  }

  rm(AD, GT_full, H_full, snp_annotated); gc()

  # -- 7. Merge per-chromosome results -----------------------------------------
  message("step 4/4: merging results")

  info_all    <- do.call(rbind, lapply(results_list, `[[`, "ad_info"))
  lhr_dup_all <- do.call(c,     lapply(results_list, `[[`, "lhr_dup"))
  Nhet_all    <- do.call(c,     lapply(results_list, `[[`, "Nhet"))
  if (deletion)
    lhr_del_all <- do.call(c,   lapply(results_list, `[[`, "lhr_del"))
  rm(results_list); gc()

  # -- 8. Early return (no permutation test) -----------------------------------
  if (!perm_test) {
    info <- if (deletion) {
      data.frame(info_all, lhr.dup = lhr_dup_all, lhr.del = lhr_del_all)
    } else {
      data.frame(info_all, lhr.dup = lhr_dup_all)
    }
    return(info)
  }

  # -- 9. Global permutation thresholds (from cross-chromosome sample pool) ----
  message("step 3/4: estimating permutation thresholds (global pool)")

  run_perm_thresh <- function(geno_acc, dep_acc, label) {
    message(sprintf("  %s...", label))
    geno_s <- sample(geno_acc, min(1e5L, length(geno_acc)))
    dep_s  <- sample(dep_acc,  min(1e5L, length(dep_acc)))
    thre <- unlist(pb_lapply(
      unique(Nhet_all), sample_and_quantile2,
      nsample = nsample, geno = geno_s, dep = dep_s,
      nrep = 10000L, quan = sig
    ))
    names(thre) <- unique(Nhet_all)
    thre
  }

  thre_dup <- run_perm_thresh(perm_acc$geno_dup, perm_acc$dep_dup,
                              "duplication thresholds")
  if (deletion)
    thre_del <- run_perm_thresh(perm_acc$geno_del, perm_acc$dep_del,
                                "deletion thresholds")
  rm(perm_acc); gc()

  # -- 10. Assemble final output ------------------------------------------------
  if (deletion) {
    info <- data.frame(
      info_all,
      lhr.dup  = lhr_dup_all,
      dup.sig  = thre_dup[as.character(Nhet_all)],
      lhr.del  = lhr_del_all,
      del.sig  = thre_del[as.character(Nhet_all)]
    )
    info$dup.stat <- ifelse(info$lhr.dup > info$dup.sig, "duplicated",  "non-duplicated")
    info$del.stat <- ifelse(info$lhr.del > info$del.sig, "deleted",     "non-deleted")
  } else {
    info <- data.frame(
      info_all,
      lhr.dup = lhr_dup_all,
      dup.sig = thre_dup[as.character(Nhet_all)]
    )
    info$dup.stat <- ifelse(info$lhr.dup > info$dup.sig, "duplicated", "non-duplicated")
  }

  info
}



sample_and_quantile2 <- function(nhet, nsample, geno, dep, nrep, quan = 0.95) {
  tem <- replicate(nrep, {
    lld <- sample(dep, nsample, replace = TRUE)
    if (nhet > 0) {
      pos <- sample.int(nsample, nhet)
      lld[pos] <- lld[pos] + sample(geno, nhet, replace = TRUE)
    }
    sum(pmax(lld, 0)) * 2
  })
  unname(stats::quantile(tem, quan,na.rm=T))
}


# sample directly from matrix entries (rows x cols) without flattening:
sample_from_mat <- function(M, n) {
  nr <- nrow(M); nc <- ncol(M); n <- min(n, nr*nc)
  r <- sample.int(nr, n, replace = TRUE)
  c <- sample.int(nc, n, replace = TRUE)
  M[cbind(r, c)]
}


cal_lld_v2.1 <- function(x,inb,nb_stat,alpha = NULL,beta = NULL,p.cor = NULL,deletion)
{

  # ------------------------------------------------------------------
  # 1. Split input row vector
  # ------------------------------------------------------------------
  nsample <- nrow(nb_stat)
  ref <- as.numeric(x[1:nsample])
  tot <- as.numeric(x[(nsample + 1):(nsample*2)])
  gt <- x[(2*nsample + 1):(nsample*3)]

  ####### allelic ratio based lld
  geno <- gt

  # ------------------------------------------------------------------
  # 2. Compute allele frequencies for N=2 quickly
  # ------------------------------------------------------------------

  p.alt2 <- sum(geno[geno != -9])/2/sum(geno != -9)
  p.ref2 <- 1 - p.alt2

  # Genotype priors under inbreeding F = inb
  pG.alt.hom <- p.alt2^2*(1-inb) + p.alt2*inb
  pG.ref.hom <- p.ref2^2*(1-inb) + p.ref2*inb
  pG.het <- 2*p.alt2*p.ref2*(1-inb)


  if(deletion){
    #### N=1
    p.alt1 <- p.alt2
    p.ref1 <- 1 - p.alt1
    geno_lld1 <- rep(1,nrow(nb_stat))
  }


  #### N=4
  ref_cor <- ref/nb_stat$mean
  tot_cor <- tot/nb_stat$mean

  p.ref4 <- sum(ref_cor)/sum(tot_cor)
  p.alt4 <- 1- p.ref4


  # ------------------------------------------------------------------
  # 3. Depth likelihoods (vectorized)
  # ------------------------------------------------------------------
  geno_lld2 <- rep(1,nrow(nb_stat))
  geno_lld4 <- rep(1,nrow(nb_stat))


  # ------------------------------------------------------------------
  # 4. Allelic ratio likelihoods
  # ------------------------------------------------------------------
  het.idx <- which(geno == 1)
  geno_lld2 <- rep(1, nsample)
  geno_lld4 <- rep(1, nsample)
  if (deletion) geno_lld1 <- rep(1, nsample)

  if (length(het.idx)) {

    r <- ref[het.idx]
    n <- tot[het.idx]

    if (is.null(p.cor)) {
      # Beta-binomial
      ab     <- alpha + beta
      a_err  <- ab * 0.999
      b_err  <- ab * 0.001

      l_ref_hom <- dbbinom(r, n, a_err, b_err)
      l_alt_hom <- dbbinom(r, n, b_err, a_err)
      l_het     <- dbbinom(r, n, alpha, beta)

      l_Aaaa <- dbbinom(r,n,
                        alpha/(alpha+3*beta)*(alpha+beta),
                        3*beta/(alpha+3*beta)*(alpha+beta))
      l_AAaa <- dbbinom(r,n, alpha, beta)
      l_AAAa <- dbbinom(r,n,
                        3*alpha/(3*alpha+beta)*(alpha+beta),
                        beta/(3*alpha+beta)*(alpha+beta))
    } else {
      # Binomial
      l_ref_hom <- dbinom(r, n, 0.99)
      l_alt_hom <- dbinom(r, n, 0.01)
      l_het     <- dbinom(r, n, p.cor)

      l_Aaaa <- dbinom(r,n, p.cor/(p.cor+3*(1-p.cor)))
      l_AAaa <- dbinom(r,n, p.cor)
      l_AAAa <- dbinom(r,n, 3*p.cor/(3*p.cor+(1-p.cor)))
    }

    # --- N=2 posterior ---
    u_ref <- pG.ref.hom * l_ref_hom
    u_alt <- pG.alt.hom * l_alt_hom
    u_het <- pG.het     * l_het
    Z     <- u_ref + u_alt + u_het

    geno_lld2[het.idx] <- (l_ref_hom*(u_ref/Z) +
                             l_het*(u_het/Z) +
                             l_alt_hom*(u_alt/Z))

    # --- N=1 posterior (if deletion) ---
    if (deletion) {
      u_a <- p.alt1 * l_alt_hom
      u_A <- p.ref1 * l_ref_hom
      Z1  <- u_a + u_A
      geno_lld1[het.idx] <- (l_alt_hom*(u_a/Z1) +
                               l_ref_hom*(u_A/Z1))
    }

    # --- N=4 prior weights ---
    pG.aaaa <- p.alt4^4*(1-inb) + p.alt4*inb
    pG.Aaaa <- 4*p.ref4*p.alt4^3*(1-inb)
    pG.AAaa <- 6*p.ref4^2*p.alt4^2*(1-inb)
    pG.AAAa <- 4*p.ref4^3*p.alt4*(1-inb)
    pG.AAAA <- p.ref4^4*(1-inb) + p.ref4*inb

    l_aaaa <- l_alt_hom
    l_AAAA <- l_ref_hom

    U4 <- (pG.AAAA*l_AAAA +
             pG.AAAa*l_AAAa +
             pG.AAaa*l_AAaa +
             pG.Aaaa*l_Aaaa +
             pG.aaaa*l_aaaa)

    geno_lld4[het.idx] <-
      (l_AAAA*(pG.AAAA*l_AAAA/U4) +
         l_AAAa*(pG.AAAa*l_AAAa/U4) +
         l_AAaa*(pG.AAaa*l_AAaa/U4) +
         l_Aaaa*(pG.Aaaa*l_Aaaa/U4) +
         l_aaaa*(pG.aaaa*l_aaaa/U4))
  }
  if (deletion) {
    return(c(geno_lld2,geno_lld1,geno_lld4))
  }
  else {
    return(c(geno_lld2,geno_lld4))
  }
}


ad_parse <- function(x, mode = c("ref_tot", "ref_alt", "tot", "tot_inlier"),
                     outlier_sd = 1.5) {
  mode <- match.arg(mode)

  if (mode == "ref_alt") {
    # vectorized: returns list(ref=int, alt=int)
    tok <- stringi::stri_split_fixed(x, ",", n = 2, simplify = TRUE)
    return(list(
      ref = as.integer(tok[, 1]),
      alt = as.integer(tok[, 2])
    ))
  }

  # everything else works for a vector OR a single string
  if (length(x) == 1L) {
    ad <- as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
    ref <- ad[1]
    alt <- ad[2]
    tot <- ref + alt

    if (mode == "ref_tot") return(data.frame(ref = ref, tot = tot))
    if (mode == "tot")     return(tot)
    if (mode == "tot_inlier") stop("mode='tot_inlier' expects a vector (e.g. column).")
  }

  # vector input
  tot <- vapply(x, function(z) {
    ad <- as.numeric(strsplit(z, ",", fixed = TRUE)[[1]])
    sum(ad)
  }, numeric(1))

  if (mode == "tot") return(tot)

  if (mode == "tot_inlier") {
    return(abs(scale(tot)) < outlier_sd)
  }

  # mode == "ref_tot" for vector input: return data.frame row-bindable
  out <- do.call(rbind, lapply(x, function(z) {
    ad <- as.numeric(strsplit(z, ",", fixed = TRUE)[[1]])
    data.frame(ref = ad[1], tot = sum(ad))
  }))
  rownames(out) <- NULL
  out
}

###### estimate negative binomial distribution parameters
db_mod <- function(x){
  q <- quantile(x, c(0.025, 0.975), na.rm=TRUE)
  xx <- x[x > q[1] & x < q[2]]

  m <- mean(xx)
  s <- sd(xx)
  size <- fitdistrplus::fitdist(xx, "nbinom")$estimate["size"]
  p <- size/(size + m)
  data.frame(mean = m, sd = s, size = size, p = p)
}



## Beta-Binomial MLE with variable n_i
## -----------------------------------
## Inputs:
##   x : integer vector of successes
##   n : integer vector of trials (same length as x)
##   start: optional named list with alpha,beta; if NULL, uses MoM start
##   fix_mu: optional number in (0,1). If provided, fit with mean fixed (estimate dispersion only)
## Returns:
##   list with estimates (alpha,beta), SEs, vcov, logLik, convergence code, etc.

betabinom_mle <- function(x, n, start = NULL, fix_mu = NULL) {
  stopifnot(length(x) == length(n), all(n >= 0), all(x >= 0), all(x <= n))
  x <- as.numeric(x); n <- as.numeric(n)
  m <- length(x)

  # ---- helper: stable log Beta-Binomial log-likelihood contribution
  # log f(x | n, alpha, beta) = lchoose(n,x) + lbeta(x+alpha, n-x+beta) - lbeta(alpha, beta)
  ll_contrib <- function(x, n, alpha, beta) {
    lchoose(n, x) + lbeta(x + alpha, n - x + beta) - lbeta(alpha, beta)
  }

  # ---- Method-of-moments starting values (robust, works for variable n)
  if (is.null(start)) {
    mu_hat <- sum(x) / sum(n)  # pooled proportion
    p_i <- ifelse(n > 0, x / n, 0)
    a <- mean(ifelse(n > 0, 1 / n, 0))
    s2 <- var(p_i)             # sample variance of per-locus proportions
    s2 <- ifelse(is.na(s2), 1e-6, s2)

    # Var(p_i) ~ mu(1-mu) * ( a + rho*(1 - a) )
    rho_hat <- (s2 / (mu_hat * (1 - mu_hat)) - a) / (1 - a)
    rho_hat <- max(rho_hat, 0) # no negative ICC
    if (is.infinite(rho_hat) || is.na(rho_hat)) rho_hat <- 0

    phi_hat <- if (rho_hat == 0) 1e4 else (1 / rho_hat - 1)  # precision = alpha+beta
    alpha0 <- max(mu_hat * phi_hat, 1e-3)
    beta0  <- max((1 - mu_hat) * phi_hat, 1e-3)

    start <- list(alpha = alpha0, beta = beta0)
  } else {
    if (is.null(start$alpha) || is.null(start$beta)) stop("start must have alpha and beta")
  }

  # ---- Two modes: (1) free alpha,beta ; (2) fix mean and estimate precision only
  if (is.null(fix_mu)) {
    # Unconstrained optimization via log-transform to keep alpha,beta > 0
    par0 <- c(log_alpha = log(start$alpha), log_beta = log(start$beta))

    nll <- function(par) {
      alpha <- exp(par[1]); beta <- exp(par[2])
      # sum negative log-likelihood
      -sum(ll_contrib(x, n, alpha, beta))
    }

    opt <- optim(par0, nll, method = "BFGS", hessian = TRUE)
    alpha_hat <- exp(opt$par[1]); beta_hat <- exp(opt$par[2])

    # Delta method for SEs: Var(theta) ~ J * Var(log-theta) * J^T with J=diag(theta)
    vcov_log <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA, 2, 2))
    J <- diag(c(alpha_hat, beta_hat))
    vcov_ab <- J %*% vcov_log %*% J
    se_alpha <- sqrt(pmax(diag(vcov_ab)[1], 0))
    se_beta  <- sqrt(pmax(diag(vcov_ab)[2], 0))

  } else {
    # Fix mean mu, estimate only precision phi => alpha=mu*phi, beta=(1-mu)*phi
    mu <- fix_mu
    stopifnot(mu > 0, mu < 1)

    par0 <- log(start$alpha + start$beta)  # log-phi

    nll_phi <- function(log_phi) {
      phi <- exp(log_phi)
      alpha <- mu * phi
      beta  <- (1 - mu) * phi
      -sum(ll_contrib(x, n, alpha, beta))
    }

    opt <- optim(par = par0, fn = nll_phi, method = "BFGS", hessian = TRUE)
    phi_hat <- exp(opt$par)
    alpha_hat <- mu * phi_hat
    beta_hat  <- (1 - mu) * phi_hat

    # SE for phi (delta from log-phi); then delta to alpha,beta if desired
    var_logphi <- tryCatch(1 / opt$hessian[1, 1], error = function(e) NA_real_)
    se_phi <- if (is.na(var_logphi)) NA_real_ else sqrt(var_logphi) * phi_hat

    # Approx SEs for alpha,beta via delta method (optional)
    se_alpha <- if (is.na(var_logphi)) NA_real_ else mu * se_phi
    se_beta  <- if (is.na(var_logphi)) NA_real_ else (1 - mu) * se_phi

    vcov_ab <- matrix(NA_real_, 2, 2)
  }

  logLik <- -opt$value
  phi_hat <- alpha_hat + beta_hat
  rho_hat <- 1 / (phi_hat + 1)     # intra-class correlation for Beta-Binomial
  mu_hat  <- alpha_hat / (alpha_hat + beta_hat)

  list(
    estimates = c(alpha = alpha_hat, beta = beta_hat, mu = mu_hat, phi = phi_hat, rho = rho_hat),
    se = c(alpha = se_alpha, beta = se_beta),
    vcov = vcov_ab,
    logLik = logLik,
    convergence = opt$convergence,
    message = opt$message
  )
}


# Compute GC content per GRanges on a FASTA
compute_gc <- function(fasta, regions, ignore_N = TRUE) {
  stopifnot(is(regions, "GRanges"))

  fa <- FaFile(fasta)
  if (!file.exists(paste0(fasta, ".fai"))) {
    message("Indexing FASTA (creating .fai)...")
    indexFa(fa)
  }

  # Harmonize seqlevels
  fa_si <- seqinfo(fa)
  common <- intersect(seqlevels(regions), seqlevels(fa_si))
  if (length(common) == 0) stop("No matching seqlevels between regions and FASTA.")
  regions <- keepSeqlevels(regions, common, pruning.mode = "coarse")

  # Pull sequences for all regions (vectorized)
  dna <- getSeq(fa, regions)

  # Count bases
  counts <- letterFrequency(dna, letters = c("A","C","G","T","N"))

  # GC numerator
  gc_num <- counts[, "C"] + counts[, "G"]

  # Denominator: either total A/C/G/T or full width
  if (isTRUE(ignore_N)) {
    denom <- rowSums(counts[, c("A", "C", "G", "T"), drop = FALSE])
  } else {
    denom <- width(dna)
  }

  # Avoid division by zero (e.g., all-N windows). Set GC to NA in such cases.
  gc <- ifelse(denom > 0, gc_num / denom, NA_real_)

  # Save results into a data.frame (and/or GRanges mcols)
  res <- data.frame(
    CHROM = as.character(seqnames(regions)),
    start    = start(regions),
    end      = end(regions),
    width    = width(regions),
    gc_percent = gc * 100,
    nonN_fraction = ifelse(width(dna) > 0,
                           rowSums(counts[, c("A","C","G","T"), drop = FALSE]) / width(dna),
                           NA_real_),
    region_idx = paste0("region",1:length(gc)),
    stringsAsFactors = FALSE
  )

  return(res)
}

# Annotate SNPs with regions



run_gc_estimation <- function(snps_df, ref, win_size = 1000L, ignore_N = TRUE) {

  if (is.null(ref) || !file.exists(ref)) {
    stop("Please provide a valid path to the reference FASTA.")
  }
  stopifnot(
    is.data.frame(snps_df),
    all(c("CHROM", "POS") %in% names(snps_df))
  )

  # Build SNP-centred windows (GRanges)
  win <- snps_to_windows(snps_df, ref, win_size)

  # Compute GC for those windows
  gc_res <- compute_gc(ref, win, ignore_N = ignore_N)

  # Bind back to original SNP coordinates
  data.frame(
    CHROM         = snps_df$CHROM,
    POS           = snps_df$POS,
    win_start     = gc_res$start,
    win_end       = gc_res$end,
    win_width     = gc_res$width,         # < win_size when clipped
    win_clipped   = win$win_clipped,
    gc_percent    = gc_res$gc_percent,
    nonN_fraction = gc_res$nonN_fraction,
    stringsAsFactors = FALSE
  )
}

snps_to_windows <- function(snps_df, ref, win_size) {

  # --- Validate ---------------------------------------------------------------
  stopifnot(
    is.data.frame(snps_df),
    all(c("CHROM", "POS") %in% names(snps_df)),
    nrow(snps_df) > 0L,
    is.character(ref), file.exists(ref),
    is.numeric(win_size), length(win_size) == 1L, win_size >= 1
  )
  win_size <- as.integer(win_size)

  # --- Open FASTA / seqinfo ---------------------------------------------------
  fa <- Rsamtools::FaFile(ref)
  if (!file.exists(paste0(ref, ".fai"))) {
    message("Indexing FASTA (creating .fai)...")
    Rsamtools::indexFa(fa)
  }
  si      <- Rsamtools::seqinfo(fa)
  seqlens <- GenomeInfoDb::seqlengths(si)

  # --- Reconcile chromosomes --------------------------------------------------
  snp_chroms <- unique(as.character(snps_df$CHROM))
  common     <- intersect(snp_chroms, names(seqlens))

  if (length(common) == 0L) {
    stop(paste0(
      "No shared seqlevels between SNP table and FASTA.\n",
      "  SNPs use  : ", paste(head(snp_chroms,      5), collapse = ", "), "\n",
      "  FASTA has : ", paste(head(names(seqlens),  5), collapse = ", "), "\n",
      "Check for chr-prefix mismatch (e.g. 'chr1' vs '1')."
    ))
  }
  if (length(common) < length(snp_chroms)) {
    dropped <- setdiff(snp_chroms, common)
    warning(sprintf(
      "%d SNP chromosome(s) absent from FASTA and dropped: %s",
      length(dropped), paste(dropped, collapse = ", ")
    ))
    snps_df <- snps_df[as.character(snps_df$CHROM) %in% common, , drop = FALSE]
  }

  # --- Build centred windows --------------------------------------------------
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(snps_df$CHROM),
    ranges   = IRanges::IRanges(start = as.integer(snps_df$POS), width = 1L),
    strand   = "*"
  )
  gr <- GenomicRanges::resize(gr, width = win_size, fix = "center")

  # Attach seqlengths so trim() clips correctly
  GenomeInfoDb::seqinfo(gr) <- si[GenomeInfoDb::seqlevels(gr)]
  gr_trimmed <- GenomicRanges::trim(gr)

  # Flag SNPs whose window was clipped by a chromosome boundary
  gr_trimmed$win_clipped <- BiocGenerics::width(gr_trimmed) < win_size
  if (any(gr_trimmed$win_clipped)) {
    message(sprintf(
      "%d SNP window(s) clipped at chromosome boundaries.",
      sum(gr_trimmed$win_clipped)
    ))
  }

  gr_trimmed
}



# ==============================================================================
# Helper: parse_region_filter
#   Converts a chroms vector or BED path into a logical row mask and a named
#   list of per-chromosome row indices (into the *original* snps_df).
#   The main function re-splits after subsetting, so the chunks here are used
#   only to build the mask.
# ==============================================================================
parse_region_filter <- function(snps_df, chroms = NULL, bed = NULL) {

  if (!is.null(chroms) && !is.null(bed))
    stop("`chroms` and `bed` are mutually exclusive; supply only one.")

  stopifnot(
    is.data.frame(snps_df),
    all(c("CHROM", "POS") %in% names(snps_df))
  )

  n       <- nrow(snps_df)
  snp_chr <- as.character(snps_df$CHROM)

  # -- No filter --------------------------------------------------------------
  if (is.null(chroms) && is.null(bed)) {
    message(sprintf("No region filter: processing all %d SNPs on %d chromosome(s).",
                    n, length(unique(snp_chr))))
    return(rep(TRUE, n))
  }

  # -- Chromosome list --------------------------------------------------------
  if (!is.null(chroms)) {
    chroms <- as.character(chroms)
    mask   <- snp_chr %in% chroms

    if (!any(mask))
      stop("None of the requested chromosomes found in VCF: ",
           paste(chroms, collapse = ", "))

    missing_chr <- setdiff(chroms, unique(snp_chr[mask]))
    if (length(missing_chr))
      warning("Chromosomes not found in VCF and skipped: ",
              paste(missing_chr, collapse = ", "))

    message(sprintf(
      "Chromosome filter: retaining %d / %d SNPs on %d / %d requested chromosome(s).",
      sum(mask), n,
      length(unique(snp_chr[mask])), length(chroms)
    ))
    return(mask)
  }

  # -- BED file ---------------------------------------------------------------
  if (!file.exists(bed))
    stop("BED file not found: ", bed)

  bed_gr <- rtracklayer::import(bed, format = "BED")
  BiocGenerics::strand(bed_gr) <- "*"

  snp_gr <- GenomicRanges::GRanges(
    seqnames = snp_chr,
    ranges   = IRanges::IRanges(start = as.integer(snps_df$POS), width = 1L)
  )

  hits <- S4Vectors::queryHits(
    GenomicRanges::findOverlaps(snp_gr, bed_gr, ignore.strand = TRUE)
  )
  mask <- logical(n)
  mask[hits] <- TRUE

  if (!any(mask))
    stop("No SNPs overlap any region in the provided BED file.")

  message(sprintf(
    "BED filter: retaining %d / %d SNPs across %d region(s) on %d chromosome(s).",
    sum(mask), n,
    length(bed_gr),
    length(unique(as.character(GenomicRanges::seqnames(bed_gr))))
  ))
  mask
}
