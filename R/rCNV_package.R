#' rCNV: Detect Copy Number Variants from SNP Data
#'
#' Functions in this package import filtered VCF files of SNP data and
#' generate datasets to detect copy number variants, visualize them, and
#' perform downstream analyses (e.g. environmental association analyses).
#'
#' @section Main WGS functions:
#' \describe{
#'   \item{\code{\link{allele.info.WGS3}}}{Duplication/deletion detection
#'     using a negative-binomial depth model; no GC correction.}
#'   \item{\code{\link{allele.info.WGS.GCcor}}}{Same framework with
#'     per-SNP GC-bias correction via SNP-centred windows.}
#' }
#'
#' @author Piyal Karunarathne \email{piyalkarumail@@yahoo.com}
#'
#' @importFrom GenomicRanges GRanges findOverlaps seqnames start end width
#'   resize trim
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevels seqlengths keepSeqlevels
#' @importFrom Rsamtools FaFile indexFa
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom BiocGenerics strand width
#' @importFrom S4Vectors queryHits
#' @importFrom rtracklayer import
#' @importFrom extraDistr dbbinom
#' @importFrom fitdistrplus fitdist
#' @importFrom pbapply pbapply pblapply
#' @importFrom stringi stri_split_fixed
#' @importFrom stats dnbinom dbinom pnorm quantile optim sd var
#' @importFrom utils head
#' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ
"_PACKAGE"
