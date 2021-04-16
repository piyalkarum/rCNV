# 1. wrapper for progressbar
apply_pb <- function(X, MARGIN, FUN, ...)
{
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total,
                       style = 3)

  wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env),
                      curVal +1)
    FUN(...)
  }
  res <- apply(X, MARGIN, wrapper, ...)
  close(pb)
  res
}
#https://ryouready.wordpress.com/2010/01/11/progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/


#' Import VCF file
#'
#' Function to import raw single/multi-sample VCF files generated from GatK or VCFtools.
#' The function required the R-package "data.table" for faster importing.
#'
#' @param vcf.file.path path to the vcf file
#'
#' @return Returns a dataframe with scaffold/chromosome positions and depth values for individual samples
#' @importFrom data.table fread
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path)
#'
#' @export
readVCF <- function(vcf.file.path){
  tt <- fread(vcf.file.path,sep="\t",skip="#CHROM")
  return(tt)
}


#' Generate heterozygote table (allele depth values)
#'
#' hetTgen extracts the read depth and coverage values for each snp for all the individuals from a vcf file generated from readVCF (or GatK VariantsToTable: see details)
#'
#' @param vcf an imported vcf file in data.frame or matrix format using "readVCF"
#'
#' @details If you generate the depth values for allele by sample using GatK VariantsToTable option, use only -F CHROM -F POS -GF AD flags to generate the table. Or keep only the CHROM POS and individual AD columns.
#'
#' @author Piyal Karunarathne
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf)
#'
#' @export
hetTgen<-function(vcf){
  xx <- vcf[,10:ncol(vcf)]
  AD<-which(strsplit(as.character(vcf[1,9]),":")[[1]]=="AD")
  h.table<-apply_pb(xx,2,function(X)do.call(rbind,lapply(X,function(x) paste(strsplit(x, ":")[[1]][AD], collapse = ':'))))
  het.table<-cbind(vcf[,1:3],h.table)
  colnames(het.table)[1]<-"CHROM"
  return(het.table)
}


