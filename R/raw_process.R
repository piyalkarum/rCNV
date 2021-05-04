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

lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

#https://ryouready.wordpress.com/2010/01/11/progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/

#2. remove non-biallelic snps
non_bi_rm<-function(vcf){
  nbal<-which(apply(vcf[,5],1,nchar)>1)
  vcf<-vcf[!nbal,]
  gtyp<-hetTgen(vcf,"GT")
  alcount<-apply(gtyp[,-c(1:3)],1,function(x){y=unique(x);y=y[y!="./."];return(length(y))})
  vcf<-vcf[!which(alcount<2),]
}




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
#' @param info.type character. "AD"=allele depth value, "GT"=genotype. Default "AD". See details.
#' @param verbose logical. whether to show the progress of the analysis
#'
#' @details If you generate the depth values for allele by sample using GatK VariantsToTable option, use only -F CHROM -F POS -GF AD flags to generate the table. Or keep only the CHROM POS and individual AD columns.
#' For info.type "GT" option is provided to extract the genotypes of individuals by snp. However to obtain the heterozygosity data necessary for getting the duplication information with the function dup.snp.info, info.type needs to be set to "AD", the default argument.
#'
#' @author Piyal Karunarathne
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf)
#'
#' @export
hetTgen<-function(vcf,info.type=c("AD","GT"),verbose=TRUE){
  if(length(which(apply(vcf[,5],1,nchar)>1))>1){
    message("vcf file contains multi-allelic variants: only bi-allelic SNPs allowed")
    ans<-readline(prompt="Remove non-biallelic SNPs (y/n) ?: ")
    if(ans=="y"){
      vcf<-non_bi_rm(vcf)
    }
    if(ans=="n"){
      stop("Non-bi-allelic variant file: 0nly bi-allelic SNPs allowed so far")
    }
  }
  xx <- vcf[,10:ncol(vcf)]
  info.type<-match.arg(info.type)
  AD<-which(strsplit(as.character(vcf[1,9]),":")[[1]]==info.type)
  if(length(AD)==0){
    AD<-which(strsplit(as.character(vcf[1,9]),":")[[1]]=="DPR")
  }
  if(verbose) {
    message("generating heterozygotes table")
    h.table<-apply_pb(xx,2,function(X)do.call(rbind,lapply(X,function(x) paste(strsplit(x, ":")[[1]][AD], collapse = ':'))))
  } else {
    h.table<-apply(xx,2,function(X)do.call(rbind,lapply(X,function(x) paste(strsplit(x, ":")[[1]][AD], collapse = ':'))))
  }
  het.table<-cbind(vcf[,1:3],h.table)
  colnames(het.table)[1]<-"CHROM"
  het.table[het.table=="NA"]<-"0,0"
  return(het.table)
}


#' Get missingness of individuals in raw vcf
#'
#' A function to get missingness of snps data on a per-individual basis similar to --missing-indiv option in vcftools
#'
#' @param vcf data frame of imported vcf file using readVCF
#' @param plot logical. Whether to plot the missingness density with a suggested cut-off.
#'
#' @author Piyal Karunarathne
#' @importFrom graphics abline polygon
#' @importFrom stats density
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' missing<-get.miss(vcf,plot=TRUE)
#'
#' @export
get.miss<-function(vcf,plot=TRUE){
  ndat<-hetTgen(vcf,"GT")
  ll<-t(apply(ndat[,-c(1:3)],2,function(x)cbind(length(which(x=="./.")),length(which(x=="./."))/length(x))))
  ll<-data.frame(indiv=colnames(vcf)[-c(1:9)],n_miss=ll[,1],f_miss=ll[,2])
  rownames(ll)<-NULL
  if(plot){
    plot(density(ll$f_miss),type="n",main="Missingness %")
    polygon(density(ll$f_miss),border="red",col="lightblue")
    abline(v=quantile(ll$f_miss,p=0.95),lty=3,col="blue")
    legend("topright",lty=3,col="blue",legend="suggested cut-off",bty="n",cex=0.8)
  }
  return(ll)
}

