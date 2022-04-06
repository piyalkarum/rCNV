#' Generate allele frequency table (from genotypes)
#'
#' Get alternative allele frequency across all individuals per SNP from the
#'  genotype table
#'
#' @param gtt a genotype table produced from \code{hetTgen} (or similar)
#' @param verbose logical. whether to show the progress of the analysis
#'
#' @details Use hetTgen function to generate the genotype table with the
#' \code{GT} option
#' @return Returns a data frame of allele frequencies calculated from genotypes
#'
#' @author Piyal Karunarathne
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf,"GT")
#' frQ<-allele.freq(het.table)
#'
#' @export
allele.freq<-function(gtt,verbose=TRUE){
  gs<-gtt[,-c(1:4)]
  if(verbose){
    tmp<-apply_pb(gs,1,function(x){
      x<-as.character(x)
      tl<-strsplit(x,"/")
      tt<-unlist(lapply(tl,function(y){length(which(y==1))/2}))
      return(tt)
    })
  } else {
    tmp<-apply_pb(gs,1,function(x){
      x<-as.character(x)
      tl<-strsplit(x,"/")
      tt<-unlist(lapply(tl,function(y){length(which(y==1))/2}))
      return(tt)
    })
  }
  tmp<-t(tmp)
  tmp[gs=="./."]<-NaN
  colnames(tmp)<-colnames(gs)
  rownames(tmp)<-paste0(gtt[,1],".",gtt[,2])
  tmp<-data.frame(gtt[,1:4],tmp)
  return(tmp)
}
