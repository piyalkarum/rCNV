## VCF statistics ###
#1. get heterozygosity per individual
## P(Homo) = F + (1-F)P(Homo by chance)
## P(Homo by chance) = p^2+q^2 for a biallelic locus.
## For an individual with N genotyped loci, we
##   1. count the total observed number of loci which are homozygous (O),
##   2. calculate the total expected number of loci homozygous by chance (E)
## Then, using the method of moments, we have
##    O = NF + (1-F)E
## Which rearranges to give
##    F = (O-E)/(N-E)

het.sity<-function(ind){
  O<-length(which(ind=="0/0"))+length(which(ind=="1/1"))
  N<-length(which(ind=="0/0"))+length(which(ind=="1/1"))+length(which(ind=="0/1"))+length(which(ind=="1/0"))
  alle<-unname(unlist(lapply(ind,strsplit,split="/")))
  p<-length(which(alle=="0"))/(length(which(alle=="0"))+length(which(alle=="1")))
  q<-length(which(alle=="1"))/(length(which(alle=="0"))+length(which(alle=="1")))
  E<-(p^2+q^2)*N
  FF<-(O-E)/(N-E)
  return(c(O,E,N,FF))
}



#' Determine per sample heterozygosity and inbreeding coefficient
#'
#' This function will calculate the heterozygosity on a per-sample basis from vcf files (snps), and most importantly inbreeding coefficient which is used to filter out the samples with bad mapping quality. See details.
#' @param vcf an imported vcf file in data.frame or matrix format using "readVCF"
#' @param plot logical. Whether to plot a boxplot of inbreeding coefficients for populations. A list of populations must be provided
#' @param pops character. A list of population names
#' @importFrom graphics boxplot
#' @return Returns a data frame of expected "E(Hom)", observed "O(Hom)" homozygotes with their inbreeding coefficients.
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' hzygots<-h.zygosity(vcf,plot=TRUE,pops=c("FR_PA","RU_PA"))
#'
#' @export
h.zygosity<-function(vcf,plot=FALSE,pops=NA){
  gtt<-hetTgen(vcf,"GT",verbose=FALSE)
  hh<-t(apply_pb(gtt[,-c(1:3)],2,het.sity))
  hh<-data.frame(rownames(hh),hh)
  colnames(hh)<-c("ind","O(Hom)","E(Hom)","total","Fis")
  rownames(hh)<-NULL
  if(plot){
    if(is.na(pops[1])){
      warning("Please provide a population list")
    } else {
      hh$pop<-"ungrouped"
      for(i in seq_along(pops)){
        hh$pop[grep(pops[i],hh$ind)]<-pops[i]
      }
      if(length(unique(hh$pop))<=1){
        warning("population names must be a part of individual names")
      }
      boxplot(hh$Fis~hh$pop,ann=F)
    }
  }
  return(hh)
}







