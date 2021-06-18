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
non_bi_rm<-function(vcf,GT.table=NULL){
  if(inherits(vcf,"list")){vcf<-vcf$vcf}
  nbal<-which(apply(vcf[,5],1,nchar)>1)
  vcf<-vcf[!nbal,]
  if(!is.null(GT.table)){
    gtyp<-GT.table
  } else {
    gtyp<-hetTgen(vcf,"GT")
  }
  alcount<-apply(gtyp[,-c(1:3)],1,function(x){y=unique(x);y=y[y!="./."];return(length(y))})
  vcf<-vcf[!which(alcount<2),]
  return(list(vcf=vcf))
}

#3. helper for genotype count in gt.format
gg<-function(x,benv=TRUE){
  tt <-unlist(strsplit(x,split = "/"))
  tt<-data.frame(table(tt))
  rownames(tt)<-tt$tt
  if(benv){
    tl <-cbind(tt["0",2],tt["1",2])
  } else {
    tl <-data.frame(z=tt["0",2],on=tt["1",2])
  }
  return(tl)
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
  return(list(vcf=tt))
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
  if(inherits(vcf,"list")){vcf<-vcf$vcf}
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
  het.table[het.table=="NA" | het.table==".,."]<-"0,0"
  het.table <- data.frame(het.table)[,colSums(het.table=="0,0")<nrow(het.table)]
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
  vcf<-vcf$vcf
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


#' Format genotype for BayEnv and BayPass
#'
#' This function generates necessary genotype count formats for BayEnv and BayPass with a subset of SNPs
#'
#' @param gt multi-vector. path to GT.FORMAT file generated from VCFTools, an imported data.frame of genotypes or genotype data frame generated by "hetTgen"
#' @param info data frame containing sample and population information. It must have "sample" and "population" columns
#' @param snp.subset logical. whether to generate a randomly sampled subset of 1/10 snps
#'
#' @return returns a list with formated genotype data: $hor - snps in horizontal format (two lines per snp); $ver - vertical format (two column per snp); $hor.chunk - a subset snps of $hor
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf,"GT")
#' info<-unique(substr(colnames(het.table)[-c(1:3)],1,8))
#' GT<-gt.format(het.table,info)
#'
#' @export
gt.format <- function(gt,info,snp.subset=TRUE) {
  if(is.character(gt)){
    gts <-as.data.frame(fread(gt))
  } else {
    gts<-gt
  }
  if(is.character(info)){
    pop.col<-NULL
    for(i in seq_along(info)){
      pop.col[grep(info[i],colnames(gts))]<-info[i]
    }
    info<-data.frame(population=pop.col)
  }
  rownames(gts)<-paste(gts$CHROM,gts$POS,sep=".")
  pp<-na.omit(unique(info$population))

  message("generating horizontal table")
  pgt.h<-lapply_pb(pp,function(x,gts,info){
    tm <- t(gts[,which(info$population==x)])
    gtt <- lapply(data.frame(tm),gg,benv=T)
    gf <- do.call(cbind,gtt)
    return(gf)
  },gts=gts,info=info)
  pgt.h<-do.call(rbind,pgt.h)

  message("generating vertical table")
  pgt.v<-lapply_pb(pp,function(x,gts,info){
    tm <- t(gts[,which(info$population==x)])
    gtt <- apply(tm,2,gg,benv=F)
    gf <- do.call(rbind,gtt)
    return(gf)
  },gts=gts,info=info)
  pgt.v<-do.call(cbind,pgt.v)

  nm <- unlist(lapply(rownames(gts),FUN=function(x)c(paste0(x,"~1"),paste0(x,"~2"))))
  colnames(pgt.h)<-nm
  rownames(pgt.h)<-pp
  nm <- unlist(lapply(pp,FUN=function(x)c(paste0(x,".1"),paste0(x,".2"))))
  colnames(pgt.v)<-nm

  GP.h <- as.matrix(pgt.h)
  GP.h[which(is.na(GP.h))] <- 0
  GP.v <- as.matrix(pgt.v)
  GP.v[which(is.na(GP.v))] <- 0
  rownames(GP.v)<-rownames(gts)

  message("generating snps subset")
  if(snp.subset){
    rn<-sample(1:10,nrow(gts),replace = T)
    snps<-unique(rownames(gts))
    subsnp<-snps[rn==2]## take 1/10 of snps randomly
    gts2<-gts[subsnp,]
    info2<-info
    pp2<-na.omit(unique(info2$population))

    pgt.chunk<-lapply_pb(pp2,function(x,gts,info){
      tm <- t(gts[,which(info$population==x)])
      gtt <- lapply(data.frame(tm),gg,benv=T)
      gf <- do.call(cbind,gtt)
      return(gf)
    },gts=gts2,info=info2)

    pgt.chunk<-do.call(rbind,pgt.chunk)
    nm <- unlist(lapply(rownames(gts2),FUN=function(x)c(paste0(x,"~1"),paste0(x,"~2"))))
    colnames(pgt.chunk)<-nm
    rownames(pgt.chunk)<-pp2
    GP.chunk <- as.matrix(pgt.chunk)
    GP.chunk[which(is.na(GP.chunk))] <- 0

  } else { GP.chunk <- NULL}

  return(list(hor=GP.h,ver=GP.v,hor.subset=GP.chunk))
}



