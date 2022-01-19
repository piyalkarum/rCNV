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

#(extra) generate colors for Rmarkdown docs [extracted from Rmarkdown guide book]
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}

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
gg<-function(x){
  tt <-unlist(strsplit(x,split = "/"))
  tt<-data.frame(table(tt))
  rownames(tt)<-tt$tt
  tl <-cbind(tt["0",2],tt["1",2])
  return(tl)
}

#' Remove MAF allele
#'
#' A function to remove the alleles with minimum allele frequency and keep only a biallelic matrix
#'
#' @param h.table a data frame of allele depth values generated from hetTgen
#' @param AD logical. If TRUE a allele depth table similar to hetTgen output will be returns; If FALSE, individual AD values per SNP will be returned in a list.
#'
#' @return a data frame or a list
#'
#' @author Piyal Karunarathne
#'
#' @export
maf<-function(h.table,AD=TRUE,verbose=TRUE){
  htab<-h.table[,-c(1:3)]
  if(verbose){
    glt<-apply_pb(htab,1,function(X){
      gg<-do.call(rbind,lapply(X,function(x){as.numeric(strsplit(x,",")[[1]])}))
      pp<-proportions(gg,1)
      if(ncol(gg)>2)gg<-gg[,-which.min(colMeans(pp,na.rm=T))]
      if(AD){return(paste0(gg[,1],",",gg[,2]))}else{return(gg)}
    })
  }else{
    glt<-apply(htab,1,function(X){
      gg<-do.call(rbind,lapply(X,function(x){as.numeric(strsplit(x,",")[[1]])}))
      pp<-proportions(gg,1)
      if(ncol(gg)>2)gg<-gg[,-which.min(colMeans(pp,na.rm=T))]
      if(AD){return(paste0(gg[,1],",",gg[,2]))}else{return(gg)}
    })
  }
  glt<-data.frame(h.table[,1:3],t(glt))
  colnames(glt)<-colnames(h.table)
  return(glt)
}



#' Import VCF file
#'
#' Function to import raw single/multi-sample VCF files generated from GatK or VCFtools.
#' The function required the R-package "data.table" for faster importing.
#'
#' @param vcf.file.path path to the vcf file
#' @param verbose logical. show progress
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
readVCF <- function(vcf.file.path,verbose=FALSE){
  tt <- fread(vcf.file.path,sep="\t",skip="#CHROM",verbose=verbose)
  return(list(vcf=tt))
}


#' Generate heterozygote table (allele depth values)
#'
#' hetTgen extracts the read depth and coverage values for each snp for all the individuals from a vcf file generated from readVCF (or GatK VariantsToTable: see details)
#'
#' @param vcf an imported vcf file in data.frame or matrix format using "readVCF"
#' @param info.type character. "AD"=allele depth value, "AD-tot"=allele total depth, "DP"=unfiltered depth (sum), "GT"=genotype, "GT-012"=genotype in 012 format, "GT-AB"=genotype in AB format. Default "AD",  See details.
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
hetTgen<-function(vcf,info.type=c("AD","AD-tot","GT","GT-012","GT-AB","DP"),verbose=TRUE){
  if(inherits(vcf,"list")){vcf<-vcf$vcf}
  if(inherits(vcf,"data.frame")){vcf<-data.table::data.table(vcf)}
  if(length(which(apply(vcf[,5],1,nchar)>1))>1){
    message("vcf file contains multi-allelic variants: only bi-allelic SNPs allowed")
    ans<-readline(prompt="Remove non-biallelic SNPs (y/n) ?: ")
    if(ans=="y"){
      vcf<-non_bi_rm(vcf)
    }
  }
  if(inherits(vcf,"list")){vcf<-vcf$vcf}

  info.type<-match.arg(info.type)
  itype<-substr(info.type,1,2)
  adn<-unname(unlist(lapply(unlist(vcf[,"FORMAT"]),function(i){which(strsplit(i,":")[[1]]==itype)})))
  xx<-data.frame(as.numeric(adn),vcf[,-c(1:9)])

  if(verbose) {
    if(itype=="AD"){message("generating allele depth table")
      h.table<-t(apply_pb(xx,1,function(X){do.call(cbind,lapply(X,function(x){paste(strsplit(x, ":")[[1]][as.numeric(X[1])], collapse = ':')}))}))
    } else if(itype=="GT"){message("generating genotypes table")
      h.table<-t(apply_pb(xx,1,function(X){do.call(cbind,lapply(X,function(x){paste(strsplit(x, ":")[[1]][as.numeric(X[1])], collapse = ':')}))}))
    }
    else if(itype=="DP"){message("generating unfiltered allele depth table")
      h.table<-t(apply_pb(xx,1,function(X){do.call(cbind,lapply(X,function(x){paste(strsplit(x, ":")[[1]][as.numeric(X[1])], collapse = ':')}))}))}
  } else {
    if(itype=="AD"){h.table<-t(apply(xx,1,function(X){do.call(cbind,lapply(X,function(x){paste(strsplit(x, ":")[[1]][as.numeric(X[1])], collapse = ':')}))}))}
    else if(itype=="DP"){h.table<-t(apply(xx,1,function(X){do.call(cbind,lapply(X,function(x){paste(strsplit(x, ":")[[1]][as.numeric(X[1])], collapse = ':')}))}))}
    else if(itype=="GT"){h.table<-t(apply(xx,1,function(X){do.call(cbind,lapply(X,function(x){paste(strsplit(x, ":")[[1]][as.numeric(X[1])], collapse = ':')}))}))}
  }
  h.table<-h.table[,-1]

  if(info.type!="DP"){h.table[is.na(h.table) | h.table==".,."]<-"./."}

  if(info.type=="AD-tot"){
    if(verbose){
      message("generating total depth values")
      h.table<-apply_pb(h.table,2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))))}))})
    } else {
      h.table<-apply(h.table,2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))))}))})
    }
  }
  if(info.type=="GT-012"){
    h.table[h.table=="0/0"]<-0
    h.table[h.table=="1/1"]<-1
    h.table[h.table=="1/0" | h.table=="0/1"] <- 2
    h.table[h.table=="./."| h.table=="."]<-NA
  }
  if(info.type=="GT-AB"){
    h.table[h.table=="0/0"]<-"AA"
    h.table[h.table=="1/1"]<-"BB"
    h.table[h.table=="1/0" | h.table=="0/1"] <- "AB"
    h.table[h.table=="./."| h.table=="."]<--9
  }
  if(info.type=="AD" ){
    h.table[h.table=="./." | h.table=="."]<-"0,0"
  }
  if(info.type=="DP"){
    h.table[is.character(h.table)]<-0
    h.table[is.na(h.table)]<-0
  }
  het.table<-as.data.frame(cbind(vcf[,1:3],h.table))
  colnames(het.table)[1]<-"CHROM"
  return(het.table)
}


#' Get missingness of individuals in raw vcf
#'
#' A function to get missingness of snps data on a per-individual basis similar to --missing-indiv option in vcftools
#'
#' @param data a list containing imported vcf file using readVCF or genotype table generated using hetTgen
#' @param type character. Get the missing percentages per sample or per SNP."samples" or "snps", default both
#' @param plot logical. Whether to plot the missingness density with a suggested cut-off
#' @param verbose logical. Whether to show progress
#'
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
get.miss<-function(data,type=c("samples","snps"),plot=TRUE,verbose=TRUE){
  if(inherits(data,"list")){
    vcf<-data$vcf
    ndat<-hetTgen(vcf,"GT",verbose=verbose)
  } else {ndat<-data}
  type<-match.arg(type,several.ok = T)
  if(any(type=="samples")){
    ll<-t(apply(ndat[,-c(1:3)],2,function(x){
      cbind(length(which(x=="./." | is.na(x) | x=="0,0" | x==".,.")),length(which(x=="./." | is.na(x) | x=="0,0" | x==".,."))/length(x))
    }))
    ll<-data.frame(indiv=colnames(ndat)[-c(1:3)],n_miss=ll[,1],f_miss=ll[,2])
    rownames(ll)<-NULL
  }
  if(any(type=="snps")){
    L<-apply(ndat[,-c(1:3)],1,function(x){
      cbind(length(which(x=="./." | is.na(x) | x=="0,0" | x==".,.")),length(which(x=="./." | is.na(x) | x=="0,0" | x==".,."))/length(x))
    })
    if(is.list(L)){
      L<-do.call(rbind,L)
    } else { L<-t(L)}
    colnames(L)<-c("n_miss","f_miss")
    L<-data.frame(ndat[,1:3],L)
  }
  if(plot){
    if(length(type)==2){par(mfrow=c(1,2))}
    #missing samples
    if(any(type=="samples")){
      plot(density(ll$f_miss),type="n",main="Missing % per sample", xlim = c(0,1))
      polygon(density(ll$f_miss),border="red",col="lightblue")
      abline(v=quantile(ll$f_miss,p=0.95),lty=3,col="blue")
      text(x=quantile(ll$f_miss,p=0.95)+0.02,y=max(density(ll$f_miss)$y)/2,round(quantile(ll$f_miss,p=0.95),3),offset=10,col=2)
      legend("topright",lty=3,col="blue",legend="95% quantile",bty="n",cex=0.8)
    }
    #missing snps
    if(any(type=="snps")){
      plot(density(L$f_miss),type="n",main="Missing % per SNP", xlim = c(0,1))
      polygon(density(L$f_miss),border="red",col="lightblue")
      abline(v=quantile(L$f_miss,p=0.95),lty=3,col="blue")
      text(x=quantile(L$f_miss,p=0.95)+0.02,y=max(density(L$f_miss)$y)/2,round(quantile(L$f_miss,p=0.95),3),offset=10,col=2)
      legend("topright",lty=3,col="blue",legend="95% quantile",bty="n",cex=0.8)
    }
    par(mfrow=c(1,1))
  }
  if(!exists("ll")){ll<-NULL}
  if(!exists("L")){L<-NULL}
  return(list(perSample=ll,perSNP=L))
}



#' Format genotype for BayEnv and BayPass
#'
#' This function generates necessary genotype count formats for BayEnv and BayPass with a subset of SNPs
#'
#' @param gt multi-vector. path to GT.FORMAT file generated from VCFTools, an imported data.frame of genotypes or genotype data frame generated by "hetTgen"
#' @param info data frame containing sample and population information. It must have "sample" and "population" columns
#' @param snp.subset logical. whether to generate a randomly sampled subset of 1/10 snps
#' @param verbose logical. If TRUE shows progress
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
gt.format <- function(gt,info,snp.subset=FALSE,verbose=FALSE) {
  if(is.character(gt)){
    gt <-as.data.frame(fread(gt))
    gts <-gt[,-c(1,2)]
  } else {
    gts<-gt[,-c(1,2)]
  }
  if(is.character(info)){
    pop.col<-NULL
    for(i in seq_along(info)){
      pop.col[grep(info[i],colnames(gts))]<-info[i]
    }
    info<-data.frame(population=pop.col)
  }
  rownames(gts)<-paste(gt$CHROM,gt$POS,sep=".")
  pp<-na.omit(unique(info$population))

  infos<-as.character(info$population)
 if(verbose){
   pgt.h<-lapply_pb(pp,function(x,gts,info){
     tm <- as.data.frame(gts[,which(info==x)])
     gtt <- lapply(1:nrow(tm),function(y,tm){gg(as.character(tm[y,]))},tm=tm)
     gf <- t(do.call(cbind,gtt))
     colnames(gf)<-x
     return(gf)
   },gts=gts,info=infos)
 } else {
   pgt.h<-lapply(pp,function(x,gts,info){
     tm <- as.data.frame(gts[,which(info==x)])
     gtt <- lapply(1:nrow(tm),function(y,tm){gg(as.character(tm[y,]))},tm=tm)
     gf <- t(do.call(cbind,gtt))
     colnames(gf)<-x
     return(gf)
   },gts=gts,info=infos)
 }
  pgt.h<-do.call(cbind,pgt.h)
  nm <- unlist(lapply(rownames(gts),FUN=function(x)c(paste0(x,"~1"),paste0(x,"~2"))))
  rownames(pgt.h)<-nm
  pgt.h <- as.matrix(pgt.h)
  pgt.h[which(is.na(pgt.h))] <- 0
  if(snp.subset){
    message("subsetting")
    rn<-sample(1:10,nrow(gts),replace = T)
    snps<-rownames(gts)
    chu<-lapply_pb(1:10,function(nn,pgt.h,snps,rn){
      sset<-snps[rn==nn]
      tmp0<-NULL
      for(k in seq_along(sset)){
        tmp<-pgt.h[grep(sset[k],rownames(pgt.h)),]
        tmp0<-rbind(tmp0,tmp)
      }
      return(tmp0)
    },pgt.h=pgt.h,snps=snps,rn=rn)
  } else { chu <- NULL}
  return(list(hor=pgt.h,ver=t(pgt.h),subsets=chu,pop=as.character(pp)))
}





