# 1. wrapper for progressbar
# adopted from https://ryouready.wordpress.com/2010/01/11/progress-bars-in-r-part-ii-a-wrapper-for-apply-functions/

apply_pb <- function(X, MARGIN, FUN, ...)
{
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total,width = 50,
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
  pb <- txtProgressBar(min = 0, max = pb_Total, width = 50,style = 3)

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

sapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, width = 50,style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- sapply(X, wrapper, ...)
  close(pb)
  res
}


combn_pb <- function(X, size, FUN, ...)
{
  env <- environment()
  n<-length(X)
  r<-size
  pb_Total <- factorial(n)/(factorial(r)*factorial(n-r))
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total,width = 50,
                       style = 3)

  wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env),
                      curVal +1)
    FUN(...)
  }
  res <- combn(X, size, wrapper, ...)
  close(pb)
  res
}

# for-loop progress bar
# for(i in seq_along(xx)) {
#   pb <- txtProgressBar(min = 0, max = length(xx), style = 3, width = 50, char = "=")
#   setTxtProgressBar(pb, i)
#   ### code
# }
# close(pb)


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
#' A function to remove the alleles with minimum allele frequency and keep only
#' a bi-allelic matrix when loci are multi-allelic
#'
#' @param h.table allele depth table generated from the function \code{hetTgen}
#' @param AD logical. If TRUE a allele depth table similar to \code{hetTgen}
#' output will be returns; If \code{FALSE}, individual AD values per SNP will be
#' returned in a list.
#' @param verbose logical. Show progress
#'
#' @return A data frame or a list of minimum allele frequency removed allele depth
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{mf<-maf(ADtable)}
#'
#' @export
maf<-function(h.table,AD=TRUE,verbose=TRUE){
  htab<-h.table[,-c(1:3)]
  if(verbose){
    glt<-apply_pb(htab,1,function(X){
      gg<-do.call(rbind,lapply(X,function(x){as.numeric(strsplit(x,",")[[1]])}))
      while(ncol(gg)>2){gg<-gg[,-which.min(colMeans(proportions(gg,1),na.rm=T))]}
      if(AD){return(paste0(gg[,1],",",gg[,2]))}else{return(gg)}
    })
  }else{
    glt<-apply(htab,1,function(X){
      gg<-do.call(rbind,lapply(X,function(x){as.numeric(strsplit(x,",")[[1]])}))
      if(ncol(gg)>2)gg<-gg[,-which.min(colMeans(proportions(gg,1),na.rm=T))]
      if(AD){return(paste0(gg[,1],",",gg[,2]))}else{return(gg)}
    })
  }
  glt<-data.frame(h.table[,1:3],t(glt))
  colnames(glt)<-colnames(h.table)
  return(glt)
}



#' Import VCF file
#'
#' Function to import raw single and multi-sample VCF files.
#' The function required the R-package \code{data.table} for faster importing.
#'
#' @param vcf.file.path path to the vcf file
#' @param verbose logical. show progress
#'
#' @return Returns a list with vcf table in a data frame, excluding meta data.
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


#' Generate allele depth or genotype table
#'
#' hetTgen extracts the read depth and coverage values for each snp for all
#' the individuals from a vcf file generated from readVCF (or GatK
#' VariantsToTable: see details)
#'
#' @param vcf an imported vcf file in a list using \code{readVCF}
#' @param info.type character. \code{AD}: allele depth value, {AD-tot}:total
#' allele depth, \code{DP}=unfiltered depth (sum), \code{GT}: genotype,
#' \code{GT-012}:genotype in 012 format, \code{GT-AB}:genotype in AB format.
#' Default \code{AD},  See details.
#' @param verbose logical. whether to show the progress of the analysis
#'
#' @details If you generate the depth values for allele by sample using GatK
#' VariantsToTable option, use only -F CHROM -F POS -GF AD flags to generate
#' the table. Or keep only the CHROM, POS, ID, ALT, and individual AD columns.
#' For info.type \code{GT} option is provided to extract the genotypes of
#' individuals by snp.
#' @return Returns a data frame of allele depth, genotype of SNPs for all the
#' individuals extracted from a VCF file
#'
#' @author Piyal Karunarathne, Klaus Schliep
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf)
#'
#' @export
hetTgen <- function(vcf,info.type=c("AD","AD-tot","GT","GT-012","GT-AB","DP"),verbose=TRUE){
  if(inherits(vcf,"list")){vcf<-vcf$vcf}
  if(inherits(vcf,"data.frame")){vcf<-data.table::data.table(vcf)}
  if(any(nchar(vcf$ALT)>1)){
    warning("vcf file contains multi-allelic variants: only bi-allelic SNPs allowed\nUse maf() to remove non-bi-allilic snps")
  }
  if(inherits(vcf,"list")){vcf<-vcf$vcf}

  info.type<-match.arg(info.type)
  itype<-substr(info.type,1,2)

  adn <- sapply(strsplit(unlist(vcf[,"FORMAT"], use.names = FALSE),":"), function(x)match(itype, x))
  max_adn <- max(adn) + 1L
  ind <- cbind(seq_along(adn), adn)
  xx<-data.frame(vcf[,-c(1:9)])

  if(info.type=="AD-tot"){
    h.table <- matrix(NA_integer_, nrow(xx), ncol(xx))
    if(verbose) message("generating total depth values")
    for(i in seq_len(ncol(xx))){
      if(verbose){
        pb <- txtProgressBar(min = 0, max = ncol(xx), style = 3, width = 50, char = "=")
        setTxtProgressBar(pb, i)
      }
      tmp <- stringr::str_split_fixed(xx[,i], ":", max_adn)[ind]
      tmp <- stringr::str_split_fixed(tmp, ",", 2L)
      h.table[, i] <- as.numeric(tmp[,1]) + as.numeric(tmp[,2]) ## as.integer???
    }
  }
  else {
    h.table <- matrix(NA_character_, nrow(xx), ncol(xx))
    if(verbose) message("generating table")
    for(i in seq_len(ncol(xx))){
      if(verbose){
        pb <- txtProgressBar(min = 0, max = ncol(xx), style = 3, width = 50, char = "=")
        setTxtProgressBar(pb, i)
      }
      h.table[, i] <- stringr::str_split_fixed(xx[,i], ":", max_adn)[ind]
    }
    if(info.type!="DP"){h.table[is.na(h.table) | h.table==".,."]<-"./."}
  }
if(verbose) {close(pb)}
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
    h.table[h.table=="./."| h.table=="."]<- -9
  }
  if(info.type=="AD" ){
    h.table[h.table=="./." | h.table=="." | is.na(h.table)]<-"0,0"
  }
  if(info.type=="DP"){
    h.table[is.character(h.table)]<-0
    h.table[is.na(h.table)]<-0
  }
  het.table<-as.data.frame(cbind(vcf[,c(1:3,5)],h.table))
  colnames(het.table)<-c("CHROM",colnames(vcf)[c(2,3,5,10:ncol(vcf))])
  return(het.table)
}

#' Get missingness of individuals in raw vcf
#'
#' A function to get the percentage of missing data of snps per SNP and per
#' sample
#'
#' @param data a list containing imported vcf file using \code{readVCF} or
#' genotype table generated using \code{hetTgen}
#' @param type character.  Missing percentages per sample
#' \dQuote{samples} or per SNP \dQuote{snps}, default both
#' @param plot logical. Whether to plot the missingness density with ninety
#' five percent quantile
#' @param verbose logical. Whether to show progress
#'
#'
#' @author Piyal Karunarathne
#' @importFrom graphics abline polygon
#' @importFrom stats density
#'
#' @returns
#' Returns a data frame of allele depth or genotypes
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
    ll<-t(apply(ndat[,-c(1:4)],2,function(x){
      cbind(sum(x=="./." | is.na(x) | x=="0,0" | x==".,."),sum(x=="./." | is.na(x) | x=="0,0" | x==".,.")/length(x))
    }))
    ll<-data.frame(indiv=colnames(ndat)[-c(1:4)],n_miss=ll[,1],f_miss=ll[,2])
    rownames(ll)<-NULL
  }
  if(any(type=="snps")){
    L<-apply(ndat[,-c(1:4)],1,function(x){
      cbind(sum(x=="./." | is.na(x) | x=="0,0" | x==".,."),sum(x=="./." | is.na(x) | x=="0,0" | x==".,.")/length(x))
    })
    if(is.list(L)){
      L<-do.call(rbind,L)
    } else { L<-t(L)}
    colnames(L)<-c("n_miss","f_miss")
    L<-data.frame(ndat[,1:3],L)
  }
  if(plot){
    opars<-par(no.readonly = TRUE)
    on.exit(par(opars))

    if(length(type)==2){par(mfrow=c(1,2))}
    cl<-makeTransparent(c(1,2),alpha = 0.6)
    #missing samples
    if(any(type=="samples")){
      plot(density(ll$f_miss),type="n",main="Missing % per sample", xlim = c(0,1))
      polygon(density(ll$f_miss),border=2,col=cl[1])
      abline(v=quantile(ll$f_miss,p=0.95),lty=3,col="blue")
      text(x=quantile(ll$f_miss,p=0.95)+0.02,y=max(density(ll$f_miss)$y)/2,round(quantile(ll$f_miss,p=0.95),3),offset=10,col=2)
      legend("topright",lty=3,col="blue",legend="95% quantile",bty="n",cex=0.8)
    }
    #missing snps
    if(any(type=="snps")){
      plot(density(L$f_miss),type="n",main="Missing % per SNP", xlim = c(0,1))
      polygon(density(L$f_miss),border=2,col=cl[1])
      abline(v=quantile(L$f_miss,p=0.95),lty=3,col="blue")
      text(x=quantile(L$f_miss,p=0.95)+0.02,y=max(density(L$f_miss)$y)/2,round(quantile(L$f_miss,p=0.95),3),offset=10,col=2)
      legend("topright",lty=3,col="blue",legend="95% quantile",bty="n",cex=0.8)
    }

  }
  if(!exists("ll")){ll<-NULL}
  if(!exists("L")){L<-NULL}
  return(list(perSample=ll,perSNP=L))
}



#' Format genotype for BayEnv and BayPass
#'
#' This function generates necessary genotype count formats for BayEnv and
#' BayPass with a subset of SNPs
#'
#' @param gt multi-vector. an imported data.frame of genotypes or genotype
#' data frame generated by \code{hetTgen} or path to GT.FORMAT
#' file generated from VCFTools
#' @param info a data frame containing sample and population information.
#' It must have \dQuote{sample} and \dQuote{population} columns
#' @param format character. output format i.e., for BayPass or BayEnv
#' @param snp.subset numerical. number of randomly selected subsets of SNPs.
#' \code{default = NULL}
#' subset
#'
#' @return Returns a list with formatted genotype data: \code{$bayenv} - snps
#' in horizontal format - for BayEnv (two lines per snp); \code{$baypass} - vertical format - for BayPass
#' (two column per snp); \code{$sub.bp} - subsets snps for BayPass \code{$sub.be} - subsets of snps for BayEnv
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' het.table<-hetTgen(vcf,"GT")
#' info<-unique(substr(colnames(het.table)[-c(1:3)],1,8))
#' GT<-gt.format(het.table,info)}
#'
#' @export
gt.format <- function(gt,info,format=c("benv","bpass"),snp.subset=NULL) {
  if(is.character(gt)){
    gt <-as.data.frame(fread(gt))
    gts <-gt[,-c(1,2)]
  } else {
    gts<-gt[,-c(1:4)]
  }
  if(is.character(info)){
    if(length(info)==ncol(gts)){
      info<-data.frame(population=info)
    } else {
      pop.col<-NULL
      for(i in seq_along(info)){
        pop.col[grep(info[i],colnames(gts))]<-info[i]
      }
      info<-data.frame(population=pop.col)
    }
  }
  rownames(gts)<-paste(gt$CHROM,gt$POS,sep=".")
  pp<-na.omit(unique(info$population))
  infos<-as.character(info$population)
  format<-match.arg(format,several.ok = TRUE)

  lgt<-split.data.frame(t(gts),f=info$population)

  if(any(format=="benv")){
    message("Formating BayEnv")
    ppe<-lapply_pb(lgt,function(X){
      out<-apply(X,2,function(x){zero<-sum((unlist(stringr::str_split(x,"/||")))==0)
      one<-sum((unlist(stringr::str_split(x,"/||")))==1)
      ot<-as.data.frame(c(zero,one),col.names=F)
      return(ot)},simplify = F)
      out<-do.call(rbind,out)
      colnames(out)<-NULL
      return(out)
    })
    ppe<-do.call(cbind,ppe)
    rownames(ppe)<-paste0(paste0(gt[,1],".",gt[,2]),"~",rep(c(1,2),nrow(gt)))
  } else {ppe<-NULL}

  if(any(format=="bpass")){
    message("Formating BayPass")
    ppp<-lapply_pb(lgt,function(X){
      out<-apply(X,2,function(x){zero<-sum((unlist(stringr::str_split(x,"/||")))==0)
      one<-sum((unlist(stringr::str_split(x,"/||")))==1)
      ot<-c(zero,one)
      return(ot)},simplify = F)
      out<-do.call(rbind,out)
      colnames(out)<-NULL
      return(out)
    })
    ppp<-do.call(cbind,ppp)
    colnames(ppp)<-paste0(rep(names(lgt),each=2),"~",rep(c(1,2),ncol(ppp)/2))
    rownames(ppp)<-paste0(gt[,1],".",gt[,2])

    if(!is.null(snp.subset)){
      rn<-sample(1:snp.subset,nrow(gts),replace = T)
      rownames(gts)<-paste0(gt[,1],".",gt[,2])
      chu.p<-split.data.frame(ppp,f=rn)
    } else { chu.p <- NULL}
  } else {ppp<-NULL}
  return(list(baypass=ppp,bayenv=ppe,sub.bp=chu.p,pop=as.character(pp)))
}



#' Correct allele depth values
#'
#' A function to correct depth values with odd number of coverage values due to
#' sequencing anomalies or miss classification where genotype is homozygous and
#' depth values indicate heterozygosity.
#' The function adds a value of one to the allele with the lowest depth value
#' for when odd number anomalies or make the depth value zero for when
#' miss-classified. The genotype table must be provided for the latter.
#'
#' @param het.table allele depth table generated from the function
#' \code{hetTgen}
#' @param gt.table genotype table generated from the function hetTgen
#' @param odd.correct logical, to correct for odd number anomalies in AD values.
#'  default \code{TRUE}
#' @param verbose logical. show progress. Default \code{TRUE}
#'
#' @return Returns the coverage corrected allele depth table similar to the
#'  output of \code{hetTgen}
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{adc<-ad.correct(ADtable)}
#'
#' @export
ad.correct<-function(het.table,gt.table=NULL,odd.correct=TRUE,verbose=TRUE){
  if(!is.null(gt.table)){
    if(verbose){
      message("correcting genotype mis-classification")
      Nw.ad<-lapply_pb(5:ncol(het.table),function(n){
        X<-het.table[,n]
        x<-gt.table[,n]
        Y<-data.frame(do.call(cbind,data.table::tstrsplit(X,",")))
        y<-which(x=="0/0" & Y$X2>0)
        Y[y,]<-0
        z<-which(x=="./.")
        Y[z,]<-0
        out<-paste0(Y$X1,",",Y$X2)
        return(out)
      })
      Nw.ad<-do.call(cbind,Nw.ad)
      Nw.ad<-data.frame(het.table[,1:4],Nw.ad)
      colnames(Nw.ad)<-colnames(het.table)
    } else {
      Nw.ad<-lapply(5:ncol(het.table),function(n){
        X<-het.table[,n]
        x<-gt.table[,n]
        Y<-data.frame(do.call(cbind,data.table::tstrsplit(X,",")))
        y<-which(x=="0/0" & Y$X2>0)
        rr<-range(Y$X2[y])#range of depth in miss classified snps
        Y[y,]<-0
        ll<-length(y)#number of miss classifications
        out<-paste0(Y$X1,",",Y$X2)
        return(out)
      })
      Nw.ad<-do.call(cbind,Nw.ad)
      Nw.ad<-data.frame(het.table[,1:4],Nw.ad)
      colnames(Nw.ad)<-colnames(het.table)
    }
    het.table<-Nw.ad
    rm(Nw.ad)
  }
  X<-data.frame(het.table[,-c(1:4)])
  if(odd.correct){
    if(verbose){
      message("correcting odd number anomalies")
      vv<-apply_pb(X,2,function(sam){
        sam<-unname(unlist(sam))
        tmp <- stringr::str_split_fixed(sam, ",", 2L)
        tmp<-cbind(as.integer(tmp[,1]),as.integer(tmp[,2]))
        tch<-which((rowSums(tmp,na.rm=TRUE)%%2)!=0)
        for(i in tch){
          tm<-tmp[i,]
          if(any(tm==0)){tm[which.max(tm)]<-tm[which.max(tm)]+1L}else{tm[which.min(tm)]<-tm[which.min(tm)]+1L}
          tmp[i,]<-tm
        }
        return(paste0(tmp[,1],",",tmp[,2]))
      })
      if(inherits(vv,"list")){
        vv<-do.call(cbind,vv)
      }
    } else {
      vv<-apply(X,2,function(sam){
        dl<-lapply(sam,function(y){
          l<-as.numeric(unlist(strsplit(y,",")))
          if((sum(l,na.rm=TRUE)%%2)!=0){
            if(any(l==0)){l}else{
              l[which.min(l)]<-l[which.min(l)]+1}
          }
          return(paste0(l,collapse = ","))
        })
      })
      if(inherits(vv,"list")){
        vv<-do.call(cbind,vv)
      }
    }
    return(cbind(het.table[,1:4],vv))
  } else { return(het.table)}
}


