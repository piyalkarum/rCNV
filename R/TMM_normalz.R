### TMM calculation of this function was adopted from edgeR package (Robinson et al. 2010)
### The methodology behind the TMM normalization was originally described by Mark Robinson (2010)

calcFactorQuantile <- function (data, lib.size, p=0.75)
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson and Gordon Smyth
  #	adopted from EdgeR package
{
  f <- rep_len(1,ncol(data))
  for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p,na.rm=TRUE)
  if(min(f)==0) warning("One or more quantiles are zero")
  f / lib.size
}

TMM <- function(obs, ref, logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10)
  #	TMM between two libraries
  #	Mark Robinson
  # adopted from edgeR package (Robinson et al. 2010)
{
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)

  nO <- sum(obs)
  nR <- sum(ref)

  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref

  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)

  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]

  if(max(abs(logR)) < 1e-6) return(1)

  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

  if(Weighting){
    f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
  } else{
    f <- mean(logR[keep], na.rm=TRUE)
  }

  if(is.na(f)) f <- 0

  2^f
}


TMMex <- function(obs,ref, logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10)
  #	Gordon Smyth
  #	adopted from edgeR package (Robinson et al. 2010)
{
  ref <- as.numeric(ref)
  obs <- as.numeric(obs)
  eps <- 1e-14 #floating zero

  pos.obs <- (obs > eps)
  pos.ref <- (ref > eps)
  npos <- 2L * pos.obs + pos.ref

  i <- which(npos==0L | is.na(npos))
  if(length(i)) {
    obs <- obs[-i]
    ref <- ref[-i]
    npos <- npos[-i]
  }
  #	Check library sizes
  libsize.obs <- sum(obs)
  libsize.ref <- sum(ref)

  zero.obs <- (npos == 1L)
  zero.ref <- (npos == 2L)
  k <- (zero.obs | zero.ref)
  n.eligible.singles <- min( sum(zero.obs), sum(zero.ref))
  if(n.eligible.singles > 0L) {
    refk <- sort(ref[k],decreasing=TRUE)[1:n.eligible.singles]
    obsk <- sort(obs[k],decreasing=TRUE)[1:n.eligible.singles]
    obs <- c(obs[!k],obsk)
    ref <- c(ref[!k],refk)
  } else {
    obs <- obs[!k]
    ref <- ref[!k]
  }

  n <- length(obs)
  if(n==0L) return(1)

  obs.p <- obs / libsize.obs
  ref.p <- ref / libsize.ref
  M <- log2( obs.p / ref.p )
  A <- 0.5 * log2( obs.p * ref.p )

  if(max(abs(M)) < 1e-6) return(1)

  obs.p.shrunk <- (obs+0.5) / (libsize.obs+0.5)
  ref.p.shrunk <- (ref+0.5) / (libsize.ref+0.5)
  M.shrunk <- log2( obs.p.shrunk / ref.p.shrunk )
  o.M <- order(M, M.shrunk)
  o.A <- order(A)

  loM <- as.integer(n * logratioTrim) + 1L
  hiM <- n + 1L - loM
  keep.M <- rep_len(FALSE,n)
  keep.M[o.M[loM:hiM]] <- TRUE
  loA <- as.integer(n * sumTrim) + 1L
  hiA <- n + 1L - loA
  keep.A <- rep_len(FALSE,n)
  keep.A[o.A[loA:hiA]] <- TRUE
  keep <- keep.M & keep.A
  M <- M[keep]

  if(Weighting) {
    obs.p <- obs.p[keep]
    ref.p <- ref.p[keep]
    v <- (1-obs.p)/obs.p/libsize.obs + (1-ref.p)/ref.p/libsize.ref
    w <- (1+1e-6) / (v+1e-6)
    TMM <- sum(w*M) / sum(w)
  } else {
    TMM <- mean(M)
  }

  2^TMM
}


quantile_normalisation <- function(df,het.table=NULL,verbose=verbose){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort,na.last=TRUE))
  df_mean <- apply(df_sorted, 1, mean,na.rm=TRUE)
  if(!is.null(het.table)){
    het.table<-het.table[,-c(1:4)]
  }
  if(verbose){
    message("\ncalculating normalized depth")
    df_final <- lapply_pb(1:ncol(df_rank), index_to_mean, my_mean=df_mean,indx=df_rank,al=het.table)
  } else {
    df_final <- lapply(1:ncol(df_rank), index_to_mean, my_mean=df_mean,indx=df_rank,al=het.table)
  }
  return(df_final)
}


index_to_mean <- function(x,indx, my_mean, al=NULL){
  my_index<-indx[,x]
  nr<-my_index[my_mean]
  if(!is.null(al)){
    alleles<-al[,x]
    sm<-do.call(cbind,lapply(data.table::tstrsplit(alleles,","),as.numeric))
    af<-proportions(sm,1)
    af<-round(af*nr,0)
    af<-paste0(af[,1],",",af[,2])
    af[af=="NaN,NaN" | af=="NA,NA"] <- NA
    return(af)
  } else {
    return(nr)
  }
}


#' Calculate normalization factor for each sample
#'
#' This function calculates the normalization factor for each sample using different methods. See details.
#'
#' @param df a data frame or matrix of allele depth values (total depth per snp per sample)
#' @param method character. method to be used (see detials). Default="TMM"
#' @param logratioTrim numeric. percentage value (0 - 1) of variation to be trimmed in log transformation
#' @param sumTrim numeric. amount of trim to use on the combined absolute levels ("A" values) for method="TMM"
#' @param Weighting logical, whether to compute (asymptotic binomial precision) weights
#' @param Acutoff numeric, cutoff on "A" values to use before trimming
#'
#' @details Originally described for normalization of RNA sequences (Robinson & Oshlack 2010), this function computes normalization (scaling) factors to convert observed library sizes into effective library sizes. It uses the method trimmed means of M-values proposed by Robinson & Oshlack (2010). See the original publication and edgeR package for more information.
#'
#' @return Returns a numerical vector of normalization factors for each sample
#'
#' @author Piyal Karunarathne
#'
#' @references Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25
#' Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path)
#' df<-hetTgen(vcf,"AD-tot",verbose=FALSE)
#' norm.fact(df)
#'
#'
#' @export
norm.fact<-function(df,method=c("TMM","TMMex"),logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10){
  if(setequal(c("CHROM","POS","ID"),colnames(df)[1:3])){df<-df[,-c(1:4)]}
  method<-match.arg(method,several.ok = TRUE)
  if(length(method)>1){ifelse(nrow(df)<10001, assign("method","TMM"),assign("method","TMMex"))}
  if(method=="TMM"){
    f75 <- suppressWarnings(calcFactorQuantile(data=as.matrix(df), lib.size=colSums(as.matrix(df)), p=0.75))
    if(median(f75) < 1e-20) {
      ref <- which.max(colSums(sqrt(as.matrix(df))))
    } else {
      ref <- which.min(abs(f75-mean(f75)))
    }
    out<-apply(df,2,FUN=TMM,ref=df[,ref],logratioTrim=logratioTrim,sumTrim=sumTrim,Weighting=Weighting,Acutoff=Acutoff)
  } else if(method=="TMMex"){
    ref<-which.max(colSums(sqrt(as.matrix(df))))
    out<-apply(df,2,FUN=TMMex,ref=df[,ref],logratioTrim=logratioTrim,sumTrim=sumTrim,Weighting=Weighting,Acutoff=Acutoff)
  }
  out<-data.frame(lib.size=colSums(df),norm.factor=out)
  return(out)
}



#' Calculate normalized depth for alleles
#'
#' This function outputs the normalized depth values separately for each allele, calculated using normalization factor with trimmed mean of M-values of sample libraries, median ratios normalization or quantile normalization, See details.
#'
#' @param het.table a data frame or matrix of coverage (output of hetTgen with "AD")
#' @param method character. method to be used (see detials). Default="TMM"
#' @param logratioTrim numeric. percentage value (0 - 1) of variation to be trimmed in log transformation
#' @param sumTrim numeric. amount of trim to use on the combined absolute levels ("A" values) for method="TMM"
#' @param Weighting logical, whether to compute (asymptotic binomial precision) weights
#' @param Acutoff numeric, cutoff on "A" values to use before trimming
#' @param verbose logical. show progress
#'
#' @details This function converts an observed depth value table to an effective depth value table using several normaliztion methods; 1. TMM normalization (See the original publication for more information). It is different from the function normz only in calculation of the counts per million is for separate alleles instead of the total depth. The "TMMex" method is an extension of the "TMM" method for large data sets containing SNPs exceeding 10000; 2. The method "MedR" is median ratio normalization; 3. QN - quantile normalization (see  Maza, Elie, et al. 2013 for a comparison of methods).
#'
#' @return Returns a list with (AD) a data frame of normalized depth values similar to the output of hetTgen function and (outliers) a list of outlier sample names
#'
#' @author Piyal Karunarathne, Qiujie Zhou
#'
#' @references Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25 \n
#' Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140 \n
#' Maza, Elie, et al. "Comparison of normalization methods for differential gene expression analysis in RNA-Seq experiments: a matter of relative size of studied transcriptomes." Communicative & integrative biology 6.6 (2013): e25849
#'
#' @examples
#' \dontrun{data(ADtable)
#' ADnormalized<-cpm.normal(ADtable)}
#'
#'
#' @export
cpm.normal<-function(het.table, method=c("TMM","TMMex","MedR","QN"),logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10,verbose=TRUE){
  method<-match.arg(method)
  if(verbose){
    message("calculating normalization factor")
    tdep<-apply_pb(het.table[,-c(1:4)],2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))),na.rm=TRUE)}))})
  } else {
    tdep<-apply(het.table[,-c(1:4)],2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))),na.rm=TRUE)}))})
  }
  #find and warn about outliers
  ot<-boxplot.stats(colSums(tdep,na.rm = T))$out
  cl<-rep("dodgerblue",ncol(tdep))
  ot.ind<-which(colnames(tdep) %in% names(ot))
  cl[ot.ind]<-2
  if(length(ot)>0){barplot(colSums(tdep,na.rm = T),col=cl,border=NA,xlab="sample",ylab="total depth")
    message("OUTLIERS DETECTED\nConsider removing the samples:")
    cat(colnames(tdep)[ot.ind])}
  if(method=="TMM" | method=="TMMex"){
    nf<-norm.fact(tdep,method = method,logratioTrim=logratioTrim,sumTrim=sumTrim,Weighting=Weighting,Acutoff=Acutoff)
    if(verbose){
      message("\ncalculating normalized depth")
      out<-apply_pb(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/(nf[,1]*nf[,2]))*1e6,2)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    } else {
      out<-apply(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/(nf[,1]*nf[,2]))*1e6,2)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    }
  } else if(method=="MedR") {
    pseudo<- apply(tdep,1,function(xx){exp(mean(log(as.numeric(xx)[as.numeric(xx)>0])))})
    nf<-  apply(tdep,2,function(xx){median(as.numeric(xx)/pseudo,na.rm=T)})
    if(verbose){
      message("\ncalculating normalized depth")
      out<-apply_pb(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/nf),0)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})
    } else {
      out<-apply(het.table[,-c(1:4)],1, function(X){y<-data.frame(do.call(rbind,strsplit(as.character(X),",")))
      y[,1]<-as.numeric(y[,1]);if(ncol(y)>1){y[,2]<-as.numeric(y[,2])}
      nt<-round((y/nf),0)
      if(ncol(nt)>1){paste0(nt[,1],",",nt[,2])}else{nt[,1]}})

    }
  } else if(method=="QN"){
    out<-t(do.call(cbind,quantile_normalisation(tdep,het.table,verbose=verbose)))
  }
  out<-data.frame(het.table[,c(1:4)],t(out))
  colnames(out)<-colnames(het.table)
  return(list(AD=out,outliers=ot.ind))
}






