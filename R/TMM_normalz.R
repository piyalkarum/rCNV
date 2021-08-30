### This set of functions were adopted from edgeR package (Robinson et al. 2010) unless otherwise stated.
### The methodology behind the normalization was originally described by Mark Robinson (2010)

calcFactorQuantile <- function (data, lib.size, p=0.75)
  #	Generalized version of upper-quartile normalization
  #	Mark Robinson and Gordon Smyth
  #	adopted from EdgeR package
{
  f <- rep_len(1,ncol(data))
  for (j in seq_len(ncol(data))) f[j] <- quantile(data[,j], probs=p)
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


#' Calculate normalization factor for each sample
#'
#' This function calculates the normalization factor for each sample using trimmed mean of M-values of sample libraries. See details.
#'
#' @param df a data frame or matrix of count values (total counts per snp per sample)
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
#' print("to be added")
#'
#' @export
norm.fact<-function(df,method=c("TMM","TMMex"),logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10){
  method<-match.arg(method)
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

#' Calculate normalized depth values
#'
#' This function outputs the normalized depth values calculated using normalization factor with trimmed mean of M-values of sample libraries. See details.
#'
#' @param df a data frame or matrix of count values (total counts per snp per sample)
#' @param method character. method to be used (see detials). Default="TMM"
#' @param logratioTrim numeric. percentage value (0 - 1) of variation to be trimmed in log transformation
#' @param sumTrim numeric. amount of trim to use on the combined absolute levels ("A" values) for method="TMM"
#' @param Weighting logical, whether to compute (asymptotic binomial precision) weights
#' @param Acutoff numeric, cutoff on "A" values to use before trimming
#' @param tot.count numeric, a total library size that depth values to be normalized to, default=1000000
#'
#' @details This function converts an observed depth value table to an effective depth vallue table using TMM normalization. See the original publication for more information.
#'
#' @return Returns a data frame of normalized depth values
#'
#' @author Piyal Karunarathne
#'
#' @references Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25
#' Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#'
#' @examples
#' print("to be added")
#'
#' @export
normz<-function(df,method=c("TMM","TMMex"),logratioTrim=.3, sumTrim=0.05, Weighting=TRUE, Acutoff=-1e10, tot.count=1e6){
  df<-as.matrix(df)
  nf<-norm.fact(df,method = method,logratioTrim=logratioTrim,sumTrim=sumTrim,Weighting=Weighting,Acutoff=Acutoff)
  cpm <- lapply(1:ncol(df), function(x,df,nf) {(df[,x]/(sum(df[,x]))*nf[x])*tot.count},df=df,nf=nf)
  cpm<-do.call(cbind,cpm)
  colnames(cpm)<-colnames(df)
  return(cpm)
}
