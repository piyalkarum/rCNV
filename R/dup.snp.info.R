#helpers
#1. TTM normalization
####### from edgeR #######

#' @importFrom methods is
#' @keywords internal
TMM<-function (object, method = c("TMM", "RLE", "upperquartile", "none"),
               refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE,
               Acutoff = -1e+10, p = 0.75)
{
  if (is(object, "DGEList")) {
    x <- as.matrix(object$counts)
    lib.size <- object$samples$lib.size
  }
  else {
    x <- as.matrix(object)
    lib.size <- colSums(x)
  }
  method <- match.arg(method)
  allzero <- rowSums(x > 0) == 0
  if (any(allzero))
    x <- x[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 1)
    method = "none"
  f <- switch(method, TMM = {
    f75 <- .calcFactorQuantile(data = x, lib.size = lib.size,
                               p = 0.75)
    if (is.null(refColumn)) refColumn <- which.min(abs(f75 -
                                                         mean(f75)))
    if (length(refColumn) == 0 | refColumn < 1 | refColumn >
        ncol(x)) refColumn <- 1
    f <- rep(NA, ncol(x))
    for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted(obs = x[,
                                                             i], ref = x[, refColumn], libsize.obs = lib.size[i],
                                                     libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                     sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)
    f
  }, RLE = .calcFactorRLE(x)/lib.size, upperquartile = .calcFactorQuantile(x,
                                                                           lib.size, p = p), none = rep(1, ncol(x)))
  f <- f/exp(mean(log(f)))
  if (is(object, "DGEList")) {
    object$samples$norm.factors <- f
    return(object)
  }
  else {
    return(f)
  }
}

.calcFactorQuantile <- function (data, lib.size, p = 0.75)
{
  y <- t(t(data)/lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
}
.calcFactorRLE<-function (data)
{
  gm <- exp(rowMeans(log(data)))
  apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorWeighted <- function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
                                 sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
{
  if (all(obs == ref))
    return(1)
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  if (is.null(libsize.obs))
    nO <- sum(obs)
  else nO <- libsize.obs
  if (is.null(libsize.ref))
    nR <- sum(ref)
  else nR <- libsize.ref
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  n <- sum(fin)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >=
                                                       loS & rank(absE) <= hiS)
  if (doWeighting)
    2^(sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],
                                                 na.rm = TRUE))
  else 2^(mean(logR[keep], na.rm = TRUE))
}
##########################
normz<-function(dat){
  DD<-apply(dat,2,function(xx){(dd<-do.call(rbind,lapply(xx,function(x){yy<-strsplit(x,",");sum(as.numeric(unlist(yy)))})))})
  colnames(DD)<-gsub(".AD","",colnames(DD))
  return(TMM(DD))### consider avoiding package dependency
}

#2. generate dupinfo for each snp
dup.info <- function(gt,nf=1){
  y<-data.frame(do.call(rbind,strsplit(as.character(gt),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2]);y<-y*nf
  hm0 <- na.omit(y[y[,1]>0 & y[,2]==0,]);hm1 <- na.omit(y[y[,1]==0 & y[,2]>0,]);ht <- na.omit(y[y[,1]>0 & y[,2]>0,]);hmm<-rbind(hm0,hm1)
  nas<-dim(y[y[,1]==0 & y[,2]==0,])[1]
  hm0.ratio<-proportions(as.matrix(hm0),margin = 1)[,1];hm1.ratio<-proportions(as.matrix(hm1),margin = 1)[,2]
  ht.ratio<-proportions(as.matrix(ht),margin = 1)[,2]
  med_ratio<-median(ht.ratio,na.rm = T);mn_ratio<-mean(ht.ratio,na.rm = T)
  all.ratio<-c(hm0.ratio,hm1.ratio,ht.ratio);med.all <- median(all.ratio,na.rm = T)
  coverage.ht<-rowSums(ht)
  ht.rts<-c(med_ratio,mn_ratio,median(coverage.ht),sum(coverage.ht))
  tot.dp.ht0<-colSums(ht)[1];tot.dp.ht1<-colSums(ht)[2]
  coverage.hm0<-rowSums(hm0);coverage.hmm <- rowSums(hmm);coverage.hm1<-rowSums(hm1)
  tot.dp.hm0<-sum(coverage.hm0);tot.dp.hm1<-sum(coverage.hm1);tot.dp.0<-tot.dp.hm0+tot.dp.ht0
  tot.dp.1<-tot.dp.hm1+tot.dp.ht1;med_depth<-median(c(coverage.ht,coverage.hmm))
  dps<-c(tot.dp.0,tot.dp.1,tot.dp.ht0,tot.dp.ht1,med_depth)
  num.rare <- length(coverage.hm1)+length(coverage.ht)
  num_sample <- length(coverage.ht)+length(coverage.hmm)
  num.het <- length(coverage.ht);num.hm <- dim(hm0)[1]
  med.hmm <- median(coverage.hmm)
  prop_ht <- num.het/num_sample;prop_hm <- num.hm/num_sample
  prop_hm1 <- length(coverage.hm1)/num_sample
  Hobs <- prop_ht;p<-prop_hm+(prop_ht/2);q<-prop_hm1+(prop_ht/2);Hexp<-2*p*q;fis<-1-Hobs/Hexp

  c(ht.rts,med.hmm, num.het, prop_hm, prop_ht, prop_hm1, num.rare,num.hm,length(coverage.hm1),fis,num_sample,dps,nas)
}


#' Generate SNPs duplication information
#'
#' This function generates duplication information of all the SNPs provided in the data set.
#' dup.snp.info takes a tab separated depth values Vs. SNPs by sample/population table and calculates numbers and proportions of hetero/homo-zygotes for each snp.
#' The resulting table is used in the downstream analysis such as CNV/duplication detection.
#'
#' @param het.table tab-delimited depth value by sample/population (output of hetTgen)
#' @param normalize logical. Whether to apply normalization to depth values using trimmed mean of M-values. See ref.
#' @importFrom stats na.omit
#'
#' @author Piyal Karunarathne
#'
#' @references
#' Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25.
#'
#' @examples
#' data(hets)
#' dup.info <- dup.snp.info(het.table=hets,normalize=FALSE)
#'
#' @export
dup.snp.info<-function(het.table,normalize=FALSE){
  gts<-het.table[,-c(1:3)]
  res<-het.table[,1:3]
  if(normalize){
    nf<-normz(gts)
    suppressWarnings(out<-t(apply_pb(gts,MARGIN = 1,dup.info,nf=nf)))
  } else {
    suppressWarnings(out<-t(apply_pb(gts,MARGIN = 1,dup.info)))
  }
  snp.dup<-data.frame(res$CHROM,res$POS,res$ID,out)#
  colnames(snp.dup)<-c("Scaffold","Position","ID","MedRatio","AvgRatio","MedCovHet","TotCovHet","MedCovHom","NumHet","PropHomFreq","PropHet","PropHomRare","NumRare","NHomFreq","NHomRare","Fis","truNsample","totDepFreq","totDepRare","totDepFreqHet","totDepRareHet","medDepth","numMiss")#
  snp.dup <- na.omit(snp.dup)
  return(snp.dup)
}








