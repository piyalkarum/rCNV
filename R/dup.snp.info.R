#helpers
#1. generate dupinfo for each snp
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
  n.ref<- length(coverage.hm0)+(length(coverage.ht)/2)
  n.alt<- length(coverage.hm1)+(length(coverage.ht)/2)
  num.rare.al <- pmin(n.ref,n.alt)
  num_sample <- length(coverage.ht)+length(coverage.hmm)
  num.het <- length(coverage.ht);num.hm <- dim(hm0)[1]
  med.hmm <- median(coverage.hmm)
  prop_ht <- num.het/num_sample;prop_hm <- num.hm/num_sample
  prop_hm1 <- length(coverage.hm1)/num_sample
  Hobs <- prop_ht;p<-prop_hm+(prop_ht/2);q<-prop_hm1+(prop_ht/2);Hexp<-2*p*q;fis<-1-Hobs/Hexp

  c(ht.rts,med.hmm, num.het, prop_hm, prop_ht, prop_hm1, num.rare.al,num.hm,length(coverage.hm1),fis,num_sample,dps,nas)
}

dup.info2 <- function(gt,nf=1,l.size=NULL){
  y<-data.frame(do.call(rbind,strsplit(as.character(gt),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2]);if(!is.null(l.size)) y<-(y/(nf*l.size))*1e6
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
  n.ref<- length(coverage.hm0)+(length(coverage.ht)/2)
  n.alt<- length(coverage.hm1)+(length(coverage.ht)/2)
  num.rare.al <- pmin(n.ref,n.alt)
  num_sample <- length(coverage.ht)+length(coverage.hmm)
  num.het <- length(coverage.ht);num.hm <- dim(hm0)[1]
  med.hmm <- median(coverage.hmm)
  prop_ht <- num.het/num_sample;prop_hm <- num.hm/num_sample
  prop_hm1 <- length(coverage.hm1)/num_sample
  Hobs <- prop_ht;p<-prop_hm+(prop_ht/2);q<-prop_hm1+(prop_ht/2);Hexp<-2*p*q;fis<-1-Hobs/Hexp

  c(ht.rts,med.hmm, num.het, prop_hm, prop_ht, prop_hm1, num.rare.al,num.hm,length(coverage.hm1),fis,num_sample,dps,nas)
}


#' Generate SNPs duplication information
#'
#' This function generates duplication information of all the SNPs provided in the data set.
#' dup.snp.info takes a tab separated depth values Vs. SNPs by sample/population table and calculates numbers and proportions of hetero/homo-zygotes for each snp.
#' The resulting table is used in the downstream analysis such as CNV/duplication detection.
#'
#' @param het.table tab-delimited depth value by sample/population (output of hetTgen)
#' @param normalize logical. Whether to apply normalization to depth values using trimmed mean of M-values. See ref.
#' @param verbose logical. if TRUE, shows the progress
#' @importFrom stats na.omit
#'
#' @author Piyal Karunarathne
#'
#' @references
#' Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25.
#'
#' @examples
#' print("deprecated")
#'
#' @export
dup.snp.info<-function(het.table,normalize=FALSE,verbose=TRUE){
  gts<-het.table[,-c(1:3)]
  res<-het.table[,1:3]
  if(normalize){
    gtst<-apply(gts,2,function(x){do.call(cbind,lapply(x,function(y){sum(as.numeric(unlist(strsplit(as.character(y),","))))}))})
    suppressWarnings(nf<-norm.fact(gtst,ifelse(nrow(gtst)>100000,"TMMex","TMM")))
    if(verbose){
      message("normalizing depth of coverage")
      suppressWarnings(out<-t(apply_pb(gts,MARGIN = 1,dup.info2,nf=nf[,2],l.size=nf[,1])))
    } else {
      suppressWarnings(out<-t(apply(gts,MARGIN = 1,dup.info2,nf=nf[,2],l.size=nf[,1])))
    }
  } else {
    if(verbose){
      suppressWarnings(out<-t(apply_pb(gts,MARGIN = 1,dup.info)))
    } else {
      suppressWarnings(out<-t(apply(gts,MARGIN = 1,dup.info)))
    }
  }
  snp.dup<-data.frame(res$CHROM,res$POS,res$ID,out)#
  colnames(snp.dup)<-c("Scaffold","Position","ID","MedRatio","AvgRatio","MedCovHet","TotCovHet","MedCovHom","NumHet","PropHomFreq","PropHet","PropHomRare","NoRareAllele","NHomFreq","NHomRare","Fis","truNsample","totDepFreq","totDepRare","totDepFreqHet","totDepRareHet","medDepth","numMiss")#
  snp.dup <- na.omit(snp.dup)
  return(snp.dup)
}








