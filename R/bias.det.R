# helpers
#1. get allele based p and chi-square values for all the samples
get.pvals0<-function(x,df,p.cal){
  snp1<-df[x,-c(1:4)]
  y<-data.frame(do.call(rbind,strsplit(as.character(unlist(snp1)),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
  rr1<-y[,2]/rowSums(y,na.rm = TRUE)
  snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
  homalt<-sum(rr1==1,na.rm=TRUE)
  homref<-sum(rr1==0,na.rm=TRUE)
  Nsamp <- nrow(snp1het)+homalt+homref
  if(nrow(snp1het)>3){
    propHet<-nrow(na.omit(snp1het))/length(na.omit(rr1))
    medRatio<-median(proportions(as.matrix(snp1het),margin = 1)[,2],na.rm = TRUE)
    p.sum<-p.cal[x,2]
    p.05<-0.5
    p.all<-p.cal[x,1]
    n<-unname(rowSums(snp1het,na.rm = TRUE))
    chi.het<-sum((n*p.sum-snp1het[,2])^2/(n*p.sum),na.rm = TRUE)
    chi.het.sum<-chi.het
    chi.05<-sum((n*p.05-snp1het[,2])^2/(n*p.05),na.rm = TRUE)
    chi.0.5.sum<-chi.05
    chi.all<-sum((n*p.all-snp1het[,2])^2/(n*p.all))
    chi.all.sum<-chi.all
    z <- (n*p.sum-snp1het[,2])/sqrt(n*p.sum*(1-p.sum))
    z.het.sum<-sum(z,na.rm = TRUE)
    z<-pnorm(sum(z,na.rm = TRUE),0,sqrt(nrow(snp1het)))
    z.05 <- (n*p.05-snp1het[,2])/sqrt(n*p.05*(1-p.05))
    z.05.sum<-sum(z.05,na.rm = TRUE)
    z.05<-pnorm(sum(z.05,na.rm = TRUE),0,sqrt(nrow(snp1het)))
    z.all<- (n*p.all-snp1het[,2])/sqrt(n*p.all*(1-p.all))
    z.all.sum<-sum(z.all,na.rm = TRUE)
    z.all<-pnorm(sum(z.all,na.rm = TRUE),0,sqrt(nrow(snp1het)))
    ll<-data.frame(NHet=nrow(snp1het),propHet,medRatio,NHomRef=homref,NHomAlt=homalt,propHomAlt=homalt/Nsamp,Nsamp,
                   pAll=p.all,pHet=p.sum,fis=1-(nrow(snp1het)/(2*(homref+(nrow(snp1het)/2))*(homalt+(nrow(snp1het)/2)))),
                   z.het=ifelse(z>0.5, (1-z)*2, z*2),
                   z.05=ifelse(z.05>0.5, (1-z.05)*2, z.05*2),
                   z.all=ifelse(z.all>0.5, (1-z.all)*2, z.all*2),
                   chi.het=pchisq(chi.het,nrow(snp1het)-1,lower.tail=F),
                   chi.05=pchisq(chi.05,nrow(snp1het)-1,lower.tail = F),
                   chi.all=pchisq(chi.all,nrow(snp1het)-1,lower.tail = F),
                   z.het.sum,z.05.sum,z.all.sum,chi.het.sum,chi.0.5.sum,chi.all.sum)
  } else {
    ll<-NA
  }
  return(ll)
}

get.pvals<-function(x,df,p.cal){
  snp1<-unlist(df[x,-c(1:4)])
  y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
  y[,1]<-as.integer(y[,1]);y[,2]<-as.integer(y[,2])
  rr1<-y[,2]/rowSums(y,na.rm = TRUE)
  snp1het<-y[-which(rr1 == 0L | rr1 == 1L | is.na(rr1)==TRUE),]
  homalt<-sum(rr1==1L,na.rm=TRUE)
  homref<-sum(rr1==0L,na.rm=TRUE)
  NHet<-nrow(na.omit(snp1het))
  Nsamp <- NHet+homalt+homref
  if(NHet>3L){
    propHet<-NHet/length(na.omit(rr1))
    medRatio<-median(proportions(as.matrix(snp1het),margin = 1)[,2],na.rm = TRUE)
    p.sum<-p.cal[x,2]
    p.05<-0.5
    p.all<-p.cal[x,1]
    n<-unname(rowSums(snp1het,na.rm = TRUE))

    chi.het<-sum((((n*p.sum)-snp1het[,2])^2L/n*p.sum)+(((n*(1L-p.sum))-snp1het[,1])^2L/n*(1L-p.sum)),na.rm = TRUE)
    chi.het.sum<-chi.het
    chi.05<-sum((((n*p.05)-snp1het[,2])^2L/n*p.05)+(((n*(1L-p.05))-snp1het[,1])^2L/n*(1L-p.05)),na.rm = TRUE)
    chi.05.sum<-chi.05
    chi.all<-sum((((n*p.all)-snp1het[,2])^2L/n*p.all)+(((n*(1L-p.all))-snp1het[,1])^2L/n*(1L-p.all)),na.rm = TRUE)
    chi.all.sum<-chi.all

    z <- (n*p.sum-snp1het[,2])/sqrt(n*p.sum*(1L-p.sum))
    z.het.sum<-sum(z,na.rm = TRUE)
    z<-pnorm(z.het.sum,0,sqrt(NHet))
    z.05 <- (n*p.05-snp1het[,2])/sqrt(n*p.05*(1L-p.05))
    z.05.sum<-sum(z.05,na.rm = TRUE)
    z.05<-pnorm(z.05.sum,0,sqrt(NHet))
    z.all<- (n*p.all-snp1het[,2])/sqrt(n*p.all*(1L-p.all))
    z.all.sum<-sum(z.all,na.rm = TRUE)
    z.all<-pnorm(z.all.sum,0,sqrt(NHet))
    ll<-data.frame(NHet=NHet,propHet,medRatio,NHomRef=homref,NHomAlt=homalt,propHomAlt=homalt/Nsamp,Nsamp,
                   pAll=p.all,pHet=p.sum,fis=1L-(NHet/(2L*(homref+(NHet/2L))*(homalt+(NHet/2L)))),
                   z.het=ifelse(z>0.5, (1L-z)*2L, z*2L),
                   z.05=ifelse(z.05>0.5, (1L-z.05)*2L, z.05*2L),
                   z.all=ifelse(z.all>0.5, (1L-z.all)*2L, z.all*2L),
                   chi.het=pchisq(chi.het,NHet-1L,lower.tail=F),
                   chi.05=pchisq(chi.05,NHet-1L,lower.tail = F),
                   chi.all=pchisq(chi.all,NHet-1L,lower.tail = F),
                   z.het.sum,z.05.sum,z.all.sum,chi.het.sum,chi.05.sum,chi.all.sum)
  } else {
    ll<-NA
  }
  return(ll)
}


#' Get allele information for duplicate detection
#'
#' The function to calculate allele median ratios, proportion of heterozygotes
#'  and allele probability values under different assumptions (see details),
#'   and their chi-square significance values for duplicate detection
#'
#' @param X allele depth table generated from the function
#' \code{hetTgen} (non-normalized)
#' @param x.norm a data frame of normalized allele coverage, output of
#' \code{cpm.normal}. If not provided, calculated using \code{X}.
#' @param method character. method to be used for normalization
#' (see \code{cpm.normal} details). Default \code{TMM}
#' @param logratioTrim numeric. percentage value (0 - 1) of variation to be
#'  trimmed in log transformation
#' @param sumTrim numeric. amount of trim to use on the combined absolute
#'  levels (\dQuote{A} values) for method \code{TMM}
#' @param Weighting logical, whether to compute (asymptotic binomial precision)
#'  weights
#' @param Acutoff numeric, cutoff on \dQuote{A} values to use before trimming
#' @param plot.allele.cov logical, plot comparative plots of allele depth
#' coverage in homozygotes and heterozygotes
#' @param verbose logical, whether to print progress
#' @param \dots further arguments to be passed to \code{plot}
#'
#' @importFrom stats pchisq pnorm na.omit
#'
#' @details
#' Allele information generated here are individual SNP based and presents the
#'  proportion of heterozygotes, number of samples, and deviation of allele
#'  detection from a 1:1 ratio of reference and alternative alleles.
#'  The significance of the deviation is tested with Z-score test
#'  \eqn{Z = \frac{ \frac{N}{2}-N_A}{ \sigma_{x}}},
#'  and chi-square test (see references for more details on the method).
#'
#' @return Returns a data frame of median allele ratio, proportion of
#' heterozygotes, number of heterozygotes, and allele probability at different
#'  assumptions with their chi-square significance
#'
#' @author Piyal Karunarathne, Pascal Milesi, Qiujie Zhou
#'
#' @references
#' \itemize{
#'  \item{McKinney, G. J., Waples, R. K., Seeb, L. W., & Seeb, J. E. (2017).
#'  Paralogs are revealed by proportion of heterozygotes and deviations in read
#'   ratios in genotyping by sequencing data from natural populations.
#'    Molecular Ecology Resources, 17(4)}
#'  \item{Karunarathne et al. 2022 (to be added)}
#' }
#'
#' @examples
#' \dontrun{data(ADtable)
#' AI<-allele.info(ADtable,x.norm=ADnorm)}
#'
#' @export
allele.info<-function(X,x.norm=NULL,method=c("TMM","TMMex","MedR","QN","pca"),logratioTrim = 0.3,sumTrim = 0.05,Weighting = TRUE,Acutoff = -1e+10,plot.allele.cov=TRUE,verbose = TRUE,...){
  method=match.arg(method)
  if(is.null(x.norm)){
    x.norm<-cpm.normal(X,method=method,logratioTrim=logratioTrim,sumTrim = sumTrim,Weighting = Weighting,Acutoff = Acutoff,verbose = verbose)
  }
  if(!inherits(x.norm,"list")){x.norm<-list(AD=x.norm)}
  if(inherits(x.norm,"list")){x.norm<-x.norm$AD}

  if(verbose){
    message("calculating probability values of alleles")
    p.cal<-apply_pb(x.norm[,-c(1:4)],1,function(snp1){
      if(is.character(unname(unlist(snp1[1])))){
        y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
        y[,1]<-as.integer(y[,1])
        y[,2]<-as.integer(y[,2])} else {y<-snp1}
      rs<-rowSums(y)
      rs[rs==0L]<-NA
      cv<-sd(unlist(y),na.rm = TRUE)/mean(unlist(y),na.rm = TRUE)
      rr1<-y[,2]/rs
      snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
      homalt<-sum(rr1==1,na.rm=TRUE)
      homref<-sum(rr1==0,na.rm=TRUE)
      covrefhomo<-sum(y[c(rr1 == 0,na.rm = TRUE),],na.rm = TRUE)
      covalthomo<-sum(y[c(rr1 == 1,na.rm = TRUE),],na.rm = TRUE)
      covalt<-sum(y[,2],na.rm = TRUE)
      covref<-sum(y[,1],na.rm = TRUE)
      NHet<-nrow(snp1het)
      if(NHet>3){
        p.all<-(covalt/(NHet+(2*homalt)))/((covalt/(NHet+(2*homalt)))+(covref/(NHet+(2*homref))))
        p.sum<-sum(snp1het[,2])/sum(snp1het)
        ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/NHet,mean.r.het=sum(snp1het[,1])/NHet,cv=cv)
      } else {
        ll<-NA
      }
      return(ll)
    })
  } else {
    p.cal<-apply(x.norm[,-c(1:4)],1,function(snp1){
      if(is.character(unname(unlist(snp1[1])))){
        y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
        y[,1]<-as.integer(y[,1])
        y[,2]<-as.integer(y[,2])} else {y<-snp1}
      rs<-rowSums(y)
      rs[rs==0]<-NA
      cv<-sd(unlist(y),na.rm = TRUE)/mean(unlist(y),na.rm = TRUE)
      rr1<-y[,2]/rs
      snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
      homalt<-sum(rr1==1,na.rm=TRUE)
      homref<-sum(rr1==0,na.rm=TRUE)
      covrefhomo<-sum(y[c(rr1 == 0,na.rm = TRUE),],na.rm = TRUE)
      covalthomo<-sum(y[c(rr1 == 1,na.rm = TRUE),],na.rm = TRUE)
      covalt<-sum(y[,2],na.rm = TRUE)
      covref<-sum(y[,1],na.rm = TRUE)
      NHet<-nrow(snp1het)
      if(NHet>3){
        p.all<-(covalt/(NHet+(2*homalt)))/((covalt/(NHet+(2*homalt)))+(covref/(NHet+(2*homref))))
        p.sum<-sum(snp1het[,2])/sum(snp1het)
        ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/NHet,mean.r.het=sum(snp1het[,1])/NHet,cv=cv)
      } else {
        ll<-NA
      }
      return(ll)
    })
  }
  if(is.list(p.cal)){
    p.cal<-do.call(rbind,p.cal)
  } else {
    p.cal<-t(p.cal)
  }
  p.cal[p.cal=="NaN"]<-0
  if(plot.allele.cov){
    p.list<-list(...)
    if(is.null(p.list$pch)) p.list$pch=19
    if(is.null(p.list$cex)) p.list$cex=0.6
    if(is.null(p.list$col)) p.list$col<-makeTransparent(colorspace::heat_hcl(12,h=c(0,-100),c=c(40,80),l=c(75,40),power=1)[11])
    if(is.null(p.list$lcol)) p.list$lcol="tomato"
    par(mfrow=c(2,2))
    par(mar=c(4,5,2,2))
    plot(p.cal$mean.a.homo,p.cal$mean.r.homo,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nalt. allele in homozygotes",
         ylab="Mean coverage of \n ref. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.a.het,p.cal$mean.r.het,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nalt. allele in heterozygotes",
         ylab="Mean coverage of \n ref. allele in heterozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.a.het,p.cal$mean.a.homo,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nalt. allele in heterozygotes",
         ylab="Mean coverage of \n alt. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.r.het,p.cal$mean.r.homo,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nref. allele in heterozygotes",
         ylab="Mean coverage of \n ref. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    par(mfrow=c(1,1))
  }

  if(verbose){
    message("calculating chi-square significance")
    pvals<-lapply_pb(1:nrow(X),get.pvals,df=X,p.cal=p.cal)
  } else {
    pvals<-lapply(1:nrow(X),get.pvals,df=X,p.cal=p.cal)
  }
  pvals<-do.call(rbind,pvals)
  pvals<-cbind(X[,1:3],pvals)
  pvals<-na.omit(pvals)
  ht<-sig.hets(pvals,plot = FALSE, verbose = verbose)
  pvals<-data.frame(pvals,eH.pval=ht[,"eH.pval"],eH.delta=ht[,"eH.delta"],cv=na.omit(p.cal[,"cv"]))

  return(pvals)
}

