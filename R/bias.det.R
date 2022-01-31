# helpers
get.pvals<-function(x,df,p.cal){
  snp1<-df[x,-c(1:4)]
  y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
  rr1<-y[,2]/rowSums(y)
  snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==T),]
  homalt<-sum(rr1==1,na.rm=T)
  homref<-sum(rr1==0,na.rm=T)
  Nsamp=nrow(snp1het)+homalt+homref
  if(nrow(snp1het)>3){
    propHet<-nrow(na.omit(snp1het))/length(na.omit(rr1))
    medRatio<-median(proportions(as.matrix(snp1het),margin = 1)[,2],na.rm = T)
    p.sum<-p.cal[x,2]
    p.05<-0.5
    p.all<-p.cal[x,1]
    n<-unname(rowSums(snp1het))
    chi.het<-sum((n*p.sum-snp1het[,2])^2/(n*p.sum))
    chi.05<-sum((n*p.05-snp1het[,2])^2/(n*p.05))
    chi.all<-sum((n*p.all-snp1het[,2])^2/(n*p.all))
    z <- (n*p.sum-snp1het[,2])/sqrt(n*p.sum*(1-p.sum))
    z<-pnorm(sum(z),0,sqrt(nrow(snp1het)))
    z.05 <- (n*p.05-snp1het[,2])/sqrt(n*p.05*(1-p.05))
    z.05<-pnorm(sum(z.05),0,sqrt(nrow(snp1het)))
    z.all<- (n*p.all-snp1het[,2])/sqrt(n*p.all*(1-p.all))
    z.all<-pnorm(sum(z.all),0,sqrt(nrow(snp1het)))
    ll<-data.frame(NHet=nrow(snp1het),propHet,medRatio,NHomRef=homref,NHomAlt=homalt,propHomAlt=homalt/Nsamp,Nsamp,
                   pAll=p.all,pHet=p.sum,fis=1-(nrow(snp1het)/(2*(homref+(nrow(snp1het)/2))*(homalt+(nrow(snp1het)/2)))),
                   z.het=ifelse(z>0.5, (1-z)*2, z*2),
                   z.05=ifelse(z.05>0.5, (1-z.05)*2, z.05*2),
                   z.all=ifelse(z.all>0.5, (1-z.all)*2, z.all*2),
                   chi.het=pchisq(chi.het,nrow(snp1het)-1,lower.tail=F),
                   chi.05=pchisq(chi.05,nrow(snp1het)-1,lower.tail = F),
                   chi.all=pchisq(chi.all,nrow(snp1het)-1,lower.tail = F))
  } else {
    ll<-NA
  }
  return(ll)
}


#' Get allele info. for duplicate detection
#'
#' The function to calculate allele median ratios, proportion of heterozygotes and allele probability values under different assumptions (see details), and their chi-square significance values for duplicate detection
#'
#' @param X a data frame of allele coverage (non-normalized). output of hetTgen()
#' @param x.norm a data frame of normalized allele coverage. output of cpm.normal(). If not provided, calculated using X.
#' @param method character. method to be used for normalization (see cpm.normal detials). Default="TMM" for data <100,000 snps and "TMMex" for >100,000
#' @param logratioTrim numeric. percentage value (0 - 1) of variation to be trimmed in log transformation
#' @param sumTrim numeric. amount of trim to use on the combined absolute levels ("A" values) for method="TMM"
#' @param Weighting logical, whether to compute (asymptotic binomial precision) weights
#' @param Acutoff numeric, cutoff on "A" values to use before trimming
#' @param plot.allele.cov logical, plot comparative plots of allele depth coverage in homozygotes and heterozygotes
#' @param verbose logical, whether to print progress
#' @param ... further arguments to be passed to plot
#'
#' @importFrom stats pchisq pnorm na.omit
#'
#' @details
#' print("to be added")
#'
#' @return Returns a data frame of median allele ratio, proportion of heterozygotes, number of heterozygotes, and allele probability at different assumptions with their chi-square significance
#'
#' @author Piyal Karunarathne, Pascal Milesi, Qiujie Zhou
#'
#' @examples
#' data(ADtable)
#' AI<-allele.info(ADtable)
#'
#' @export
allele.info<-function(X,x.norm=NULL,method=c("TMM", "TMMex"),logratioTrim = 0.3,sumTrim = 0.05,Weighting = TRUE,Acutoff = -1e+10,plot.allele.cov=TRUE,verbose = TRUE,...){
  method=match.arg(method)
  if(is.null(x.norm)){
    x.norm<-cpm.normal(X,method=method,logratioTrim=logratioTrim,sumTrim = sumTrim,Weighting = Weighting,Acutoff = Acutoff,verbose = verbose)
  }
  if(verbose){
    message("calculating probability values of alleles")
    p.cal<-apply_pb(x.norm[,-c(1:4)],1,function(snp1){
      y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
      rr1<-y[,2]/rowSums(y)
      snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==T),]
      homalt<-sum(rr1==1,na.rm=T)
      homref<-sum(rr1==0,na.rm=T)
      covrefhomo<-sum(y[c(rr1 == 0,na.rm = T),],na.rm = T)
      covalthomo<-sum(y[c(rr1 == 1,na.rm = T),],na.rm = T)
      covalt<-sum(y[,2])
      covref<-sum(y[,1])
      if(nrow(snp1het)>3){
        p.all<-(colSums(y)[2]/(nrow(snp1het)+(2*homalt)))/((colSums(y)[2]/(nrow(snp1het)+(2*homalt)))+(colSums(y)[1]/(nrow(snp1het)+(2*homref))))
        p.sum<-sum(snp1het[,2])/sum(snp1het)
        ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/nrow(snp1het),mean.r.het=sum(snp1het[,1])/nrow(snp1het))
      } else {
        ll<-NA
      }
      return(ll)
    })
  } else {
    p.cal<-apply(x.norm[,-c(1:4)],1,function(snp1){
      y<-data.frame(do.call(rbind,strsplit(as.character(snp1),",")));y[,1]<-as.numeric(y[,1]);y[,2]<-as.numeric(y[,2])
      rr1<-y[,2]/rowSums(y)
      snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==T),]
      homalt<-sum(rr1==1,na.rm=T)
      homref<-sum(rr1==0,na.rm=T)
      covrefhomo<-sum(y[c(rr1 == 0,na.rm = T),],na.rm = T)
      covalthomo<-sum(y[c(rr1 == 1,na.rm = T),],na.rm = T)
      covalt<-sum(y[,2])
      covref<-sum(y[,1])
      if(nrow(snp1het)>3){
        p.all<-(colSums(y)[2]/(nrow(snp1het)+(2*homalt)))/((colSums(y)[2]/(nrow(snp1het)+(2*homalt)))+(colSums(y)[1]/(nrow(snp1het)+(2*homref))))
        p.sum<-sum(snp1het[,2])/sum(snp1het)
        ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/nrow(snp1het),mean.r.het=sum(snp1het[,1])/nrow(snp1het))
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

  if(plot.allele.cov){
    p.list<-list(...)
    if(is.null(p.list$pch)) p.list$pch=19
    if(is.null(p.list$cex)) p.list$cex=0.6
    if(is.null(p.list$col)) p.list$col<-makeTransparent(colorspace::heat_hcl(12,h=c(0,-100),c=c(40,80),l=c(75,40),power=1)[11])
    if(is.null(p.list$lcol)) p.list$lcol="tomato"
    par(mfrow=c(2,2))
    par(mar=c(4,5,2,2))
    plot(p.cal$mean.a.homo,p.cal$mean.r.homo,pch=p.list$pch,cex=p.list$cex,col=p.list$col,xlab="Mean coverage of \nalt. allele in homozygotes",ylab="Mean coverage of \n ref. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.a.het,p.cal$mean.r.het,pch=p.list$pch,cex=p.list$cex,col=p.list$col,xlab="Mean coverage of \nalt. allele in heterozygotes",ylab="Mean coverage of \n ref. allele in heterozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.a.het,p.cal$mean.a.homo,pch=p.list$pch,cex=p.list$cex,col=p.list$col,xlab="Mean coverage of \nalt. allele in heterozygotes",ylab="Mean coverage of \n alt. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.r.het,p.cal$mean.r.homo,pch=p.list$pch,cex=p.list$cex,col=p.list$col,xlab="Mean coverage of \nref. allele in heterozygotes",ylab="Mean coverage of \n ref. allele in homozygotes",cex.lab=0.8)
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
  return(pvals)
}

