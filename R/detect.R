#1 calculate expected proportions of heterozygotes from real data
ex.prop<-function(rs,method=c("fisher","chi.sq")){
  n<-unlist(c(rs[4]))
  p<-(rs[1]+rs[2]/2)/n
  ob<-unlist(c(rs[1:3]))
  eX<-unlist(c((p^2) * n,2*p*(1-p) * n,((1-p)^2) * n))
  delta <- rs[2]-(2*p*(1-p) * n)
  stat<-match.arg(method)
  pval<-switch(stat,fisher=suppressWarnings(fisher.test(cbind(ob,eX),workspace = 2e8))$p.value,chi.sq=suppressWarnings(chisq.test(cbind(ob,eX)))$p.value)
  return(c(eX/n,pval,delta))
}

#' Identify significantly different heterozygotes from SNPs data
#'
#' This function will recognize the SNPs with a proportion of heterozygotes significantly higher than expected under HWE and plot putatively duplicated snps
#'
#' @param a.info allele info table generated from filtered vcfs using the function '\link[rCNV]{allele.info}
#' @param method character. Method for testing significance. Fisher exact test ("fisher") or Chi squre test ("chi.sq")
#' @param plot logical. Whether to plot the identified duplicated snps with the expected values
#' @param verbose logical, if TRUE, the progress is shown
#' @param ... other arguments passed to \link[graphics]{plot}
#'
#' @return A matrix of expected heterozygote proportions from the observed data with p-value indicating significance of deviation.
#'
#' @author Piyal Karunarathne, Pascal Milesi
#'
#' @examples
#' \dontrun{data(alleleINF)
#' AI <- alleleINF
#' duplicates<-sig.hets(AI,plot=TRUE)}
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats fisher.test median quantile rbinom sd smooth.spline chisq.test
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom colorspace rainbow_hcl
#' @export
sig.hets<-function(a.info,method=c("fisher","chi.sq"),plot=TRUE,verbose=TRUE,...){
  d<-a.info[,c("NHomRef","NHet","NHomAlt","Nsamp")]
  colnames(d)<-c("h1","het","h2","truNsample")
  method<-match.arg(method)
  if(verbose){
    message("assessing excess of heterozygotes")
    df<-data.frame(t(apply_pb(d,1,ex.prop,method=method)))
  } else {
    df<-data.frame(t(apply(d,1,ex.prop,method=method)))
  }
  colnames(df)<-c("p2","het","q2","pval","delta")
  df$dup.stats<-"singlet";df$dup.stats[which(df$pval < 0.05 & df$delta > 0 )]<-"duplicated"
  if(plot){
    l<-list(...)
    if(is.null(l$cex)) l$cex=0.2
    if(is.null(l$pch)) l$pch=19
    if(is.null(l$xlim)) l$xlim=c(0,1)
    if(is.null(l$ylim)) l$ylim=c(0,1)
    if(is.null(l$col)) cols<-makeTransparent(rainbow_hcl(2),alpha=0.3) else cols<-makeTransparent(l$col,alpha=0.3)

    d$Color <- cols[1]
    d$Color [which(df$dup.stats=="duplicated")]<- cols[2]#& df$delta > 0
    plot(a.info$propHet~a.info$propHomAlt, pch=l$pch, cex=l$cex,col=d$Color,xlim=l$xlim,ylim=l$ylim,
         xlab="Proportion of Alternate Homozygotes",ylab="Proportion of Heterozygotes")
    lines((smm<-smooth.spline(df$het~df$q2)),col="blue")
    legend("bottomright", c("singlet","duplicate","expected"), col = c(cols,"blue"), lty = c(0, 0, 1), lwd = c(0, 0, 1),pch = c(l$pch, l$pch, NA),
           cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
  }
  return(data.frame(cbind(a.info[,c(1:3)],df[,c(4:6)]),row.names = NULL))
}


#' Plot duplicates
#'
#' The function plots detected duplicates from functions "sig.snps", and "dupGet"
#'
#' @param ds a data frame of detected duplicates
#' @param ... other graphical parameters to be passed to the function \link[graphics]{plot}
#'
#' @importFrom colorspace rainbow_hcl
#'
#' @return Plots duplicates on proportion of heterozygotes vs allele median ratio
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{data(alleleINF)
#' DD<-dupGet(alleleINF,plot=FALSE)
#' dup.plot(DD)}
#'
#' @export
dup.plot<-function(ds,...){
  l<-list(...)
  if(is.null(l$cex)) l$cex=0.2
  if(is.null(l$pch)) l$pch=19
  if(is.null(l$xlim)) l$xlim=c(0,1)
  if(is.null(l$ylim)) l$ylim=c(0,1)
  if(is.null(l$alpha)) l$alpha=0.3
  if(is.null(l$col)) l$col<-makeTransparent(c("tomato","#2297E6FF"))#colorspace::terrain_hcl(12,c=c(65,0),l=c(45,90),power=c(1/2,1.5))[2]
  ds$Color <- l$col[2]
  ds$Color [ds$dup.stat=="duplicated"]<- l$col[1]
  plot(ds$medRatio~ds$propHet, pch=l$pch, cex=l$cex,col=ds$Color,xlim=l$xlim,ylim=l$ylim,frame=F,
       ylab="Allele Median Ratio",xlab="Proportion of Heterozygotes")
  legend("bottomright", c("duplicates","singlets"), col = makeTransparent(l$col,alpha=1), pch=l$pch,
         cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")

}


#' Detect duplicates from SNPs
#'
#' Detect duplicated snps using excess of heterozygotes (alleles that do not follow HWE) and snp deviates (alleles that do not follow a normal/chi-square distribution). See details.
#'
#' @param data data frame of the output of \link[rCNV]{allele.info}
#' @param test character. type of test to be used for significance. See detials
#' @param intersection logical, whether to use the intersection of the methods specified in test (if more than one)
#' @param method character. method for testing excess of heterozygotes. Fisher exact test ("fisher") or Chi squre test ("chi.sq")
#' @param plot logical. whether to plot the detected singlets and duplicates on allele ratio vs. proportion of heterozygotes plot.
#' @param verbose logical. show progress
#' @param ... additional parameters to be passed on to \link[graphics]{plot}
#'
#' @importFrom colorspace terrain_hcl
#'
#' @return Returns a data frame of snps/alleles with their duplication status
#'
#' @details Duplicates are detected with both excess of heterozygosity according to HWE and deviant SNPs where deviants are detected using the following methods \cr
#' 1. Z-score test \deqn{Z =  \frac{ \frac{N}{2} -  N_{A}  }{    \sigma _{x}   }} \cr
#' 2. chi-square test (see references for more details on the method) \cr
#' Users can pick among Z-score for heterozygotes (z.het, chi.het), all allele combinations (z.all, chi.all) and the assumption of no probe bias p=0.5 (z.05, chi.05)
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{data(alleleINF)
#' DD<-dupGet(alleleINF)}
#'
#' @export
dupGet<-function(data,test=c("z.het","z.05","z.all","chi.het","chi.05","chi.all"),intersection=FALSE,method=c("fisher","chi.sq"),plot=TRUE,verbose=TRUE,...){
  #data check
  data<-as.data.frame(data)
  if(!any(colnames(data)=="propHet")){
    stop("please provide the data with the output of allele.info()")
  } else {
    ht<-sig.hets(data,plot=F,verbose=verbose)
    test<-match.arg(test,several.ok = TRUE)
    if(length(test)==6){
      pp<-data[,c("z.all","chi.all")]
    } else {
      pp<-data.frame(data[,test])
    }
    if(intersection){
      df<-matrix(NA,nrow = nrow(pp),ncol = ncol(pp))
      for(i in 1:ncol(pp)){
        df[which(pp[,i]<0.05/nrow(pp)),i]<-1
      }
      pp$dup.stat<-"singlet"
      pp$dup.stat[which(rowSums(df)==(ncol(pp)-1))]<-"duplicated"
    } else {
      pp$dup.stat<-"singlet"
      for(i in 1:ncol(pp)){
        pp$dup.stat[pp[,i]<0.05/nrow(pp)]<-"duplicated"
      }
    }
    pp$dup.stat[which(ht$dup.stats=="duplicated")]<-"duplicated"
    pp<-data.frame(data[,1:10],dup.stat=pp$dup.stat)

    if(plot){
      l<-list(...)
      if(is.null(l$cex)) l$cex=0.2
      if(is.null(l$pch)) l$pch=19
      if(is.null(l$xlim)) l$xlim=c(0,1)
      if(is.null(l$ylim)) l$ylim=c(0,1)
      if(is.null(l$alpha)) l$alpha=0.3
      if(is.null(l$col)) l$col<-makeTransparent(c("tomato","#2297E6FF"))#colorspace::terrain_hcl(12,c=c(65,0),l=c(45,90),power=c(1/2,1.5))[2]
      Color <- rep(l$col[2],nrow(pp))
      Color[pp$dup.stat=="duplicated"]<- l$col[1]
      plot(pp$medRatio~pp$propHet, pch=l$pch, cex=l$cex,col=Color,xlim=l$xlim,ylim=l$ylim,frame=F,
           ylab="Allele Median Ratio",xlab="Proportion of Heterozygotes")
      legend("bottomright", c("duplicates","singlets"), col = makeTransparent(l$col,alpha=1), pch=l$pch,
             cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
    }
  }
  return(pp)
}







