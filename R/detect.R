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

#1 pvals for sig.hets when AD table is provided
get.eHpvals<-function(df,method=c("fisher","chi.sq")){
  snp1<-df[-c(1:4)]
  y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
  y[,1]<-as.integer(y[,1]);y[,2]<-as.integer(y[,2])
  rr1<-y[,2]/rowSums(y,na.rm = T)
  snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==T),]
  NHet<-nrow(snp1het)
  if(NHet>3){
    propHet<-NHet/length(na.omit(rr1))
    medRatio<-median(proportions(as.matrix(snp1het),margin = 1)[,2],na.rm = T)
    homalt<-sum(rr1==1,na.rm=T)
    homref<-sum(rr1==0,na.rm=T)
    Nsamp<-NHet+homalt+homref
    propHomAlt<-homalt/Nsamp
    rs<-c(homref,NHet,homalt,Nsamp)
    n<-Nsamp
    p<-(rs[1]+rs[2]/2)/n
    ob<-c(homref,NHet,homalt)
    eX<-unlist(c((p^2) * n,2*p*(1-p) * n,((1-p)^2) * n))
    delta <- rs[2]-(2*p*(1-p) * n)
    stat<-match.arg(method)
    pval<-switch(stat,fisher=suppressWarnings(fisher.test(cbind(ob,eX),workspace = 2e8))$p.value,chi.sq=suppressWarnings(chisq.test(cbind(ob,eX)))$p.value)
    ll<-c(medRatio,propHomAlt,propHet,eX/n,pval,delta)
  } else {
    ll<-NA
  }
  return(ll)
}

#' Identify significantly different heterozygotes from SNPs data
#'
#' This function will recognize the SNPs with a proportion of heterozygotes
#'  significantly higher than expected under HWE and plot deviant snps based
#'  only on the excess of heterozygotes.
#'
#' @param a.info allele info table generated from filtered vcfs using the
#' function \code{allele.info} or allele depth table generated from \code{hetTgen}
#' @param method character. Method for testing significance. Fisher exact test
#'  (\code{fisher}) or Chi-square test (\code{chi.sq})
#' @param plot logical. Whether to plot the identified duplicated snps with
#' the expected values
#' @param verbose logical, if TRUE, the progress is shown
#' @param ... other arguments passed to \code{plot}
#'
#' @return A matrix of expected heterozygote proportions from the observed
#' data with p-value indicating significance of deviation.
#'
#' @author Piyal Karunarathne, Pascal Milesi
#'
#' @examples
#' \dontrun{data(alleleINF)
#' AI <- alleleINF
#' duplicates<-sig.hets(AI,plot=TRUE)}
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats fisher.test median quantile rbinom sd smooth.spline
#' chisq.test
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom colorspace rainbow_hcl
#' @export
sig.hets<-function(a.info,method=c("fisher","chi.sq"),plot=TRUE,verbose=TRUE,...){
  if(!any(colnames(a.info)=="NHomRef")){
    if(verbose){message("assessing excess of heterozygotes")
      df<-apply_pb(a.info,1,get.eHpvals,method=method)
    } else {df<-apply(a.info,1,get.eHpvals,method=method)}
    df<-data.frame(do.call(rbind,df))
    colnames(df)<-c("medRatio","propHomAlt","propHet","p2","het","q2","eH.pval","eH.delta")
    #df<-na.omit(df)
  } else {
    d<-a.info[,c("NHomRef","NHet","NHomAlt","Nsamp")]
    colnames(d)<-c("h1","het","h2","truNsample")
    method<-match.arg(method)
    if(verbose){
      message("assessing excess of heterozygotes")
      df<-data.frame(t(apply_pb(d,1,ex.prop,method=method)))
    } else {
      df<-data.frame(t(apply(d,1,ex.prop,method=method)))
    }
    colnames(df)<-c("p2","het","q2","eH.pval","eH.delta")
    df$propHomAlt <- a.info$propHomAlt
    df$propHet<-a.info$propHet
    df$medRatio<-a.info$medRatio
  }

  df$dup.stat<-"non-deviant";df$dup.stat[which(df$eH.pval < 0.05/nrow(df) & df$eH.delta > 0 )]<-"deviant"
  d<-na.omit(df[,c("p2","het","q2")])
  df<-na.omit(data.frame(cbind(a.info[,c(1:4)],df[,c("medRatio","propHomAlt","propHet","eH.pval","eH.delta","dup.stat")]),row.names = NULL))

  if(plot){
    l<-list(...)
    if(is.null(l$cex)) l$cex=0.2
    if(is.null(l$pch)) l$pch=19
    if(is.null(l$xlim)) l$xlim=c(0,1)
    if(is.null(l$ylim)) l$ylim=c(0,1)
    if(is.null(l$col)) cols<-makeTransparent(rainbow_hcl(2),alpha=0.3) else cols<-makeTransparent(l$col,alpha=0.3)

    d$Color <- cols[1]
    d$Color [which(df$dup.stat=="deviant")]<- cols[2]#& df$delta > 0
    plot(df$propHet~df$propHomAlt, pch=l$pch, cex=l$cex,col=d$Color,xlim=l$xlim,ylim=l$ylim,
         xlab="Proportion of Alternate Homozygotes",ylab="Proportion of Heterozygotes")
    lines((smm<-smooth.spline(d$het~d$q2)),col="blue")
    legend("bottomright", c("non-deviants","deviants","expected"), col = c(cols,"blue"), lty = c(0, 0, 1), lwd = c(0, 0, 1),pch = c(l$pch, l$pch, NA),
           cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
  }

  return(df)
}


#' Plot duplicates
#'
#' The function plots detected duplicates from functions \code{sig.snps}, and
#' \code{dupGet}
#'
#' @param ds a data frame of detected duplicates
#' @param \dots other graphical parameters to be passed to the function
#' \code{plot}
#'
#' @importFrom colorspace rainbow_hcl
#'
#' @return Returns no value, only plots proportion of heterozygotes vs allele
#' median ratio seperated by duplication status
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
  st<-sort(unique(ds$dup.stat))
  ds$Color [ds$dup.stat==st[1]]<- l$col[1]
  plot(ds$medRatio~ds$propHet, pch=l$pch, cex=l$cex,col=ds$Color,xlim=l$xlim,ylim=l$ylim,frame=F,
       ylab="Allele Median Ratio",xlab="Proportion of Heterozygotes")
  legend("bottomright", st, col = makeTransparent(l$col,alpha=1), pch=l$pch,
         cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")

}


#' Detect duplicates from SNPs
#'
#' Detect duplicated snps using excess of heterozygotes
#' (alleles that do not follow HWE) and snp deviates
#' (alleles that do not follow a normal or chi-square distribution).
#'  See details.
#'
#' @param data data frame of the output of \code{allele.info}
#' @param test character. type of test to be used for significance. See details
#' @param intersection logical, whether to use the intersection of the methods
#'  specified in \code{test} (if more than one)
#' @param method character. method for testing excess of heterozygotes.
#'  Fisher exact test (\code{fisher}) or Chi-square test (\code{chi.sq})
#' @param plot logical. whether to plot the detected singlets and duplicates
#'  on allele ratio vs. proportion of heterozygotes plot.
#' @param verbose logical. show progress
#' @param \dots additional parameters passed on to \code{plot}
#'
#' @importFrom colorspace terrain_hcl
#'
#' @return Returns a data frame of snps/alleles with their duplication status
#'
#' @details SNP deviants are detected with both excess of heterozygosity
#' according to HWE and deviant SNPs where depth values fall outside of the
#' normal distribution are detected using the
#'  following methods:
#' 1. Z-score test \eqn{Z = \frac{\frac{N}{2} -  N_{A}}{\sigma _{x}}}
#' 2. chi-square test (see references for more details on the method)
#'
#' Users can pick among Z-score for heterozygotes (\code{z.het, chi.het}),
#' all allele combinations (\code{z.all, chi.all}) and the assumption of no
#' probe bias p=0.5 (\code{z.05, chi.05})
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
    if(!(any(colnames(data)=="eH.pval"))) {
      ht<-sig.hets(data,plot=F,verbose=verbose)
    } else {
      ht <-data[,c("eH.pval","eH.delta")]
      ht$dup.stat<-"non-deviant"
      ht$dup.stat[which(ht$eH.pval < 0.05/nrow(ht) & ht$eH.delta > 0 )]<-"deviant"
    }

    test<-match.arg(test,several.ok = TRUE)
    if(length(test)==6){
      if(verbose){cat(paste0("categorizing deviant SNPs with \n excess of heterozygotes ","z.all & chi.all"))}
      pp<-data[,c("z.all","chi.all")]
    } else {
      pp<-data.frame(data[,test])
      if(verbose){cat(paste0("categorizing deviant SNPs with \n", " excess of heterozygotes \n ",paste0(unlist(test),collapse = "\n ")))}
    }
    if(intersection){
      df<-matrix(NA,nrow = nrow(pp),ncol = ncol(pp))
      for(i in 1:ncol(pp)){
        df[which(pp[,i]<0.05/nrow(pp)),i]<-1
      }
      pp$dup.stat<-"non-deviant"
      pp$dup.stat[which(rowSums(df)==(ncol(pp)-1))]<-"deviant"
    } else {
      pp$dup.stat<-"non-deviant"
      for(i in 1:ncol(pp)){
        pp$dup.stat[pp[,i]<0.05/nrow(pp)]<-"deviant"
      }
    }
    pp$dup.stat[which(ht$dup.stat=="deviant")]<-"deviant"
    pp<-data.frame(data[,1:10],eH.pval=ht[,"eH.pval"],eH.delta=ht[,"eH.delta"],dup.stat=pp$dup.stat)

    if(plot){
      l<-list(...)
      if(is.null(l$cex)) l$cex=0.2
      if(is.null(l$pch)) l$pch=19
      if(is.null(l$xlim)) l$xlim=c(0,1)
      if(is.null(l$ylim)) l$ylim=c(0,1)
      if(is.null(l$alpha)) l$alpha=0.3
      if(is.null(l$col)) l$col<-makeTransparent(c("tomato","#2297E6FF"))#colorspace::terrain_hcl(12,c=c(65,0),l=c(45,90),power=c(1/2,1.5))[2]
      Color <- rep(l$col[2],nrow(pp))
      Color[pp$dup.stat=="deviant"]<- l$col[1]
      plot(pp$medRatio~pp$propHet, pch=l$pch, cex=l$cex,col=Color,xlim=l$xlim,ylim=l$ylim,frame=F,
           ylab="Allele Median Ratio",xlab="Proportion of Heterozygotes")
      legend("bottomright", c("deviants","non-deviants"), col = makeTransparent(l$col,alpha=1), pch=l$pch,
             cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
    }
  }
  return(pp)
}

#' Find CNVs from deviants
#'
#' Categorize deviant and non-deviant into "singlets" and "duplicates" based on the statistical approaches specified by the user.
#' The intersection of all the stats provided will be used in the categorization. If one would like to use the intersection of at least two stats, this can be specified in the \code{n.ints}
#'
#' @param data A data frame of allele information generated with the function
#' \code{allele.info}
#' @param test vector of characters. Type of test to be used for significance.
#' See details
#' @param filter character. Type of filter to be used for filtering CNVs.
#' default \code{kmeans}. See details.
#' @param ft.threshold confidence interval for filtering \code{default = 0.05}
#' @param plot logical. Plot the detection of duplicates. default \code{TRUE}
#' @param verbose logical. show progress
#' @param ... other arguments to be passed to \code{plot}
#'
#' @return Returns a data frame of SNPs with their detected duplication status
#'
#' @importFrom colorspace terrain_hcl
#' @importFrom stats kmeans
#'
#' @details SNP deviants are detected with both excess of heterozygosity
#' according to HWE and deviant SNPs where depth values fall outside of the
#' normal distribution are detected using the
#'  following methods:
#' 1. Z-score test \eqn{Z = \frac{\frac{N}{2} -  N_{A}}{\sigma _{x}}}
#' 2. chi-square test (see references for more details on the method)
#'
#' Users can pick among Z-score for heterozygotes (\code{z.het, chi.het}),
#' all allele combinations (\code{z.all, chi.all}) and the assumption of no
#' probe bias p=0.5 (\code{z.05, chi.05})
#'
#' \code{filter} will determine whether the \code{intersection} or \code{kmeans}
#' clustering of the provided \code{test}s should be used in filtering CNVs.
#' The intersection uses threshold values for filtering and kmeans use
#' unsupervised clustering. Kmeans clustering is recommended if one is uncertain
#' about the threshold values.
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{data(alleleINF)
#' DD<-cnv(alleleINF)}
#'
#' @export
cnv<-function(data,test=c("z.het","z.05","z.all","chi.het","chi.05","chi.all"),filter=c("intersection","kmeans"),ft.threshold=0.05,plot=TRUE,verbose=TRUE,...){
  #data check
  data<-as.data.frame(data)
  data$z.het.sum<-abs(data$z.het.sum)
  data$z.05.sum<-abs(data$z.05.sum)
  data$z.all.sum<-abs(data$z.all.sum)
  if(!any(colnames(data)=="propHet")){
    stop("please provide the output of allele.info()")
  } else {
    ht<-data[,c("eH.pval","eH.delta")]
    filter<-match.arg(filter,several.ok = TRUE)
    if(length(filter)>1){filter="kmeans"}

    test<-match.arg(test,several.ok = TRUE)
    if(length(test)==6){
      if(verbose){cat(paste0("categorizing putative duplicates with \n excess of heterozygotes ","z.all & chi.all"))}
      pp<-data[,c("z.all","chi.all")]
    } else {
      pp<-data.frame(data[,test])
      if(verbose){cat(paste0("categorizing putative duplicates with \n", "excess of heterozygotes \n",paste0(unlist(test),collapse = "\n")))}
    }
    if(filter=="intersection"){

      df<-as.data.frame(apply(pp,2,function(x){ifelse(x<ft.threshold/length(x),1,0)}))
      df$eH<-0
      df$eH[which(ht$eH.pval < 0.05/nrow(ht) & ht$eH.delta > 0 )]<-1
      pp$dup.stat<-"singlet"
      pp$dup.stat[which(rowSums(df)>=2)]<-"duplicated"
      pp$dup.stat[which(df$eH==1)]<-"duplicated"
      pp<-data.frame(data[,1:10],dup.stat=pp$dup.stat)
    }

    if(filter=="kmeans"){
      test<-paste0(test,".sum")
      candidate<-data.frame(data[,test]/data[,"NHet"],data[,c("eH.delta","cv")])
      candidate<-scale(candidate)
      cls<-kmeans(candidate, centers=2, nstart = 50)
      pp$dup.stat<-rep("singlet",nrow(data))
      mn<-which.min(table(cls$cluster))
      pp$dup.stat[which(cls$cluster==mn)]<-"duplicated"
      pp$dup.stat[which(ht$eH.pval < 0.05/nrow(ht) & ht$eH.delta > 0 )]<-"duplicated"
      pp<-data.frame(data[,1:10],dup.stat=pp$dup.stat)
    }

    if(plot){
      l<-list(...)
      if(is.null(l$cex)) l$cex=0.2
      if(is.null(l$pch)) l$pch=19
      if(is.null(l$xlim)) l$xlim=c(0,1)
      if(is.null(l$ylim)) l$ylim=c(0,1)
      if(is.null(l$alpha)) l$alpha=0.3
      if(is.null(l$col)) l$col<-makeTransparent(c("tomato","#2297E6FF"))
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





