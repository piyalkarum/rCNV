## VCF statistics ###

#1. get heterozygosity per individual
## P(Homo) = F + (1-F)P(Homo by chance)
## P(Homo by chance) = p^2+q^2 for a biallelic locus.
## For an individual with N genotyped loci, we
##   1. count the total observed number of loci which are homozygous (O),
##   2. calculate the total expected number of loci homozygous by chance (E)
## Then, using the method of moments, we have
##    O = NF + (1-F)E
## Which rearranges to give
##    F = (O-E)/(N-E)
het.sity1 <- function(ind){  ### calculate expected number of hets
  tab <- table(factor(ind, levels=c("0/0", "1/1", "1/0", "0/1", "./.", ".")))
  N <- tab["0/0"] + tab["1/1"] + tab["0/1"] + tab["1/0"]
  p <- (2 * tab["0/0"] + tab["0/1"] + tab["1/0"])/ (2*N)
  q <- 1 - p
  E <-(p^2+q^2)
  return(E)
}

het.sity2 <- function(ind,eh){
  tab <- table(factor(ind, levels=c("0/0", "1/1", "1/0", "0/1", "./.", ".")))
  O <- tab["0/0"] + tab["1/1"]
  N <- tab["0/0"] + tab["1/1"] + tab["0/1"] + tab["1/0"]
  E <- sum(eh[which(ind == "0/0" | ind == "1/1" | ind == "1/0" | ind == "0/1" )])
  FF <-(O-E)/(N-E)
  return(c(O,E,N,FF))
}
#2 relatedness (according to Yang et al. 2010 equation no. 6 in the paper)
#1. get alt allele freq for all snps
#2. get genotype at each snp per individual  gt=genotype1+genotype2 [1/0]  >> aa=0,Aa=1, AA=2
#     denominator div = 1.0/(2.0*freq*(1.0-freq))
#3. calculate Ajk pairwise
#     if j=k >> Ajk[ui][ui] += (x[ui]*x[ui] - (1 + 2.0*freq)*x[ui] + 2.0*freq*freq) * div
#     if j!=k >> Ajk[ui][uj] += (x[ui] - 2.0*freq) * (x[uj] - 2.0*freq) * div
##gt2<-function(x,gg,freq){two<-gg[,x[1]];one<-gg[,x[2]]
##vp<-NULL
##if(x[1]==x[2]){
##  for(i in seq_along(one)){
##    vp[i]<-((one[i]*one[i])-((1+(2*freq[i]))*one[i])+(2*freq[i]*freq[i]))/(2*freq[i]*(1-freq[i]))
##  }
##  vp<-na.omit(vp)
##  Ajk<-((sum(vp))/length(vp))+1
##} else {
##  for(i in seq_along(one)){
##    vp[i]<-((one[i]-(2*freq[i]))*(two[i]-(2*freq[i])))/(2*freq[i]*(1-freq[i]))
##  }
##  vp<-na.omit(vp)
##  Ajk<-sum(vp)/length(vp)
##}
##V<-c(colnames(gg)[x[2]],colnames(gg)[x[1]],Ajk)
##return(V)}

gt2 <- function(x,gg,freq){
  two<-gg[,x[1]]
  one<-gg[,x[2]]
  if(x[1]==x[2]){
    vp<-((one*one)-((1+(2*freq))*one)+(2*freq*freq))/(2*freq*(1-freq))
    vp<-na.omit(vp)
    Ajk<-((sum(vp))/length(vp))+1
  } else {
    vp<-((one-(2*freq))*(two-(2*freq)))/(2*freq*(1-freq))
    vp<-na.omit(vp)
    Ajk<-sum(vp)/length(vp)
  }
  V<-c(colnames(gg)[x[2]],colnames(gg)[x[1]],Ajk)
  return(V)
}

#' Determine per sample heterozygosity and inbreeding coefficient
#'
#' This function will calculate the heterozygosity on a per-sample basis from
#' vcf files (snps), and most importantly inbreeding coefficient which is used
#'  to filter out the samples with bad mapping quality.
#' @param vcf an imported vcf file in in a list using
#' \code{readVCF} or a data frame of genotypes generated using
#' \code{hetTgen}
#' @param plot logical. Whether to plot a boxplot of inbreeding coefficients
#' for populations. A list of populations must be provided
#' @param pops character. A list of population names with the same length and
#' order as the number of samples in the vcf
#' @param verbose logical. Show progress
#' @param parallel logical. Parallelize the process
#' @importFrom graphics boxplot
#' @importFrom parallel parApply detectCores parLapply stopCluster makeCluster
#' @return Returns a data frame of expected \dQuote{E(Hom)} and observed
#' \dQuote{O(Hom)} homozygotes with their inbreeding coefficients.
#'
#' @author Piyal Karunarathne, Pascal Milesi, Klaus Schliep
#'
#' @examples
#' \dontrun{vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' pp<-substr(colnames(vcf$vcf)[-c(1:9)],1,8)
#' hzygots<-h.zygosity(vcf,plot=TRUE,pops=pp)}
#'
#' @export
h.zygosity<-function(vcf,plot=FALSE,pops=NA,verbose=TRUE,parallel=FALSE){
  if(inherits(vcf,"list")) {
    vcf<-vcf$vcf
    gtt<-hetTgen(vcf,"GT",verbose=verbose)
  } else {
    if(any(colnames(vcf)=="REF")){gtt<-hetTgen(vcf,"GT",verbose=verbose)}
    else {gtt <-vcf}
  }
  if(parallel){
    numCores<-detectCores()-1
    cl<-makeCluster(numCores)
    clusterExport(cl, c("het.sity1","het.sity2"))
    eh<-unlist(t(parApply(cl=cl,gtt[,-c(1:4)],1,het.sity1)))
    hh<-t(parApply(cl=cl,gtt[,-c(1:4)],2,het.sity2,eh=eh))
    stopCluster(cl)
  } else {
    if(verbose){
      message("assessing per sample homozygosity")
      eh<-unlist(apply_pb(gtt[,-c(1:4)],1,het.sity1))
      hh<-t(apply_pb(gtt[,-c(1:4)],2,het.sity2,eh=eh))
    } else {
      eh<-unlist(t(apply(gtt[,-c(1:4)],1,het.sity1)))
      hh<-t(apply(gtt[,-c(1:4)],2,het.sity2,eh=eh))
    }
  }
  hh<-data.frame(rownames(hh),hh)
  colnames(hh)<-c("ind","O(Hom)","E(Hom)","total","Fis")
  rownames(hh)<-NULL
  if(plot){
    if(is.na(pops[1])){
      warning("Please provide a population list")
    } else {
      hh$pop<-pops
      if(length(unique(hh$pop))<=1){
        warning("individuals must come from at least two populations")
      }
      boxplot(hh$Fis~hh$pop,ann=F)
    }
  }
  return(hh)
}


#' Determine pairwise relatedness
#'
#' Relatedness is determined according to genome-wide relationship assessment
#' of Yang et al. 2010 equation 6, on a per sample basis (with itself and
#' others), using SNPs.
#'
#' @param vcf an imported vcf file in a list using \code{readVCF}
#'  or a data frame of genotypes generated using \code{hetTgen}
#' @param plot logical. Whether to plot relatedness of samples against
#' themselves, among themselves and outliers
#' @param threshold numerical. A value indicating to filter the individuals of
#' relatedness among themselves. Default \code{0.5} (siblings)
#' @param verbose logical. Show progress.
#' @param parallel logical. Parallelize the process
#' @importFrom graphics hist
#' @importFrom parallel parApply detectCores parLapply stopCluster makeCluster
#'
#' @return
#' A data frame of individuals and relatedness score \eqn{A_{jk}}
#'
#' @details
#' According to Yang et al. (2010), out breeding non-related pairs should have a
#' relatedness value of zero while the individual with itself will have a
#' relatedness value of one. Relatedness value of ~0.5 indicates siblings.
#'
#' @author Piyal Karunarathne, Klaus Schliep
#'
#' @references Yang, J., Benyamin, B., McEvoy, B. et al. Common SNPs explain a
#' large proportion of the heritability for human height. Nat Genet 42, 565569
#'  (2010).
#'
#' @examples
#' \dontrun{vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' relate<-relatedness(vcf)}
#'
#' @export
relatedness<-function(vcf,plot=TRUE,threshold=0.5,verbose=TRUE,parallel=FALSE){
  if(inherits(vcf,"data.frame") & any(colnames(vcf)=="INDV1")){
    T2<-vcf
  } else {
    if(!inherits(vcf,"list")) {
      if(any(colnames(vcf)=="REF")){gtt<-hetTgen(vcf,"GT",verbose=verbose)}
      else {gtt <-vcf}
    } else {
      vcf<-vcf$vcf
      gtt<-hetTgen(vcf,"GT",verbose=verbose,parallel=parallel)
    }
    gt<-gtt[,-c(1:4)]
    freq<-apply(gt,1,function(xx){aal<-stringr::str_split(xx,"/",simplify = T)
    return(sum(aal=="1")/(sum(aal=="0")+sum(aal=="1")))})

    gg<-apply(gt,2,function(x){XX<-rep(NA,length(x))
    XX[which(x=="0/0")]<-0
    XX[which(x=="1/1")]<-2
    XX[which(x=="1/0" | x=="0/1")]<-1
    XX})

    comb<-expand.grid(1:ncol(gg),1:ncol(gg))
   if(parallel){
     numCores<-detectCores()-1
     cl<-makeCluster(numCores)
     T2<-parApply(cl=cl,comb,1,gt2,gg=gg,freq=freq)
     stopCluster(cl)
   } else {
     if(verbose){
       message("assessing pairwise relatedness")
       T2<-apply_pb(comb,1,gt2,gg=gg,freq=freq)
     } else {
       T2<-apply(comb,1,gt2,gg=gg,freq=freq)
     }
   }
    T2<-data.frame(t(T2))
    T2[,3]<-as.numeric(T2[,3])
    colnames(T2)<-c("indv1","indv2","relatedness_Ajk")
  }
  if(plot){
    same = T2[T2[,1] == T2[,2], ]
    diff = T2[T2[,1] != T2[,2], ]
    outliers = diff[diff[,3] > threshold, ]

    opars<-par(no.readonly = TRUE)
    on.exit(par(opars))

    par(mfrow=c(3, 1))
    hist(same[,3],
         main="Samples against themselves",
         col="grey",
         breaks=seq(-1000, 1000, by=0.05),
         xlim=c(-1, 2))
    hist(diff[,3],
         main="Samples among themselves",
         col="grey",
         breaks=seq(-100, 100, by=0.05),
         xlim=c(-1, 2))
    hist(outliers[,3],
         main="Outlier samples",
         col="grey",
         breaks=seq(-100, 100, by=0.05),
         xlim=c(-1, 2))
  }
  return(T2)
}

#' Get sequencing quality statistics of raw VCF files
#' (with GatK generated vcf files only)
#'
#' This function will generate a table similar to VariantsToTable option in
#'  GatK from raw vcf files for filtering purposes. The function will also
#'  plot all the parameters (see details & values).
#'
#' @param vcf an imported vcf file in data.frame or matrix format using
#' \code{readVCF}
#' @param plot logical. Whether to plot the (12) parameters
#' @param ... other arguments passed on to \code{plot}
#' (e.g. col,border)
#'
#' @importFrom stats density
#' @importFrom graphics polygon
#'
#' @return Returns a data frame with quality parameters from the INFO. field of
#'  the vcf
#' \itemize{
#'  \item{QUAL: The Phred-scaled probability that a REF/ALT polymorphism exists
#'   at this site given sequencing data}
#'  \item{AC: Allele count}
#'  \item{AF: Allele frequency}
#'  \item{DP: unfiltered depth}
#'  \item{QD: QualByDepth - This is the variant confidence (from the QUAL
#'  field) divided by the unfiltered depth of non-hom-ref samples}
#'  \item{FS: FisherStrand - This is the Phred scaled probability that there is
#'   strand bias at the site}
#'  \item{SOR: StrandOddsRatio - This is another way to estimate strand bias
#'  using a test similar to the symmetric odds ratio test}
#'  \item{MQ: RMSMappingQuality - This is the root mean square mapping quality
#'   over all the reads at the site}
#'  \item{MQRankSum: MappingQualityRankSumTest - This is the u-based
#'  z-approximation from the Rank Sum Test for mapping qualities}
#'  \item{ReadPosRankSum: ReadPosRankSumTest: This is the u-based
#'  z-approximation from the Rank Sum Test for site position within reads}
#' }
#' @details
#' For more details see instructions of GatK
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' statistics<-vcf.stat(vcf,plot=TRUE)}
#'
#' @export
vcf.stat<-function(vcf,plot=TRUE,...){
  vcf<-vcf$vcf
  tbb<-apply(vcf[,8],1,function(xx){
    tmp<-unlist(strsplit(as.character(xx),";"))
    AC<-gsub(".*=","\\1",tmp[grep("^AC\\=",tmp)])
    AF<-gsub(".*=","\\1",tmp[grep("^AF\\=",tmp)])
    DP<-gsub(".*=","\\1",tmp[grep("^DP\\=",tmp)])
    MQ<-gsub(".*=","\\1",tmp[grep("^MQ\\=",tmp)])
    QD<-gsub(".*=","\\1",tmp[grep("^QD\\=",tmp)])
    FS<-gsub(".*=","\\1",tmp[grep("^FS\\=",tmp)])
    SOR<-gsub(".*=","\\1",tmp[grep("^SOR\\=",tmp)])
    HET<-gsub(".*=","\\1",tmp[grep("^ExcessHet\\=",tmp)])
    MQRankSum<-gsub(".*=","\\1",tmp[grep("^MQRankSum\\=",tmp)])
    ReadPosRankSum<-gsub(".*=","\\1",tmp[grep("^ReadPosRankSum\\=",tmp)])
    InbrCo<-gsub(".*=","\\1",tmp[grep("^InbreedingCoeff\\=",tmp)])
    ll<-list(AC,AF,DP,MQ,QD,FS,SOR,HET,MQRankSum,ReadPosRankSum,InbrCo)
    ll<-unlist(lapply(ll,function(x){if(length(x)==0)x<-NA else x<-x}))
    return(ll)
  })
  tbb<-t(tbb)
  colnames(tbb)<-c("AC","AF","DP","MQ","QD","FS","SOR","HET","MQRankSum","ReadPosRankSum","InbrCo")
  tbb<-data.frame(cbind(vcf[,c(1:3,6)],tbb))
  if(plot){
    opars<-par(no.readonly = TRUE)
    on.exit(par(opars))
    #par(mfrow=c(4,3))
    pl<-list(...)
    if(is.null(pl$col)) pl$col<-"lightblue"
    if(is.null(pl$border)) pl$border<-"firebrick"
    pp<-sapply(colnames(tbb[,-c(1:3)]),function(x,pl){plot((dd<-density(as.numeric(tbb[,x]),na.rm = T)),main = x,typ="n",xlab=NA)
      polygon(dd,col=pl$col,border=pl$border)},pl=pl)
    #par(mfrow=c(1,1))
  }
  return(tbb)
}






