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

het.sity<-function(ind){
  O<-length(which(ind=="0/0"))+length(which(ind=="1/1"))
  N<-length(which(ind=="0/0"))+length(which(ind=="1/1"))+length(which(ind=="0/1"))+length(which(ind=="1/0"))
  alle<-unname(unlist(lapply(ind,strsplit,split="/")))
  p<-length(which(alle=="0"))/(length(which(alle=="0"))+length(which(alle=="1")))
  q<-length(which(alle=="1"))/(length(which(alle=="0"))+length(which(alle=="1")))
  E<-(p^2+q^2)*N
  FF<-(O-E)/(N-E)
  return(c(O,E,N,FF))
}

#2 relatedness (according to Yang et al. 2010 equation no. 6 in the paper)
#1. get alt allele freq for all snps
#2. get genotype at each snp per individual  gt=genotype1+genotype2 [1/0]  >> aa=0,Aa=1, AA=2
#     denominator div = 1.0/(2.0*freq*(1.0-freq))
#3. calculate Ajk pairwise
#     if j=k >> Ajk[ui][ui] += (x[ui]*x[ui] - (1 + 2.0*freq)*x[ui] + 2.0*freq*freq) * div
#     if j!=k >> Ajk[ui][uj] += (x[ui] - 2.0*freq) * (x[uj] - 2.0*freq) * div
gt2<-function(x,gg,freq){two<-gg[,x[1]];one<-gg[,x[2]]
vp<-NULL
if(x[1]==x[2]){
  for(i in seq_along(one)){
    vp[i]<-((one[i]*one[i])-((1+(2*freq[i]))*one[i])+(2*freq[i]*freq[i]))/(2*freq[i]*(1-freq[i]))
  }
  vp<-na.omit(vp)
  Ajk<-((sum(vp))/length(vp))+1
} else {
  for(i in seq_along(one)){
    vp[i]<-((one[i]-(2*freq[i]))*(two[i]-(2*freq[i])))/(2*freq[i]*(1-freq[i]))
  }
  vp<-na.omit(vp)
  Ajk<-sum(vp)/length(vp)
}
V<-c(colnames(gg)[x[2]],colnames(gg)[x[1]],Ajk)
return(V)}



#' Determine per sample heterozygosity and inbreeding coefficient
#'
#' This function will calculate the heterozygosity on a per-sample basis from vcf files (snps), and most importantly inbreeding coefficient which is used to filter out the samples with bad mapping quality. See details.
#' @param vcf an imported vcf file in data.frame or matrix format using "readVCF"
#' @param plot logical. Whether to plot a boxplot of inbreeding coefficients for populations. A list of populations must be provided
#' @param pops character. A list of population names
#' @importFrom graphics boxplot
#' @return Returns a data frame of expected "E(Hom)", observed "O(Hom)" homozygotes with their inbreeding coefficients.
#'
#' @author Piyal Karunarathne Pascal Milesi
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' hzygots<-h.zygosity(vcf,plot=TRUE,pops=c("FR_PA","RU_PA"))
#'
#' @export
h.zygosity<-function(vcf,plot=FALSE,pops=NA){
  if(inherits(vcf,"list")) {
    vcf<-vcf$vcf
    message("reading genotypes")
    gtt<-hetTgen(vcf,"GT",verbose=TRUE)
  } else {
    gtt <-vcf
  }
  message("assessing per sample homozygosity")
  hh<-t(apply_pb(gtt[,-c(1:3)],2,het.sity))
  hh<-data.frame(rownames(hh),hh)
  colnames(hh)<-c("ind","O(Hom)","E(Hom)","total","Fis")
  rownames(hh)<-NULL
  if(plot){
    if(is.na(pops[1])){
      warning("Please provide a population list")
    } else {
      hh$pop<-"ungrouped"
      for(i in seq_along(pops)){
        hh$pop[grep(pops[i],hh$ind)]<-pops[i]
      }
      if(length(unique(hh$pop))<=1){
        warning("population names must be a part of individual names")
      }
      boxplot(hh$Fis~hh$pop,ann=F)
    }
  }
  return(hh)
}


#' Determine pairwise relatedness
#'
#' Relatedness is determined according to genome-wide relationship assessment of Yang et al. 2010 (doi:10.1038/ng.608) equation 6, on a per sample basis (with itself and others), using SNPs.
#'
#' @param vcf an imported vcf file in data.frame or matrix format using "readVCF"
#' @param plot logical. Whether to plot relatedness of samples against themselves, among themselves and outliers
#' @param threshold numerical. A value indicating to filter the individuals of relatedness among themselves. Default=0.5 (siblings)
#' @importFrom graphics hist
#'
#' @return
#' A data frame of individuals and relatedness score Ajk
#'
#' @details
#' According to Yang et al. (2010), outbreeding non-related pairs should have a relatedness value of 0 while the individual with itself will have a relatedness value of 1. Relatedness value of ~0.5 indicates siblings.
#'
#' @author Piyal Karunarathne
#'
#' @references Yang, J., Benyamin, B., McEvoy, B. et al. Common SNPs explain a large proportion of the heritability for human height. Nat Genet 42, 565–569 (2010).
#' https://doi.org/10.1038/ng.608
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' relate<-relatedness(vcf)
#'
#' @export
relatedness<-function(vcf,plot=TRUE,threshold=0.5){
  if(!inherits(vcf,"list")) {
    gtt <-vcf
  } else {
    vcf<-vcf$vcf
    message("reading genotypes")
    gtt<-hetTgen(vcf,"GT",verbose=TRUE)
  }
  gt<-gtt[,-c(1:3)]
  freq<-apply(gt,1,function(xx){aal<-unlist(strsplit(as.character(xx),"/"))
  return(length(which(aal=="1"))/(length(which(aal=="0"))+length(which(aal=="1"))))})

  gg<-apply(gt,2,function(x){XX<-rep(NA,length(x))
  XX[which(x=="0/0")]<-0
  XX[which(x=="1/1")]<-2
  XX[which(x=="1/0" | x=="0/1")]<-1
  XX})

  comb<-expand.grid(1:ncol(gg),1:ncol(gg))
  message("assessing pairwise relatedness")
  T2<-apply_pb(comb,1,gt2,gg=gg,freq=freq)
  T2<-data.frame(t(T2))
  T2[,3]<-as.numeric(T2[,3])
  colnames(T2)<-c("indv1","indv2","relatedness_Ajk")
  if(plot){
    same = T2[T2[,1] == T2[,2], ]
    diff = T2[T2[,1] != T2[,2], ]
    outliers = diff[diff[,3] > threshold, ]

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

#' Get sequencing quality statistics of raw VCF files (with GatK generated vcf files only)
#'
#' This function will generate a table similar to VariantsToTable option in GatK from raw vcf files for filtering purposes. The fucntion will aslo plot all the parameters (see details).
#'
#' @param vcf an imported vcf file in data.frame or matrix format using "readVCF"
#' @param plot logical. Whether to plot the (12) parameters
#' @param ... other arguments passed on to plot (e.g. col,border)
#'
#' @importFrom stats density
#' @importFrom graphics polygon
#'
#' @return returns a data frame with quality parameters from the info. field of the vcf
#' QUAL   The Phred-scaled probability that a REF/ALT polymorphism exists at this site given sequencing data
#' AC     Allele count
#' AF     Allele frequencey
#' DP     unfiltered depth
#' MQ     .....
#' QD     QualByDepth - This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-hom-ref samples
#' FS     FisherStrand - This is the Phred-scaled probability that there is strand bias at the site.
#' SOR    StrandOddsRatio - This is another way to estimate strand bias using a test similar to the symmetric odds ratio test
#' MQ     RMSMappingQuality - This is the root mean square mapping quality over all the reads at the site
#' MQRankSum  MappingQualityRankSumTest - This is the u-based z-approximation from the Rank Sum Test for mapping qualities
#' ReadPosRankSum ReadPosRankSumTest - This is the u-based z-approximation from the Rank Sum Test for site position within reads
#' ...
#'
#' @details
#' For more details see https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path=vcf.file.path)
#' statistics<-vcf.stat(vcf,plot=TRUE)
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
    par(mfrow=c(4,3))
    pl<-list(...)
    if(is.null(pl$col)) pl$col<-"lightblue"
    if(is.null(pl$border)) pl$border<-"firebrick"
    pp<-sapply(colnames(tbb[,-c(1:3)]),function(x,pl){plot((dd<-density(as.numeric(tbb[,x]),na.rm = T)),main = x,typ="n",xlab=NA)
      polygon(dd,col=pl$col,border=pl$border)},pl=pl)
  }
  return(tbb)
}






