% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bias.det.R
\name{allele.info.WGS}
\alias{allele.info.WGS}
\title{ClrCNV: multicopy detection for WGS}
\usage{
allele.info.WGS(
  ad,
  gt,
  fis = NULL,
  vcf = NULL,
  parallel = FALSE,
  numCores = NULL,
  ...
)
}
\arguments{
\item{ad}{allele depth table generated with \code{hetTgen}}

\item{gt}{genotype table generated with \code{hetTgen}}

\item{fis}{global inbreeding coefficient calculated with \code{h.zygosity}}

\item{vcf}{if fis is not provided, vcf file imported using \code{readVCF}}

\item{parallel}{logical. to parallelize over multiple cores}

\item{numCores}{numeric. if parallel TRUE, number of cores to use; if NULL, use all cores}

\item{...}{other arguments passed to makeCluster}
}
\value{
Returns a data frame with likelihood ratios per SNP and duplication status
}
\description{
This function is catered for detecting multicopy regions from whole genome sequences (WGS)
using Likelihood Ratios
}
\examples{
\dontrun{
vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path=vcf.file.path)
ad<-hetTgen(vcf,"AD")
gt<-hetTgen(vcf,"GT")
hz<-h.zygosity(vcf,verbose=FALSE)
Fis<-mean(hz$Fis,na.rm = TRUE)
AI_WGS<-allele.info.WGS(ad,gt,Fis=Fis)}

}
\author{
Qiujie Zhou, Pascal Milesi, Piyal Karunarathne
}
