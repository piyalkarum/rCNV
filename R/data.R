#' Duplication info example data
#'
#' Semi-randomly generated data from the function dup.snp.info
#' Data contains depth and proportion values of 19800 snps
#'
#' @docType data
#'
#' @usage data(dup.info)
#'
#' @keywords datasets
#'
#' @references Karunarathne et al. (2021)
#'
#' @source \href{url}{dataset name}
#'
#' @examples
#' data(dup.info)
#' props<-dup.info[,c("PropHomRare","PropHet")]
#' \donttest{plot(props, pch=19,cex=0.2)}
'dup.info'


#' Heterozygosity example data
#'
#' example SNPs data of Picea abies from Chen et al. 2019
#' The data contains only a partial snps data set of exome-capture data after filtering
#'
#' @docType data
#' @usage data(hets)
#' @references Chen et al. 2019 (add the correct reference)
#'
'hets'
