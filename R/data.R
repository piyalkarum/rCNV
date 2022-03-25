#' Allele info example data
#'
#' Semi-randomly generated data from the function dup.snp.info
#' Data contains depth and proportion values of 19800 snps
#'
#' @docType data
#'
#' @usage data(alleleINF)
#'
#' @keywords datasets
#'
#' @references Karunarathne et al. (2021)
#'
#' @source \href{https://zenodo.org/record/5025423#.Yj2XKRDMLyQ}{Chinook Salmon sequence reads McKinney et al. 2017}
#'
#' @examples
#' data(alleleINF)
#' with(alleleINF,plot(medRatio~propHet))
'alleleINF'


#' Allele Depth (AD) example data
#'
#' example SNPs data of Picea abies from Chen et al. 2019
#' The data contains only a partial snps data set of exome-capture data after filtering
#'
#' @docType data
#' @usage data(ADtable)
#' @references McKinney et al. 2017
#'
'ADtable'


#' Normalized allele depth example data
#'
#' normalized example SNPs data of Picea abies from Chen et al. 2019
#' The data has been normalized with TMM
#'
#' @docType data
#' @usage  data(ADnorm)
#'
'ADnorm'


