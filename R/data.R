#' Allele info example data
#'
#' Semi-randomly generated data from the function dup.snp.info.
#' Data contains depth and proportion values of 2857 snps
#'
#' @docType data
#'
#' @usage data(alleleINF)
#'
#' @keywords datasets
#'
#' @references
#' \itemize{
#'  \item{Larson, W. A., Seeb, L. W., Everett, M. V., Waples, R. K., Templin, W. D., & Seeb, J. E. (2014). Genotyping by sequencing resolves #' shallow population structure to inform conservation of Chinook salmon (Oncorhynchus tshawytscha). Evolutionary Applications, 7(3), 355-369}
#'  \item{McKinney, G. J., Waples, R. K., Seeb, L. W., & Seeb, J. E. (2017). Paralogs are revealed by proportion of heterozygotes and deviations in read ratios in genotyping‐by‐sequencing data from natural populations. Molecular Ecology Resources, 17(4), 656-669.}
#' }
#'
#' @source \href{https://zenodo.org/record/5025423#.Yj2XKRDMLyQ}{Chinook Salmon sequence reads McKinney et al. 2017}
#'
#'
#' @examples
#' data(alleleINF)
#' with(alleleINF,plot(medRatio~propHet))
'alleleINF'


#' Allele Depth (AD) example data
#'
#' Example SNPs data of Chinook Salmon from Larson et al. et al. 2014.
#' The data contains only a partial snps data set of RadSeq data after filtering.
#'
#' @docType data
#' @usage data(ADtable)
#' @references
#' \itemize{
#'  \item{Larson, W. A., Seeb, L. W., Everett, M. V., Waples, R. K., Templin, W. D., & Seeb, J. E. (2014). Genotyping by sequencing resolves shallow population structure to inform conservation of Chinook salmon (Oncorhynchus tshawytscha). Evolutionary Applications, 7(3), 355-369.}
#'  \item{McKinney, G. J., Waples, R. K., Seeb, L. W., & Seeb, J. E. (2017). Paralogs are revealed by proportion of heterozygotes and deviations in read ratios in genotyping‐by‐sequencing data from natural populations. Molecular Ecology Resources, 17(4), 656-669.}
#' }
#'
'ADtable'


#' Normalized allele depth example data
#'
#' Normalized example SNPs data of Chinook Salmon from Larson et al. 2014\cr
#' The data has been normalized with TMM
#'
#' @docType data
#' @usage  data(ADnorm)
#'
#' @references Larson, W. A., Seeb, L. W., Everett, M. V., Waples, R. K., Templin, W. D., & Seeb, J. E. (2014). Genotyping by sequencing resolves shallow population structure to inform conservation of Chinook salmon (Oncorhynchus tshawytscha). Evolutionary Applications, 7(3), 355-369
#'
'ADnorm'


