
#' Export VCF files
#'
#' A function to export tables/matrices in VCF format to VCF files
#'
#' @param out.vcf a matrix or data frame in vcf file format to be exported
#' @param out.path a character string of output path for the vcf file;
#' should end in the name as the vcf file and .vcf. See examples
#' @param compress logical. whether to compress the output file. If \code{TRUE},
#'  the file will be \code{.gz} compressed
#'
#' @importFrom utils packageVersion
#' @importFrom R.utils gzip
#'
#' @return export a vcf file
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' \dontrun{vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
#' vcf <- readVCF(vcf.file.path)
#' exportVCF(vcf,"../exVcf.vcf")}
#'
#' @export
exportVCF<-function(out.vcf, out.path, compress=TRUE){
  if(inherits(out.vcf,"list")){out.vcf<-out.vcf$vcf}
  out.vcf<-as.matrix(out.vcf)
  fcon<-out.path
  header<-colnames(out.vcf)
  hr<-header[1:9]
  if(!setequal(c("#CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT"),hr)){
    stop("input vcf header incorrect\n
         must have #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT, followed by samples")
  } else {
    cat(paste0('##fileformat=VCFv4.0\n',
               '##fileDate=',gsub("-","",Sys.Date()),'\n',
               '##source=','"R rCNV',packageVersion('rCNV'),'"\n',
               '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n',
               '##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">\n',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
               '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n',
               '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Depth">\n',
               '##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype Likelihood">\n'),file=fcon,sep="\t")
    cat(colnames(out.vcf),file=fcon,append=T,sep="\t")
    cat("\n",file=fcon,append=T)
    for(i in 1:nrow(out.vcf)){
      cat(out.vcf[i,],file=fcon,append=T,sep="\t")
      cat("\n",file=fcon,append=T)
    }
    if(compress){
      message("compressing file")
      gzip(fcon,ext="gz",FUN=gzfile)
    }
  }
}










