# ClrCNV: multicopy detection for WGS

This function is catered for detecting multicopy regions from whole
genome sequences (WGS) using Likelihood Ratios

## Usage

``` r
allele.info.WGS(
  ad,
  gt,
  fis = NULL,
  vcf = NULL,
  parallel = FALSE,
  numCores = NULL,
  ...
)
```

## Arguments

- ad:

  allele depth table generated with `hetTgen`

- gt:

  genotype table generated with `hetTgen`

- fis:

  global inbreeding coefficient calculated with `h.zygosity`

- vcf:

  if fis is not provided, vcf file imported using `readVCF`

- parallel:

  logical. to parallelize over multiple cores

- numCores:

  numeric. if parallel TRUE, number of cores to use; if NULL, use all
  cores

- ...:

  other arguments passed to makeCluster

## Value

Returns a data frame with likelihood ratios per SNP and duplication
status

## Author

Qiujie Zhou, Pascal Milesi, Piyal Karunarathne

## Examples

``` r
if (FALSE) { # \dontrun{
vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path=vcf.file.path)
ad<-hetTgen(vcf,"AD")
gt<-hetTgen(vcf,"GT")
hz<-h.zygosity(vcf,verbose=FALSE)
Fis<-mean(hz$Fis,na.rm = TRUE)
AI_WGS<-allele.info.WGS(ad,gt,Fis=Fis)} # }
```
