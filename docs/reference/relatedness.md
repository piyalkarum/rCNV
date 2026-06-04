# Determine pairwise relatedness

Relatedness is determined according to genome-wide relationship
assessment of Yang et al. 2010 equation 6, on a per sample basis (with
itself and others), using SNPs.

## Usage

``` r
relatedness(
  vcf,
  plot = TRUE,
  threshold = 0.5,
  verbose = TRUE,
  parallel = FALSE
)
```

## Arguments

- vcf:

  an imported vcf file in a list using `readVCF` or a data frame of
  genotypes generated using `hetTgen`

- plot:

  logical. Whether to plot relatedness of samples against themselves,
  among themselves and outliers

- threshold:

  numerical. A value indicating to filter the individuals of relatedness
  among themselves. Default `0.5` (siblings)

- verbose:

  logical. Show progress.

- parallel:

  logical. Parallelize the process

## Value

A data frame of individuals and relatedness score \\A\_{jk}\\

## Details

According to Yang et al. (2010), out breeding non-related pairs should
have a relatedness value of zero while the individual with itself will
have a relatedness value of one. Relatedness value of ~0.5 indicates
siblings.

## References

Yang, J., Benyamin, B., McEvoy, B. et al. Common SNPs explain a large
proportion of the heritability for human height. Nat Genet 42, 565569
(2010).

## Author

Piyal Karunarathne, Klaus Schliep

## Examples

``` r
if (FALSE) vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path=vcf.file.path)
relate<-relatedness(vcf) # \dontrun{}
#> generating table
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
#> assessing pairwise relatedness
#>   |                                                          |                                                  |   0%  |                                                          |=                                                 |   3%  |                                                          |===                                               |   6%  |                                                          |====                                              |   8%  |                                                          |======                                            |  11%  |                                                          |=======                                           |  14%  |                                                          |========                                          |  17%  |                                                          |==========                                        |  19%  |                                                          |===========                                       |  22%  |                                                          |============                                      |  25%  |                                                          |==============                                    |  28%  |                                                          |===============                                   |  31%  |                                                          |=================                                 |  33%  |                                                          |==================                                |  36%  |                                                          |===================                               |  39%  |                                                          |=====================                             |  42%  |                                                          |======================                            |  44%  |                                                          |========================                          |  47%  |                                                          |=========================                         |  50%  |                                                          |==========================                        |  53%  |                                                          |============================                      |  56%  |                                                          |=============================                     |  58%  |                                                          |===============================                   |  61%  |                                                          |================================                  |  64%  |                                                          |=================================                 |  67%  |                                                          |===================================               |  69%  |                                                          |====================================              |  72%  |                                                          |======================================            |  75%  |                                                          |=======================================           |  78%  |                                                          |========================================          |  81%  |                                                          |==========================================        |  83%  |                                                          |===========================================       |  86%  |                                                          |============================================      |  89%  |                                                          |==============================================    |  92%  |                                                          |===============================================   |  94%  |                                                          |================================================= |  97%  |                                                          |==================================================| 100%

```
