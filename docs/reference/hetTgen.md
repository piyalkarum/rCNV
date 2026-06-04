# Generate allele depth or genotype table

hetTgen extracts the read depth and coverage values for each snp for all
the individuals from a vcf file generated from readVCF (or GatK
VariantsToTable: see details)

## Usage

``` r
hetTgen(
  vcf,
  info.type = c("AD", "AD-tot", "GT", "GT-012", "GT-AB", "DP"),
  verbose = TRUE,
  parallel = FALSE
)
```

## Arguments

- vcf:

  an imported vcf file in a list using `readVCF`

- info.type:

  character. `AD`: allele depth value, `AD-tot`:total allele depth,
  `DP`=unfiltered depth (sum), `GT`: genotype, `GT-012`:genotype in 012
  format, `GT-AB`:genotype in AB format. Default `AD`, See details.

- verbose:

  logical. whether to show the progress of the analysis

- parallel:

  logical. whether to parallelize the process

## Value

Returns a data frame of allele depth, genotype of SNPs for all the
individuals extracted from a VCF file

## Details

If you generate the depth values for allele by sample using GatK
VariantsToTable option, use only -F CHROM -F POS -GF AD flags to
generate the table. Or keep only the CHROM, POS, ID, ALT, and individual
AD columns. For info.type `GT` option is provided to extract the
genotypes of individuals by snp.

## Author

Piyal Karunarathne, Klaus Schliep

## Examples

``` r
if (FALSE) vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path=vcf.file.path)
het.table<-hetTgen(vcf) # \dontrun{}
#> generating table
#>   |                                                          |                                                  |   0%  |                                                          |=====                                             |  10%  |                                                          |==========                                        |  20%  |                                                          |===============                                   |  30%  |                                                          |====================                              |  40%  |                                                          |=========================                         |  50%  |                                                          |==============================                    |  60%  |                                                          |===================================               |  70%  |                                                          |========================================          |  80%  |                                                          |=============================================     |  90%  |                                                          |==================================================| 100%
```
