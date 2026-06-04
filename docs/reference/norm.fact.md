# Calculate normalization factor for each sample

This function calculates the normalization factor for each sample using
different methods. See details.

## Usage

``` r
norm.fact(
  df,
  method = c("TMM", "TMMex", "MedR", "QN"),
  logratioTrim = 0.3,
  sumTrim = 0.05,
  Weighting = TRUE,
  Acutoff = -1e+10
)
```

## Arguments

- df:

  a data frame or matrix of allele depth values (total depth per snp per
  sample)

- method:

  character. method to be used (see details). Default `TMM`

- logratioTrim:

  numeric. percentage value (0 - 1) of variation to be trimmed in log
  transformation

- sumTrim:

  numeric. amount of trim to use on the combined absolute levels (“A”
  values) for method `TMM`

- Weighting:

  logical, whether to compute (asymptotic binomial precision) weights

- Acutoff:

  numeric, cutoff on “A” values to use before trimming

## Value

Returns a numerical vector of normalization factors for each sample

## Details

Originally described for normalization of RNA sequences (Robinson &
Oshlack 2010), this function computes normalization (scaling) factors to
convert observed library sizes into effective library sizes. It uses the
method trimmed means of M-values proposed by Robinson & Oshlack (2010).
See the original publication and `edgeR` package for more information.
The method `MedR` is median ratio normalization; QN - quantile
normalization (see Maza, Elie, et al. 2013 for a comparison of methods).

## References

- Robinson MD, and Oshlack A (2010). A scaling normalization method for
  differential expression analysis of RNA-seq data. Genome Biology 11,
  R25

- Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor
  package for differential expression analysis of digital gene
  expression data. Bioinformatics 26

## Author

Piyal Karunarathne

## Examples

``` r
if (FALSE) vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path)
df<-hetTgen(vcf,"AD-tot",verbose=FALSE)
norm.fact(df) # \dontrun{}
#>             lib.size norm.factor
#> AT_PA_06_12    31164   1.0000000
#> CH_PA_02_03    90526   0.9381430
#> DE_PA_10_12   120888   0.9891654
#> FI_PA_17_02   136291   0.9478968
#> FR_PA_21_01   138419   0.9425358
#> FR_PA_21_17   122046   0.9582529
#> NO_PA_13_22   116385   0.9546818
#> RU_PA_19_06   154211   0.9571496
#> RU_PA_26_01    85270   0.9873622
#> SE_PA_16_18   113656   0.9596083

```
