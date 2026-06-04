# Allele info example data

Semi-randomly generated data from the function dup.snp.info. Data
contains depth and proportion values of 2857 snps

## Usage

``` r
data(alleleINF)
```

## Format

An object of class `data.frame` with 2857 rows and 28 columns.

## Source

[Chinook Salmon sequence reads McKinney et al.
2017](https://zenodo.org/record/5025423#.Yj2XKRDMLyQ)

## References

- Larson, W. A., Seeb, L. W., Everett, M. V., Waples, R. K., Templin, W.
  D., & Seeb, J. E. (2014). Genotyping by sequencing resolves \#'
  shallow population structure to inform conservation of Chinook salmon
  (Oncorhynchus tshawytscha). Evolutionary Applications, 7(3)

- McKinney, G. J., Waples, R. K., Seeb, L. W., & Seeb, J. E. (2017).
  Paralogs are revealed by proportion of heterozygotes and deviations in
  read ratios in genotyping by sequencing data from natural populations.
  Molecular Ecology Resources, 17(4)

## Examples

``` r
data(alleleINF)
with(alleleINF,plot(medRatio~propHet))
```
