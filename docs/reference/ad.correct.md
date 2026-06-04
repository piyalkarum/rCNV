# Correct allele depth values

A function to correct depth values with odd number of coverage values
due to sequencing anomalies or miss classification where genotype is
homozygous and depth values indicate heterozygosity. The function adds a
value of one to the allele with the lowest depth value for when odd
number anomalies or make the depth value zero for when miss-classified.
The genotype table must be provided for the latter.

## Usage

``` r
ad.correct(
  het.table,
  gt.table = NULL,
  odd.correct = TRUE,
  verbose = TRUE,
  parallel = FALSE
)
```

## Arguments

- het.table:

  allele depth table generated from the function `hetTgen`

- gt.table:

  genotype table generated from the function hetTgen

- odd.correct:

  logical, to correct for odd number anomalies in AD values. default
  `TRUE`

- verbose:

  logical. show progress. Default `TRUE`

- parallel:

  logical. whether to parallelize the process

## Value

Returns the coverage corrected allele depth table similar to the output
of `hetTgen`

## Author

Piyal Karunarathne

## Examples

``` r
if (FALSE) adc<-ad.correct(ADtable) # \dontrun{}
```
