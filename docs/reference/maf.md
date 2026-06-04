# Remove MAF allele

A function to remove the alleles with minimum allele frequency and keep
only a bi-allelic matrix when loci are multi-allelic

## Usage

``` r
maf(h.table, AD = TRUE, verbose = TRUE, parallel = FALSE)
```

## Arguments

- h.table:

  allele depth table generated from the function `hetTgen`

- AD:

  logical. If TRUE a allele depth table similar to `hetTgen` output will
  be returns; If `FALSE`, individual AD values per SNP will be returned
  in a list.

- verbose:

  logical. Show progress

- parallel:

  logical. whether to parallelize the process

## Value

A data frame or a list of minimum allele frequency removed allele depth

## Author

Piyal Karunarathne

## Examples

``` r
if (FALSE) mf<-maf(ADtable) # \dontrun{}
```
