# Import VCF file

Function to import raw single and multi-sample VCF files. The function
required the R-package `data.table` for faster importing.

## Usage

``` r
readVCF(vcf.file.path, verbose = FALSE)
```

## Arguments

- vcf.file.path:

  path to the vcf file

- verbose:

  logical. show progress

## Value

Returns a list with vcf table in a data frame, excluding meta data.

## Author

Piyal Karunarathne

## Examples

``` r
if (FALSE) vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path) # \dontrun{}
```
