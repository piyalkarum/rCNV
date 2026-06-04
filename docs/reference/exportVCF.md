# Export VCF files

A function to export tables/matrices in VCF format to VCF files

## Usage

``` r
exportVCF(out.vcf, out.path, compress = TRUE)
```

## Arguments

- out.vcf:

  a matrix or data frame in vcf file format to be exported

- out.path:

  a character string of output path for the vcf file; should end in the
  name as the vcf file and .vcf. See examples

- compress:

  logical. whether to compress the output file. If `TRUE`, the file will
  be `.gz` compressed

## Value

Exports a vcf file to a given destination

## Author

Piyal Karunarathne

## Examples

``` r
if (FALSE) vcf.file.path <- paste0(path.package("rCNV"), "/example.raw.vcf.gz")
vcf <- readVCF(vcf.file.path)
exportVCF(vcf,"../exVcf.vcf") # \dontrun{}
#> compressing file
```
