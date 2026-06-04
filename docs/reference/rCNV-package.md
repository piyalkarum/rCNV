# rCNV: Detect Copy Number Variants from SNP Data

Functions in this package import filtered VCF files of SNP data and
generate datasets to detect copy number variants, visualize them, and
perform downstream analyses (e.g. environmental association analyses).

## Main WGS functions

- [`allele.info.WGS3`](allele.info.WGS3.md):

  Duplication/deletion detection using a negative-binomial depth model;
  no GC correction.

- [`allele.info.WGS.GCcor`](allele.info.WGS.GCcor.md):

  Same framework with per-SNP GC-bias correction via SNP-centred
  windows.

## See also

Useful links:

- <https://piyalkarum.github.io/rCNV/>

- <https://cran.r-project.org/package=rCNV>

- Report bugs at <https://github.com/piyalkarum/rCNV/issues>

## Author

Piyal Karunarathne <piyalkarumail@yahoo.com>
