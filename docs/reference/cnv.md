# Find CNVs from deviants

Categorize deviant and non-deviant into "singlets" and "duplicates"
based on the statistical approaches specified by the user. The
intersection of all the stats provided will be used in the
categorization. If one would like to use the intersection of at least
two stats, this can be specified in the `n.ints`

## Usage

``` r
cnv(
  data,
  test = c("z.het", "z.05", "z.all", "chi.het", "chi.05", "chi.all"),
  filter = c("intersection", "kmeans"),
  WGS = TRUE,
  ft.threshold = 0.05,
  plot = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame of allele information generated with the function
  `allele.info`

- test:

  vector of characters. Type of test to be used for significance. See
  details

- filter:

  character. Type of filter to be used for filtering CNVs. default
  `kmeans`. See details.

- WGS:

  logical. test parameter. See details

- ft.threshold:

  confidence interval for filtering `default = 0.05`

- plot:

  logical. Plot the detection of duplicates. default `TRUE`

- verbose:

  logical. show progress

- ...:

  other arguments to be passed to `plot`

## Value

Returns a data frame of SNPs with their detected duplication status

## Details

SNP deviants are detected with both excess of heterozygosity according
to HWE and deviant SNPs where depth values fall outside of the normal
distribution are detected using the following methods:

- Z-score test \\Z\_{x} = \sum\_{i=1}^{n} Z\_{i}\\; \\Z\_{i} =
  \frac{\left ( (N\_{i}\times p)- N\_{Ai} \right )}{\sqrt{N\_{i}\times
  p(1-p)}}\\

- chi-square test \\X\_{x}^{2} = \sum\_{i-1}^{n} X\_{i}^{2}\\;
  \\X\_{i}^{2} = (\frac{(N\_{i}\times p - N\_{Ai})^2}{N\_{i}\times p} +
  \frac{(N\_{i}\times (1 - p)- (N\_{i} - N\_{Ai}))^2}{N\_{i}\times
  (1-p)})\\

See references for more details on the methods

Users can pick among Z-score for heterozygotes (`z.het, chi.het`), all
allele combinations (`z.all, chi.all`) and the assumption of no probe
bias p=0.5 (`z.05, chi.05`)

`filter` will determine whether the `intersection` or `kmeans`
clustering of the provided `test`s should be used in filtering CNVs. The
intersection uses threshold values for filtering and kmeans use
unsupervised clustering. Kmeans clustering is recommended if one is
uncertain about the threshold values.

`WGS` is a test parameter to include or exclude coefficient of variance
(cv) in kmeans. For data sets with more homogeneous depth distribution,
excluding cv improves CNV detection. If you're not certain about this,
use `TRUE` which is the default.

## Author

Piyal Karunarathne Qiujie Zhou

## Examples

``` r
if (FALSE) data(alleleINF)
DD<-cnv(alleleINF) # \dontrun{}
#> categorizing putative duplicates with 
#>  excess of heterozygotes z.all & chi.all

```
