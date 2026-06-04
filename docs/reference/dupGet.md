# Detect deviants from SNPs; classify SNPs

Detect deviant SNPs using excess of heterozygotes (alleles that do not
follow HWE) and allelic-ratio deviations (alleles with ratios that do
not follow a normal Z-score or chi-square distribution). See details.

## Usage

``` r
dupGet(
  data,
  Fis,
  test = c("z.het", "z.05", "z.all", "chi.het", "chi.05", "chi.all"),
  intersection = FALSE,
  method = c("fisher", "chi.sq"),
  plot = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  data frame of the output of `allele.info`

- Fis:

  numeric. Inbreeding coefficient calculated using
  [`h.zygosity()`](h.zygosity.md) function

- test:

  character. type of test to be used for significance. See details

- intersection:

  logical, whether to use the intersection of the methods specified in
  `test` (if more than one)

- method:

  character. method for testing excess of heterozygotes. Fisher exact
  test (`fisher`) or Chi-square test (`chi.sq`)

- plot:

  logical. whether to plot the detected singlets and duplicates on
  allele ratio vs. proportion of heterozygotes plot.

- verbose:

  logical. show progress

- ...:

  additional parameters passed on to `plot`

## Value

Returns a data frame of snps/alleles with their duplication status

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

## Author

Piyal Karunarathne Qiujie Zhou

## Examples

``` r
if (FALSE) data(alleleINF)
DD<-dupGet(alleleINF,Fis=0.1,test=c("z.05","chi.05")) # \dontrun{}
#> categorizing deviant SNPs with 
#>  excess of heterozygotes 
#>  z.05
#>  chi.05

```
