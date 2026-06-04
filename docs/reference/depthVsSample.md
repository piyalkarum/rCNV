# Simulate median allele ratios for varying number of samples and depth values

This function will simulate the expected median allele ratios under HWE
for given ranges of no. of samples and depth coverage values. This is
useful if you need to find the cutoff values of allele ratios for
different no. of samples and depth of coverage values in your data set.

## Usage

``` r
depthVsSample(
  cov.len = 100,
  sam.len = 100,
  nsims = 1000,
  plot = TRUE,
  col = c("#1C86EE", "#00BFFF", "#DAA520", "#FF0000")
)
```

## Arguments

- cov.len:

  max value of depth of coverage to be simulated

- sam.len:

  maximum no. of samples to be simulated

- nsims:

  numerical. no. of simulations to be done for each combination of
  samples and depth depth and no. samples ranges

- plot:

  logical. Whether to plot the output (a plot of no. samples vs median
  depth of coverage colored by median allele ratios)

- col:

  character. Two colors to add to the gradient

## Value

A matrix of median allele ratios where rows are the number of samples
and columns are depth of coverage values

## Author

Pascal Milesi, Piyal Karunarathne

## Examples

``` r
if (FALSE) depthVsSample(cov.len=100,sam.len=100) # \dontrun{}
```
