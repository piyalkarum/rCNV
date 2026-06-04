# Simulate and plot detection power of bias in allele ratios

This function simulates 95% confidence level Z-score based detection
power of allele biases for a given number of samples and a range of
depths

## Usage

``` r
power.bias(
  Dlist = c(2, 4, 8, 16),
  sam = 100,
  intensity = 0.005,
  nsims = 1000,
  p = 0.5,
  plot = TRUE
)
```

## Arguments

- Dlist:

  numerical. vector of depths values to be tested

- sam:

  numerical. number of samples

- intensity:

  numerical. frequency of bias

- nsims:

  numerical. number of simulations to be done for each sample

- p:

  numerical. expected allele ratio (0.5 for data with known sequencing
  biases)

- plot:

  logical. plot the output

## Value

Returns a list of detection probability values for the given range of
samples and depth

## Author

Pascal Milesi, Piyal Karunarathne
