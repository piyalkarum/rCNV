# Simulate Allele Frequencies

This function simulates allele frequencies of a desired population size
under HWE

## Usage

``` r
sim.als(n = 500, nrun = 10000, res = 0.001, plot = TRUE)
```

## Arguments

- n:

  desired populations size (set this value same as your actual
  population size for an accurate simulation)

- nrun:

  number of simulations to run on each allele frequency. The higher this
  number, the closer the simulations will be to the theoretical values
  (at the cost of computer power); 10000 is an optimal value.

- res:

  desired resolution of the theoretical allele frequency

- plot:

  logical. whether to plot the simulation

## Value

A list of two matrices:

1.  allele_freqs: theoretical allele frequency

2.  simulated_freqs: simulated frequencies at different confidence
    intervals

## Author

Piyal Karunarathne, Pascal Milesi

## Examples

``` r
if (FALSE) alleles <- sim.als(n=200,nrun=1000,res=0.001,plot=TRUE) # \dontrun{}
```
