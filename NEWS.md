# rCNV 1.1.0 (first update)

## Changes
* sig.het function is 10x faster
* dupGet function categorizes SNPs into deviants and non-deviants not duplicates and singlets
* The new function cnv categorizes deviant SNPs to CNVs using the given statistical approaches or using unsupervised clustering
* minor issues of chi-square and Z-score calculations are fixed
* allele.info calculates excess of heterozygotes in addition to Z-score and chi-square values
* sig.hets function also accepts allele depth table for input
* new function added 'power.bias' to plot detection bias from simulations
* webpage updated to match the changes

# rCNV 1.0.0 (Release date: 29/03/2022)

## Changes
* All the major functions completed
* bugs fixed
* Website updated
* this is the first release of the package
