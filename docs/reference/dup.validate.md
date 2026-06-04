# Validate detected deviants/cnvs

This function will validate the detected duplicated-SNPs (deviants/cnvs)
using a moving window approach (see details)

## Usage

``` r
dup.validate(d.detect, window.size = 100, scaf.size = 10000)
```

## Arguments

- d.detect:

  a data frame of detected duplicates or deviants from the outputs of
  `dupGet` or `cnv`

- window.size:

  numerical. a single value of the desired moving window size (default
  `100` bp)

- scaf.size:

  numerical. scaffold size to be checked. i.e. the chromosome/scaffolds
  will be split into equal pieces of this size default=10000

## Value

A data frame of deviant/cnv ratios (column `cnv.ratio`) for a split of
the chromosome/scaffold given by the `scaf.size`; this ratio is an
average value of the percentage of deviants/cnvs present within the
given `window.size` for each split
(`chromosome/scaffold length/sacf.size`); the start and the end
positions of each split is given in the `start` and `end` columns

## Details

Loci/SNP positions correctly ordered according to a reference sequence
is necessary for this function to work properly. The list of
deviants/cnvs provided in the `d.detect` will be split into pices of
`scaf.size` and the number of deviants/cnvs will be counted along each
split with a moving window of `window.size`. The resulting percentages
of deviants/cnvs will be averaged for each scaf.size split; this is the
`cnv.ratio` column in the output. Thus, ideally, the `cnv.ratio` is a
measure of how confident the detected deviants/cnvs are in an actual
putative duplicated region withing the given `scaf.size`. This ratio is
sensitive to the picked window size and the scaf.size; as a rule of
thumb, it is always good to use a known gene length as the scaf.size, if
you need to check a specific gene for the validity of the detected
duplicates. Please also note that this function is still in its
`beta-testing` phase and also under development for non-mapped reference
sequences. Therefore, your feedback and suggestions will be highly
appreciated.

## Author

Piyal Karunarathne

## Examples

``` r
if (FALSE) { # \dontrun{
# suggestion to visualize dup.validate output

library(ggplot2)
library(dplyr)

dvs<-dupGet(alleleINF,test=c("z.05","chi.05"))
dvd<-dup.validate(dvs,window.size = 1000)

# Example data frame
df <- data.frame(dvd[,3:5])
df$cnv.ratio<-as.numeric(df$cnv.ratio)

# Calculate midpoints
df <- df %>%
  mutate(midpoint = (start + end) / 2)

ggplot() +
  # Horizontal segments for each start-end range
  geom_segment(data = df, aes(x = start, xend = end,
  y = cnv.ratio, yend = cnv.ratio), color = "blue") +
  # Midpoints line connecting midpoints of each range
  geom_path(data = df, aes(x = midpoint, y = cnv.ratio), color = "red") +
  geom_point(data = df, aes(x = midpoint, y = cnv.ratio), color = "red") +
  # Aesthetic adjustments
  theme_minimal() +
  labs(title = "CNV Ratio along a Continuous Axis with Midpoint Fluctuation",
      x = "Genomic Position",
       y = "CNV Ratio")
} # }

```
