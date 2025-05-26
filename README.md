<!-- badges: start -->

[![R-CMD-check](https://github.com/piyalkarum/rCNV/workflows/R-CMD-check/badge.svg)](https://github.com/piyalkarum/rCNV/actions)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov test
coverage](https://codecov.io/gh/piyalkarum/rCNV/branch/master/graph/badge.svg)](https://app.codecov.io/gh/piyalkarum/rCNV?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/rCNV)](https://CRAN.R-project.org/package=rCNV)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/rCNV)](https://cran.r-project.org/package=rCNV)
<!-- badges: end -->

# rCNV <img src='man/figures/logo.png' align='right' height='139' />

# rCNV: An R package for detecting copy number variants from SNPs data

<span style="color: dodgerblue;">Piyal Karunarathne, Qiujie Zhou, Klaus Schliep, and
Pascal Milesi</span>

**rCNV** was designed to identify duplicates (CNV) from SNPs data with
ease.

<img src="vignettes/dup.plot.parrotfish.png" width="400" />

***For a comprehensive tutorial on the package, go to
<https://piyalkarum.github.io/rCNV/> and navigate to “Get started” where
all the functions and usage are explained with ample examples.***

## NEWS 
rCNV is currently developing methods to detect multicopy regions from whole genome sequencing (WGS) data using maximum likelihood ratios. 
*For a sneak peek, go to the section 2.4 at <https://piyalkarum.github.io/rCNV/>*

## Installation

-   CRAN link <https://cran.r-project.org/package=rCNV>

<!-- -->

    install.packages("rCNV")

-   You can install the development version of rCNV from
    [GitHub](https://github.com/) with:

<!-- -->

        if (!requireNamespace("devtools", quietly = TRUE)) 
            install.packages("devtools") 
        devtools::install_github("piyalkarum/rCNV", build_vignettes = TRUE)

Please don’t forget to cite us if you use the package.

## How to cite

-   Karunarathne P, Zhou Q, Schliep K, Milesi P. A comprehensive framework for detecting copy number variants from single nucleotide polymorphism data: 'rCNV', a versatile r package for paralogue and CNV detection. Mol Ecol Resour. 2023 Jul 29. doi:<http://doi.org/10.1111/1755-0998.13843>

-   Karunarathne, P., Zhou, Q., Schliep, K., & Milesi, P. (2022). A new framework for detecting copy number variants from single nucleotide polymorphism data: ‘rCNV’, a versatile R package for paralogs and CNVs detection. BioRxiv, 2022.10.14.512217. doi:<http://doi.org/10.1101/2022.10.14.512217>
