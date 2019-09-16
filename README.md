
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FluxGapsR

<!-- badges: start -->

<!-- badges: end -->

This is a package including four gap-filling methods for soil
respiration data investigated in the study of Zhao et al. 2019 (to be
published). The four methods are referred to as non-linear least squares
(NLS), artificial neural networks (ANN), singular spectrum analysis
(SSA) and expectation-maximization (EM).

**Package installation**

First, make sure the package `devtools` is installed in R. If not,
install the package by:

``` r
install.packages("devtools")
```

Then, install the `FluxGapsR` package in R by:

``` r
devtools::install_github("junbinzhao/FluxGapsR")
```

The functioning of the package is based on other R packages:
`tidyverse`,`spectral.methods`,`minpack.lm`,`mtsdi`,`neuralnet` and they
must be installed before using the functions in the `FluxGapsR` package.

**Please cite the package as:** (TBD)
