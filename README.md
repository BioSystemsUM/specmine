
<!-- README.md is generated from README.Rmd. Please edit that file -->

# specmine

<!-- badges: start -->

[![R build
status](https://github.com/BioSystemsUM/specmine/workflows/R-CMD-check/badge.svg)](https://github.com/BioSystemsUM/specmine/actions)
[![R build
status](https://github.com/BioSystemsUM/specmine/workflows/CRAN/badge.svg)](https://github.com/BioSystemsUM/specmine/actions)
<!-- badges: end -->

The goal of *specmine* is to provide a set of methods for metabolomics
data analysis, including data loading in different formats,
pre-processing, metabolite identification, univariate and multivariate
data analysis, machine learning, feature selection and pathway analysis.
Case studies can be found on the website:
<http://bio.di.uminho.pt/metabolomicspackage/index.html>. This package
suggests ‘rcytoscapejs’, a package not in mainstream repositories. If
you need to install it, use:
`devtools::install_github('cytoscape/r-cytoscape.js@v0.0.7')`.

## Installation

You can install the released version of *specmine* from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("specmine")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BioSystemsUM/specmine")
```

## Example

This is a basic example which shows you how to load the namespace of
*specmine* and add it to your search list:

``` r
library(specmine)
```
