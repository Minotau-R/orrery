
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/Minotau-R/orrery)](https://github.com/Minotau-R/orrery/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/Minotau-R/orrery)](https://github.com/Minotau-R/orrery/pulls)
[![R
BiocCheck](https://github.com/Minotau-R/orrery/actions/workflows/test.yml/badge.svg)](https://github.com/Minotau-R/orrery/actions/workflows/test.yml)
<!-- badges: end -->

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install the released version of
`orrery`:

``` r
# Not yet on CRAN!  install.packages('orrery')
```

Or install the development version from this repository:

``` r
install.packages("remotes")
remotes::install_github("minotau-R/orrery")
```

## Getting started using orrery

### Generate random data

``` r
library(orrery)
library(dirmult)
```

``` r
# Create a table with 20 features and 2 samples
x <- rdirichlet(2, runif(20))

# CLR transform
x.clr <- dclr(x)

# compute distance
dist.x <- dist(x.clr, method = "euclidean")

# Copy x into y
y <- rbind(x, perturb_by_relab(x[1L, ], by = 10))

y.clr <- dclr(y)

dist.y <- dist(y.clr, method = "euclidean")

dist.x
#>          1
#> 2 20.92217
dist.y
#>           1         2
#> 2 20.922171          
#> 3  9.847272 22.894288
```

[See the vignettes on the package
site.](https://minotau-r.github.io/orrery/articles/orrery.html)
