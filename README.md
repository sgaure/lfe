
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lfe

<!-- badges: start -->

[![R-CMD-check](https://github.com/MatthieuStigler/lfe/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MatthieuStigler/lfe/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of lfe is to speed up the estimation of linear models with
large fixed effects. It includes support for instrumental variables,
conditional F statistics for weak instruments, robust and multi-way
clustered standard errors, as well as limited mobility bias correction.
See Gaure (2013) <doi:10.1016/j.csda.2013.03.024> and Gaure 2014
<doi:10.1002/sta4.68>.

## Installation

You can install the development version of lfe like so:

``` r
remotes::install_github("MatthieuStigler/lfe")
```

## Example

This is a basic example which shows you the speed improvement over base
R for fixed effects estimation.

``` r
library(lfe) # fixed effects estimation
library(tradepolicy) # intl trade data
library(dplyr) # data cleaning/transforming

training_data <- agtpa_applications %>% 
  mutate(
    log_trade = log(trade),
    log_dist = log(dist),
    exp_year = paste(exporter, year, sep = "_"),
    imp_year = paste(importer, year, sep = "_")
  ) %>% 
  filter(trade > 0, exporter != importer, year %in% seq(1986, 2006, 4)) %>% 
  select(year, log_trade, log_dist, cntg, lang, clny, rta, exp_year, imp_year)

# note the difference with the | operator to indicate the FEs
# this is just an example, here I am not estimating a PPML model or anything
# in the state of the art
fml1 <- 0 + log_trade ~ 
  log_dist + cntg + lang + clny + rta + exp_year + imp_year # base

fml2 <- log_trade ~ 
  log_dist + cntg + lang + clny + rta | exp_year + imp_year # lfe

lm(fml1, data = training_data)

felm(fml2, data = training_data)
```

## Testing

For a complete test with `devtools::check()`, you need to run
`sudo apt-get devtools` or similar before. The package is written in C,
in the future I shall try to rewrite it in C++ to ease long term
maintenance.

The package also needs additional testing. At the present time,the tests
it cover around 30% of the written lines.
