# covidseir

> Bayesian SEIR model to estimate physical-distancing effects

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/covidseir)](https://CRAN.R-project.org/package=covidseir)
<!-- badges: end -->

## Overview

__covidseir__ fits a Bayesian SEIR (Susceptible, Exposed, Infectious, Recovered)
model to daily COVID-19 case data. The package focuses on estimating the
fraction of usual contacts encountered for individuals participating in physical
distancing (social distancing). The model is coded in
[__Stan__](https://mc-stan.org/). The model can accommodate multiple types of
case data at once (e.g., reported cases, hospitalizations, ICU admissions) and
accounts for delays between symptom onset and case appearance.

The model is a continuation of the model described in the preprint [Estimating
the impact of COVID-19 control measures using a Bayesian model of physical
distancing](https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1).

### The model is very much a work in progress. Arguments and output format may still change and not all functionality has been tested yet. ###

## Installation

Before installation, you will need a [C++ compiler installed to compile the Stan model](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

In particular, the following must return `TRUE` before continuing:

```r
pkgbuild::has_build_tools(debug = TRUE)
```

Then, install the package with:

```r
# install.packages("remotes")
remotes::install_github("seananderson/covidseir")
```

If you also want to read the compiled vignette (slower):

```r
remotes::install_github("seananderson/covidseir", build_vignettes = TRUE)
browseVignettes("covidseir")
```

Alternatively, read the source code for the vignette [here](https://github.com/seananderson/covidseir/tree/master/vignettes).

Or if you have a local copy of the repository:

```r
# with vignettes:
devtools::install("covidseir", build_vignettes = TRUE)
browseVignettes("covidseir")

# without vignettes (faster):
devtools::install("covidseir")
```
