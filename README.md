# covidseir

> Bayesian SEIR Modelling for Multivariate COVID-19 Case Data

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/seananderson/covidseir.svg?branch=master)](https://travis-ci.org/seananderson/covidseir)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

## Overview

**covidseir** fits a Bayesian SEIR (Susceptible, Exposed, Infectious, Recovered)
model to daily COVID-19 case data. The package focuses on estimating the
fraction of the usual contact rate for individuals participating in physical
distancing (social distancing). The model is coded in
[**Stan**](https://mc-stan.org/). The model can accommodate multiple types of
case data at once (e.g., reported cases, hospitalizations, ICU admissions) and
accounts for delays between symptom onset and case appearance.

The model is a continuation of the one described in the preprint:\
Anderson, S. C., Edwards, A. M., Yerlanov, M., Mulberry, N., Stockdale, J., Iyaniwura, S. A., Falcao, R. C., Otterstatter, M. C., Irvine, M. A., Janjua, N. Z., Coombs, D., & Colijn, C. (2020). Estimating the impact of COVID-19 control measures using a Bayesian model of physical distancing. MedRxiv, 2020.04.17.20070086. https://doi.org/10.1101/2020.04.17.20070086

### This package is a work in progress. Arguments and output format may still change and not all functionality has been tested yet. ###

## Installation

Before installation, you will need a [C++ compiler installed to compile the Stan model](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). **Note that if you are on Windows, you will need to follow the 'Configuration' section [here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows) not the one from the previous sentence**. However, you likely do not need to follow the rest of the instructions from that page to install rstan itself from source.

In particular, the following must return `TRUE` before continuing:

```r
# install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```

Then, install the package with:

```r
# install.packages("remotes")
remotes::install_github("seananderson/covidseir")
```

See the examples in `?fit_seir` and `?project_seir`.
