# Bayesian SEIR model to estimate physical-distancing effects

This is a work in progress based on the model in the preprint [Estimating the impact of COVID-19 control measures using a Bayesian model of physical distancing](https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1) with extensions to include multiple data types and user-friendly helper functions.

### The model is very much a work in progress. Arguments and output format may still change.

Before installation, you will need a [C++ compiler installed to compile the Stan model](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

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
