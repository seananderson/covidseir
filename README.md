# Bayesian SEIR model to estimate physical-distancing effects

This is a work in progress based on the model in the preprint [Estimating the impact of COVID-19 control measures using a Bayesian model of physical distancing](https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1).

Install the package with:

```r
# install.packages("remotes")
remotes::install_github("seananderson/covidseir")
```

If you also want to read the vignette:

```r
remotes::install_github("seananderson/covidseir", build_vignettes = TRUE)
browseVignettes("covidseir")
```

Or if you have a local copy of the repository:

```r
# with vignettes:
devtools::install("covidseir", build_vignettes = TRUE)
browseVignettes("covidseir")

# without vignettes (faster):
devtools::install("covidseir")
```
