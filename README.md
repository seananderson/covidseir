# Bayesian SEIR model to estimate physical-distancing effects

The model includes the following components: susceptible (S), exposed to the virus, asymptomatic and not infectious (E1), exposed, asymptomatic and infectious (E2), infectious (I), quarantined (Q), and recovered or deceased (R).

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
