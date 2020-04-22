# Bayesian SIR model to estimate physical-distancing effects

The model includes the following components: susceptible (S), exposed to the virus, asymptomatic and not infectious (E1), exposed, asymptomatic and infectious (E2), infectious (I), quarantined (Q), and recovered and deceased (R).

Install the package with:

```r
# install.packages("remotes")
remotes::install_github("seananderson/sirbayes")
```

If you also want to read the vignette:

```r
remotes::install_github("seananderson/sirbayes", build_vignettes = TRUE)
browseVignettes("sirbayes")
```

Or if you have a local copy of the repository:

```r
# with vignettes:
devtools::install("sirbayes", build_vignettes = TRUE)
browseVignettes("sirbayes")

# without (faster):
devtools::install("sirbayes")
```
