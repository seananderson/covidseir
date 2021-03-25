# remove.packages(c("StanHeaders", "rstan"))
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() / 2)
library(covidseir)

# stan_mod <- rstan::stan_model(system.file("stan", "seir.stan", package = "covidseir"))
stan_mod <- stan_model("inst/stan/seir.stan")

cases <- c(
  0, 0, 1, 3, 1, 8, 0, 6, 5, 0, 7, 7, 18, 9, 22, 38, 53, 45, 40,
  77, 76, 48, 67, 78, 42, 66, 67, 92, 16, 70, 43, 53, 55, 53, 29,
  26, 37, 25, 45, 34, 40, 35
)
s1 <- c(rep(0.1, 13), rep(0.2, length(cases) - 13))
m <- fit_seir(
  cases,
  chains = 1,
  stan_model = stan_mod,
  iter = 100,
  fit_type = "optimizing",
  samp_frac_fixed = s1
)
print(m)

p <- project_seir(m, stan_model = stan_mod)


obs_dat <- data.frame(day = seq_along(cases), value = cases)

tidy_seir(p) %>% plot_projection(obs_dat = obs_dat)
