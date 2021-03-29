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

cases <- c(1, 1, 4, 6, 4, 6, 1, 3, 5, 11, 20, 12,
  25, 12, 27, 46, 41, 60, 47, 62, 46, 43, 88, 59, 91, 62, 62, 37,
  30, 82, 56, 47, 55, 31, 21, 13, 51, 40, 33, 34, 41, 28, 15, 30,
  44, 14, 51, 25, 30, 15, 34, 63, 33, 37, 72, 61, 15, 58, 41, 26,
  29, 23, 33, 21, 16, 24, 29, 21, 14, 8, 17, 10, 13, 16, 14, 21,
  8, 11, 5, 19, 11, 22, 12, 4, 9, 8, 7, 10, 5, 11, 9, 16, 4, 23,
  2, 4, 12, 7, 9, 10, 14, 12, 17, 14, 14, 9, 11, 19, 9, 5, 9, 6,
  17, 11, 12, 21, 12, 9, 13, 3, 12, 17, 10, 9, 12, 16, 7, 10, 19,
  18, 22, 25, 23, 20, 13, 19, 24, 33, 44, 19, 29, 33, 34, 32, 30,
  30, 21, 21, 24, 47, 26, 39, 52, 29, 48, 28, 43, 46, 47, 51, 38,
  42, 46, 86, 73, 90, 100, 86, 51, 77, 65, 84, 90, 115, 79, 81)

# s1 <- c(rep(0.1, 13), rep(0.2, length(cases) - 13))
s1 <- create_segments_vector(
  start_date = "2020-01-01",
  total_days = length(cases),
  segments = c("2020-01-13"),
  values = c(0.1, 0.2)
)

transmission_vec <- create_ramp_vector(
  start_date = "2020-01-01",
  total_days = length(cases),
  start_ramp = "2020-04-01",
  ramp_length = 50,
  ramp_max = 1.001
)
print(transmission_vec)
plot(cases)
plot(transmission_vec)

# fake the data so we can fit it:
cases_fake <- cases
cases_fake[91:length(cases_fake)] <- round(cases_fake[91:length(cases_fake)] + exp(0.06 * 1:86))
plot(cases_fake)
points(cases, col = "green")
abline(v = 91)

m <- fit_seir(
  cases_fake,
  stan_model = stan_mod,
  iter = 100,
  fit_type = "optimizing",
  samp_frac_fixed = s1,
  transmission_vec = transmission_vec
  # algorithm = "BFGS"
)
print(m)

transmission_vec <- create_ramp_vector(
  start_date = "2020-01-01",
  total_days = length(cases) + 30,
  start_ramp = "2020-02-12",
  ramp_length = 7 * 4,
  ramp_max = 10.0
)
print(transmission_vec)

p <- project_seir(m,
  stan_model = stan_mod,
  forecast_days = 30,
  transmission_vec = transmission_vec
)

obs_dat <- data.frame(day = seq_along(cases), value = cases)

tidy_seir(p) %>% plot_projection(obs_dat = obs_dat)
