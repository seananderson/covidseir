## You *must* have rstan 2.26 installed
## Currently, that means doing this:

# remove.packages(c("StanHeaders", "rstan"))
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

## If you need to undo this and revert to the current CRAN rstan:
# remove.packages(c("StanHeaders", "rstan"))
# install.packages("StanHeaders")
# install.packages("rstan")

library(rstan)
library(tidyverse)
library(lubridate)
rstan_options(auto_write = TRUE) # to cache the model building below
options(mc.cores = parallel::detectCores() / 2)
library(covidseir)

## sometimes the model caching breaks for me, in that case run:
# file.remove(system.file("stan", "seir.rds", package = "covidseir"))

## You *must* build the Stan model first to use rstan 2.26
## The following assumes you have installed the latest `stan-2.26` branch covidseir
## remotes::install_github("seananderson/covidseir", ref = "stan-2.26")
stan_mod <- rstan::stan_model(system.file("stan", "seir.stan", package = "covidseir"))

## or we've been developing with this, but you must be in the covidseir working directory locally:
file.remove("inst/stan/seir.rds")
stan_mod <- stan_model("inst/stan/seir.stan")

# vignette data:
dat <- structure(list(value = c(1, 1, 4, 6, 4, 6, 1, 3, 5, 11, 20, 12,
  25, 12, 27, 46, 41, 60, 47, 62, 46, 43, 88, 59, 91, 62, 62, 37,
  30, 82, 56, 47, 55, 31, 21, 13, 51, 40, 33, 34, 41, 28, 15, 30,
  44, 14, 51, 25, 30, 15, 34, 63, 33, 37, 72, 61, 15, 58, 41, 26,
  29, 23, 33, 21, 16, 24, 29, 21, 14, 8, 17, 10, 13, 16, 14, 21,
  8, 11, 5, 19, 11, 22, 12, 4, 9, 8, 7, 10, 5, 11, 9, 16, 4, 23,
  2, 4, 12, 7, 9, 10, 14, 12, 17, 14, 14, 9, 11, 19, 9, 5, 9, 6,
  17, 11, 12, 21, 12, 9, 13, 3, 12, 17, 10, 9, 12, 16, 7, 10, 19,
  18, 22, 25, 23, 20, 13, 19, 24, 33, 44, 19, 29, 33, 34, 32, 30,
  30, 21, 21, 24, 47, 26, 39, 52, 29, 48, 28, 43, 46, 47, 51, 38,
  42, 46, 86, 73, 90, 100, 86, 51, 77, 65, 84, 90, 115, 79, 81),
  day = 1:176), row.names = c(NA, -176L),
  class = "data.frame")
dat <- dplyr::as_tibble(dat)
dat$date <- ymd("2020-03-01") + dat$day - 1

# the usual fit setup:
samp_frac <- c(rep(0.14, 13), rep(0.21, 38))
samp_frac <- c(samp_frac, rep(0.37, nrow(dat) - length(samp_frac)))

f_seg <- c(0, rep(1, nrow(dat) - 1))
day_new_f <- which(dat$date == ymd("2020-05-01"))
f_seg[seq(day_new_f, length(f_seg))] <- 2
day_ch <- which(dat$date == ymd("2020-06-01"))
f_seg[seq(day_ch, length(f_seg))] <- 3
f_seg

first_day <- min(dat$date)
last_day <- 300 # how many days to create dates for
lut <- dplyr::tibble(
  day = seq_len(last_day),
  date = seq(first_day, first_day + length(day) - 1, by = "1 day")
)

## VOC example:

# fake some data so we can fit it
# On July 1 2020, we will add a step change in R0:
cases <- dat$value
replacement_mu <- cases[121:length(cases)] + exp(0.12 * 1:56)
set.seed(12038)
replacement <- MASS::rnegbin(length(replacement_mu), mu = replacement_mu, theta = 5)
plot(replacement)
dat$value_voc <- dat$value
dat$value_voc[121:length(dat$value_voc)] <- replacement

plot(dat$value_voc)
points(dat$value, col = "green")
abline(v = 91)
abline(v = 121)

plot(dat$value_voc, log = "y")
points(dat$value, col = "green")
abline(v = 91)
abline(v = 121)

## create a vector of multipliers on R0
## you can do it by hand or use this helper function:
transmission_vec <- covidseir::create_ramp_vector(
  start_date = "2020-03-01",
  total_days = length(dat$value),
  start_ramp = "2020-07-01",
  ramp_length = 40,
  ramp_max = 1.4
)
print(transmission_vec)
plot(dat$value)
plot(transmission_vec)

# fit once, ignoring VoCs:
fit_voc_ignore <- covidseir::fit_seir(
  daily_cases = dat$value_voc,
  stan_model = stan_mod, # Note that for now you must pass the rstan model here!
  samp_frac_fixed = samp_frac,
  f_seg = f_seg,
  i0_prior = c(log(8), 1),
  e_prior = c(0.8, 0.05),
  start_decline_prior = c(log(15), 0.1),
  end_decline_prior = c(log(22), 0.1),
  f_prior = cbind(c(0.4, 0.5, 0.6), c(0.2, 0.2, 0.2)),
  R0_prior = c(log(2.6), 0.2),
  N_pop = 5.1e6,
  iter = 100,
  fit_type = "optimizing",
  ode_control = c(1e-08, 1e-07, 1e+07)
)

# fit with the VoC ramp:
fit_voc <- covidseir::fit_seir(
  daily_cases = dat$value_voc,
  stan_model = stan_mod, # Note that for now you must pass the rstan model here!
  samp_frac_fixed = samp_frac,
  f_seg = f_seg,
  i0_prior = c(log(8), 1),
  e_prior = c(0.8, 0.05),
  start_decline_prior = c(log(15), 0.1),
  end_decline_prior = c(log(22), 0.1),
  f_prior = cbind(c(0.4, 0.5, 0.6), c(0.2, 0.2, 0.2)),
  R0_prior = c(log(2.6), 0.2),
  N_pop = 5.1e6,
  iter = 100,
  fit_type = "optimizing",
  transmission_vec = transmission_vec,
  ode_control = c(1e-08, 1e-07, 1e+07)
)

# plot
p_voc <- project_seir(fit_voc,
  stan_model = stan_mod, # Note that for now you must pass the rstan model here!
  forecast_days = 0, iter = 1:30
)
p_voc_tidy <- tidy_seir(p_voc)
p_voc_tidy <- dplyr::left_join(p_voc_tidy, lut, by = "day")
g1 <- p_voc_tidy %>%
  plot_projection(obs_dat = dat, date_column = "date", value_column = "value_voc") +
  geom_vline(xintercept = ymd("2020-07-01"), lty = 2)

p_voc_ignore <- project_seir(fit_voc_ignore,
  stan_model = stan_mod, # Note that for now you must pass the rstan model here!
  forecast_days = 0, iter = 1:30
)
p_voc_ignore_tidy <- tidy_seir(p_voc_ignore)
p_voc_ignore_tidy <- dplyr::left_join(p_voc_ignore_tidy, lut, by = "day")
g2 <- p_voc_ignore_tidy %>%
  plot_projection(obs_dat = dat, date_column = "date", value_column = "value_voc") +
  geom_vline(xintercept = ymd("2020-07-01"), lty = 2)

cowplot::plot_grid(g1 + ggtitle("Account for VoC"), g2 + ggtitle("Ignore VoC"), ncol = 1)

# see the final f posterior has changed too:
par(mfrow = c(1, 2))
plot(density(fit_voc$post$f_s[,3]), main = "Last f: Account for VoC")
plot(density(fit_voc_ignore$post$f_s[,3]), main = "Last f: Ignore VoC")

# with forecast_days > 0
p_voc2 <- project_seir(fit_voc,
  stan_model = stan_mod,
  forecast_days = 10, iter = 1:15
)
p_voc_tidy <- tidy_seir(p_voc2)
p_voc_tidy <- dplyr::left_join(p_voc_tidy, lut, by = "day")
p_voc_tidy %>%
  plot_projection(obs_dat = dat, date_column = "date", value_column = "value_voc") +
  geom_vline(xintercept = ymd("2020-07-01"), lty = 2)
