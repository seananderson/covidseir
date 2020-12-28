library(covidseir)
library(dplyr)
library(ggplot2)
ymd <- lubridate::ymd
options(future.rng.onMisuse="ignore")

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
ggplot(dat, aes(date, value)) + geom_line()

# Based on estimation with hospital data in other model:
samp_frac <- c(rep(0.14, 13), rep(0.21, 38))
samp_frac <- c(samp_frac, rep(0.37, nrow(dat) - length(samp_frac)))
f_seg <- c(0, rep(1, nrow(dat) - 1))
day_new_f <- which(dat$date == ymd("2020-05-01"))
f_seg[seq(day_new_f, length(f_seg))] <- 2
day_ch <- which(dat$date == ymd("2020-06-01"))
f_seg[seq(day_ch, length(f_seg))] <- 3

fit <- covidseir::fit_seir(
  daily_cases = dat$value,
  samp_frac_fixed = samp_frac,
  f_seg = f_seg,
  i0_prior = c(log(8), 1),
  e_prior = c(0.8, 0.05),
  start_decline_prior = c(log(15), 0.1),
  end_decline_prior = c(log(22), 0.1),
  f_prior = cbind(c(0.4, 0.5, 0.6), c(0.2, 0.2, 0.2)),
  R0_prior = c(log(2.6), 0.2),
  N_pop = 5.1e6, # BC population
  iter = 500, # number of posterior samples
  fit_type = "optimizing" # for speed only
)

f_s <- fit$post$f_s
mean(f_s[,ncol(f_s)])

p <- project_seir(fit,
  forecast_days = 21,
  f_fixed_start = nrow(dat) + 1,
  f_fixed = rep(0.5, 21),
  iter = 1:30
)
g1 <- plot_projection(tidy_seir(p), dat) + coord_cartesian(ylim = c(0, 300), expand = FALSE)

# roughly same with f_multi:

# get final f posterior:
final_f <- f_s[,ncol(f_s)]
hist(final_f)

# desired mean f = 0.5

ratio <- 0.5 / final_f
hist(ratio)
(f_multi <- mean(ratio))

p1 <- project_seir(fit,
  forecast_days = 21,
  f_fixed_start = nrow(dat) + 1,
  f_multi = rep(f_multi, 21),
  f_multi_seg = ncol(fit$post$f_s),
  iter = 1:30
)
g2 <- plot_projection(tidy_seir(p1), dat) + coord_cartesian(ylim = c(0, 300), expand = FALSE)

cowplot::plot_grid(g1, g2)
