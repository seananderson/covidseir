library(dplyr)
library(ggplot2)
ymd <- lubridate::ymd

dat <- structure(list(value = c(1, 1, 4, 6, 4, 6, 1, 3, 5, 11, 20, 12,
  25, 12, 27, 46, 41, 60, 47, 62, 46, 43, 88, 59, 91, 62, 62, 37,
  30, 82, 56, 47, 55, 31, 21, 13, 51, 40, 33, 34, 41, 28, 15, 30,
  44, 14, 51, 25, 30, 15, 34, 63, 33, 37, 72, 61, 15, 58, 41, 26,
  29, 23, 33, 21, 16, 24, 29, 21, 14, 8, 17, 10, 13, 16, 14, 21,
  8, 11, 5, 19, 11, 22, 12, 4, 9, 8, 7, 10, 5, 11, 9, 16, 4, 23,
  2, 4, 12, 7, 9, 10, 14, 12, 17, 14, 14, 9, 11, 19, 9, 5, 9, 6,
  17, 11, 12, 21, 12, 9, 13, 3, 12, 17, 10, 9, 12, 16, 7, 10, 19,
  18, 22, 25, 23, 20, 13, 19, 24, 33, 44, 19, 29, 33, 34, 32, 30,
  30, 21, 21, 24, 47, 26, 39, 52, 29, 48, 28, 43, 46, 47, 52, 38,
  43, 46, 87, 75, 91, 102, 89, 51, 77, 63, 86, 91, 117, 80, 81,
  54, 58, 65, 125, 98, 114, 93, 58, 97, 89, 114, 136, 126, 108,
  83, 89, 155, 120, 154, 118, 83, 59, 126, 179, 143, 122, 116,
  130, 100, 83), day = 1:206), row.names = c(NA, -206L), class = c("tbl_df",
    "tbl", "data.frame"))
dat$date <- ymd("2020-03-01") + dat$day - 1

ggplot(dat, aes(date, value)) + geom_line()

samp_frac <- c(rep(0.14, 13), rep(0.21, 38))
samp_frac <- c(samp_frac, rep(0.37, nrow(dat) - length(samp_frac)))

f_seg <- c(0, rep(1, nrow(dat) - 1))
day_new_f <- which(dat$date == ymd("2020-05-01"))
f_seg[seq(day_new_f, length(f_seg))] <- 2
day_ch <- which(dat$date == ymd("2020-06-01"))
f_seg[seq(day_ch, length(f_seg))] <- 3
f_seg

dow <- data.frame(day_of_week = rep(gl(7, 1), 999)[-c(1:6)])[1:nrow(dat),,drop=FALSE]
dow$fake <- rep(0, nrow(dow))

# factors:
X <- model.matrix(fake ~ 0 + day_of_week, dow)

# library(mgcv)
# dow$day_of_week <- as.numeric(dow$day_of_week)
# dow$cases <- dat$value
# mgam <- gam(cases ~ 0 + s(day_of_week, bs = "cc", k = 5), data = dow, family = nb)
# plot(mgam)
# X <- model.matrix(mgam)
X[which(dat$date < ymd("2020-07-01")), ] <- 1/7

# sine wave

# xs <- sin(2 * pi * dat$day/7)
# X <- matrix(xs, ncol = 1)
# plot(dat$day, xs, type = "l")


# d2 <- dat
# d2$day_of_week <- dow$day_of_week
# filter(d2, date >= ymd("2020-07-01")) %>%
# # d2 %>%
#   group_by(day_of_week) %>%
#   summarise(m = mean(value))

fit <- fit_seir(
  daily_cases = dat$value,
  samp_frac_fixed = samp_frac,
  f_seg = f_seg,
  i0_prior = c(log(8), 1),
  e_prior = c(0.8, 0.05),
  start_decline_prior = c(log(15), 0.1),
  end_decline_prior = c(log(22), 0.1),
  f_prior = cbind(c(0.4, 0.5, 0.6), c(0.2, 0.2, 0.2)),
  R0_prior = c(log(2.6), 0.2),
  N_pop = 5.1e6,
  iter = 500,
  fit_type = "optimizing",
  X = X
)

print(fit)

purrr::map_dfr(1:7, ~{
  tibble(dow = .x, b = fit$post[[paste0("beta[", .x, "]")]])
}) %>%
  ggplot(aes(dow, exp(b), group = dow)) + geom_violin() +
  scale_x_continuous(breaks = 1:7, labels = c("M", "T", "W", "T", "F", "S", "S"))
#
# dow2 <- data.frame(day_of_week = rep(gl(7, 1), 999)[-c(1:6)])[1:nrow(dat),,drop=FALSE]
# dow2$fake <- rep(0, nrow(dow))
# # dow$day_of_week[dat$date < ymd("2020-07-01")] <- sort(unique(dow$day_of_week))[4]

# X <- model.matrix(fake ~ 0 + day_of_week, dow)
# X[which(dat$date < ymd("2020-08-01")), ] <- 1/7

p <- project_seir(fit, iter = 1:200)
p2 <- tidy_seir(p)
plot_projection(p2, obs_dat = dat)

p2 <- left_join(p2, dat)

plot_residuals(p2, dat, fit, date_column = "date")


xx <- left_join(p, dat)

xx$log_prob <- stats::dnbinom(xx$value, mu = xx$mu, size = xx$phi, log = TRUE)

group_by(xx, day, date) %>%
  # mutate(log_prob = exp(log_prob)) %>%
  summarise(med = median(log_prob), lwr = quantile(log_prob, probs = 0.05),
    upr =  quantile(log_prob, probs = 0.95)) %>%
  ggplot(aes(date, med)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
  ylab("Log probability of observed data") +
  theme(axis.title.x = element_blank()) +
  ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")
