sir <- function(.t, state, theta, x_r, x_i) {

  S     = state[1]
  E1    = state[2]
  E2    = state[3]
  I     = state[4]
  Q     = state[5]
  R     = state[6]
  Sd    = state[7]
  E1d   = state[8]
  E2d   = state[9]
  Id    = state[10]
  Qd    = state[11]
  Rd    = state[12]

  N     = x_r[1]
  D     = x_r[2]
  k1    = x_r[3]
  k2    = x_r[4]
  q     = x_r[5]
  ud     = x_r[6]
  ur    = x_r[7]
  f0    = x_r[8]
  f_ramp_rate = x_r[9]
  imported_cases = x_r[10]
  imported_window = x_r[11]

  last_day_obs = x_i[1]

  R0 = theta[1]
  start_decline = theta[2]
  end_decline = theta[3]
  f1 = theta[4]

  n_f = x_i[2] #// the number of f parameters in time (after f0)
  f_seg_id <- rep(NA, n_f) #// a lookup vector to grab the appropriate f parameter in time

  f <- NA #// will store the f value for this time point
  introduced <- NA #// will store the introduced cases

  dydt <- rep(NA, 12)

  # X
  # // integer version of the day for this time point:
  #   // (must use this workaround since floor() in Stan produces a real number
  #     // that can't be used to index an array)
  # day
  day = 1
  while ((day + 1) < floor(.t)) day = day + 1


  for (i in 1:n_f) {
    # // `i + 3` because of number of x_i before `f_seg_id`
    # // `+ 3` at end because of number of thetas before f_s thetas
    f_seg_id[i] = x_i[i + 3] + 3 # TO CHANGE: CHANGED i + 2 to i + 3!!
  }

  f = f0 #// business as usual before physical distancing
  if (.t < start_decline) {
    f = f0
  }
  #// allow a ramp-in of physical distancing:
  if (f_ramp_rate == 0.0) {
    if (.t >= start_decline && .t < end_decline) {
      f = f1 + (end_decline - .t) *
        (f0 - f1) / (end_decline - start_decline)
    }
  } else {
    if (.t >= start_decline && .t < end_decline) {
      X = (f0 - f1) / (exp(f_ramp_rate * (end_decline - start_decline)) - 1)
      f = f0 - X * (exp(f_ramp_rate * (.t - start_decline)) - 1)
    }
  }
  if (.t >= end_decline) {
    f = f1 # // default until a break gets possibly triggered...
    n_breaks <- x_i[3] # NEED TO ADD THIS TO x_i
    for (s in 1:n_breaks) { # EDITTED
      if (.t >= theta[4 + n_breaks + s]) { # EDITTED
        f = theta[4 + s] # EDITTED
      }
    }
  }
  cat("t =", .t, "\n")
  cat("f =", f, "\n\n")

  if (.t > last_day_obs) { # // THIS HAS MOVED AND CHANGED
    f = theta[f_seg_id[day]] # // the respective f segment
  }

  if (.t > last_day_obs && .t <= (last_day_obs + imported_window)) {
    introduced = imported_cases / imported_window
  } else {
    introduced = 0.0
  }

  dydt[1]  = -(R0/(D+1/k2)) * (I + E2 + f*(Id+E2d)) * S/N - ud*S + ur*Sd
  dydt[2]  = (R0/(D+1/k2)) * (I + E2 + f*(Id+E2d)) * S/N - k1*E1 -ud*E1 + ur*E1d
  dydt[3]  = k1*E1 - k2*E2 - ud*E2 + ur*E2d + introduced
  dydt[4]  = k2*E2 - q*I - I/D - ud*I + ur*Id
  dydt[5]  = q*I - Q/D - ud*Q + ur*Qd
  dydt[6]  = I/D + Q/D - ud*R + ur*Rd

  dydt[7]  = -(f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N + ud*S - ur*Sd
  dydt[8]  = (f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N - k1*E1d +ud*E1 - ur*E1d
  dydt[9]  = k1*E1d - k2*E2d + ud*E2 - ur*E2d
  dydt[10] = k2*E2d - q*Id - Id/D + ud*I - ur*Id
  dydt[11] = q*Id - Qd/D + ud*Q - ur*Qd
  dydt[12] = Id/D + Qd/D + ud*R - ur*Rd

  dydt
}


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
dat$date <- lubridate::ymd("2020-03-01") + dat$day - 1
samp_frac <- c(rep(0.14, 13), rep(0.21, 40 - 13), rep(0.21, 11))
samp_frac <- c(samp_frac, rep(0.37, nrow(dat) - length(samp_frac)))

f_seg <- c(0, rep(1, nrow(dat) - 1))
day_new_f <- which(dat$date == lubridate::ymd("2020-05-01"))
f_seg[seq(day_new_f, length(f_seg))] <- 2
day_ch <- which(dat$date == lubridate::ymd("2020-06-01"))
f_seg[seq(day_ch, length(f_seg))] <- 3
f_seg
fit <- covidseir::fit_seir(
  daily_cases = dat$value,
  samp_frac_fixed = samp_frac,
  f_seg = f_seg,
  i0_prior = c(log(8), 1),
  start_decline_prior = c(log(15), 0.1),
  end_decline_prior = c(log(22), 0.1),
  f_prior = cbind(c(0.4, 0.5, 0.6), c(0.2, 0.2, 0.2)),
  R0_prior = c(log(2.6), 0.2),
  N_pop = 5.1e6, # BC population
  iter = 500, # number of samples
  fit_type = "optimizing" # for speed only
)

d <- fit$stan_data
theta <- rep(NA, d$S + 3)
theta[1] = 2.6
theta[2] = 4
theta[3] = 10
f_s <- c(0.66, 0.77, 0.88)
for (s in 1:d$S) {
  # // `s + 3` because of number of thetas before f_s
  theta[s + 3] = f_s[s];
}
theta
theta <- c(theta, c(22, 44))

d$x_i <- c(d$x_i[1:2], c("n_breaks" = max(d$x_i[-c(1:2)]) - 1), d$x_i[-c(1:2)])

out <- sapply(15:50, function(x) sir(.t = x, state = c(1000, 1000, d$y0_vars), theta = theta, x_r = d$x_r, x_i = d$x_i))
