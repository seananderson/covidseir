#' Get expected cases
#'
#' Returns expected number of cases on day d
#'
#' @param out Data from the ODE solver
#' @param pars Parameters
#' @param day The current date
#' @param sampFrac A vector of sample fractions
#' @param delayShape Delay shape parameter
#' @param delayScale Delay scale parameter
#'
#' @author Caroline Colijn, Jessica Stockdale
#' @export
#' @return The expected number of cases as a numeric value
getlambd <- function(out,
                     pars,
                     day,
                     sampFrac,
                     delayShape = 1.73,
                     delayScale = 9.85) {
  # relevant times to identify new cases
  ii <- which(out$time > day - 45 & out$time <= day)
  dx <- out$time[ii[2]] - out$time[ii[1]]
  # all new cases arising at each of those times
  incoming <- pars$k2 * (out$E2[ii] + out$E2d[ii])
  if (day <= length(sampFrac)) {
    thisSamp <- sampFrac[day]
  } else {
    thisSamp <- sampFrac[length(sampFrac)]
  }
  # each of the past times' contribution to this day's case count
  ft <- thisSamp * incoming * stats::dweibull(
    x = max(out$time[ii]) - out$time[ii],
    shape = delayShape,
    scale = delayScale
  )
  # return numerical integral of ft
  return(0.5 * (dx) * (ft[1] + 2 * sum(ft[2:(length(ft) - 1)]) + ft[length(ft)]))
}

#' Social distancing model
#'
#' The main ODEs for evaluating in R
#'
#' SEIR-type model with time-dependent social distancing. Social distancing
#' reduces frequency of contact. Individuals can move between distanced and not
#' distanced compartments.
#'
#' @param t time
#' @param state (S, E1, E2, I, Q, R, Sd, E1d, E2d, Id, Qd, Rd) S: Susceptible,
#'   E1: Exposed but not infectious, E2: Exposed and Infectious, I: Infectious,
#'   can be quarantined, R: Removed. The d compartments denote socially
#'   distanced individuals.
#' @param parms (N, D, R0, k1, k2, q, r, ur, f) f: strength of social distancing,
#'   r/(r+ur): frac of population who are distancing
#' @param sdtiming_function timing of social distancing function
#'
#' @return Time derivatives for input to ODE solver
#' @export
#' @author Caroline Colijn
socdistmodel <- function(t,
                         state,
                         parms,
                         sdtiming_function) {
  with(as.list(c(
    state,
    parms
  )), {
    f <- sdtiming_function(t, f1 = parms$f1, f2 = parms$f2, last_day_obs = parms$last_day_obs)
    dSdt <- -(R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) *
      S / N - r * S + ur * Sd
    dE1dt <- (R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) *
      S / N - k1 * E1 - r * E1 + ur * E1d
    dE2dt <- k1 * E1 - k2 * E2 - r * E2 + ur * E2d
    dIdt <- k2 * E2 - q * I - I / D - r * I + ur * Id
    dQdt <- q * I - Q / D - r * Q + ur * Qd
    dRdt <- I / D + Q / D - r * R + ur * Rd

    dSddt <- -(f * R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) *
      Sd / N + r * S - ur * Sd
    dE1ddt <- (f * R0 / (D + 1 / k2)) * (I + E2 + f * (Id + E2d)) *
      Sd / N - k1 * E1d + r * E1 - ur * E1d
    dE2ddt <- k1 * E1d - k2 * E2d + r * E2 - ur * E2d
    dIddt <- k2 * E2d - q * Id - Id / D + r * I - ur * Id
    dQddt <- q * Id - Qd / D + r * Q - ur * Qd
    dRddt <- Id / D + Qd / D + r * R - ur * Rd
    list(c(
      dSdt,
      dE1dt,
      dE2dt,
      dIdt,
      dQdt,
      dRdt,
      dSddt,
      dE1ddt,
      dE2ddt,
      dIddt,
      dQddt,
      dRddt
    ))
  })
}

# #' Linear decrease in f between two time points
# #'
# #' @param t time
# #' @param start_decline Day to start transition from f1 to f2
# #' @param end_decline Day to end transition from f1 to f2
# #' @param f1 Value before decline
# #' @param f2 Value after decline
# #'
# #' @author Andrew Edwards
# sdtiming_gradual <- function(t,
#   start_decline = 15, # start the decline the next day
#   end_decline = 22, # end decline at f2
#   f1 = pars$f1, # f value before decline
#   f2 = pars$f2) { # f value after decline
#   if (t < start_decline) {
#     return(f1)
#   }
#   if (t >= start_decline & t < end_decline) {
#     return(f2 + (end_decline - t) * (f1 - f2) / (end_decline - start_decline))
#   }
#   if (t >= end_decline) {
#     return(f2)
#   }
# }

# The social-distance timing function requires three arguments:
# t for time, will come from the model
# last_day_obs which will come from the model
# f1, which will come from pars$f1
# f2, which will come from pars$f2
# can use any other arguments it wants:
# It must match what was done in the model up to the last observation
# After that the sky is the limit.
return_f <- function(t, last_day_obs, f1, f2,
                     start_decline = 15, end_decline = 22, f_vec) {
  if (t < start_decline) {
    return(f1)
  }
  if (t >= start_decline && t < end_decline) {
    return(f2 + (end_decline - t) * (f1 - f2) / (end_decline - start_decline))
  }
  floor_t <- floor(t)
  if (t >= end_decline && floor_t <= last_day_obs) {
    return(f2)
  }
  if (t >= end_decline && floor_t > last_day_obs) {
    if (!is.null(f_vec)) {
      return(f_vec[floor_t])
    } else {
      return(f2)
    }
  }
}

project_fit_i <- function(obj, max_day = max(obj$time),
                          .i, sdfunc, data_col = 1L,
                          start_time = min(obj$time)) {
  .pars <- as.list(obj$pars)
  .pars$R0 <- obj$post$R0[.i]
  .pars$f2 <- obj$post$f2[.i]
  .pars$phi <- obj$post$phi[.i, data_col]
  .pars$last_day_obs <- obj$last_day_obs
  time <- seq(start_time, max_day, by = obj$time[2] - obj$time[1])
  .d <- as.data.frame(ode(
    y = obj$state_0,
    times = time,
    func = socdistmodel,
    parms = .pars,
    method = "rk4",
    sdtiming_func = sdfunc
  ))
  mu <- vapply(seq_len(max_day), function(x) {
    getlambd(
      .d,
      pars = .pars,
      sampFrac = obj$sampFrac[, data_col, drop = TRUE],
      day = x)
  }, FUN.VALUE = numeric(1L))
  out <- data.frame(
    day = seq(1, max_day),
    lambda_d = mu,
    y_rep = MASS::rnegbin(max_day, mu, theta = .pars$phi),
    iterations = .i,
    R0 = .pars$R0, f2 = .pars$f2, phi = .pars$phi
  )
  list(states = dplyr::mutate(.d, iterations = .i), cases = out)
}

#' Make projections with a fitted object
#'
#' @param obj An object from [fit_seir()].
#' @param proj_days The number of days to project beyond the last day of fitted
#'   data.
#' @param i A vector of posterior samples to use.
#' @param f_vec An optional vector of `f` values to use after the last day of
#'   fitted data. If left as the default of `NULL`, then the estimated `f2`
#'   value will be used in the projection
#'
#' @return A named list with states and cases as data frames.
#' @importFrom stats median
#' @export
project_fit <- function(obj,
                        proj_days = 0L,
                        i = seq_len(50L),
                        f_vec = NULL,
                        data_col = 1L,
                        start_time = min(obj$time)) {
  max_day <- obj$last_day_obs + proj_days
  # out <- future.apply::future_lapply(i, function(x) {
  # out <- furrr::future_map(i, function(x) {

  do_projection <- function(x) {
    project_fit_i(
      obj = obj,
      sdfunc = function(t, last_day_obs, f1, f2) {
        return_f(
          t = t,
          last_day_obs = last_day_obs,
          f1 = f1,
          f2 = f2,
          f_vec = f_vec,
          start_decline = obj$pars[["start_decline"]],
          end_decline = obj$pars[["end_decline"]]
        )
      },
      max_day = max_day,
      .i = x,
      data_col = data_col,
      start_time = start_time
    )
  }
  # out <- purrr::map(i, do_projection)
  out <- furrr::future_map(i, do_projection)
  states <- purrr::map_dfr(out, "states")
  cases <- purrr::map_dfr(out, "cases")
  list(states = states, cases = cases)
}

getu <- function(f, r) (r - f * r) / f

# J <- ncol(obj$post$phi)
# N <- max_day
# time_day_id0 <- obj$stan_data$time_day_id0
# time_day_id <- obj$stan_data$time_day_id
# y_hat <- .d
# x_r <- obj$stan_data$x_r
# sampFrac <- obj$stan_data$sampFrac
# .T <- obj$stan_data$T
# ft <- rep(0, .T)
# lambda_d <- matrix(nrow = N, ncol = J)
# delayShape <- obj$stan_data$delayShape
# delayScale <- obj$stan_data$delayScale
# dx <- obj$time[2] - obj$time[1]
#
# for (j in 1:J) {
#   for (n in 1:N) {
#     this_samp = sampFrac[n,j]
#     for (t in 1:.T) {
#       ft[t] = 0 # initialize at 0 across the full 1:T
#     }
#     # a fancy way of moving across a window of time:
#       for (t in time_day_id0[n]:time_day_id[n]) { # t is an increment here
#         k2 = x_r[4]
#         E2 = y_hat[t,3]
#         E2d = y_hat[t,9]
#         ft[t] = this_samp * k2 * (E2 + E2d) *
#           dweibull(time[time_day_id[n]] - time[t], delayShape[j], delayScale[j])
#       }
#     sum_ft_inner = 0 # initialize at 0
#     for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
#       sum_ft_inner = sum_ft_inner + ft[t]
#     }
#     lambda_d[n,j] = 0.5 * dx *
#       (ft[time_day_id0[n]] + 2 * sum_ft_inner + ft[time_day_id[n]])
#   }
# }
# out <- data.frame(
#   day = seq(1, max_day),
#   lambda_d = lambda_d[,data_col],
#   y_rep = MASS::rnegbin(max_day, lambda_d[,data_col], theta = .pars$phi),
#   iterations = .i,
#   R0 = .pars$R0, f2 = .pars$f2, phi = .pars$phi
# )
