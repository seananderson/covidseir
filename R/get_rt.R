#' Get equivalent Rt values
#'
#' @param obj Output from [fit_seir()].
#' @param lag Time lag for the inputs for Rt.
#' @param forecast_days Number of projection days.
#' @param iter MCMC iterations to include. Defaults to all.
#'
#' @return A data frame of state values with a column of Rt values. Number of
#'   rows is 1 per time step. An `NA` for Rt where insufficient lagged data are
#'   available.
#' @export
#'
#' @details
#' See [project_seir()] for an example.

get_rt <- function(obj, lag = 30, forecast_days = 0,
                   iter = seq_along(obj$post$R0)) {
  p <- project_seir(obj,
    forecast_days = forecast_days,
    iter = iter, return_states = TRUE
  )
  states <- tidyr::pivot_wider(p, names_from = variable, values_from = value)
  temp <- dplyr::group_by(states, .iteration)
  temp <- dplyr::group_split(temp)
  temp <- purrr::map_dfr(temp,
    get_rt_one_iteration,
    pars = obj$pars, lag = lag
  )
  states$Rt <- temp$Rt
  states
}

get_rt_one_iteration <- function(states, pars, lag = 30) {
  lag_avail <- states$time - min(states$time)
  inds <- which(lag_avail >= lag)
  dt <- states$time[2] - states$time[1]
  numback <- lag / (dt) - 1
  denominator <- vector(length = nrow(states))

  inc_e2 <- pars["k1"] * (states$E1 + states$E1d)

  # have experienced that it doesn't matter which incident compartment we use
  # and this doesn't require reading in posterior time-dep f vals, and R0s
  for (i in inds) {
    pastinds <- i - (1:numback)
    pastcases <- inc_e2[pastinds]
    s <- seq(from = 0, by = dt, length.out = length(pastinds))
    integrand <- pastcases * infkernel(
      rate1 = pars[["k1"]], rate2 = pars[["k2"]],
      rate3 = pars[["q"]] + 1 / pars[["D"]], .t = s
    )
    trapezoid_weights <- rep(1, length(integrand))
    trapezoid_weights[1] <- 0.5
    trapezoid_weights[length(trapezoid_weights)] <- 0.5
    denominator[i] <- dt * sum(integrand * trapezoid_weights)
  }
  denominator[which(denominator == 0)] <- NA
  dplyr::tibble(Rt = inc_e2 / denominator)
}

# @param rate1 rate of first exponential
# @param rate2 rate of second exponential
# @param t time
# @return density for the convolution of two exponential processes f1 * f2 (t)
convodens <- function(rate1, rate2, .t) {
  if (rate1 == rate2) {
    rate1 * rate2 * .t * exp(-rate1 * .t)
  } else {
    rate1 * rate2 * exp(-rate2 * .t) *
      (exp((rate2 - rate1) * .t) - 1) / (rate2 - rate1)
  }
}

# @param rate1 rate of first exponential (here, k1)
# @param rate2 rate of second exponential (here, k2)
# @param rate3 rate of third exponential (here, rate leaving infectious class)
# @param t time
# @return infectivity kernel based density
infkernel <- function(rate1, rate2, rate3, .t) {
  if (rate1 == rate2 || rate2 == rate3 || rate3 == rate1) {
    stop("Rates cannot be identical.", call. = FALSE)
  }
  (rate1 * rate2 * rate3 / (rate2 + rate3)) *
    (exp(-rate1 * .t) * (rate2 + rate3 - rate1) / ((rate1 - rate2) * (rate1 - rate3)) +
      exp(-rate2 * .t) * rate3 / ((rate1 - rate2) * (rate3 - rate2)) +
      exp(-rate3 * .t) * rate2 / ((rate1 - rate3) * (rate2 - rate3)))
}
