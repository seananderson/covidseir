#' Project a SEIR fit
#'
#' Project a fit from [fit_seir()], possibly with a forecast. By default, the
#' forecast uses the estimated f values (fraction of normal
#' contacts encountered for those physical distancing). The function also includes
#' functionality to specify a vector of fixed f values starting on a given
#' future date.
#'
#' @param obj Output from [fit_seir()].
#' @param forecast_days Number of forecast days.
#' @param f_s_fixed_start Optional day to start changing f. Must be set if
#'   `f_s_fixed` is set.
#' @param f_s_fixed Optional fixed f values for forecast. Should be length
#'   `forecast_days - (f_s_fixed_start - nrow(daily_cases) - 1)`. I.e. one value
#'   per day after `f_s_fixed_start` day.
#' @param iter MCMC iterations to include. Defaults to all.
#' @param ... Other arguments to pass to [rstan::sampling()].
#'
#' @details
#' Set a [future::plan()] and this function will operate in parallel
#' across MCMC iterations using \pkg{future.apply}.
#'
#' @return
#' A data frame:
#' \describe{
#'   \item{day}{Day}
#'   \item{data_type}{Data-type column from the case data}
#'   \item{mu}{Expected number of cases}
#'   \item{y_rep}{Posterior predictive replicate observation}
#'   \item{phi}{Posterior draw of phi, the NB2 dispersion parameter, if included}
#'   \item{.iteration}{The MCMC iteration}
#' }
#' @export
#' @examples
#' cases <- c(
#'   0, 0, 1, 3, 1, 8, 0, 6, 5, 0, 7, 7, 18, 9, 22, 38, 53, 45, 40,
#'   77, 76, 48, 67, 78, 42, 66, 67, 92, 16, 70, 43, 53, 55, 53, 29,
#'   26, 37, 25, 45, 34, 40, 35
#' )
#' # Example fixed sample fractions:
#' s1 <- c(rep(0.1, 13), rep(0.2, length(cases) - 13))
#'
#' # To use parallel processing:
#' # options(mc.cores = parallel::detectCores() / 2)
#'
#' m <- fit_seir(
#'   cases,
#'   iter = 150,
#'   chains = 1,
#'   samp_frac_fixed = s1
#' )
#' print(m)
#'
#' p <- project_seir(m)
#' p
#'
#'
#' plot_ts <- function(p) {
#'   if (require("ggplot2")) { # just for R CMD check in this example
#'     ggplot(p, aes(day, mu, group = .iteration)) +
#'       geom_line(alpha = 0.4) +
#'       geom_point(aes(y = y_rep), alpha = 0.2, pch = 21) +
#'       geom_point(aes(x = day, y = cases), data.frame(day = seq_along(cases), cases),
#'         colour = "red", inherit.aes = FALSE
#'       ) +
#'       labs(y = "Reported cases", x = "Day")
#'   }
#' }
#' plot_ts(p)
#'
#' p <- project_seir(m,
#'   forecast_days = 100,
#'   f_s_fixed_start = 53,
#'   f_s_fixed = c(rep(0.7, 60), rep(0.2, 30)),
#'   iter = 1:25
#' )
#' p
#' plot_ts(p)

project_seir <- function(
                    obj,
                    forecast_days = 30,
                    f_s_fixed_start = NULL,
                    f_s_fixed = NULL,
                    iter = seq_along(obj$post$R0),
                    ...) {
  if (!identical(class(obj), "covidseir")) {
    stop("`obj` must be of class `covidseir`.")
  }

  d <- obj$stan_data
  p <- obj$post

  stopifnot(
    (is.null(f_s_fixed_start) && is.null(f_s_fixed_start)) ||
      (!is.null(f_s_fixed_start) && !is.null(f_s_fixed_start))
  )

  stopifnot(length(f_s_fixed) == forecast_days - (f_s_fixed_start - nrow(d$daily_cases) - 1))

  days <- seq(1L, nrow(d$daily_cases) + forecast_days)
  time_increment <- d$time[2] - d$time[1]
  time <- seq(-30, max(days), time_increment)

  # find the equivalent time of each day (end):
  time_day_id <- vapply(days, get_time_id, numeric(1), time = time)
  # find the equivalent time of each day (start):
  time_day_id0 <- vapply(days, get_time_day_id0, numeric(1),
    time = time, days_back = obj$days_back
  )

  d$time_day_id0 <- time_day_id0
  d$time_day_id <- time_day_id
  d$time <- time
  d$days <- days
  d$T <- length(time)
  d$N <- length(days)

  .s <- d$samp_frac_seg[length(d$samp_frac_seg)]
  added_length <- nrow(d$daily_cases) + forecast_days - length(d$samp_frac_seg)
  d$samp_frac_seg <- c(d$samp_frac_seg, rep(.s, added_length))

  .s <- d$samp_frac_fixed[nrow(d$samp_frac_fixed), ]
  d$samp_frac_fixed <- rbind(
    d$samp_frac_fixed,
    do.call(rbind, replicate(added_length, .s, simplify = FALSE))
  )

  .f_id <- d$x_i[length(d$x_i)]

  max_f_seg_id <- max(d$x_i[-c(1:2)]) # 1:2 is the non-f_s x_i values

  if (!is.null(f_s_fixed_start)) {
    d$S <- d$S + length(f_s_fixed)
    d$x_i <- c(d$x_i, rep(
      d$x_i[length(d$x_i)],
      f_s_fixed_start - nrow(d$daily_cases) - 1
    ))
    d$x_i <- c(d$x_i, seq(max_f_seg_id + 1, max_f_seg_id + 1 + length(f_s_fixed)))
  } else {
    d$x_i <- c(d$x_i, rep(d$x_i[length(d$x_i)], forecast_days))
  }
  d$x_i[["n_f_s"]] <- length(d$x_i) - 2 # 2 is number of non-f_s x_i values
  d$n_x_i <- length(d$x_i)

  initf_project <- function(post, i, stan_data) {
    R0 <- post$R0[i]
    f_s <- array(c(post$f_s[i, ], f_s_fixed))
    if ("phi" %in% names(post)) {
      phi <- array(post$phi[i, ])
    } else {
      phi <- 1 # fake
    }
    if ("samp_frac" %in% names(post)) {
      samp_frac <- array(post$samp_frac[i, ])
    } else {
      samp_frac <- rep(1, stan_data$n_samp_frac) # fake
    }
    list(R0 = R0, f_s = f_s, phi = phi, samp_frac = samp_frac)
  }

  # out <- furrr::future_map_dfr(iter, function(i) {
  # out <- purrr::map_dfr(iter, function(i) {
  # out <- lapply(iter, function(i) {
  out <- future.apply::future_lapply(iter, function(i) {
    fit <- rstan::sampling(
      stanmodels$seir,
      data = d,
      iter = 1L,
      chains = 1L,
      init = function() initf_project(p, i, d),
      pars = c("R0", "f_s", "phi", "mu", "y_rep"),
      refresh = 0L,
      algorithm = "Fixed_param",
      ... = ...
    )
    post <- rstan::extract(fit)
    df <- data.frame(day = rep(d$days, d$J))
    df$data_type <- rep(seq_len(d$J), each = d$N)
    df$mu <- as.vector(post$mu[1, , ])
    df$y_rep <- as.vector(post$y_rep[1, , ])
    df$phi <- rep(as.vector(post$phi[1, ]), each = d$N)
    df$.iteration <- i
    df
  })
  do.call(rbind, out)
}
