#' Make a forecast with a SEIR fit
#'
#' Project a fit from [fit_seir()], possibly with a forecast. By default, the
#' forecast uses the estimated f values (fraction of normal
#' contacts encountered for those physical distancing). The function also includes
#' functionality to specify a vector of fixed f values starting on a given
#' future date.
#'
#' @param obj Output from [fit_seir()].
#' @param forecast_days Number of forecast days.
#' @param f_fixed_start Optional day to start changing f. Must be set if
#'   `f_fixed` is set.
#' @param f_fixed An optional vector of fixed f values for the forecast.
#'   Should be length `forecast_days - (f_fixed_start - nrow(daily_cases) -
#'   1)`. I.e. one value per day after `f_fixed_start` day.
#' @param iter MCMC iterations to include. Defaults to all.
#' @param return_states Logical for whether to return the ODE states.
#' @param ... Other arguments to pass to [rstan::sampling()].
#'
#' @importFrom dplyr bind_rows
#'
#' @details
#' Set a [future::plan()] and this function will operate in parallel
#' across MCMC iterations using \pkg{furrr}.
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
#' # To use parallel processing with multiple chains:
#' # options(mc.cores = parallel::detectCores() / 2)
#'
#' # Only using 100 iterations and 1 chain for a quick example:
#' m <- fit_seir(
#'   cases,
#'   iter = 100,
#'   chains = 1,
#'   samp_frac_fixed = s1
#' )
#' print(m)
#'
#' # For parallel processing (more important for more iterations):
#' # library(future)
#' # plan(multisession)
#' p <- forecast_seir(m)
#' p
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
#' p <- forecast_seir(m,
#'   forecast_days = 100,
#'   f_fixed_start = 53,
#'   f_fixed = c(rep(0.7, 60), rep(0.2, 30)),
#'   iter = 1:25
#' )
#' p
#' plot_ts(p)
forecast_seir <- function(
                         obj,
                         forecast_days = 30,
                         f_fixed_start = NULL,
                         f_fixed = NULL,
                         iter = seq_along(obj$post$R0),
                         return_states = FALSE,
                         ...) {
  if (!identical(class(obj), "covidseir")) {
    stop("`obj` must be of class `covidseir`.")
  }

  d <- obj$stan_data
  p <- obj$post

  stopifnot(
    (is.null(f_fixed_start) && is.null(f_fixed_start)) ||
      (!is.null(f_fixed_start) && !is.null(f_fixed_start))
  )

  stopifnot(length(f_fixed) == forecast_days - (f_fixed_start - nrow(d$daily_cases) - 1))

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

  if (!is.null(f_fixed_start)) {
    d$S <- d$S + length(f_fixed)
    est_f_forecast_days <- f_fixed_start - nrow(d$daily_cases) - 1
    d$x_i <- c(d$x_i, rep(d$x_i[length(d$x_i)], est_f_forecast_days))
    fixed_f_forecast_ids <- seq(max_f_seg_id + 1, max_f_seg_id + length(f_fixed))
    d$x_i <- c(d$x_i, fixed_f_forecast_ids)
  } else {
    d$x_i <- c(d$x_i, rep(d$x_i[length(d$x_i)], forecast_days))
  }
  d$x_i[["n_f_s"]] <- length(d$x_i) - 2 # 2 is number of non-f_s x_i values
  d$n_x_i <- length(d$x_i)

  initf_project <- function(post, i, stan_data) {
    R0 <- post$R0[i]
    f_s <- array(c(post$f_s[i, ], f_fixed))
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

  # incl_progress <- (d$last_day_obs + forecast_days) > 100

  pars <- c("R0", "f_s", "phi", "mu", "y_rep")
  if (return_states) pars <- c("y_hat")

  out <- furrr::future_map_dfr(iter, function(i) {
  # out <- purrr::map_dfr(iter, function(i) {
  # out <- lapply(iter, function(i) {
  # out <- future.apply::future_lapply(iter, function(i) {
    fit <- rstan::sampling(
      stanmodels$seir,
      data = d,
      iter = 1L,
      chains = 1L,
      init = function() initf_project(p, i, d),
      pars = pars,
      refresh = 0L,
      algorithm = "Fixed_param",
      ... = ...
    )
    post <- rstan::extract(fit)
    if (!return_states) {
      df <- data.frame(day = rep(d$days, d$J))
      df$data_type <- rep(seq_len(d$J), each = d$N)
      df$mu <- as.vector(post$mu[1, , ])
      df$y_rep <- as.vector(post$y_rep[1, , ])
      df$phi <- rep(as.vector(post$phi[1, ]), each = d$N)
    } else {
      df <- reshape2::melt(post$y_hat)
      df$iterations <- NULL
    }
    df$.iteration <- i
    df
  })
  out <- bind_rows(out)

  if (return_states) {
    variables_df <- dplyr::tibble(
      variable = names(obj$state_0),
      variable_num = seq_along(obj$state_0)
    )
    ts_df <- dplyr::tibble(time = d$time, time_num = seq_along(d$time))
    out <- dplyr::rename(out, time_num = Var2, variable_num = Var3)
    out <- dplyr::left_join(out, variables_df, by = "variable_num")
    out <- dplyr::left_join(out, ts_df, by = "time_num")
  } else {
    .forecast <- c(rep(FALSE, d$last_day_obs), rep(TRUE, forecast_days))
    out$forecast <- rep(.forecast, length(iter) * d$J)
  }

  if (!is.null(f_fixed_start) && !return_states) {
    .f_fixed <- c(
      rep(FALSE, d$last_day_obs),
      rep(FALSE, est_f_forecast_days),
      rep(TRUE, length(fixed_f_forecast_ids))
    )
    out$f_fixed <- rep(.f_fixed, length(iter) * d$J)
  }
  tibble::as_tibble(out)
}
