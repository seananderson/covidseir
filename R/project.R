#' Project a SEIR fit
#'
#' Project a fit from [fit_seir()], possibly with a forecast.
#'
#' @param obj Output from [fit_seir()].
#' @param forecast_days Number of forecast days.
#' @param fixed_f_forecast Optional fixed f value for forecast.
#' @param ... Other arguments to pass to [rstan::sampling()].
#'
#' @details
#' Set a [future::plan()] and this function will operate in parallel
#' across MCMC iterations using \pkg{furrr}.
#'
#' @return
#' A [dplyr::tibble()] from [tidybayes::spread_draws()]:
#' \describe{
#'   \item{day}{Day}
#'   \item{data_type}{Data-type column from the case data}
#'   \item{mu}{Expected number of cases}
#'   \item{y_rep}{Posterior predictive replicate observation}
#'   \item{phi}{Posterior draw of phi, the NB2 dispersion parameter, if included}
#'   \item{R0}{Posterior draw of R0}
#'   \item{f_seq}{f parameter segment; segment 0 is the initial no physical distancing}
#'   \item{f_s}{f values; estimated for f_seq > 0}
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
#' # To use parallel processing, set a 'future' plan:
#' # library(future)
#' # plan(multisession, workers = parallel::detectCores() / 2)
#' p <- project_fit(m)
#' p
#'
#' library(ggplot2)
#' ggplot(p, aes(day, mu, group = .iteration)) +
#'   geom_line(alpha = 0.4) +
#'   geom_point(aes(y = y_rep), alpha = 0.2, pch = 21) +
#'   geom_point(aes(x = day, y = cases), tibble(day = seq_along(cases), cases),
#'     colour = "red", inherit.aes = FALSE) +
#'   labs(y = "Reported cases", x = "Day")

project_fit <- function(obj, forecast_days = 30, fixed_f_forecast = NULL, ...) {

  if (!identical(class(obj), "covidseir"))
    stop("`obj` must be of class `covidseir`.")

  d <- obj$stan_data
  p <- obj$post

  if (is.null(fixed_f_forecast)) fixed_f_forecast <- 0

  d$x_r[["fixed_f_forecast"]] <- fixed_f_forecast
  d$x_r

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
  d$x_i <- c(d$x_i, rep(.f_id, added_length))
  d$x_i[1] <- length(d$x_i) - 1L
  d$n_x_i <- length(d$x_i)

  initf_project <- function(post, i, stan_data) {
    R0 <- post$R0[i]
    f_s <- array(post$f_s[i, ])
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

  out <- furrr::future_map_dfr(seq_along(p$R0), function(i) {
  # out <- purrr::map_dfr(seq_along(p$R0), function(i) {
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

    if ("phi" %in% names(post)) {
      df <- tidybayes::spread_draws(
        fit, mu[day, data_type], y_rep[day, data_type], phi[data_type], R0)
    } else {
      df <- tidybayes::spread_draws(
        fit, mu[day, data_type], y_rep[day, data_type], R0)
    }
    df$.chain <- NULL
    df$.draw <- NULL
    df$.iteration <- NULL
    df2 <- tidybayes::spread_draws(fit, f_s[f_seg])
    df2$.chain <- NULL
    df2$.draw <- NULL
    df2$.iteration <- NULL
    df2 <- dplyr::bind_rows(df2, tibble(f_seg = 0, f_s = 1))
    df2 <- dplyr::left_join(tibble(f_seg = d$x_i[-1]), df2, by = "f_seg")
    df2$day <- seq_len(nrow(df2))
    df <- dplyr::left_join(df, df2, by = "day")
    df$.iteration <- i
    df
  })
  out
}

