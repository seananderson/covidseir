#' Make a projection with a SEIR fit
#'
#' Predict from a [fit_seir()] object, possibly with a future prediction. By default, the
#' projection uses the estimated f values (fraction of normal
#' contacts encountered for those physical distancing). The function also includes
#' functionality to specify a vector of fixed f values starting on a given
#' future date.
#'
#' @param obj Output from [fit_seir()].
#' @param forecast_days Number of projection days.
#' @param f_fixed_start Optional day to start changing f. Must be set if
#'   `f_fixed` is set.
#' @param f_fixed An optional vector of fixed f values for the projection.
#'   Should be length `forecast_days - (f_fixed_start - nrow(daily_cases) -
#'   1)`. I.e. one value per day after `f_fixed_start` day.
#' @param f_multi Multiplicative vector of f values. Same structure as `f_fixed`.
#' @param f_multi_seg Which `f` to use for `f_multi`.
#' @param imported_cases Number of cases to import starting on first projection day.
#' @param imported_window Number of days over which to distribute imported cases.
#' @param iter MCMC iterations to include. Defaults to all.
#' @param return_states Logical for whether to return the ODE states.
#' @param parallel Use parallel processing via \pkg{future} and \pkg{furrr}?
#' @param X An optional model matrix that acts additively on log expected cases.
#' @param ... Other arguments to pass to [rstan::sampling()].
#'
#' @importFrom dplyr bind_rows
#'
#' @details
#' Set a [future::plan()] and this function will operate in parallel
#' across MCMC iterations using \pkg{future} and \pkg{furrr}.
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
#' # Using the MAP estimate for speed in this example:
#' m <- fit_seir(
#'   cases,
#'   iter = 100,
#'   i0_prior = c(log(8), 0.2),
#'   samp_frac_fixed = s1,
#'   fit_type = "optimizing"
#' )
#' print(m)
#'
#' p <- project_seir(m)
#' p
#'
#' obs_dat <- data.frame(day = seq_along(cases), value = cases)
#'
#' library(magrittr) # for %>%
#' tidy_seir(p) %>% plot_projection(obs_dat = obs_dat)
#'
#' tidy_seir(p) %>%
#'   plot_residuals(obs_dat = obs_dat, obj = m)
#'
#' # for parallel processing (optional)
#' # future::plan(future::multisession)
#' p <- project_seir(m,
#'   forecast_days = 100,
#'   f_fixed_start = 53,
#'   f_fixed = c(rep(0.7, 60), rep(0.2, 30)),
#'   iter = 1:25
#' )
#' p
#' tidy_seir(p) %>% plot_projection(obs_dat = obs_dat)
#'
#' # fake example to show optional Rt colouring:
#' proj <- tidy_seir(p)
#' plot_projection(proj, obs_dat = obs_dat, Rt = rlnorm(nrow(proj)), col = "#00000050")
#'
#'
#' # Get threshold for increase:
#' # future::plan(future::multisession) # for parallel processing (optional)
#' # (only using 40 iterations for a fast example)
#' thresh <- get_threshold(m, iter = 1:40, show_plot = TRUE)
#' mean(thresh)
#'
#' # Get doubling time (prevalence is declining in this example)
#' get_doubling_time(m, iter = 1:40, show_plot = TRUE)
#'
#' states_with_Rt <- get_rt(m)
#' states_with_Rt
#'
#' library(ggplot2)
#' ggplot(states_with_Rt, aes(time, Rt, group = .iteration)) +
#'   geom_line(alpha = 0.5)
project_seir <- function(
                         obj,
                         forecast_days = 0,
                         f_fixed_start = NULL,
                         f_fixed = NULL,
                         f_multi = NULL,
                         f_multi_seg = NULL,
                         iter = seq_along(obj$post$R0),
                         return_states = FALSE,
                         imported_cases = 0,
                         imported_window = 1,
                         parallel = TRUE,
                         X = obj$stan_data$X,
                         ...) {
  if (!identical(class(obj), "covidseir")) {
    stop("`obj` must be of class `covidseir`.")
  }
  d <- obj$stan_data
  p <- obj$post

  # stopifnot(
  #   (is.null(f_fixed_start) && is.null(f_fixed_start)) ||
  #     (!is.null(f_fixed_start) && !is.null(f_fixed_start))
  # )

  if (!is.null(f_fixed) && !is.null(f_multi)) {
    stop("!is.null(f_fixed) && !is.null(f_multi)", call. = FALSE)
  }
  if (!is.null(f_fixed)) {
    stopifnot(length(f_fixed) == forecast_days - (f_fixed_start - nrow(d$daily_cases) - 1))
  }

  if (!is.null(f_multi)) {
    stopifnot(length(f_multi) == forecast_days - (f_fixed_start - nrow(d$daily_cases) - 1))
  }

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

  if (ncol(X) == 0L) {
    X <- matrix(0, ncol = 0L, nrow = d$N)
    d$X <- X
  }
  if (nrow(X) < d$N) {
    stop("`nrow(X)` < `length(days)`", call. = FALSE)
  }

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
  beta_sd_ext <- d$f_prior[d$S, 2]
  beta_mean_ext <- d$f_prior[d$S, 1]
  if (!is.null(f_fixed)) {
    d$f_prior <- rbind(
      d$f_prior,
      matrix(c(
        rep(beta_mean_ext, length(f_fixed)),
        rep(beta_sd_ext, length(f_fixed))
      ), ncol = 2, nrow = length(f_fixed))
    )
    d$S <- d$S + length(f_fixed)
    est_f_forecast_days <- f_fixed_start - nrow(d$daily_cases) - 1
    d$x_i <- c(d$x_i, rep(d$x_i[length(d$x_i)], est_f_forecast_days))
    fixed_f_forecast_ids <- seq(max_f_seg_id + 1, max_f_seg_id + length(f_fixed))
    d$x_i <- c(d$x_i, fixed_f_forecast_ids)
  } else if (!is.null(f_multi)) {
    # FIXME: DRY
    d$f_prior <- rbind(
      d$f_prior,
      matrix(c(
        rep(beta_mean_ext, length(f_multi)),
        rep(beta_sd_ext, length(f_multi))
      ), ncol = 2, nrow = length(f_multi))
    )
    d$S <- d$S + length(f_multi)
    est_f_forecast_days <- f_fixed_start - nrow(d$daily_cases) - 1
    d$x_i <- c(d$x_i, rep(d$x_i[length(d$x_i)], est_f_forecast_days))
    fixed_f_forecast_ids <- seq(max_f_seg_id + 1, max_f_seg_id + length(f_multi))
    d$x_i <- c(d$x_i, fixed_f_forecast_ids)
  } else {
    d$x_i <- c(d$x_i, rep(d$x_i[length(d$x_i)], forecast_days))
  }

  d$x_i[["n_f_s"]] <- length(d$x_i) - 2 # 2 is number of non-f_s x_i values
  d$n_x_i <- length(d$x_i)

  if (!"imported_cases" %in% names(d$x_r)) {
    warning("Adding to `d$x_r`. Looks like an old model.", call. = FALSE)
    d$x_r <- c(d$x_r, "imported_cases" = 0)
  }
  if (!"imported_window" %in% names(d$x_r)) {
    warning("Adding to `d$x_r`. Looks like an old model.", call. = FALSE)
    d$x_r <- c(d$x_r, "imported_window" = 1)
  }
  stopifnot(identical(names(d$x_r), c(
    "N", "D", "k1", "k2", "q", "ud", "ur", "f0", "use_ramp",
    "imported_cases", "imported_window"
  )))
  d$x_r[names(d$x_r) == "imported_cases"] <- imported_cases
  d$x_r[names(d$x_r) == "imported_window"] <- imported_window
  d$n_x_r <- length(d$x_r)

  initf_project <- function(post, i, stan_data) {
    R0 <- post$R0[i]
    i0 <- post$i0[i]
    ur <- post$ur[i]
    if ("beta[1]" %in% names(post)) {
      beta_n <- length(grep("beta\\[", names(post)))
      .beta <- rep(0, beta_n)
      for (k in seq_len(beta_n)) {
        .beta[k] <- post[[paste0("beta[", k, "]")]][i]
      }
      beta <- array(.beta)
    } else {
      beta <- array(numeric(0))
    }
    if (!is.null(f_fixed)) {
      f_s <- array(c(post$f_s[i, ], f_fixed))
    }
    if (!is.null(f_multi)) {
      fs <- post$f_s[i, ]
      if (f_multi_seg > length(fs)) {
        stop("`f_multi_seg` to large.", call. = FALSE)
      }
      proj_f <- fs[f_multi_seg]
      f_s <- array(c(post$f_s[i, ], f_multi * proj_f))
      if (max(f_s) >= 0.999) {
        warning("f_s >= 0.999! Setting to 0.999.", call. = FALSE)
        f_s <- ifelse(f_s >= 0.999, 0.999, f_s)
      }
    }
    if (is.null(f_fixed) && is.null(f_multi)) {
      f_s <- array(c(post$f_s[i, ], f_fixed)) # TODO: what is this?
    }
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
    start_decline <- post$start_decline[i]
    end_decline <- post$end_decline[i]
    list(
      R0 = R0, i0 = i0, f_s = f_s, ur = ur, phi = phi, samp_frac = samp_frac,
      start_decline = start_decline, end_decline = end_decline, beta = beta
    )
  }

  pars <- c("R0", "i0", "f_s", "e", "ur", "phi", "mu", "y_rep")
  if (return_states) pars <- c("y_hat")

  proj_func <- function(i) {
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
  }
  if (parallel) {
    out <- furrr::future_map_dfr(iter, proj_func, .options = furrr::furrr_options(seed = TRUE))
  } else {
    out <- purrr::map_dfr(iter, proj_func)
  }

  if (return_states) {
    states <- c("S", "E1", "E2", "I", "Q", "R", "Sd", "E1d", "E2d", "Id", "Qd", "Rd")
    variables_df <- dplyr::tibble(
      variable = states,
      variable_num = seq_along(states)
    )
    ts_df <- dplyr::tibble(time = d$time, time_num = seq_along(d$time))
    out <- dplyr::rename(out, time_num = Var2, variable_num = Var3)
    out <- dplyr::left_join(out, variables_df, by = "variable_num")
    out <- dplyr::left_join(out, ts_df, by = "time_num")
    out$variable_num <- NULL
    out$time_num <- NULL
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
