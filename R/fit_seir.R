#' Fit a Stan SEIR model
#'
#' @param daily_cases Either a vector of daily new cases if fit into a single
#'   data type or a matrix of case data is fitting to multiple data types. Each
#'   data type should be in its own column. Can have NA values (will be ignored
#'   in likelihood). A vector will be turned into a one column matrix.
#' @param obs_model Type of observation model
#' @param forecast_days Number of days into the future to forecast. The model
#'   will run faster with fewer forecasted days.
#' @param time_increment Time increment for ODEs and Weibull delay-model
#'   integration
#' @param days_back Number of days to go back for Weibull delay-model
#'   integration
#' @param R0_prior Lognormal log mean and SD for R0 prior
#' @param phi_prior SD of `1/sqrt(phi) ~ Normal(0, SD)` prior, where NB2(mu,
#'   phi) and `Var(Y) = mu + mu^2 / phi`.
#'   <https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations>
#' @param f2_prior Beta mean and SD for `f2` parameter
#' @param samp_frac_prior `samp_frac` prior if `samp_frac_type` is "estimated"
#'   or "rw" or "segmented". In the case of the random walk, this specifies the
#'   initial state prior. The two values correspond to the mean and SD of a Beta
#'   distribution. Only applies to first time series. (Currently disabled!)
#' @param samp_frac_type How to treat the sample fraction. Fixed, estimated, or
#'   a constrained random walk. Also segmented. Only applies to the first data
#'   type column if there are multiple data types. The other data types always
#'   have a fixed sample fraction.
#' @param samp_frac_seg A vector of sample fraction segment indexes of length
#'   `daily_cases` + `forecast_days`. Should start at 1.
#' @param rw_sigma The standard deviation on the optional `samp_frac` random
#'   walk on the first data type.
#' @param f_seg A vector of segment indexes of length `daily_cases` +
#'   `forecast_days`. The the segment index values should start at 0 to
#'   represent the fixed "no social distancing" value `f1` from the `pars`
#'   argument.
#' @param seed MCMC seed
#' @param chains Number of MCMC chains
#' @param iter MCMC iterations per chain
#' @param samp_frac_fixed A vector (or matrix) of sampled fractions. Should be
#'   of dimensions: `nrow(daily_cases) + forecast_days` (rows) by
#'   `ncol(daily_cases` (columns). A vector will be turned into a one column
#'   matrix.
#' @param fixed_f_forecast Optional fixed `f` (fraction of normal distancing)
#'   for the forecast.
#' @param day_start_fixed_f_forecast Day to start using `fixed_f_forecast`.
#' @param pars A named numeric vector of fixed parameter values
#' @param i0 Infected people infected at initial point in time.
#' @param fsi Fraction socially distancing. Derived parameter.
#' @param nsi Fraction not socially distancing. Derived parameter.
#' @param state_0 Initial state: a named numeric vector.
#' @param save_state_predictions Include the state predictions? `y_hat` Will
#'   make the resulting model object much larger.
#' @param delay_scale Weibull scale parameter for the delay in reporting. If
#'   there are multiple datatypes then this should be a vector of the same
#'   length as the data types.
#' @param delay_shape Weibull shape parameter for the delay in reporting. Same
#'   format as for `delay_scale`.
#' @param ode_control Control options for the Stan ODE solver. First is relative
#'   difference, that absolute difference, and then maximum iterations. You
#'   probably don't need to touch these.
#' @param ... Other arguments to pass to [rstan::sampling()].
#' @export
#' @return A named list object
#' @examples
#' # Example daily case data:
#' cases <- c(
#' 0, 0, 1, 3, 1, 8, 0, 6, 5, 0, 7, 7, 18, 9, 22, 38, 53, 45, 40,
#' 77, 76, 48, 67, 78, 42, 66, 67, 92, 16, 70, 43, 53, 55, 53, 29,
#' 26, 37, 25, 45, 34, 40, 35
#' )
#'
#' # Example assume sampling fractions of positive cases:
#' s1 <- c(rep(0.1, 13), rep(0.2, length(cases) - 13))
#'
#' # Using only 100 iterations and 1 chain for a quick example:
#' m <- fit_seir(
#'   cases,
#'   iter = 100,
#'   chains = 1,
#'   samp_frac_fixed = s1
#' )
#' print(m)
#' names(m)
#' names(m$post)
#' post_mu <- m$fit %>% tidybayes::spread_draws(mu[day, data_type])
#' post_mu
#' post_y_rep <- m$fit %>% tidybayes::spread_draws(mu[day, data_type])
#' post_y_rep
#'
#' # Add hospitalization data and estimate 2 sample-fraction blocks
#' # for the reported cases:
#' samp_frac_seg <- c(rep(1, 13), rep(2, length(cases) - 13))
#'
#' s2 <- rep(0.7, 42) # Assuming 7% of positive individuals are hospitalized
#' hosp <- c(
#'   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 11, 9, 4,
#'   7, 19, 23, 13, 11, 3, 13, 21, 14, 17, 29, 21, 19, 19, 10, 6,
#'   11, 8, 11, 8, 7, 6
#' )
#' \donttest{
#' m2 <- fit_seir(
#'   daily_cases = cbind(cases, hosp),
#'   iter = 100,
#'   chains = 1,
#'   samp_frac_type = "segmented", # `samp_frac_type` affects the first data type only
#'   samp_frac_seg = samp_frac_seg,
#'   samp_frac_fixed = cbind(s1, s2), # s1 is ignored and could be anything
#'   delay_scale = c(9.8, 11.2),
#'   delay_shape = c(1.7, 1.9),
#' )
#' print(m2)
#' }
#'
#' # Estimate a second estimated f_s block:
#' f_seg <- c(rep(0, 14), rep(1, 20), rep(2, length(cases) - 20 - 14))
#' m3 <- fit_seir(
#'   cases,
#'   iter = 100,
#'   chains = 1,
#'   f_seg = f_seg,
#'   samp_frac_fixed = s1
#' )
#' print(m3)

fit_seir <- function(daily_cases,
                     obs_model = c("NB2", "Poisson"),
                     forecast_days = 0,
                     time_increment = 0.1,
                     days_back = 45,
                     R0_prior = c(log(2.6), 0.2),
                     phi_prior = 1,
                     f2_prior = c(0.4, 0.2),
                     samp_frac_prior = c(0.4, 0.2),
                     samp_frac_type = c("fixed", "estimated", "rw", "segmented"),
                     samp_frac_seg = NULL,
                     rw_sigma = 0.1,
                     f_seg = c(rep(0, 14), rep(1, nrow(daily_cases) + forecast_days - 14)),
                     seed = 42,
                     chains = 4,
                     iter = 2000,
                     samp_frac_fixed = NULL,
                     fixed_f_forecast = NULL,
                     day_start_fixed_f_forecast = nrow(daily_cases) + 1,
                     pars = c(
                       N = 5.1e6, D = 5, k1 = 1 / 5,
                       k2 = 1, q = 0.05,
                       r = 0.1, ur = 0.02, f1 = 1.0,
                       start_decline = 15,
                       end_decline = 22
                     ),
                     i0 = 8,
                     fsi = pars[["r"]] / (pars[["r"]] + pars[["ur"]]),
                     nsi = 1 - fsi,
                     state_0 = c(
                       S = nsi * (pars[["N"]] - i0),
                       E1 = 0.4 * nsi * i0,
                       E2 = 0.1 * nsi * i0,
                       I = 0.5 * nsi * i0,
                       Q = 0,
                       R = 0,
                       Sd = fsi * (pars[["N"]] - i0),
                       E1d = 0.4 * fsi * i0,
                       E2d = 0.1 * fsi * i0,
                       Id = 0.5 * fsi * i0,
                       Qd = 0,
                       Rd = 0
                     ),
                     save_state_predictions = FALSE,
                     delay_scale = 9.85,
                     delay_shape = 1.73,
                     ode_control = c(1e-6, 1e-5, 1e5),
                     ...) {
  obs_model <- match.arg(obs_model)
  obs_model <-
    if (obs_model == "Poisson") {
      0L
    } else if (obs_model == "NB2") {
      1L
    }
  x_r <- pars

  samp_frac_type <- match.arg(samp_frac_type)
  n_samp_frac <-
    if (samp_frac_type == "fixed") {
      0L
    } else if (samp_frac_type == "estimated") {
      1L
    } else if (samp_frac_type == "segmented") {
      max(samp_frac_seg)
    } else { # random walk:
      nrow(daily_cases)
    }

  stopifnot(
    names(x_r) ==
      c("N", "D", "k1", "k2", "q", "r", "ur", "f1", "start_decline", "end_decline")
  )
  stopifnot(
    names(state_0) == c("S", "E1", "E2", "I", "Q", "R", "Sd", "E1d", "E2d", "Id", "Qd", "Rd")
  )

  # Checks and type conversions:
  if (!is.matrix(daily_cases)) daily_cases <- matrix(daily_cases, ncol = 1)
  if (!is.matrix(samp_frac_prior)) samp_frac_prior <- matrix(samp_frac_prior, ncol = 1)
  if (!is.matrix(samp_frac_fixed)) samp_frac_fixed <- matrix(samp_frac_fixed, ncol = 1)
  stopifnot(length(delay_scale) == ncol(daily_cases))
  stopifnot(length(delay_shape) == ncol(daily_cases))

  days <- seq(1, nrow(daily_cases) + forecast_days)
  last_day_obs <- nrow(daily_cases)
  time <- seq(-30, max(days), time_increment)
  x_r <- c(x_r, if (!is.null(fixed_f_forecast)) fixed_f_forecast else 0)
  names(x_r)[length(x_r)] <- "fixed_f_forecast"
  x_r <- c(x_r, c("last_day_obs" = last_day_obs))
  x_r <- c(x_r, c("day_start_fixed_f_forecast" = day_start_fixed_f_forecast))

  # find the equivalent time of each day (end):
  time_day_id <- vapply(days, get_time_id, numeric(1), time = time)
  # find the equivalent time of each day (start):
  time_day_id0 <- vapply(days, get_time_day_id0, numeric(1),
    time = time, days_back = days_back
  )

  stopifnot(nrow(samp_frac_fixed) == length(days))
  stopifnot(ncol(samp_frac_fixed) == ncol(daily_cases))

  beta_sd <- f2_prior[2]
  beta_mean <- f2_prior[1]
  beta_shape1 <- get_beta_params(beta_mean, beta_sd)$alpha
  beta_shape2 <- get_beta_params(beta_mean, beta_sd)$beta

  if (samp_frac_type == "fixed") {
    samp_frac_prior <- c(1, 1) # fake
  }
  # samp_frac_prior_trans <- matrix(nrow = 2, ncol = ncol(daily_cases))
  # for (i in seq_len(ncol(daily_cases))) {
  #   samp_frac_prior_trans[1,i] <-
  #     get_beta_params(samp_frac_prior[1,i], samp_frac_prior[2,i])$alpha
  #   samp_frac_prior_trans[2,i] <-
  #     get_beta_params(samp_frac_prior[1,i], samp_frac_prior[2,i])$beta
  # }

  samp_frac_prior_trans <- c(
    get_beta_params(samp_frac_prior[1], samp_frac_prior[2])$alpha,
    get_beta_params(samp_frac_prior[1], samp_frac_prior[2])$beta
  )

  if (9999999L %in% daily_cases) {
    stop("covidseir uses `9999999` as a 'magic' number for `NA`.", call. = FALSE)
  }

  daily_cases_stan <- daily_cases
  daily_cases_stan[is.na(daily_cases_stan)] <- 9999999L # magic number for NA

  if (is.null(samp_frac_seg)) {
    samp_frac_seg <- rep(1, length(days))
  }

  stan_data <- list(
    T = length(time),
    days = days,
    daily_cases = daily_cases_stan,
    J = ncol(daily_cases),
    N = length(days),
    S = length(unique(f_seg)) - 1, # - 1 because of 0 for fixed f1 before soc. dist.
    y0 = state_0,
    t0 = min(time) - 0.000001,
    time = time,
    x_r = x_r,
    n_x_i = length(f_seg) + 1L,
    x_i = c(length(f_seg), f_seg),
    delay_shape = array(delay_shape),
    delay_scale = array(delay_scale),
    samp_frac_fixed = samp_frac_fixed,
    samp_frac_seg = samp_frac_seg,
    samp_frac_type = if (samp_frac_type == "segmented") 4L else 1L, # FIXME
    time_day_id = time_day_id,
    time_day_id0 = time_day_id0,
    R0_prior = R0_prior,
    phi_prior = phi_prior,
    f2_prior = c(beta_shape1, beta_shape2),
    samp_frac_prior = samp_frac_prior_trans,
    n_samp_frac = n_samp_frac,
    rw_sigma = rw_sigma,
    priors_only = 0L,
    last_day_obs = last_day_obs,
    obs_model = obs_model,
    ode_control = ode_control,
    est_phi = if (obs_model %in% 1L) ncol(daily_cases) else 0L
  )
  initf <- function(stan_data) {
    R0 <- stats::rlnorm(1, log(R0_prior[1]), R0_prior[2])
    f <- stats::rbeta(
      1,
      get_beta_params(f2_prior[1], f2_prior[2])$alpha,
      get_beta_params(f2_prior[1], f2_prior[2])$beta
    )
    f_s <- array(f, dim = stan_data$S)
    init <- list(R0 = R0, f_s = f_s)
    init
  }
  pars_save <- c("R0", "f_s", "phi", "mu", "y_rep", "samp_frac")
  if (save_state_predictions) pars_save <- c(pars_save, "y_hat")
  fit <- rstan::sampling(
    stanmodels$seir,
    data = stan_data,
    iter = iter,
    chains = chains,
    init = function() initf(stan_data),
    seed = seed, # https://xkcd.com/221/
    pars = pars_save,
    ... = ...
  )
  post <- rstan::extract(fit)
  structure(list(
    fit = fit, post = post, phi_prior = phi_prior, R0_prior = R0_prior,
    f2_prior = f2_prior, obs_model = obs_model,
    samp_frac_fixed = samp_frac_fixed, state_0 = state_0,
    daily_cases = daily_cases, days = days, time = time,
    last_day_obs = last_day_obs, pars = x_r, f2_prior_beta_shape1 = beta_shape1,
    f2_prior_beta_shape2 = beta_shape2, stan_data = stan_data, days_back = days_back
  ), class = "covidseir")
}

get_beta_params <- function(mu, sd) {
  var <- sd^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  list(alpha = alpha, beta = beta)
}

get_time_id <- function(day, time) max(which(time <= day))

get_time_day_id0 <- function(day, time, days_back) {
  # go back `days_back` or to beginning if that's negative time:
  check <- time <= (day - days_back)
  if (sum(check) == 0L) {
    1L
  } else {
    max(which(check))
  }
}

getu <- function(f, r) (r - f * r) / f
