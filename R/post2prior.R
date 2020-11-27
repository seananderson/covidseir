#' Extract new priors from a previous fit
#'
#' This function takes a model fit with [fit_seir()] and converts the posterior
#' distribution into priors that can be fed into a second model. This allows
#' splitting a time series into multiple blocks.
#'
#' @param obj A model fit with [fit_seir()].
#' @param iter Iterations to use went creating projections to extract the state
#'   posterior.
#' @param time_slice When to take the state priors at. Defaults to the
#'   end of the fit.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom stats sd
#' @examples
#' # See the vignette 'Fitting case data in multiple blocks'
#' # for a complete example
#'
#' cases <- c(
#'   0, 0, 1, 3, 1, 8, 0, 6, 5, 0, 7, 7, 18, 9, 22, 38, 53, 45, 40,
#'   77, 76, 48, 67, 78, 42, 66, 67, 92, 16, 70, 43, 53, 55, 53, 29,
#'   26, 37, 25, 45, 34, 40, 35
#' )
#' s1 <- c(rep(0.1, 13), rep(0.2, length(cases) - 13))
#' m <- fit_seir(
#'   cases,
#'   iter = 200,
#'   fit_type = "optimizing",
#'   samp_frac_fixed = s1
#' )
#' post2prior(m, iter = seq_len(10)) # just 10 for example speed

post2prior <- function(obj, iter = seq_len(100), time_slice = max(proj$time)) {
  p <- obj$post
  R0 <- p$R0
  R0_hat <- unname(MASS::fitdistr(R0, "lognormal")$estimate)

  proj <- covidseir::project_seir(obj, iter = iter, return_states = TRUE)
  i30 <- dplyr::filter(
    proj,
    variable %in% c("I", "Id", "E1", "E2", "E1d", "E2d") & time == time_slice - 30
  ) %>%
    dplyr::group_by(.iteration) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    dplyr::pull(value)
  i30_hat <- unname(MASS::fitdistr(i30, "lognormal")$estimate)
  e_hat <- c(mean(p$e), stats::sd(p$e))
  start_decline_fake <- c(log(1), 0.1)
  end_decline_fake <- c(log(1), 0.1)
  phi_hat <- unname(MASS::fitdistr(p$phi, "lognormal")$estimate)
  f <- p$f_s[, ncol(p$f_s)]
  f_hat <- c(mean(f), sd(f))

  # state0 at time -30
  S0 <- proj %>%
    dplyr::filter(time == time_slice - 30) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(value = mean(value), .groups = "drop")
  state_order <- c("S", "E1", "E2", "I", "Q", "R", "Sd", "E1d", "E2d", "Id", "Qd", "Rd")
  S0_hat <- S0$value
  names(S0_hat) <- S0$variable
  S0_hat <- S0_hat[state_order]

  i0 <- mean(i30)
  state_0 <- c(
    E1_frac = S0_hat[["E1"]] / ((1 - e_hat[1]) * i0),
    E2_frac = S0_hat[["E2"]] / ((1 - e_hat[1]) * i0),
    I_frac = S0_hat[["I"]] / ((1 - e_hat[1]) * i0),
    Q_num = S0_hat[["Q"]],
    R_num = S0_hat[["R"]],
    E1d_frac = S0_hat[["E1d"]] / (e_hat[1] * i0),
    E2d_frac = S0_hat[["E2d"]] / (e_hat[1] * i0),
    Id_frac = S0_hat[["Id"]] / (e_hat[1] * i0),
    Qd_num = S0_hat[["Qd"]],
    Rd_num = S0_hat[["Rd"]]
  )

  list(
    R0 = R0_hat,
    i0 = i30_hat,
    e = e_hat,
    phi = phi_hat,
    f = f_hat,
    start_decline = start_decline_fake,
    end_decline = end_decline_fake,
    state_0 = state_0
  )
}
