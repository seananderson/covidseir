#' Get the contact rate fraction threshold for increases
#'
#' @param obj Output from [fit_seir()].
#' @param iter Vector of MCMC iterations to work with. More iterations will be
#'   slower and not necessarily render a different result. ~200 total iterations
#'   may be plenty.
#' @param forecast_days Days to use in forecast.
#' @param fs Contact rate fractions to test.
#' @param show_plot Make a diagnostic plot?
#' @param window_check The window of days to use from the last day forecasted.
#' @param ... Other arguments for [project_seir()].
#'
#' @return
#' The threshold value.
#' @export
#' @importFrom tibble tibble
#' @importFrom stats lm predict coef
#' @examples
#' # See ?project_seir
get_threshold <- function(obj, iter = seq_along(obj$post$R0),
                          forecast_days = 25,
                          fs = seq(0.3, 0.8, length.out = 4),
                          show_plot = TRUE,
                          window_check = 25,
                          ...) {

  m_fs <- purrr::map(fs, function(.f) {
    cat("Projecting", round(.f, 2), "\n")
    project_seir(obj,
      forecast_days = forecast_days,
      iter = iter,
      f_fixed_start = nrow(obj$daily_cases) + 1,
      f_fixed = rep(.f, forecast_days),
      return_states = TRUE, ...
    )
  })
  slopes <- purrr::map2_df(m_fs, fs, function(x, y) {
    temp <- dplyr::filter(
      x, time > max(x$time) - window_check,
      variable %in% c("I", "Id")
    )
    temp <- dplyr::group_by(temp, .iteration, time)
    temp <- dplyr::summarize(temp,
      I = value[variable == "I"], Id = value[variable == "Id"],
      prevalence = I + Id, .groups = "drop_last"
    )
    iters <- dplyr::summarise(dplyr::group_by(temp, .iteration),
      iter = .iteration[[1]], .groups = "drop_last"
    )

    temp <- dplyr::group_by(temp, .iteration)
    temp <- dplyr::group_split(temp)
    temp <- purrr::map(temp, ~ lm(log(prevalence) ~ time, data = .x))
    temp <- purrr::map_df(temp, ~ tibble(slope = coef(.x)[[2]]))
    temp <- dplyr::mutate(temp, f = y)
    temp <- dplyr::ungroup(temp)
    dplyr::mutate(temp, .iteration = iters$iter)
  })
  if (show_plot) plot(slopes$f, slopes$slope)

  nd <- tibble(f = seq(0.01, 0.99, length.out = 2000))
  out <- dplyr::group_by(slopes, .iteration)
  out <- dplyr::group_split(out)
  out <- purrr::map(out, ~ lm(slope ~ f, data = .x))
  purrr::map_dbl(out, function(.x) {
    nd$predicted_slope <- predict(.x, newdata = nd)
    out2 <- dplyr::filter(nd, predicted_slope > 0)
    out2[1, "f", drop = TRUE]
  })
}

#' Get the doubling time
#'
#' @param obj Output from [fit_seir()].
#' @param iter Vector of MCMC iterations to work with.
#' @param forecast_days Days to use in forecast.
#' @param window_check The window of days to use from the last day forecasted.
#' @param show_plot Make a diagnostic plot? Each line is a posterior sample of
#'   the log prevalence time series. Each line should be well approximated by a
#'   linear model
#' @param ... Other arguments for [project_seir()].
#'
#' @details Assumes that prevalence is increasing; negative values imply halving
#'   times. Use `show_plot` if you want to check the time series.
#' @return A vector of doubling times across iterations.
#' @export
#' @importFrom graphics lines plot
#' @examples
#' # See ?project_seir
get_doubling_time <- function(obj, iter = seq_along(obj$post$R0),
                              forecast_days = 25, window_check = 25,
                              show_plot = TRUE, ...) {
  m_proj <- project_seir(obj,
    forecast_days = forecast_days,
    iter = iter,
    return_states = TRUE, ...
  )
  temp <- dplyr::filter(
    m_proj, time > max(m_proj$time) - window_check,
    variable %in% c("I", "Id")
  )
  temp <- dplyr::group_by(temp, .iteration, time)
  temp <- dplyr::summarize(temp,
    I = value[variable == "I"], Id = value[variable == "Id"],
    prevalence = I + Id, .groups = "drop_last"
  )
  if (show_plot) {
    plot(temp$time, log(temp$prevalence),
      type = "n", xlab = "Time",
      ylab = "log(prevalence)"
    )
    for (i in unique(temp$.iteration)) {
      lines(temp[temp$.iteration == i, "time", drop = TRUE],
        log(temp[temp$.iteration == i, "prevalence", drop = TRUE]),
        col = "#00000060"
      )
    }
  }
  temp <- dplyr::group_by(temp, .iteration)
  temp <- dplyr::group_split(temp)
  temp <- purrr::map(temp, ~ lm(log(prevalence) ~ time, data = .x))
  temp <- purrr::map_df(temp, ~ tibble(slope = coef(.x)[[2]]))
  log(2) / temp$slope
}

# get_doubling_ts <- function(obj, iter = seq_along(obj$post$R0), ...) {
#   m_proj <- project_seir(obj,
#     forecast_days = 0,
#     iter = iter,
#     return_states = TRUE, ...
#   )
#   temp <- dplyr::filter(
#     m_proj, time > 0,
#     variable %in% c("I", "Id")
#   )
#   temp <- dplyr::group_by(temp, .iteration, time)
#   temp <- dplyr::summarize(temp,
#     I = value[variable == "I"], Id = value[variable == "Id"],
#     prevalence = I + Id, .groups = "drop_last"
#   )
#
#   browser()
#
#
#
#
#
#   temp <- dplyr::group_by(temp, .iteration)
#   temp <- dplyr::group_split(temp)
#   temp <- purrr::map(temp, ~ lm(log(prevalence) ~ time, data = .x))
#   temp <- purrr::map_df(temp, ~ tibble(slope = coef(.x)[[2]]))
#   log(2) / temp$slope
# }
