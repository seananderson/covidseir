#' Make a projection plot
#'
#' @param models A list of models from [fit_seir()]. If there are multiple elements than the names of the list elements will be used for the facet labels.
#' @param cumulative Logical for whether or not it should be a cumulative plot.
#' @param first_date The first date in the plot.
#' @param ylim The y-axis limits.
#' @param outer_quantile A vector representing the lower and upper outer quantiles of the credible interval.
#' @param facet Logical for whether or not to facet the panels. If false, then will put multiple ribbons on the same axes.
#' @param linetype Whether or not the line should represent the mean of the expectation or the median of the observation distribution.
#' @param omitted_days An optional vector of days to omit from the plot.
#' @param y_rep_dat An optional posterior predictive replicate data frame to override the input from `models`.
#' @param mu_dat An optional expected value data frame to override the input from `models`.
#' @param points_size The point size of dots.
#' @param sc_order An optional vector of scenario names to force the order of the panels. This character vector should match the named elements of `models`.
#'
#' @return A ggplot
#' @export
#' @importFrom ggplot2 theme facet_wrap guides ylab xlab scale_color_manual
#' scale_fill_manual geom_point labs aes_string ggplot annotate geom_ribbon coord_cartesian
#' element_blank geom_line facet_grid
#' @importFrom dplyr summarise group_by mutate rename filter summarize ungroup tibble as_tibble
#' @importFrom deSolve ode
make_projection_plot <- function(models, cumulative = FALSE,
  first_date = "2020-03-01", ylim = c(0, max(out$upr) * 1.03), outer_quantile = c(0.05, 0.95),
  facet = TRUE, cols = NULL, linetype = c("mu", "obs"),
  omitted_days = NULL, y_rep_dat = NULL, mu_dat = NULL, points_size = 1.25,
  sc_order = NULL
  ) {

  linetype <- match.arg(linetype)
  obj <- models[[1]]
  actual_dates <- seq(lubridate::ymd(first_date),
    lubridate::ymd(first_date) + max(obj$days), by = "1 day")

  if (is.null(y_rep_dat) || is.null(mu_dat)) {
    out <- purrr::map_df(models, function(.x) {

      temp <- tidybayes::spread_draws(obj$fit, y_rep[day, data_type])
      if (cumulative) {
        temp <- temp %>%
          group_by(iterations, data_type) %>%
          mutate(y_rep = cumsum(y_rep)) %>%
          ungroup()
      }

      temp %>%
        group_by(day, data_type) %>%
        summarise(
          lwr = stats::quantile(y_rep , probs = outer_quantile[1]),
          lwr2 = stats::quantile(y_rep , probs = 0.25),
          upr = stats::quantile(y_rep , probs = outer_quantile[2]),
          upr2 = stats::quantile(y_rep , probs = 0.75),
          med = stats::median(y_rep )
        ) %>%
        ungroup() %>%
        mutate(day = actual_dates[day])
    }, .id = "Scenario")

    lambdas <- purrr::map_df(models, function(.x) {

      temp <- tidybayes::spread_draws(obj$fit, lambda_d[day, data_type])
      if (cumulative) {
        temp <- temp %>%
          group_by(iterations, data_type) %>%
          mutate(lambda_d = cumsum(lambda_d)) %>%
          ungroup()
      }

      temp %>%
        group_by(day, data_type) %>%
        summarise(
          med = median(lambda_d)
        ) %>%
        ungroup() %>%
        mutate(day = actual_dates[day])
    }, .id = "Scenario")
  } else {
    out <- y_rep_dat
    lambdas <- mu_dat
  }

  .d <- as.data.frame(obj$daily_cases)
  names(.d) <- seq(1, ncol(.d))
  .d$day <- actual_dates[1:obj$last_day_obs]
  dat <- tidyr::pivot_longer(.d, -day, names_to = "data_type")

  if (cumulative) {
    dat <- dat %>% group_by(data_type) %>%
      mutate(value = cumsum(value)) %>%
      ungroup()
  }
  if (is.null(cols)) {
    cols <- rep("#3182BD", 99)
  }

  if (!is.null(sc_order)) {
    out$Scenario <- factor(out$Scenario, levels = sc_order)
    lambdas$Scenario <- factor(lambdas$Scenario, levels = sc_order)
  }
  g <- ggplot(out, aes_string(x = "day", y = "med", ymin = "lwr", ymax = "upr", colour = "Scenario",
    fill = "Scenario")) +
    annotate("rect", xmin = actual_dates[obj$last_day_obs], xmax = max(out$day), ymin = 0, ymax = ylim[2], fill = "grey95") +
    coord_cartesian(expand = FALSE, ylim = ylim, xlim = range(out$day)) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    geom_ribbon(alpha = 0.2, mapping = aes_string(ymin = "lwr2", ymax = "upr2"), colour = NA)

  if (linetype == "obs")
    g <- g + geom_line(alpha = 1, lwd = 1)
  if (linetype == "mu")
    g <- g + geom_line(data = lambdas, aes_string(x = "day", y = "med", colour = "Scenario"), alpha = 1, lwd = 1, inherit.aes = FALSE)

  g <- g +
    geom_line(
      data = dat,
      col = "black", inherit.aes = FALSE, aes_string(x = "day", y = "value"), lwd = 0.35,
      alpha = 0.9
    )

  if (!is.null(omitted_days)) {
    g <- g +
      geom_point(
        data = dat[omitted_days,,drop=FALSE],
        col = "grey30", inherit.aes = FALSE, aes_string(x = "day", y = "value"), pch = 4, fill = "grey95",  size = points_size
      ) +
      geom_point(
        data = dat[-omitted_days, ,drop=FALSE],
        col = "grey30", inherit.aes = FALSE, aes_string(x = "day", y = "value"), pch = 21, fill = "grey95",  size = points_size
      )
  } else {
    g <- g + geom_point(
      data = dat,
      col = "grey30", inherit.aes = FALSE, aes_string(x = "day", y = "value"), pch = 21, fill = "grey95", size = points_size
    )
  }

  g <- g +
    ylab(if (!cumulative) "Reported cases" else "Cumulative reported cases") +
    xlab("Day") +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(colour = "Projection scenario", fill = "Projection scenario") +
    theme(axis.title.x = element_blank())

  if (facet && (length(unique(out$Scenario)) > 1 || length(unique(out$data_type)) > 1))
    g <- g + facet_grid(Scenario~data_type) + theme(legend.position = "none")

  if (length(unique(out$Scenario)) == 1) {
    g <- g + guides(fill = FALSE, colour = FALSE)
  }
  g
}
