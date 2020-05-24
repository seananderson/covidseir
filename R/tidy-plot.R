#' Generate a 'tidy' data frame of summarized projection output
#'
#' Returns a variety of quantiles and means of the posterior predictive distribution and the expected case numbers.
#'
#' @param x Output from [project_seir()].
#' @param resample_y_rep The number of times to resample from the observation distribution. Currently only set up for NB2 models.
#' @param data_type_names An optional named character vector to translate the data type numbers into names. E.g. `c("Reported cases" = 1, "Hospitalizations" = 2)`.
#'
#' @details
#' See [project_seir()] for an example.
#'
#' @return A data frame
#' @export
tidy_seir <- function(x, resample_y_rep = 10, data_type_names = NULL) {
  # FIXME: check if NB2 or Poisson and stop or adjust
  if (resample_y_rep > 0) {
    x <- purrr::map_dfr(seq_len(resample_y_rep), function(i) {
      temp <- x
      temp$y_rep <- MASS::rnegbin(length(temp$y_rep),
        temp$mu,
        theta = temp$phi
      )
      temp
    })
  }
  out <- dplyr::group_by(x, data_type, day)
  out <- dplyr::summarise(out,
      y_rep_0.05 = stats::quantile(y_rep, probs = 0.05),
      y_rep_0.25 = stats::quantile(y_rep, probs = 0.25),
      y_rep_mean = mean(y_rep),
      y_rep_0.50 = stats::quantile(y_rep, probs = 0.50),
      y_rep_0.75 = stats::quantile(y_rep, probs = 0.75),
      y_rep_0.95 = stats::quantile(y_rep, probs = 0.95),
      mu_0.05 = stats::quantile(mu, probs = 0.05),
      mu_0.25 = stats::quantile(mu, probs = 0.25),
      mu_mean = mean(mu),
      mu_0.50 = stats::quantile(mu, probs = 0.50),
      mu_0.75 = stats::quantile(mu, probs = 0.75),
      mu_0.95 = stats::quantile(mu, probs = 0.95),
      mu_0.5 = stats::quantile(mu, probs = 0.50)
    )
  if (!is.null(data_type_names)) {
    out$data_type <- names(data_type_names[as.numeric(out$data_type)])
  }
  out
}

#' Plot a SEIR model projection
#'
#' @param pred_dat Output from [tidy_seir()].
#' @param obs_dat A data frame of observed data. Should have a
#'   column of `day` and `value` that matches `pred_dat`. If
#'   using multiple data types, should also have `data_type`.
#' @param col Colour for the line and ribbon.
#'
#' @details
#' See [project_seir()] for an example.
#'
#' @export
#' @return A ggplot object
#'
#' @importFrom ggplot2 geom_ribbon ggplot facet_wrap coord_cartesian
#'   ylab theme geom_line geom_point element_blank aes_string
#' @importFrom dplyr tibble
plot_projection <- function(pred_dat, obs_dat, col = "#377EB8") {
  if (!"value" %in% names(obs_dat)) {
    stop("`obs_dat` must contain a column named `value` that contains the case counts.", call. = FALSE)
  }
  if (!"day" %in% names(obs_dat)) {
    stop("`obs_dat` must contain a column named `day` that contains the numeric day (or date).", call. = FALSE)
  }
  g <- ggplot(pred_dat, aes_string(x = "day")) +
    geom_ribbon(aes_string(ymin = "y_rep_0.05", ymax = "y_rep_0.95"),
      alpha = 0.2, fill = col) +
    geom_ribbon(aes_string(ymin = "y_rep_0.25", ymax = "y_rep_0.75"),
      alpha = 0.2, fill = col) +
    geom_line(aes_string(y = "mu_0.50"), lwd = 0.9, col = col) +
    facet_wrap(~data_type) +
    coord_cartesian(expand = FALSE, xlim = range(pred_dat$day)) +
    ylab("Cases") +
    theme(axis.title.x = element_blank())
  g <- g +
    geom_line(
      data = obs_dat,
      col = "black", inherit.aes = FALSE,
      aes_string(x = "day", y = "value"),
      lwd = 0.35, alpha = 0.9
    ) +
    geom_point(
      data = obs_dat,
      col = "grey30", inherit.aes = FALSE,
      aes_string(x = "day", y = "value"),
      pch = 21, fill = "grey95", size = 1.25
    )
  g
}
