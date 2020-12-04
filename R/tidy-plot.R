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
#' @importFrom stats quantile
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
    y_rep_0.05 = quantile(y_rep, probs = 0.05),
    y_rep_0.25 = quantile(y_rep, probs = 0.25),
    y_rep_mean = mean(y_rep),
    y_rep_0.50 = quantile(y_rep, probs = 0.50),
    y_rep_0.75 = quantile(y_rep, probs = 0.75),
    y_rep_0.95 = quantile(y_rep, probs = 0.95),
    mu_0.05 = quantile(mu, probs = 0.05),
    mu_0.25 = quantile(mu, probs = 0.25),
    mu_mean = mean(mu),
    mu_0.50 = quantile(mu, probs = 0.50),
    mu_0.75 = quantile(mu, probs = 0.75),
    mu_0.95 = quantile(mu, probs = 0.95),
    mu_0.5 = quantile(mu, probs = 0.50),
    .groups = "drop_last"
  )
  if (!is.null(data_type_names)) {
    out$data_type <- names(data_type_names[as.numeric(out$data_type)])
  }
  out
}

#' Plot a SEIR model projection
#'
#' @param pred_dat Output from [tidy_seir()].
#' @param obs_dat A data frame of observed data. Should have a column of `day`
#'   and `value` that matches `pred_dat`. If using multiple data types, should
#'   also have `data_type`.
#' @param col Colour for the line and ribbon.
#' @param value_column Column in `obs_dat` that contains the reported cases.
#' @param date_column Date or day column name. A column by this name must be in
#'   `obs_dat` and `pred_dat`. Note that [tidy_seir()] will return a data frame
#'   with `day` in it. If you want dates then you will need to add such a
#'   column.
#' @param ylab Y axis label.
#' @param Rt An optional vector of Rt values to colour the plot by.
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
#' @importFrom glue glue
#' @examples
#' # See ?project_seir
plot_projection <- function(pred_dat, obs_dat, col = "#377EB8",
                            value_column = "value", date_column = "day",
                            ylab = "Reported cases", Rt = NULL) {
  if (!value_column %in% names(obs_dat)) {
    stop(glue("`obs_dat` must contain a column `{value_column}` that contains the reported case counts."), call. = FALSE)
  }
  if (!date_column %in% names(obs_dat)) {
    stop(glue("`obs_dat` must contain a column named `{date_column}` that contains the numeric day (or date)."), call. = FALSE)
  }
  if (!date_column %in% names(pred_dat)) {
    stop(glue("`pred_dat` must contain a column named `{date_column}` that contains the numeric day (or date)."), call. = FALSE)
  }

  if (!is.null(Rt)) {
    pred_dat$Rt <- Rt
    make_poly <- function(df, lwr = "y_rep_0.05", upr = "y_rep_0.95") {
      .l <- list()
      for (i in seq(1, nrow(pred_dat) - 1)) {
        .l[[i]] <- dplyr::bind_rows(
          tibble(x = pred_dat[[date_column]][i], y = pred_dat[[lwr]][i]),
          tibble(x = pred_dat[[date_column]][i], y = pred_dat[[upr]][i]),
          tibble(x = pred_dat[[date_column]][i + 1], y = pred_dat[[upr]][i + 1]),
          tibble(x = pred_dat[[date_column]][i + 1], y = pred_dat[[lwr]][i + 1]),
        )
        .l[[i]][["Rt"]] <- pred_dat[["Rt"]][i]
        .l[[i]][["group"]] <- i
      }
      dplyr::bind_rows(.l)
    }
    ci1 <- make_poly(pred_dat, lwr = "y_rep_0.05", upr = "y_rep_0.95")
    ci2 <- make_poly(pred_dat, lwr = "y_rep_0.25", upr = "y_rep_0.75")
    pal <- rev(c("#EF8A62", "#F7F7F7", "#67A9CF")) # RColorBrewer::brewer.pal(3, "RdBu")
  }

  g <- ggplot(pred_dat, aes_string(x = date_column))

  if (!is.null(Rt)) {
    g <- g + ggplot2::geom_polygon(aes_string(x = "x", y = "y", fill = "Rt", group = "group"),
      data = ci1, alpha = 0.55
    ) +
      ggplot2::geom_polygon(aes_string(x = "x", y = "y", fill = "Rt", group = "group"),
        data = ci2, alpha = 0.9
      ) +
      ggplot2::scale_fill_gradient2(
        midpoint = 0, high = pal[3], mid = pal[2],
        low = pal[1], trans = "log10"
      )
  } else {
    g <- g + geom_ribbon(aes_string(ymin = "y_rep_0.05", ymax = "y_rep_0.95"),
      alpha = 0.2, , fill = col
    ) +
      geom_ribbon(aes_string(ymin = "y_rep_0.25", ymax = "y_rep_0.75"),
        alpha = 0.2, fill = col
      )
  }
  g <- g + geom_line(aes_string(y = "mu_0.50"), lwd = 0.9, col = col) +
    coord_cartesian(expand = FALSE, xlim = range(pred_dat[[date_column]])) +
    ylab(ylab) +
    theme(axis.title.x = element_blank())
  g <- g +
    geom_line(
      data = obs_dat,
      col = "black", inherit.aes = FALSE,
      aes_string(x = date_column, y = value_column),
      lwd = 0.35, alpha = 0.9
    ) +
    geom_point(
      data = obs_dat,
      col = "grey30", inherit.aes = FALSE,
      aes_string(x = date_column, y = value_column),
      pch = 21, fill = "grey95", size = 1.25
    )

  if (max(pred_dat[["data_type"]]) > 1) g <- g + facet_wrap(~data_type)
  g
}

#' @rdname plot_projection
#' @param obj Outut from [fit_seir()].
#' @param date_function A function to translate the character representation of
#'   the date if `date_column` is a date.
#' @param type Raw (observed - expected) or quantile residuals (`covidseir:::qres_nbinom2`)?
#' @param date_function A function to translate the character representation of
#' @param return_residuals Return residuals instead of plot?
#' @export

plot_residuals <- function(pred_dat, obs_dat, obj,
                           value_column = "value", date_column = "day",
                           type = c("raw", "quantile"),
                           ylab = if (type == "raw") "Residual" else "Quantile residual",
                           date_function = lubridate::ymd,
                           return_residuals = FALSE) {
  type <- match.arg(type)
  temp <- obs_dat
  temp$mu_0.50 <- pred_dat$mu_0.50
  if (type == "raw") {
    temp$resid <- temp[[value_column]] - pred_dat$mu_0.50
  } else {
    if (obj$stan_data$est_phi) {
      temp$resid <- qres_nbinom2(
        temp[[value_column]], pred_dat$mu_0.50,
        stats::median(obj$post$phi)
      )
    } else {
      temp$resid <- qres_pois(temp[[value_column]], pred_dat$mu_0.50)
    }
  }
  f_breaks <- as.numeric(obj$stan_data$x_i[grep("^f_seg*", names(obj$stan_data$x_i))])
  f_breaks <- obs_dat[[date_column]][diff(f_breaks) == 1][-1] + 1
  f_breaks <- f_breaks[-length(f_breaks)]
  .s <- min(obs_dat[[date_column]]) + stats::median(obj$post$start_decline)
  .e <- min(obs_dat[[date_column]]) + stats::median(obj$post$end_decline)
  f_breaks <- c(f_breaks, .s, .e)
  if (identical(class(obs_dat[[date_column]]), "Date")) {
    f_breaks <- date_function(f_breaks)
  }
  if (!return_residuals) {
    g <- ggplot(temp, aes_string(date_column, "resid")) +
      geom_point() +
      ggplot2::geom_smooth(
        se = TRUE, col = "red", method = "loess",
        formula = "y ~ x"
      ) +
      ylab(ylab) +
      theme(axis.title.x = element_blank())
    if (identical(class(obs_dat[[date_column]]), "Date")) {
      g <- g + ggplot2::scale_x_date(date_breaks = "1 month", date_labels = "%b")
    }
    for (i in seq_along(f_breaks)) {
      g <- g + ggplot2::geom_vline(xintercept = f_breaks, lty = 2, col = "grey50")
    }
    g
  } else {
    temp$resid
  }
}

qres_nbinom2 <- function(y, mu, phi) {
  a <- stats::pnbinom(y - 1, size = phi, mu = mu)
  b <- stats::pnbinom(y, size = phi, mu = mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}

qres_pois <- function(y, mu) {
  a <- stats::ppois(y - 1, mu)
  b <- stats::ppois(y, mu)
  u <- stats::runif(n = length(y), min = a, max = b)
  stats::qnorm(u)
}
