#' Create segments vector
#'
#' @param start_date date at start of data (in ymd format)
#' @param total_days int Number of days in data or projection. Gives total
#'   length of array
#' @param segments vector of dates (in ymd format)
#' @param values vector length of segments + 1
#' @returns vector of length total_days
#' @examples
#' start_date <- "2020-01-01"
#' segments <- c("2020-01-05", "2020-01-20")
#' total_days <- 30
#' values <- c(0, 0.5, 1)
#' create_segments_vector(start_date, segments, total_days, values = values)
#' @importFrom lubridate ymd days
#' @export
create_segments_vector <- function(start_date, total_days, segments, values = NULL) {
  if (!is.null(values) & (length(values) != (length(segments) + 1))) {
    stop("Length values must be length segments + 1.", call. = FALSE)
  }

  if (is.null(values)) {
    values <- seq(1, length(segments) + 1)
  }

  vec <- rep(1, total_days)

  for (i in seq_len(length(values))) {
    if (i == 1) {
      begin_segment <- ymd(start_date)
    } else {
      begin_segment <- ymd(segments[[i - 1]])
    }

    if (i == length(values)) {
      days_sofar <- as.numeric(ymd(segments[[i - 1]]) - ymd(start_date))
      end_segment <- ymd(segments[[i - 1]]) + days(total_days - days_sofar - 1)
    } else {
      end_segment <- ymd(segments[[i]])
    }

    # start of time point of segment
    s1 <- as.numeric(begin_segment - ymd(start_date)) + 1
    s2 <- as.numeric(end_segment - ymd(start_date)) + 1
    vec[s1:s2] <- values[i]
  }
  vec
}

#' Create ramp vector
#'
#' @description Creates a vector of length `total_days` starting at 1 and
#'   increasing to `ramp_max` between `start_ramp` and `end_ramp`. Useful for
#'   characterizing the change in transmission due to Variants of Concern.
#' @param start_date date at start of data (in ymd format)
#' @param total_days int Number of days in data or projection. Gives total
#'   length of array
#' @param start_ramp date at start of ramp (in ymd format)
#' @param ramp_length length of ramp (in days)
#' @param ramp_max final value at end of ramp
#' @returns vector of length total_days
#' @examples
#' start_date <- "2020-01-01"
#' start_ramp <- "2020-01-05"
#' total_days <- 30
#' ramp_length <- 10
#' ramp_max <- 2
#' create_ramp_vector(
#'   start_date, total_days, start_ramp, ramp_length,
#'   ramp_max
#' )
#' @importFrom lubridate ymd days
#' @export
create_ramp_vector <- function(start_date, total_days, start_ramp,
                               ramp_length, ramp_max) {
  vec <- rep(1, total_days)

  # start of ramp
  s1 <- as.numeric(ymd(start_ramp) - ymd(start_date)) + 1

  if (s1 > total_days) {
    warning("Ramp starts after `total_days`")
  } else {
    # end of ramp (ensure doesn't go past total_days)
    s2 <- min(s1 + ramp_length, total_days)

    vec[seq(s1, s2)] <- 1 + (ramp_max - 1) * (seq(0, s2 - s1)) / ramp_length

    # only include if total_days greathr than end of ramp
    if (total_days > s2) {
      vec[seq(s2 + 1, total_days)] <- ramp_max
    }
  }
  vec
}
