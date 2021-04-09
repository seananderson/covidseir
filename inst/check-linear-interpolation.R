linear_interpolation <- function(x, x_pred, y_pred) {
  K <- length(x_pred)
  if (length(y_pred) != K) {
    stop("x_pred and y_pred aren't of the same size")
  }
  if (x <= min(x_pred)) {
    ans <- min(y_pred)
  } else if (x >= max(x_pred)) {
    ans <- max(y_pred)
  } else {
    i <- 0
    for (k in 1:K) {
      if (x_pred[k] - x <= 0)  i <- i + 1
    }
    x0 <- x_pred[i]
    x1 <- x_pred[i + 1]
    y0 <- y_pred[i]
    y1 <- y_pred[i + 1]
    # https://en.wikipedia.org/wiki/Linear_interpolation
    ans <- (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
  }
  ans
}

out_approx <- sapply(seq(1, 10, length.out = 100), function(x)
  approx(1:10, seq(1, 2, length.out = 10), xout = x)$y)
out_approx
plot(out_approx)

out_custom <- sapply(seq(1, 10, length.out = 100), function(x)
  linear_interpolation(x, 1:10, seq(1, 2, length.out = 10)))
out_custom <- unlist(out_custom)
points(out_custom, col = "red")

assertthat::are_equal(out_approx, out_custom)
