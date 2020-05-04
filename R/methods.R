#' @export
#' @import methods
print.covidseir <- function(x, pars = c("R0", "phi", "f_s"), ...) {
  print(x$fit, pars = pars, ...)
}
