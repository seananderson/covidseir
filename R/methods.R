#' @export
#' @import methods
print.covidseir <- function(x, pars = c("R0", "phi", "f_s", "start_decline", "end_decline"), ...) {
  print(x$fit, pars = pars, ...)
}
