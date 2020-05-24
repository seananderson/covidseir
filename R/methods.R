#' @export
#' @import methods
print.covidseir <- function(x,
  pars = c("R0", "i0", "e", "f_s", "start_decline", "end_decline", "phi"), ...) {
  print(x$fit, pars = pars, ...)
}
