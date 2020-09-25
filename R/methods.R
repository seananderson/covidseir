#' @export
#' @import methods
print.covidseir <- function(x,
  pars = c("R0", "i0", "e", "f_s", "start_decline", "end_decline", "phi", "beta"), ...) {
  if ("fit_type" %in% names(x)) {
    if (x$fit_type == "optimizing") {
      beta_n <- grep("beta\\[", colnames(x$fit$theta_tilde))
      phi_n <- grep("phi\\[", colnames(x$fit$theta_tilde))
      e_n <- grep("^e$", colnames(x$fit$theta_tilde))
      pars_n <- c(seq_len(phi_n), e_n, beta_n)
      cat("MAP estimate:\n")
      print(round(x$fit$par[pars_n], 2), ...)
      cat("Mean in constrained space of MVN samples:\n")
      print(apply(x$fit$theta_tilde[,pars_n], 2, function(y) round(mean(y), 2)), ...)
      cat("SD in constrained space of MVN samples:\n")
      print(apply(x$fit$theta_tilde[,pars_n], 2, function(y) round(stats::sd(y), 2)), ...)
    } else {
      print(x$fit, pars = pars, ...)
    }
  } else {
    warning("This model was fit with an old version of covidseir.", call. = FALSE)
    print(x$fit, pars = pars, ...)
  }
}
