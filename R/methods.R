#' @export
#' @import methods
print.covidseir <- function(x,
                            pars = c("R0", "i0", "f_s", "start_decline", "end_decline", "phi"), ...) {
  print(x$fit, pars = pars, ...)
}

#' @export
loo.covidseir <- function(x,
                        pars = "log_lik",
                        ...,
                        save_psis = FALSE,
                        cores = getOption("mc.cores", 1L)) {
  stopifnot(length(pars) == 1L)
  LLarray <- loo::extract_log_lik(
    stanfit = x$fit,
    parameter_name = pars,
    merge_chains = FALSE
  )
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::loo.array(LLarray,
    r_eff = r_eff,
    cores = cores,
    save_psis = save_psis
  )
}
