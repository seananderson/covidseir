#' The 'covidseir' package.
#'
#' @description Bayesian SEIR Modelling for Multivariate COVID-19 Case Data
#'
#' @docType package
#' @name covidseir-package
#' @aliases covidseir
#' @useDynLib covidseir, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "value", "Var2", "iterations", "date", "day", "data_type",
    "lambda_d", "y_rep", "variable", "R0", "f_s", "f_seg", "mu", "phi",
    "Var3", "Var2", ".iteration"
  ))
}
