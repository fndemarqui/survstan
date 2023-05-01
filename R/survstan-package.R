#' The 'survstan' package.
#'
#' @description The aim of the R package survstan is to provide a toolkit for fitting survival models using Stan. The R package survstan can be used to fit right-censored survival data under independent censoring. The implemented models allow the fitting of survival data in the presence/absence of covariates. All inferential procedures are currently based on the maximum likelihood (ML) approach.
#'
#' @docType package
#' @name survstan-package
#' @aliases survstan
#' @useDynLib survstan, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2022). RStan: the R interface to Stan. R package version 2.21.5. https://mc-stan.org
#'
NULL
