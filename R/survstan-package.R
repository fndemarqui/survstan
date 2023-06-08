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
#' @import survival
#' @importFrom dplyr %>% all_of bind_cols desc select
#' @importFrom ggplot2 aes ggplot geom_abline geom_point geom_jitter geom_smooth position_jitter labs xlab xlim ylim
#' @importFrom Rdpack reprompt
#' @importFrom rlang .data
#' @importFrom rstan sampling
#'
#' @references
#'
#' \insertRef{rstan}{survstan}
#'
#' \insertRef{bookLawless}{survstan}
#'
#' \insertRef{1983Bennett}{survstan}
#'
#' \insertRef{2000Chen}{survstan}
#'
#' \insertRef{2021_Demarqui}{survstan}
#'
NULL

#'
NULL
