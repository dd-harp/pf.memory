#' Compute log10 of the red blood cell population
#'
#' @param a the age of a host cohort
#' @param par a [list] of parameters
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
log10RBC = function(a,par){
  UseMethod("log10RBC", par)
}

#' Compute log10 of the red blood cell population
#' @description
#' A static model
#'
#' @inheritParams log10RBC
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
log10RBC.static = function(a,par){par$lRBCmax}

#' Set up parameters for [log10RBC.static]
#'
#' @param lRBCmax the maximum
#'
#' @return a [list]
#' @export
#'
par_lRBC_static = function(lRBCmax=13){
  par = list()
  class(par) <- "static"
  par$lRBCmax = lRBCmax
  par
}

