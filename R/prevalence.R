
#' Compute the true PR in a cohort as a function of age and exposure
#'
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
truePRa = function(a, FoIpar, tau=0, hhat=1, r=1/200){
  1-exp(-meanMoI(a, FoIpar, tau, hhat, r))
}


