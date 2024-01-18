
#' Compute the density of parasite infections in a host cohort
#' @description
#' Solutions to the time and age version of the Ross-Macdonald model.
#' @param alpha a parasites age of infection
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#'
#' @return
#' @export
#'
#' @examples
zda = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  FoI(a-alpha, FoIpar,tau,hhat)*exp(-r*alpha)
}
