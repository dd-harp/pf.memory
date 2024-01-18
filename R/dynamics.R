
#' Title
#' @description
#'
#' @param alpha a parasites age of infection
#' @param a the host cohort age
#' @param FoIpar
#' @param hhat
#' @param tau
#' @param r
#'
#' @return
#' @export
#'
#' @examples
zda = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  FoI(a-alpha, FoIpar,tau,hhat)*exp(-r*alpha)
}
