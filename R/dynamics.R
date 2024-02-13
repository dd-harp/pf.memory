
#' Compute infection density in a cohort of humans, \eqn{z_\tau(\alpha, a)}
#'
#' @description
#' Given a function describing the FoI (\eqn{h_\tau(a)}), and a parameter
#' describing the clearance rate of infections (\eqn{r}),
#' the density of parasites of age \eqn{\alpha} in a cohort of humans of
#' age \eqn{a} is \deqn{z_\tau(\alpha, a) = e^{-r \alpha} h_\tau(a-\alpha)}
#'
#' @param alpha the age of an infection, \eqn{\alpha}
#' @inheritParams FoI
#' @param hhat scaling parameter for [FoI]
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#' @seealso [pf.memory::FoI()]
#' @export
#'
zda = function(alpha, a, FoIpar, tau=0, hhat=1, r=1/200){
  hhat*FoI(a-alpha, FoIpar, tau)*exp(-r*alpha)
}
