
#' Compute infection density in a cohort of humans, \eqn{z_\tau(\alpha, a |h)}
#'
#' @description
#' Given a function describing the FoI (\eqn{h_\tau(a)}), and a parameter
#' describing the clearance rate of infections (\eqn{r}),
#' the density of parasite clones of age \eqn{\alpha} distributed among a cohort of humans of
#' age \eqn{a} is \deqn{z_\tau(\alpha, a) = e^{-r \alpha} h_\tau(a-\alpha)}
#'
#' @param alpha the age of an infection, \eqn{\alpha}
#' @inheritParams FoI
#' @param hhat scaling parameter for [FoI]
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#' @seealso [ramp.pf.infection::FoI()]
#' @export
#'
zda = function(alpha, a, FoIpar, tau=0, hhat=1, r=1/200){
  hhat*FoI(a-alpha, FoIpar, tau)*exp(-r*alpha)
}

#' The mean MoI in a host cohort of age \eqn{a}
#'
#' @description
#' The mean multiplicity of infection (MoI) is \deqn{m_\tau(a|h) = \int_0^a z_\tau(\alpha, a|h) d \alpha}
#'
#' @inheritParams FoI
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
meanMoI = function(a, FoIpar, tau=0, hhat=1, r=1/200){
  moif = function(a, FoIpar, tau, hhat,r){
    stats::integrate(zda, 0, a, a=a, FoIpar=FoIpar, tau=tau, hhat=hhat, r=r)$value
  }
  if(length(a)==1){return(moif(a,FoIpar,tau,hhat,r))} else{
    (return(sapply(a,moif,FoIpar=FoIpar,tau=tau,hhat=hhat,r=r)))}
}


#' Compute the true PR in a cohort as a function of age and exposure
#'
#' @description
#' The true PR is \deqn{p_\tau(a|h) = 1 - e^{-m_\tau(a|h)}}
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
truePR = function(a, FoIpar, tau=0, hhat=1, r=1/200){
  1-exp(-meanMoI(a, FoIpar, tau, hhat, r))
}

