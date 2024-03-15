
#' Compute the density function for AoI
#'
#' @description
#' The density of the AoI is given by
#' \deqn{f_A(\alpha | a, \tau, h) = \frac{\int_0^a z(\alpha, a)}{m_\tau(a)}}
#'
#' @inheritParams zda
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
dAoI = function(alpha, a, FoIpar, tau=0, hhat=1, r=1/200){
  dAoIcompute = function(alpha, a, FoIpar, tau, hhat, r){
    zda(alpha, a, FoIpar, tau, hhat, r)/meanMoI(a, FoIpar,  tau, hhat, r)
  }

  if(length(alpha)==1){
    return(dAoIcompute(alpha, a, FoIpar, tau, hhat, r))
  }else{
    return(sapply(alpha, dAoIcompute, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r))
  }
}

#' Compute the distribution function for AoI
#'
#' The distribution function for the AoI is given by
#' \deqn{F_A(\alpha | a, \tau, h) = \int_0^\alpha f_A(\alpha, a, \tau | h) d \alpha}
#'
#' @inheritParams zda
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
pAoI = function(alpha, a, FoIpar, tau=0, hhat=1, r=1/200){
  pAoIfunction = function(alpha, a, FoIpar, tau, hhat, r){
    stats::integrate(dAoI,0,alpha,a=a,FoIpar=FoIpar,tau=tau, hhat=hhat,r=r)$value
  }
  if(length(alpha)==1) {return(pAoIfunction(alpha, a, FoIpar, tau, hhat, r))} else{
    return(sapply(alpha, pAoIfunction, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r))}
}

#' Random numbers for the AoI
#'
#' @description
#' Draw random numbers for the AoI from a cohort, \eqn{\hat A_\tau(a)}
#'
#' @param N the number of observations
#' @inheritParams FoI
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of the AoI to return
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
rAoI = function(N, a, FoIpar, tau=0, hhat=1, r=1/200, alphamin=0){
  stopifnot(N>0)
  alpha = alphamin:a
  scdf = pAoI(alpha, a, FoIpar, tau, hhat, r)
  pdf = diff(scdf)
  pdf / sum(pdf)
  sample(alpha[-length(alpha)], N, replace=T, prob=pdf)
}

#' Compute the moments for the AoI density function for a cohort of age a
#'
#' @inheritParams meanMoI
#' @param n the moment desired
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
momentAoI = function(a, FoIpar, tau=0, hhat=1, r=1/200, n=1){

  fAda = function(a, FoIpar, tau, hhat, r0, n){
    ff = function(alpha, a, FoIpar, tau, hhat, r, n){
      alpha^n*zda(alpha, a, FoIpar, tau, hhat, r)
    }
    m =  meanMoI(a,FoIpar,tau,hhat,r)
    stats::integrate(ff, 0, a,a=a,FoIpar=FoIpar, tau=tau, hhat=hhat, r=r, n=n)$value/m
  }

  if(length(a)==1){return(fAda(a, FoIpar, tau, hhat, r, n))} else{
    sapply(a, fAda,FoIpar=FoIpar,tau=tau,hhat=hhat,r=r, n=n)}
}
