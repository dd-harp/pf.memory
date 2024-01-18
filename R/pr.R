
#' Compute the true PR in a cohort as a function of age and exposure
#'
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
#' @examples
truePRa = function(a, FoIpar, hhat=NULL, tau=0, r=1/200){
  1-exp(-meanMoI(a, FoIpar, hhat, tau, r))
}



#' Compute the derivatives for MoI and true PR
#'
#' @param a the host age
#' @param M state variables
#' @param par the model parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return the derivatives, as a [list]
#' @export
#'
#' @examples
dpda = function(a,M,par,FoIpar){with(as.list(c(M,par)),{
  foi = FoI(a,FoIpar,tau,h)
  R = function(m){ ifelse(m==0,r,r*m/(exp(m)-1))}
  dp  = foi*(1-p) - R(m)*p
  dm  = foi - r*m
  list(c(dp, dm))
})}

#' Solve a system of differential equations to compute the true PR and the MoI
#'
#' @param h
#' @param FoIpar parameters that define an FoI function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Tmax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return
#' @export
#'
#' @examples
solve_dpda = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1){
  tms = seq(0, Tmax, by = dt)
  prms = c(h=h,r=r,FoI=FoI,tau=tau)
  inits = c(p=0,m=0)
  data.frame(lsode(inits, times=tms, dpda, prms, FoIpar=FoIpar))
}
