
#' The youngest of N infections, density function
#'
#' @param N the number of infections
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length a + 1
#' @export
#'
dAoYN = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  cdf = pAoYN(N, a, FoIpar, hhat, tau, r)
  pdf = diff(cdf)
  pdf/sum(pdf)
}

#' The youngest of N infections, distribution function
#'
#' @param N the MoI
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length a + 1
#' @export
#'
pAoYN = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  alpha = 0:a
  py = pAoI(alpha, a, FoIpar, hhat, tau, r)
  1-(1-py)^N
}

#' The youngest of N infections, random numbers
#'
#' @param R the number of observations
#' @param N the number of infections per person
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of the AoI to return
#'
#' @return a [numeric] vector of length R
#' @export
#'
rAoYN = function(R, N, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=0){
  matrix(rAoI(R*N,a,FoIpar,hhat,tau,r,alphamin), nrow=N, ncol=R)
}
