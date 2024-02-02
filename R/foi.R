
#' Trace Function - Force of Infection
#'
#' @param a host age
#' @param par a formatted list
#' @param tau cohort birthday
#' @param hhat overrides the value of hbar in par
#'
#' @return [numeric]
#' @export
FoI = function(a, par, tau=0, hhat=NULL){with(par,{
  ifelse(is.null(hhat), hbar, hhat)*
    ageFoI(a,agePar)*
    seasonalFoI(a+tau,seasonPar)*
    trendFoI(a+tau,trendPar)
})}

#' Trace Function - Force of Infection
#'
#' @param par a formatted list
#' @param maxAge age in years to normalize over
#'
#' @return [numeric]
#' @export
make_FoI = function(par, maxAge = 50){
  norm = integrate(FoI, 0, maxAge*365, tau=0, par=par)$value
  par$hbar = maxAge*par$hbar/norm
  ff = function(a, tau=0){
    FoI(a, par, tau)
  }
  return(ff)
}
