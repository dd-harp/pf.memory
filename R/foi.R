
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
