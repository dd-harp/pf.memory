
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

make_FoI = function(par){
  par1 = par
  par1$agePar = par_flatAge()
  norm1 = integrate(FoI, 0, 365, par=par1)$value
  par2 = par
  par2$seasonalFoI = par_flatSeason()
  par2$trendFoI = par_flatTrend()
  norm2 = integrate(FoI, 0, 365, par=par2)$value
  par$hbar = 1/norm1/norm2
  ff = function(a, tau=0){
    FoI(a, par, tau)
  }
  return(ff)
}

make_FoI = function(par){
  norm = integrate(FoI, 0, 55*365, par=par1)$value
  par$hbar = 1/norm/55
  ff = function(a, tau=0){
    FoI(a, par, tau)
  }
  return(ff)
}


f1 = make_FoI(foiP3)
integrate(f1, 0, 365)$value
