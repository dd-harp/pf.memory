
#' Compute immune tracking variables as a function of host age and exposure
#' @description
#' The function dispatches on `class(par)`
#'
#' @param a age of a host cohort
#' @param FoIpar parameters that define an FoI function
#' @param tau cohort birthday
#' @param hhat overrides the value of hbar in par
#' @param par parameters in a [list]
#'
#' @return a [numeric] vector of length(a)
#' @export
Wda = function(a, FoIpar, tau=0, hhat=1, par=par_Wda_none()){
  UseMethod("Wda", par)
}

#' Compute immune tracking variables as a function of host age and exposure
#'
#' @inheritParams Wda
#'
#' @return a [numeric] vector of 0's of length(a)
#' @export
Wda.none = function(a, FoIpar, tau=0, hhat=1, par=par_Wda_none()){
  0*a
}

#' Make a parameter set for [Wda.none]
#'
#' @return a [list]
#' @export
par_Wda_none = function(){
  par = list()
  class(par) <- "none"
  return(par)
}


#' Compute immune tracking variables as a function of host age and exposure
#'
#' @inheritParams Wda
#'
#' @return a [numeric] vector of length(a)
#' @export
Wda.delta = function(a, FoIpar, tau=0, hhat=1, par=par_Wda_delta()){with(par,{
  Wd = function(a,FoIpar,tau,hhat,delta){
    ff = function(s,a,FoIpar,tau,hhat,delta){
      hhat*FoI(a-s,FoIpar,tau)*exp(-delta*(a-s))
    }
    integrate(ff,0,a,a=a,FoIpar=FoIpar,tau=tau,hhat=hhat,delta=delta)$value
  }
  if(length(a)==1){
    return(Wd(a,FoIpar,tau,hhat,delta))
  } else {
    return(sapply(a,Wd,FoIpar=FoIpar,tau=tau,hhat=hhat,delta=delta))
  }
})}


#' Make a parameter set for [Wda.none]
#'
#' @param delta a decay rate
#'
#' @return a [list]
#' @export
par_Wda_delta = function(delta=0.001){
  par = list()
  class(par) <- "delta"
  par$delta=delta
  return(par)
}
