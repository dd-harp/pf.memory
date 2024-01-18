
#' Compute immune tracking variables as a function of host age and exposure
#' @description
#' The function dispatches on `class(par)`
#'
#' @param a age of a host cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat overrides the value of hbar in par
#' @param tau cohort birthday
#' @param par parameters in a [list]
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
#' @examples
Wda = function(a, FoIpar, hhat=NULL, tau=0, par=par_Wda.0()){
  UseMethod("Wda", par)
}

#' Compute immune tracking variables as a function of host age and exposure
#'
#' @inheritParams Wda
#'
#' @return a [numeric] vector of 0's of length(a)
#' @export
#'
#' @examples
Wda.0 = function(a, FoIpar, hhat=NULL, tau=0,  par=par_Wda.0()){
  0*a
}

#' Make a parameter set for [Wda.0]
#'
#' @return a [list]
#' @export
#'
#' @examples
par_Wda.0 = function(){
  par = list()
  class(par) <- "0"
  return(par)
}


#' Compute immune tracking variables as a function of host age and exposure
#'
#' @inheritParams Wda
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
#' @examples
Wda.delta = function(a, FoIpar, hhat=NULL, tau=0, par=par_Wda.0()){with(par,{
  Wd = function(a,FoIpar,hhat,tau,delta){
    ff = function(s,a,FoIpar,hhat,tau,delta)  FoI(a-s,FoIpar,tau,hhat)*exp(-delta*(a-s))
    integrate(ff,0,a,a=a,FoIpar=FoIpar,hhat=hhat,tau=tau,delta=delta)$value
  }
  if(length(a)==1){
    return(Wd(a,FoIpar,hhat,tau,delta))
  } else {
    return(sapply(a,Wd,FoIpar=FoIpar,hhat=hhat,tau=tau,delta=delta))
  }
})}


#' Make a parameter set for [Wda.delta]
#'
#' @param delta a decay rate
#'
#' @return a [list]
#' @export
#'
#' @examples
par_Wda.delta = function(delta=0.001){
  par = list()
  class(par) <- "delta"
  par$delta=delta
  return(par)
}
