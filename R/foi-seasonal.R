
#' Add a seasonal pattern to the FoI function
#'
#' @param t the time
#' @param par a [list]
#'
#' @return [numeric]
#' @export
seasonalFoI = function(t, par){
  UseMethod("seasonalFoI", par)
}

#' The "no seasonality" pattern
#'
#' @inheritParams seasonalFoI
#'
#' @return [numeric]
#' @export
seasonalFoI.flat= function(t, par){1}

#' Return a list that dispatches seasonalFoI.flat
#'
#' @return a [list]
#' @export
par_flatSeason = function(){
  pars = list()
  class(pars) <- "flat"
  pars
}

#' A generalized sinusoidal seasonal pattern
#'
#' @inheritParams seasonalFoI
#'
#' @return [numeric]
#' @export
seasonalFoI.sin = function(t, par=par_sinSeason()){with(par,{
  foi = shift + (1+lift+sin(2*pi*(t+phase)/365))^pwr
  pmax(0, foi)/norm
})}

#' Return a list to configure and dispatch seasonalFoI.sin
#'
#' @param phase shift the phase / time of year
#' @param pwr raise the seasonal pattern to a power
#' @param lift the amplitude of the
#' @param shift shift everything up or down
#'
#' @return a [list]
#' @export
par_sinSeason = function(phase=0, pwr=1, lift=0, shift=0){
  pars = list()
  class(pars) <- "sin"
  pars$phase=phase
  pars$pwr=pwr
  pars$lift=lift
  pars$shift=shift
  pars$norm = 1
  tot = stats::integrate(seasonalFoI.sin, 0, 365, par=pars)$value
  pars$norm = tot/365
  pars
}

#' A generalized sinusoidal seasonal pattern, exponentiated
#'
#' @param t the time
#' @param par a list
#'
#' @return [numeric]
#' @export
seasonalFoI.exp= function(t, par=par_expSeason()){with(par,{
  exp(shift + (1+lift+sin(2*pi*(t+phase)/365))^pwr)/norm
})}

#' Return a list to dispatch seasonalFoI.exp
#'
#' @param phase shift the phase / time of year
#' @param pwr raise the seasonal pattern to a power
#' @param lift the amplitude of the
#' @param shift shift everything up or down
#'
#' @return a [list]
#' @export
par_expSeason = function(phase=0, pwr=1, lift=0, shift=0){
  pars = par_sinSeason(phase,pwr,lift,shift)
  pars$norm = 1
  tot = stats::integrate(seasonalFoI.exp, 0, 365, par=pars)$value
  pars$norm = tot/365
  class(pars) <- "exp"
  pars
}
