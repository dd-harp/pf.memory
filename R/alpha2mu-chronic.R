#' Compute expected log10 parasite densities, `mu`, as a function of the age of infection `alpha`
#' @description Compute mean log10 parasites for a model with no immunity
#'
#' @inheritParams alpha2mu
#'
#' @return mean log10 parasite densities
#' @export
#'
alpha2mu.chronic = function(alpha, W, par_Fmu){with(par_Fmu,{
  B = tildel + (tildeb-tildel)*exp(-(Sa*(alpha)))
})}

#' Set up parameters for [alpha2mu.chronic]
#'
#' @param tildeb The maximum expected log10 parasite densities
#' @param tildel The minimum expected log10 parasite densities
#' @param Sa The decline in mu with respect to alpha
#'
#' @return a [list]
#' @export
#'
par_Fmu_chronic = function(tildeb=10.3, tildel=2, Sa=0.0033){
  par = list()
  class(par) <- "chronic"
  par$tildeb=tildeb
  par$tildel=tildel
  par$Sa=Sa
  par
}
