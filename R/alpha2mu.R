#' Compute mean, expected parasite densities `mu` as a function of the age of infection `alpha`
#'
#' @param alpha a parasite's age of infection
#' @param W the immune tracking variables
#' @param par a [list] that defines a model
#'
#' @return mean log10 parasite densities
#' @export
#'
#' @examples
alpha2mu = function(alpha, W, par){
  UseMethod("alpha2mu", par)
}

#' Compute mean, expected parasite densities `mu` for multiple values of the age of infection `alpha`
#'
#' @param alpha a parasite's age of infection
#' @param W immune tracking variables
#' @param par a [list] that defines a model
#'
#' @return mean log10 parasite densities
#' @export
#'
#' @examples
salpha2mu = function(alpha, W, par){
  sapply(alpha, alpha2mu, W=W, par=par)
}

#' Compute expected log10 parasite densities, `mu`, as a function of the age of infection `alpha`
#' @description Compute mean log10 parasites for a model with no immunity
#'
#' @inheritParams alpha2mu
#'
#' @return mean log10 parasite densities
#' @export
#'
#' @examples
alpha2mu.0 = function(alpha, W, par){with(par,{
  B = tildel + (tildeb-tildel)*exp(-(Sa*(alpha-D)))
  ix = which(alpha<D)
  if(length(ix>0)) B[ix] = tildel+(tildeb-tildel)*alpha[ix]/D
  #  ix = which(alpha<liver)
  #  if(length(ix>0)) B[ix] = sqrt(-1)
  B
})}

#' Set up parameters for alpha2mu.0
#'
#' @param D
#' @param liver
#' @param tildeb
#' @param tildel
#' @param Sa
#'
#' @return mean log10 parasite densities
#' @export
#'
#' @examples
par_alpha2mu.0 = function(D=20, liver=7, tildeb=10.3, tildel=2, Sa=0.0033){
  par = list()
  class(par) <- "0"
  par$D=D
  par$liver=liver
  par$tildeb=tildeb
  par$tildel=tildel
  par$Sa=Sa
  par
}


#' Compute expected log10 parasite densities, `mu`, as a function of the age of infection `alpha`
#' @description Compute mean log10 parasites for a model with no immunity
#'
#' @inheritParams alpha2mu
#'
#' @return mean log10 parasite densities
#' @export
#'
#' @examples
alpha2mu.W = function(alpha, W, par){with(par,{
  B = tildel + (tildeb-tildel)*exp(-Sa*(alpha-D)-Sw*W)
  ix = which(alpha<=D)
  if(length(ix>0)) B[ix] = tildel+(tildeb-tildel)*exp(-Sw*W)*alpha[ix]/D
  # ix = which(alpha<=liver)
  #  if(length(ix>0)) B[ix] = sqrt(-1)
  B
})}

#' Set up parameters for alpha2mu.0
#'
#' @param D The age of infection (in days) when parasite densities peak
#' @param liver The age of infection (in days) when parasites emerge from the liver
#' @param tildeb The maximum expected log10 parasite densities
#' @param tildel The minimum expected log10 parasite densities
#' @param Sa The decline in mu with respect to alpha
#' @param Sw The decline in mu with respect to immunity
#'
#' @return mean log10 parasite densities
#' @export
#'
#' @examples
par_alpha2mu.W = function(D=20, liver=7, tildeb=10.3, tildel=2, Sa=0.0033, Sw=0.001){
  par = list()
  class(par) <- "W"
  par$D=D
  par$liver=liver
  par$tildeb=tildeb
  par$tildel=tildel
  par$Sa=Sa
  par$Sw=Sw
  par
}
