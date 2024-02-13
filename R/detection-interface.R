#' Detection of infection given parasitemia
#'
#' @param xi mean \eqn{\log_{10}} parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
Detect = function(xi, a=0, Cpar=par_nbCounts()){
  UseMethod("Detect", Cpar)
}

#' density function of counts
#'
#' @param hatxi mean \eqn{\log_{10}} parasite counts
#' @param xi mean \eqn{\log_{10}} parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
dCounts = function(hatxi, xi, a=0, Cpar=par_nbCounts()){
  UseMethod("dCounts", Cpar)
}

#' cumulative density function of counts
#'
#' @param hatxi mean log10 parasite densities
#' @param xi mean \eqn{\log_{10}} parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
pCounts = function(hatxi, xi, a=0, Cpar=par_nbCounts()){
  UseMethod("pCounts", Cpar)
}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param bins breakpoints for summarizing outputs
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
binnedCounts = function(xi, a, bins, Cpar=par_nbCounts()){
  UseMethod("binnedCounts", Cpar)
}

