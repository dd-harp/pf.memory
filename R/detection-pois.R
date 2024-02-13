
#' Detection of infection given parasitemia

#' @param q the fraction of blood volume sampled
#' @param bvm blood volume as log10 red blood cells
#'
#' @return par a [list]
#' @export
par_poisCounts = function(q=6, bvm = par_lRBC_static()){
  par = list()
  class(par) <- "pois"
  par$q=q
  par$bvm = bvm
  par
}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
Detect.pois = function(xi, a=0, Cpar=par_poisCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  1-stats::dpois(0, 10^(xi-lRBC+q))
})}



#' Detection of infection given parasitemia
#'
#' @param hatxi mean log10 parasite densities
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
dCounts.pois = function(hatxi, xi, a=0, Cpar=par_poisCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  stats::dpois(round(10^hatxi), 10^(xi-lRBC+q))
})}



#' Detection of infection given parasitemia
#'
#' @param hatxi mean log10 parasite counts
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
pCounts.pois = function(hatxi, xi, a=0, Cpar=par_poisCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  ppois(10^hatxi, 10^(xi-lRBC+q))
})}



#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param bins breakpoints for binning data
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
binnedCounts.pois = function(xi, a, bins, Cpar=par_poisCounts()){
  if(is.null(bins)) bins=c(1:5,13)
  p0 = Detect(xi,a,Cpar)
  p1 = pCounts(bins, xi, a, Cpar)
  diff(c(1-p0,p1))/p0
}
