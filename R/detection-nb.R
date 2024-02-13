
#' Detection of infection given parasitemia
#'
#' @param q the fraction of blood examined
#' @param sz the negative binomial size parameter
#' @param bvm blood volume as log10 red blood cells
#'
#' @return binary detection result
#' @export
#'
par_nbCounts = function(q=6, sz=0.31, bvm = par_lRBC_static()){
  par = list()
  class(par) <- "nb"
  par$q=q
  par$sz=sz
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
Detect.nb = function(xi, a=0, Cpar=par_nbCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  1-stats::dnbinom(0, mu=10^(xi-lRBC+q), size=sz)
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
dCounts.nb = function(hatxi, xi, a=0, Cpar=par_nbCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  stats::dnbinom(round(10^hatxi), mu=10^(xi-lRBC+q), size=sz)
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
pCounts.nb = function(hatxi, xi, a=0, Cpar=par_nbCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  pnbinom(10^hatxi, mu=10^(xi-lRBC+q), size=sz)
})}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param bins break points for binning counts data
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
binnedCounts.nb = function(xi, a, bins, Cpar=par_nbCounts()){
  if(is.null(bins)) bins=c(1:5,13)
  p0 = Detect(xi,a,Cpar)
  p1 = pCounts(bins, xi, a, Cpar)
  diff(c(1-p0,p1))/p0
}
