#' The probability of detecting parasites
#'
#' @param xi \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_detect = function(xi, bvm=13, pars=par_nb()){
  UseMethod("d_detect", pars)
}

#' The PDF for counting \eqn{x} parasites
#'
#' @param x raw parasite counts
#' @param xi mean \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_counts = function(x, xi, bvm=13, pars=par_nb()){
  UseMethod("d_counts", pars)
}

#' The CDF for counting \eqn{x} parasites
#'
#' @param x raw parasite counts
#' @param xi \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
p_counts = function(x, xi, bvm=13, pars=par_nb()){
  UseMethod("p_counts", pars)
}

#' The PDF for non-zero counts, \eqn{\log_{10}} transformed
#'
#' @param hatxi \eqn{\log_{10}} parasite counts
#' @param xi \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_nz_counts_log = function(hatxi, xi, bvm=13, pars=par_nb()){
  nz = 1-d_detect(xi, bvm, pars)
  dc = d_counts(10^hatxi, xi, bvm, pars)
  return(dc/nz)
}

#' The CDF for non-zero counts, \eqn{\log_{10}} transformed
#'
#' @param hatxi \eqn{\log_{10}} parasite counts
#' @param xi \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
p_nz_counts_log = function(hatxi, xi, bvm=13, pars=par_nb()){
  nz = 1-d_detect(xi, bvm, pars)
  dc = p_counts(10^hatxi, xi, bvm, pars)
  return(dc/nz)
}

#' The binned PDF for non-zero counts, \eqn{\log_{10}} transformed
#'
#' @param xi \eqn{\log_{10}} parasite densities
#' @param bins breakpoints for summarizing outputs as \eqn{\log_{10}} counts
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_nz_counts_log_binned = function(xi, bins=NULL, bvm=13, pars=par_nb()){
  if(is.null(bins)) bins = c(1:5,13)
  p0 = d_detect(xi, bvm, pars)
  p1 = p_counts(10^bins, xi, bvm, pars)
  diff(c(1-p0,p1))/p0
}

#' The binned PDF for non-zero counts, \eqn{\log_{10}} transformed
#'
#' @param xi \eqn{\log_{10}} parasite densities
#' @param bins breakpoints for summarizing outputs as \eqn{\log_{10}} counts
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param pars parameters that define a detection function
#'
#' @return binary detection result
#' @export
p_nz_counts_log_binned = function(xi, bins=NULL, bvm=13, pars=par_nb()){
  if(is.null(bins)) bins = c(1:5,13)
  p0 = d_detect(xi, bvm, pars)
  p1 = p_counts(10^bins, xi, bvm, pars)-p0
  p1/p0
}
