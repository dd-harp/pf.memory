
#' The density function for parasite densities in a simple malaria infection
#'
#' @param xi vector of quantiles for \eqn{log_{10}} parasite densities
#' @param mu the expected value for \eqn{log_{10}} parasite densities
#' @param par_mu2dens parameters to compute parasite densities as a function of mu
#' @param bvm blood volume as \eqn{log_{10}} red blood cells
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
d_mu2density = function(xi, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  UseMethod("d_mu2density", par_mu2dens)
}

#' Random generation of parasite densities from a simple malaria infection
#'
#' @param n the number of observations
#' @param mu the expected value for \eqn{log_{10}} parasite densities
#' @param par_mu2dens parameters defining parasite densities as a function of the mu
#' @param bvm blood volume as \eqn{log_{10}} red blood cells
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
r_mu2density = function(n, mu, par_mu2dens=par_mu2dens_beta() , bvm=13){
  UseMethod("r_mu2density", par_mu2dens)
}

#' The density function for parasite densities as a function of the mean
#'
#' @param xi a vector of probabilities for \eqn{log_{10}} parasite densities
#' @param mu the expected value for \eqn{log_{10}} parasite densities
#' @param par_mu2dens parameters defining parasite densities as a function of the mu
#' @param bvm blood volume as \eqn{log_{10}} red blood cells
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
p_mu2density = function(xi, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  UseMethod("p_mu2density", par_mu2dens)
}

#' The quantile function for parasite densities in a simple malaria infection
#'
#' @param xi a vector of quantiles
#' @param mu the expected value for \eqn{log_{10}} parasite densities
#' @param par_mu2dens parameters defining parasite densities as a function of the mu
#' @param bvm blood volume as \eqn{log_{10}} red blood cells
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
q_mu2density = function(xi, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  UseMethod("q_mu2density", par_mu2dens)
}
